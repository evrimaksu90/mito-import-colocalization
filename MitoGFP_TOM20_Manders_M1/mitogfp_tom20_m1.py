#!/usr/bin/env python3
"""
MitoGFP-TOM20 colocalization (Manders' M1) in transfected cells, single-cell 3D analysis.

Input : Leica .lif file with 3-channel z-stacks (DAPI, MitoGFP, TOM20).
Output: per-cell table, table for GraphPad, QC images.

M1 = sum(GFP in TOM20+ voxels) / sum(GFP in the cell)    (no threshold on GFP)

See README.md for the full description of every step.
"""
import argparse, os, json
import numpy as np, pandas as pd
from scipy import ndimage as ndi
from skimage import filters, morphology, measure
from skimage.segmentation import watershed
from readlif.reader import LifFile
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt


def parse_groups(s):
    """'CTRL:Control,PAT1:Patient 1' -> ordered list of (prefix, group)."""
    return [tuple(x.split(":")) for x in s.split(",")]


def group_of(name, rules):
    for prefix, g in rules:
        if name.startswith(prefix):
            return g
    return None


def load_series(im, nch):
    """Return float32 array (C, Z, Y, X)."""
    nz = im.dims.z
    return np.stack([np.stack([np.asarray(im.get_frame(z=z, c=c)) for z in range(nz)])
                     for c in range(nch)]).astype(np.float32)


def analyze_series(st, px, ch_dapi, ch_gfp, ch_tom, a):
    D, G, T = st[ch_dapi], st[ch_gfp], st[ch_tom]
    # --- background: channel median (cells are sparse); noise SD from MAD ---
    gm = np.median(G); gs = max(1.4826 * np.median(np.abs(G - gm)), 0.5); tm = np.median(T)
    Gs = np.clip(ndi.gaussian_filter(G, (0, 1, 1)) - gm, 0, None)
    Ts = np.clip(ndi.gaussian_filter(T, (0, 1, 1)) - tm, 0, None)

    # --- 1. GFP-positive regions (2D, max projection) ---
    gmip = ndi.gaussian_filter(G.max(0) - gm, 1 / px)
    thr = max(filters.threshold_triangle(gmip), 4 * gs)
    m = ndi.binary_fill_holes(morphology.binary_closing(gmip > thr, morphology.disk(int(round(1.5 / px)))))
    m = morphology.remove_small_objects(m, int(a.min_area / px ** 2))
    lab0 = measure.label(m)

    # --- 2. split regions with >1 nucleus by nucleus-seeded watershed ---
    Dm = ndi.gaussian_filter(D.max(0), 1.5 / px); Dm = Dm - ndi.uniform_filter(Dm, int(40 / px))
    nm = Dm > filters.threshold_otsu(Dm)
    nm = ndi.binary_opening(morphology.remove_small_objects(nm, int(40 / px ** 2)), iterations=2)
    nl = measure.label(nm)
    sig = ndi.gaussian_filter((G.max(0) - gm) / max(np.percentile(G.max(0), 99.5), 1)
                              + (T.max(0) - tm) / max(np.percentile(T.max(0), 99.5), 1), 2 / px)
    lab = np.zeros_like(lab0); nxt = 1
    for rp in measure.regionprops(lab0):
        cm = lab0 == rp.label
        ids = [q.label for q in measure.regionprops(nl) if cm[q.coords[:, 0], q.coords[:, 1]].mean() > 0.5]
        if len(ids) >= 2:
            markers = np.zeros_like(lab0)
            for j, q in enumerate(ids, 1):
                markers[(nl == q) & cm] = j
            ws = watershed(-sig, markers, mask=cm)
            for j in range(1, len(ids) + 1):
                part = ws == j
                if part.sum() * px ** 2 >= a.min_area:
                    lab[part] = nxt; nxt += 1
        else:
            lab[cm] = nxt; nxt += 1

    # --- 3. GFP autofluorescence of non-transfected (TOM20+, GFP-) cells in the same field ---
    excl = morphology.binary_dilation(m, morphology.disk(int(3 / px)))
    Tm = ndi.gaussian_filter(T.max(0) - tm, 1 / px)
    tcell = (Tm > filters.threshold_otsu(Tm)) & ~excl
    if tcell.sum() < 300:
        tcell = (Tm > filters.threshold_triangle(Tm)) & ~excl
    AF = np.percentile(G.max(0)[tcell] - gm, 95) if tcell.sum() > 300 else np.nan

    # --- 4. mitochondrial compartment: image-wise Otsu threshold on TOM20 ---
    TTHR = filters.threshold_otsu(Ts[Ts > 2])

    rows = []
    for rp in measure.regionprops(lab):
        cm = lab == rp.label
        g3, t3 = Gs[:, cm], Ts[:, cm]
        zs = g3.sum(1); keep = zs >= a.focus_frac * zs.max()          # in-focus planes
        g, t = g3[keep].ravel(), t3[keep].ravel()
        if g.sum() <= 0:
            continue
        tp = t > TTHR
        M1 = g[tp].sum() / g.sum()
        rows.append(dict(cell=rp.label, area_um2=rp.area * px ** 2, GFP_mean=float(g.mean()),
                         GFP_p95=float(np.percentile(G.max(0)[cm] - gm, 95)), AF_p95=AF, M1=M1))
    return rows, lab


def qc_image(st, ch_gfp, ch_tom, lab, df_img, path, title):
    G, T = st[ch_gfp].max(0), st[ch_tom].max(0)
    rgb = np.zeros(G.shape + (3,))
    rgb[..., 0] = np.clip((T - np.median(T)) / np.percentile(T, 99.8), 0, 1)
    rgb[..., 1] = np.clip((G - np.median(G)) / max(np.percentile(G, 99.8), 1), 0, 1)
    fig, ax = plt.subplots(figsize=(10, 10)); ax.imshow(rgb)
    inc = dict(zip(df_img.cell, df_img.included)); m1 = dict(zip(df_img.cell, df_img.M1))
    for rp in measure.regionprops(lab):
        ok = inc.get(rp.label, False)
        y0, x0, y1, x1 = rp.bbox; y0, x0 = max(y0 - 2, 0), max(x0 - 2, 0); y1, x1 = min(y1 + 2, lab.shape[0]), min(x1 + 2, lab.shape[1])
        ax.contour(np.arange(x0, x1), np.arange(y0, y1), (lab[y0:y1, x0:x1] == rp.label).astype(float), [0.5],
                   colors="w" if ok else "#888", linewidths=1.0 if ok else 0.6, linestyles="-" if ok else ":")
        txt = f"#{rp.label}" + (f"\n{100 * m1[rp.label]:.0f}%" if ok and rp.label in m1 else "")
        ax.text(rp.centroid[1], rp.centroid[0], txt, color="yellow" if ok else "#aaa", fontsize=10, weight="bold", ha="center")
    ax.set_title(f"{title} - white: included (M1 %), grey dotted: excluded"); ax.axis("off")
    plt.tight_layout(); plt.savefig(path, dpi=70); plt.close()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("lif", help="Leica .lif file")
    p.add_argument("-o", "--outdir", default="results_mitogfp_tom20")
    p.add_argument("--groups", required=True,
                   help='series-name prefix -> group, comma separated, e.g. "CTRL:Control,PAT1:Patient 1"')
    p.add_argument("--exclude", default="", help="comma-separated series names to skip")
    p.add_argument("--ch-dapi", type=int, default=0); p.add_argument("--ch-gfp", type=int, default=1)
    p.add_argument("--ch-tom20", type=int, default=2)
    p.add_argument("--sar", type=float, default=3.0, help="min. GFP p95 / autofluorescence p95 to call a cell transfected")
    p.add_argument("--min-area", type=float, default=150.0, help="min. region area (um2)")
    p.add_argument("--max-area", type=float, default=3000.0, help="max. cell area (um2); larger = merged clusters")
    p.add_argument("--focus-frac", type=float, default=0.4, help="in-focus planes: GFP sum >= this fraction of max")
    p.add_argument("--no-qc", action="store_true", help="skip QC images")
    a = p.parse_args()

    rules = parse_groups(a.groups); skip = set(x for x in a.exclude.split(",") if x)
    os.makedirs(a.outdir, exist_ok=True); qcdir = os.path.join(a.outdir, "qc"); os.makedirs(qcdir, exist_ok=True)
    lif = LifFile(a.lif); all_rows = []; cache = {}
    for im in lif.get_iter_image():
        n = im.name; g = group_of(n, rules)
        if g is None or n in skip:
            print(f"skip  {n}"); continue
        px = 1.0 / im.scale[0]                                  # um / pixel
        st = load_series(im, im.channels)
        rows, lab = analyze_series(st, px, a.ch_dapi, a.ch_gfp, a.ch_tom20, a)
        for r in rows:
            r.update(group=g, series=n)
        all_rows += rows
        if not a.no_qc:
            cache[n] = (im, lab)
        print(f"done  {n}: {len(rows)} GFP+ regions", flush=True)

    df = pd.DataFrame(all_rows)
    df["AF_p95"] = df.groupby("group").AF_p95.transform(lambda s: s.fillna(s.median()))
    df["GFP_to_AF"] = df.GFP_p95 / df.AF_p95
    df["included"] = (df.GFP_to_AF >= a.sar) & (df.area_um2 <= a.max_area)
    df["exclusion_reason"] = np.where(df.GFP_to_AF < a.sar, f"GFP < {a.sar}x autofluorescence",
                                      np.where(df.area_um2 > a.max_area, f"> {a.max_area} um2 (merged)", ""))
    cols = ["group", "series", "cell", "included", "exclusion_reason", "M1", "GFP_mean", "GFP_to_AF", "area_um2"]
    df[cols].to_csv(os.path.join(a.outdir, "cells.csv"), index=False)

    s = df[df.included]; groups = [g for _, g in rules if g in set(s.group)]
    groups = list(dict.fromkeys(groups))
    pd.concat([pd.Series((100 * s[s.group == g].M1).round(2).values, name=g) for g in groups], axis=1) \
      .to_csv(os.path.join(a.outdir, "M1_for_GraphPad.csv"), index=False)
    print("included cells per group:", s.group.value_counts().reindex(groups).to_dict())

    if not a.no_qc:
        for n, (im, lab) in cache.items():
            st = load_series(im, im.channels)
            qc_image(st, a.ch_gfp, a.ch_tom20, lab, df[df.series == n], os.path.join(qcdir, f"{n}_QC.png"), n)
    json.dump(vars(a), open(os.path.join(a.outdir, "parameters.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
