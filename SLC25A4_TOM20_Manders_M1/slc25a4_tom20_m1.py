#!/usr/bin/env python3
"""
SLC25A4-TOM20 colocalization (Manders' M1), single-cell 3D analysis of endogenous staining.

Input : folder with Zeiss .czi z-stacks (SLC25A4, DAPI, TOM20).
Output: per-cell table (with QC flags), table for GraphPad, QC images.

M1 = sum(SLC25A4 in TOM20+ voxels) / sum(SLC25A4 in the cell)   (no threshold on SLC25A4;
     nuclear voxels and artefact masks excluded)

See README.md for the full description of every step.
"""
import argparse, os, glob, json, re
import numpy as np, pandas as pd
import czifile
from scipy import ndimage as ndi
from skimage import filters, morphology, measure
from skimage.segmentation import watershed
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt


def parse_groups(s):
    return [tuple(x.split(":")) for x in s.split(",")]


def group_of(name, rules):
    for prefix, g in rules:
        if name.startswith(prefix):
            return g
    return None


def read_czi(path, binning):
    """Return (C,Z,Y,X) float32 after binning, and pixel size in um."""
    with czifile.CziFile(path) as c:
        x = c.asarray()[0, 0, :, :, :, :, 0].astype(np.float32)
        md = c.metadata()
    m = re.search(r'<Distance Id="X">\s*<Value>([^<]+)', md)
    px = float(m.group(1)) * 1e6 if m else 0.0988
    if binning > 1:
        b = binning; Y, X = x.shape[2] // b * b, x.shape[3] // b * b
        x = x[:, :, :Y, :X].reshape(x.shape[0], x.shape[1], Y // b, b, X // b, b).mean((3, 5))
    return x, px * binning


def bgcorr(img):
    """Remove additive, column-/row-wise background (e.g. vertical stripe of the scan field) per z-plane."""
    out = np.empty_like(img)
    for z in range(img.shape[0]):
        p = ndi.gaussian_filter1d(np.percentile(img[z], 5, axis=0), 15)
        r = img[z] - p[None, :]
        q = ndi.gaussian_filter1d(np.percentile(r, 5, axis=1), 15)
        out[z] = r - q[:, None]
    return out


def artefact_mask(S, T, cellmask, nucmask, px, q=98.5, open_um=0.7, min_area=1.5, max_area=200, min_sol=0.75):
    """Bright, thick/round objects inconsistent with mitochondrial morphology (aggregates, debris, dead cells):
    pixels above the q-th percentile of the summed normalised signal (or the 99th of TOM20) -> morphological
    opening with a ~1.4 um disk (removes thin mitochondrial tubules) -> compact objects of 1.5-200 um2,
    outside nuclei -> dilated by 0.8 um."""
    ref = cellmask if cellmask.sum() > 1000 else np.ones_like(cellmask)
    Sm = ndi.gaussian_filter(S.max(0), 0.3 / px); Tm = ndi.gaussian_filter(T.max(0), 0.3 / px)
    I = Sm / np.percentile(Sm[ref], 99.5) + Tm / np.percentile(Tm[ref], 99.5)
    Tn = Tm / np.percentile(Tm[ref], 99.5)
    b = ((I > np.percentile(I[ref], q)) | (Tn > np.percentile(Tn[ref], 99.0))) & ~nucmask
    b = morphology.binary_opening(b, morphology.disk(int(round(open_um / px))))
    lab = measure.label(b); keep = np.zeros_like(b)
    for rp in measure.regionprops(lab):
        if min_area <= rp.area * px ** 2 <= max_area and rp.solidity >= min_sol:
            keep[lab == rp.label] = True
    return morphology.binary_dilation(keep, morphology.disk(int(round(0.8 / px))))


def analyze_image(x, px, a):
    S = bgcorr(x[a.ch_slc]); Dp = bgcorr(x[a.ch_dapi]); T = bgcorr(x[a.ch_tom20])
    Ss = np.clip(ndi.gaussian_filter(S, (0, 1, 1)), 0, None); Ts = np.clip(ndi.gaussian_filter(T, (0, 1, 1)), 0, None)
    # nuclei (2D) and nuclear voxels (3D)
    dm = ndi.gaussian_filter(Dp.max(0), 1.0 / px)
    nucl = dm > filters.threshold_triangle(dm)
    nucl = morphology.remove_small_objects(ndi.binary_fill_holes(ndi.binary_opening(nucl, iterations=2)), int(30 / px ** 2))
    nl = measure.label(nucl)
    D3 = ndi.gaussian_filter(Dp, (0, 1, 1)); nuc3 = (D3 > filters.threshold_otsu(D3[:, nucl]) * 0.5) & nucl[None]
    # cell area from SLC25A4 + TOM20, split into cells by nucleus-seeded watershed
    sig = ndi.gaussian_filter(S.max(0) / np.percentile(S.max(0), 99.5) + T.max(0) / np.percentile(T.max(0), 99.5), 1.0 / px)
    cell = sig > filters.threshold_otsu(sig) * 0.8
    cell = morphology.remove_small_objects(ndi.binary_fill_holes(morphology.binary_closing(cell, morphology.disk(5))), int(100 / px ** 2))
    lab = watershed(-sig, nl * (cell | nucl), mask=cell | nucl)
    BM = artefact_mask(S, T, cell, nucl, px)
    vals = Ts[:, cell & ~BM]; TTHR = filters.threshold_otsu(vals[vals > 0])     # image-wise TOM20 threshold
    rows = []
    for rp in measure.regionprops(lab):
        if rp.area * px ** 2 < 150:
            continue
        cm = lab == rp.label
        border = rp.bbox[0] == 0 or rp.bbox[1] == 0 or rp.bbox[2] == lab.shape[0] or rp.bbox[3] == lab.shape[1]
        ids = np.unique(nl[cm & nucl]); ids = ids[ids > 0]
        nuc_area = max([((nl == i) & cm).sum() for i in ids], default=0) * px ** 2
        art = (BM & cm).sum() / cm.sum(); cmb = cm & ~BM
        s3, t3, nv = Ss[:, cmb], Ts[:, cmb], nuc3[:, cmb]
        zz = (s3 + t3).sum(1); k = zz >= a.focus_frac * zz.max()
        s, t = s3[k][~nv[k]], t3[k][~nv[k]]
        if s.sum() <= 0:
            continue
        tp = t > TTHR
        rows.append(dict(cell=rp.label, area_um2=rp.area * px ** 2, nucleus_area_um2=nuc_area, touches_border=border,
                         artefact_area_fraction=art, M1=s[tp].sum() / s.sum()))
    return rows, lab, BM, (S, Dp, T)


def qc_image(SDT, lab, BM, df_img, path, title):
    S, D, T = [c.max(0) for c in SDT]
    rgb = np.dstack([np.clip(T / np.percentile(T, 99.7), 0, 1), np.clip(S / np.percentile(S, 99.7), 0, 1),
                     np.clip(D / np.percentile(D, 99.7), 0, 1) * 0.7])
    fig, ax = plt.subplots(figsize=(10, 10)); ax.imshow(rgb); ax.contour(BM, [0.5], colors="cyan", linewidths=1.2)
    d = df_img.set_index("cell")
    for rp in measure.regionprops(lab):
        if rp.label not in d.index:
            continue
        ok = bool(d.loc[rp.label, "included"])
        y0, x0, y1, x1 = rp.bbox
        ax.contour(np.arange(x0, x1), np.arange(y0, y1), (lab[y0:y1, x0:x1] == rp.label).astype(float), [0.5],
                   colors="w" if ok else "#ff9d00", linewidths=0.8)
        ax.text(rp.centroid[1], rp.centroid[0], f"#{rp.label}\n{100 * d.loc[rp.label, 'M1']:.0f}%",
                color="yellow" if ok else "#ff9d00", fontsize=10, weight="bold", ha="center", va="center")
    ax.set_title(f"{title} - white: included, orange: excluded (QC), cyan: artefact mask"); ax.axis("off")
    plt.tight_layout(); plt.savefig(path, dpi=70); plt.close()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("indir", help="folder with .czi files")
    p.add_argument("-o", "--outdir", default="results_slc25a4_tom20")
    p.add_argument("--groups", required=True,
                   help='file-name prefix -> group, comma separated, e.g. "CTRL:Control,PAT1:Patient 1"')
    p.add_argument("--omit-groups", default="", help="groups kept in cells.csv but left out of the GraphPad table")
    p.add_argument("--ch-slc", type=int, default=0); p.add_argument("--ch-dapi", type=int, default=1)
    p.add_argument("--ch-tom20", type=int, default=2)
    p.add_argument("--binning", type=int, default=2)
    p.add_argument("--focus-frac", type=float, default=0.4)
    p.add_argument("--min-cell-area", type=float, default=300.0, help="smaller objects = debris/dead cells (um2)")
    p.add_argument("--max-cell-area", type=float, default=5000.0); p.add_argument("--max-nucleus-area", type=float, default=700.0)
    p.add_argument("--min-nucleus-area", type=float, default=40.0); p.add_argument("--max-artefact", type=float, default=0.05)
    p.add_argument("--no-qc", action="store_true")
    a = p.parse_args()

    rules = parse_groups(a.groups)
    os.makedirs(os.path.join(a.outdir, "qc"), exist_ok=True)
    rows = []
    for f in sorted(glob.glob(os.path.join(a.indir, "*.czi"))):
        n = os.path.basename(f)[:-4]; g = group_of(n, rules)
        x, px = read_czi(f, a.binning)
        if g is None or x.shape[1] < 2:
            print(f"skip  {n}"); continue
        r, lab, BM, SDT = analyze_image(x, px, a)
        for d in r:
            d.update(group=g, image=n)
        df_img = pd.DataFrame(r)
        reasons = []
        for _, c in df_img.iterrows():
            rr = []
            if c.area_um2 < a.min_cell_area: rr.append("debris/dead cell (<%g um2)" % a.min_cell_area)
            if c.touches_border: rr.append("touches image border")
            if c.area_um2 > a.max_cell_area or c.nucleus_area_um2 > a.max_nucleus_area: rr.append("merged cells")
            if c.nucleus_area_um2 < a.min_nucleus_area: rr.append("fragment without nucleus")
            if c.artefact_area_fraction > a.max_artefact: rr.append("artefact > %g%%" % (100 * a.max_artefact))
            reasons.append("; ".join(rr))
        df_img["exclusion_reason"] = reasons; df_img["included"] = df_img.exclusion_reason == ""
        rows.append(df_img)
        if not a.no_qc:
            qc_image(SDT, lab, BM, df_img, os.path.join(a.outdir, "qc", f"{n}_QC.png"), n)
        print(f"done  {n}: {len(df_img)} cells, {int(df_img.included.sum())} included", flush=True)

    df = pd.concat(rows, ignore_index=True)
    cols = ["group", "image", "cell", "included", "exclusion_reason", "M1",
            "area_um2", "nucleus_area_um2", "artefact_area_fraction", "touches_border"]
    df["exclusion_reason"] = df.exclusion_reason.fillna("")
    df[cols].to_csv(os.path.join(a.outdir, "cells.csv"), index=False)
    s = df[df.included & ~df.group.isin(a.omit_groups.split(","))]
    groups = list(dict.fromkeys(g for _, g in rules if g in set(s.group)))
    pd.concat([pd.Series((100 * s[s.group == g].M1).round(2).values, name=g) for g in groups], axis=1) \
      .to_csv(os.path.join(a.outdir, "M1_for_GraphPad.csv"), index=False)
    json.dump(vars(a), open(os.path.join(a.outdir, "parameters.json"), "w"), indent=2)
    print("included cells per group:", s.group.value_counts().reindex(groups).to_dict())


if __name__ == "__main__":
    main()
