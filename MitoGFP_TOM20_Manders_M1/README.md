# MitoGFP–TOM20 colocalization (Manders' M1), single-cell 3D analysis

This script measures how much of a transfected mitochondria-targeted GFP (MitoGFP) reaches mitochondria, as an indirect read-out of mitochondrial protein import. Mitochondria are labelled by TOM20 immunofluorescence. Only transfected cells are analysed, and each cell is analysed separately in 3D.

**Main read-out:** Manders' coefficient M1

```
M1 = Σ GFP in TOM20-positive voxels  /  Σ GFP in the cell
```

Cells are first selected as transfected if their GFP signal is at least 3 times the autofluorescence (see Pipeline, step 4). Within each selected cell, no intensity threshold is applied to the GFP voxels: the total GFP signal of the cell is in the denominator. M1 = 1 means that all MitoGFP is located in the TOM20-positive compartment.

## Requirements

Python ≥ 3.9.

```bash
pip install -r requirements.txt
```

## Input

A Leica `.lif` file containing 3-channel z-stacks. The default channel order is 0 = DAPI, 1 = MitoGFP, 2 = TOM20. Groups are assigned from the beginning of the series name with `--groups`, for example `CTRL_1` → Control and `PAT1_3` → Patient 1.

## Usage

Example:

```bash
python mitogfp_tom20_m1.py experiment.lif -o results_mitogfp_tom20 \
    --groups "CTRL:Control,PAT1:Patient 1,PAT2:Patient 2"
```

Use `python mitogfp_tom20_m1.py -h` to list all options.

| Option | Default | Meaning |
|---|---|---|
| `--groups` | required | series-name prefix → group |
| `--exclude` | – | series to skip (e.g. repeated acquisitions of the same field) |
| `--ch-dapi / --ch-gfp / --ch-tom20` | 0 / 1 / 2 | channel indices |
| `--sar` | 3.0 | minimum GFP / autofluorescence ratio for a cell to count as transfected |
| `--min-area` | 150 | minimum region area (µm²); removes debris and small fragments |
| `--max-area` | 3000 | maximum cell area (µm²); larger regions are treated as merged clusters |
| `--focus-frac` | 0.4 | in-focus planes: per-cell GFP sum ≥ 40 % of the maximum plane |
| `--no-qc` | – | skip QC images |

## Pipeline

1. **Background.** The median of each channel is subtracted, and the image is smoothed with a 1-pixel Gaussian. The noise SD is estimated from the median absolute deviation (MAD).
2. **GFP-positive regions.** The GFP maximum projection is smoothed (σ = 1 µm) and thresholded with the triangle method (at least 4 × noise SD). The mask is then closed (1.5 µm), holes are filled, and regions smaller than 150 µm² are removed to exclude debris and small GFP-positive fragments (e.g. detached processes), which are much smaller than a whole myoblast.
3. **Single cells.** Regions containing more than one DAPI nucleus are split by nucleus-seeded watershed on the combined GFP + TOM20 signal.
4. **Transfection criterion.** GFP autofluorescence is measured in non-transfected, TOM20-positive cells of the same field (95th percentile). A region is kept only if its GFP 95th percentile is at least `--sar` times this value. Laser power and detector gain affect both values equally, so the ratio can be compared between images. Regions larger than `--max-area` are excluded as merged clusters.
5. **Mitochondrial compartment.** An Otsu threshold is applied to TOM20 in each image, using voxels with TOM20 signal.
6. **Per-cell 3D measurement.** For each cell, the in-focus z-planes are selected and M1 is computed over all voxels of the cell.

## Output (`-o` folder)

| File | Content |
|---|---|
| `cells.csv` | one row per GFP+ region: group, series, cell no., `included`, `exclusion_reason`, M1, mean GFP, GFP / autofluorescence ratio, area |
| `M1_for_GraphPad.csv` | one column per group with M1 (%) of the included cells, ready to paste into a GraphPad "Column" table |
| `qc/<series>_QC.png` | max projection (green GFP, red TOM20) with cell outlines: white = included (with M1 %), grey = excluded |
| `parameters.json` | all parameters used |

## Reference

Manders EMM, Verbeek FJ, Aten JA. Measurement of co-localization of objects in dual-colour confocal images. *J Microsc.* 1993;169(3):375–382. doi:10.1111/j.1365-2818.1993.tb03313.x
