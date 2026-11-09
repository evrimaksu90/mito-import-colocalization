# SLC25A4–TOM20 colocalization (Manders' M1), single-cell 3D analysis

This script measures what fraction of endogenous SLC25A4 (ANT1) is located in the TOM20-positive mitochondrial compartment. Each cell is analysed separately in 3D, from immunofluorescence z-stacks.

**Main read-out:** Manders' coefficient M1

```
M1 = Σ SLC25A4 in TOM20-positive voxels  /  Σ SLC25A4 in the cell
```

No threshold is applied to SLC25A4. Nuclear voxels and artefact masks are excluded from both sums.

## Requirements

Python ≥ 3.9.

```bash
pip install -r requirements.txt
```

## Input

A folder of Zeiss `.czi` z-stacks. The default channel order is 0 = SLC25A4 (AF488), 1 = DAPI, 2 = TOM20 (AF568). Single-plane files (for example `*_Snap.czi`) are skipped automatically. Groups are assigned from the beginning of the file name with `--groups`, for example `PAT1_2.czi` → Patient 1.

## Usage

Example:

```bash
python slc25a4_tom20_m1.py /path/to/czi_folder -o results_slc25a4_tom20 \
    --groups "CTRL:Control,PAT1:Patient 1,PAT2:Patient 2"
```

Use `python slc25a4_tom20_m1.py -h` to list all options.

| Option | Default | Meaning |
|---|---|---|
| `--groups` | required | file-name prefix → group |
| `--omit-groups` | – | groups kept in `cells.csv` but left out of the GraphPad table |
| `--ch-slc / --ch-dapi / --ch-tom20` | 0 / 1 / 2 | channel indices |
| `--binning` | 2 | xy binning (0.099 µm → 0.198 µm pixels) |
| `--min-cell-area` | 300 | smaller objects are treated as debris or dead cells (µm²) |
| `--max-cell-area` / `--max-nucleus-area` | 5000 / 700 | above these: merged cells (µm²) |
| `--min-nucleus-area` | 40 | below this: fragment without a nucleus (µm²) |
| `--max-artefact` | 0.05 | maximum fraction of the cell area covered by the artefact mask |
| `--focus-frac` | 0.4 | in-focus planes: per-cell SLC25A4 + TOM20 sum ≥ 40 % of the maximum plane |
| `--no-qc` | – | skip QC images |

## Pipeline

1. **Reading.** Each `.czi` is read and binned 2 × 2 in xy.
2. **Background correction.** For every z-plane, an additive column-wise and row-wise background is subtracted (smoothed 5th-percentile profiles). This removes the vertical stripe present in the scan field.
3. **Nuclei.** DAPI maximum projection with a triangle threshold. Nuclear voxels (3D) are excluded from the measurement.
4. **Cell segmentation.** The cell area is defined from the combined SLC25A4 + TOM20 signal (0.8 × Otsu). It is split into single cells by watershed seeded with the nuclei.
5. **Artefact mask.** Bright objects that are thick or round, and therefore inconsistent with mitochondrial morphology (aggregates, debris, dead cells), are masked. They are found as high-intensity pixels that survive a morphological opening with a ~1.4 µm disk; this opening removes thin mitochondrial tubules. Objects must be compact, 1.5–200 µm², and outside nuclei. The mask is dilated by 0.8 µm.
6. **Mitochondrial compartment.** An Otsu threshold is applied to TOM20 in each image, within the cell area and outside the artefact mask.
7. **Per-cell 3D measurement.** For each cell, the in-focus planes (SLC25A4 + TOM20 sum ≥ 40 % of the maximum plane) are selected. M1 is computed after removing nuclear and artefact voxels.
8. **Predefined QC exclusion.** The same criteria are applied to all groups. A cell is excluded if it:
   - is smaller than 300 µm²,
   - touches the image border,
   - is a merged cell (area > 5000 µm² or nucleus > 700 µm²),
   - is a fragment without a nucleus (nucleus < 40 µm²), or
   - has more than 5 % of its area covered by artefacts.

## Output (`-o` folder)

| File | Content |
|---|---|
| `cells.csv` | one row per cell: group, image, cell no., `included`, `exclusion_reason`, M1, area, nucleus area, artefact fraction |
| `M1_for_GraphPad.csv` | one column per group with M1 (%) of the included cells, ready to paste into a GraphPad "Column" table |
| `qc/<image>_QC.png` | max projection (green SLC25A4, red TOM20, blue DAPI): white = included cell (M1 %), orange = excluded cell, cyan = artefact mask |
| `parameters.json` | all parameters used |

## Reference

Manders EMM, Verbeek FJ, Aten JA. Measurement of co-localization of objects in dual-colour confocal images. *J Microsc.* 1993;169(3):375–382. doi:10.1111/j.1365-2818.1993.tb03313.x
