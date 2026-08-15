# Scheib_2026

Code and data repository for the Scheib et al., 2026 manuscript in eLife from
the Kerlin Lab at the University of Minnesota, Twin Cities.

The repository contains the per-animal behavioral and two-photon calcium-imaging data used in the
study, the pre-computed trace alignments and example-ROI exports the figures are built from, and the
six Jupyter notebooks (plus one shared Python library) that generate every manuscript figure panel.

The notebooks are **self-contained**: they import only `js_manuscript_final.py` and standard
scientific-Python packages. No lab-internal packages are required.

SLeD refers to the lab specific environment (Synapses Learning and Dendrites) however most of required
tools should be available in the single package file.

NOTE: preparation of this readme was performed with assistance from Claude claude-opus-5 however most
of the actual code was written by Zachary Newman, Jackson Scheib and Aaron Kerlin with minimal AI
assistance.

Please send any questions to newmanza@umn.edu.


---

## Repository structure

```
Scheib_2026/
├── Data/                              # ~8.7 GB, 301 files
│   ├── <animal ID>/                   # 22 animal folders (e.g. B00002213999)
│   │   ├── behavior/                  # Bpod behavioral session data (.pkl)
│   │   ├── somas/                     # Per-session somatic ROI calcium data (.pkl)
│   │   └── dendrites/                 # Per-session dendritic ROI calcium data (.pkl, not all animals)
│   ├── examples/                      # Per-figure example ROIs / images used by individual panels
│   │   ├── fig1_behav_examples/       # Side + bottom camera .tiff frames for the Fig1 lick examples
│   │   ├── fig2_hist_examples/        # Histology crop + image data (example animal)
│   │   ├── fig2_indiv_examples/       # Movie-preview frames + traces for two example ROIs
│   │   ├── fig2_zstack_examples/      # Structural z-stack .tif projections (example animal)
│   │   ├── fig3_indiv_examples/       # Example ROI traces, NMF footprints, ROI manifest
│   │   ├── fig4_indiv_examples/       # Example ROI traces (go-cue aligned)
│   │   ├── fig5_indiv_examples/       # Example ROI traces (contact-time aligned)
│   │   └── fig6_indiv_examples/       # Example ROI traces (go-cue aligned)
│   └── traceAlignments/               # Pre-computed, split/compressed cluster-level trace alignments
│       ├── fig3/  fig3_licking/
│       └── fig5/  fig5_licking/
├── Figures/
│   ├── js_manuscript_final.py         # Shared analysis/plotting library (~15.4k lines, 152 functions)
│   ├── Fig1.ipynb … Fig6.ipynb        # One notebook per manuscript figure
├── requirements.txt                   # pip dependency pins
├── Scheib2026_env.yml                 # conda environment definition
├── LICENSE                            # MIT License
└── README.md
```

---

## Quick start

### 1. Create the environment

The reference environment is **Python 3.8** with a NumPy 1.x / SciPy 1.10 stack. Either route works:

```bash
# conda (recommended — reproduces the reference environment)
conda env create -f Scheib2026_env.yml
conda activate Scheib2026_env

# or pip, into an existing Python 3.8 environment
pip install -r requirements.txt
```

Register the kernel so the notebooks can find it:

```bash
python -m ipykernel install --user --name Scheib2026_env --display-name "Scheib2026_env"
```

> The notebooks were last saved against a kernel named `sled_env` (Python 3.8.20). After creating
> `Scheib2026_env`, select it from the Jupyter kernel picker — no notebook edits are needed.

### 2. Point the notebooks at your copy of the data

Every notebook's first cell defines three paths. Edit them to match your machine:

```python
PARENT_DIR = "/home/range6-raid17/kerlinlab/sled/sled-klab/manuscripts/Scheib2026/VOR"
DATA_DIR   = "/export/general/pooledData/Jackson/repo_data/"
figSaveDir = os.path.join(PARENT_DIR, "figurePanels", "Fig1")
PDF_export_active = False
```

- **`DATA_DIR`** must point at this repository's `Data/` directory.
- **`figSaveDir`** is where PDF panels are written. It must already exist if you enable export.
- **`PDF_export_active`** gates all PDF writing. It ships as `False`, so figures render inline and
  nothing is written to disk. Set it to `True` once `figSaveDir` exists to export the panels.
  (The `maybe_pdf()` / `_NoPdf` helper defined in the same cell implements this; note that a handful
  of cells still call `PdfPages(...)` directly and will write regardless — see
  [Known limitations](#known-limitations).)

The same cell also offers an alternative to `import js_manuscript_final as jsm`: a commented-out
`importlib.util.spec_from_file_location` block for loading the library from an explicit path, and a
commented-out `ipynbname` block for deriving `PARENT_DIR` from the notebook's own location. Both are
optional.

### 3. Run a notebook

All figures are produced by running notebooks interactively, cell by cell. None take command-line
arguments. Start with the `init` cell, then the data-loading cell, then whichever panel section you
want.

Fonts: panels are exported with `font.family='sans-serif'`, `font.sans-serif=['Arial']` and
`pdf.fonttype=42` (editable text in the PDF). **Arial must be installed locally** for the panels to
render and export as intended.

---

## Data/

### Per-animal raw data

22 subjects, folders named with barcode-style IDs (e.g. `B00002213999`). Each animal folder has up to
three subfolders; not every animal has all three, and session counts vary.

| File | Contents |
|---|---|
| `behavior/<ID>_behavior.pkl` | One per animal. Bpod state-machine data from the directional-lickport, delayed-response task with a sensorimotor "shift" manipulation. Includes `trajs` (lick trajectories), `mats` (session info), `port_locs` (pre/post-shift lickport positions), `raw`. |
| `somas/<ID>_somas.pkl` | Small per-animal metadata file: `daysRelativeToShift`, `roiScores`. |
| `somas/<ID>_somas_session_N.pkl` | One per imaging session (N typically 0–5). Per-ROI deconvolved calcium "impulse trace" activity aligned to trial frames, somatic ROIs. |
| `dendrites/<ID>_dendrites*.pkl` | Same structure, dendritic ROIs. Session counts run up to N=7; several animals have no `dendrites/` folder (`'NOT RECORDED'` in the master dictionary). |

`jsm.masterData_loading(DATA_DIR)` walks this tree and assembles the in-memory `masterData`
dictionary the notebooks operate on. The per-animal shift metadata (`shiftDay`, `shiftTrial`,
`shiftDir`) and the list of sessions to load are hard-coded inside that function, so the animal set is
fixed — you do not need a separate consolidated dataset file.

**Memory note:** `masterData_loading()` loads every session of every animal at once. Budget well
above the ~7 GB on-disk size of the animal folders.

### `Data/examples/`

Small, pre-extracted slices of the raw data that let individual example panels be reproduced without
re-deriving them from the full dataset. Loaded via `jsm.unpickler()` and
`jsm.load_all_grouped_aligned_traces_split()`. The trace arrays are stored trimmed to only the
ROIs/frames a panel actually reads and are scattered back into full-shape arrays at load time.

### `Data/traceAlignments/`

Cluster-level trace alignments for Fig3 and Fig5, written by
`save_ROI_clustering_reduced_split` / `save_all_grouped_aligned_traces_split` as a `_meta` file plus
one part file per cluster, and read back with `jsm.load_ROI_clustering_reduced_split()` and
`jsm.load_all_grouped_aligned_traces_split()`. Both loaders take filter arguments
(`ROI_type_keys`, `ROI_scores`, `clusterList`, `trialGroupings`, `trace_types`, …) so a notebook only
reads the parts it needs.

---

## Figures/

### `js_manuscript_final.py`

A single shared library (imported everywhere as `jsm`) holding 152 functions. Its only third-party
imports are `numpy`, `scipy`, `matplotlib`, `seaborn`, `statsmodels`, `scikit-learn` and `tqdm`; it
also sets the manuscript-wide matplotlib rcParams at import time. Major functional groups:

- **Data loading / IO** — `masterData_loading()`, `pickler()`, `unpickler()`,
  `load_ROI_clustering_reduced_split()`, `load_all_grouped_aligned_traces_split()`,
  `resolve_align_trace_type()`: assemble `masterData` from `Data/`, and read back the split,
  optionally gzip/bz2/lzma-compressed alignment and example files.
- **Lick geometry & behavior** — `get_contacts()`, `get_lickAngle()`, `get_aziEle()`,
  `get_relDistances()`, `get_errors_LDA()`, `get_secondLickSide()`, `get_MER_bySecondLick()`,
  `get_sl_errors()`, `get_relDistFromLDA()`, `plot_lick_boundary_3d()`: 3D lick azimuth/elevation
  relative to an LDA left/right decision boundary, and classification of motor errors (ME) vs.
  directional errors (DE).
- **Trial structure & timing** — `get_tracker()`, `get_cTimes()/get_cTimes2()`, `get_lickTimes()`,
  `get_goCues()`, `get_lastLick()`, `get_autoWaterTrials()`, `get_behaviorShiftMask()`,
  `load_trial_structure_times()`, `trial_epoch_frameLabels()`.
- **Trace alignment** — `tsep_traces()`, `get_tsep()`, `psTsep_func()`, `downSampleTsep()`,
  `symmetric_orthogonalization()`, `extract_ROI_trace()`, `extract_cluster_traces()`: peri-event
  time-series extraction, downsampling, and vector orthogonalization for population activity.
- **Population / coding-direction (CD) analysis** — `pullAnmHH()`, `pullAnmTsep()` (and their `_ALT` /
  `_SMO` variants), `gather_anmShuff_megaTsepHH()`, `gather_anmShuff_megaTsepFull()`,
  `get_SMO_fast()`, `get_SMO_CD2_fast()`, `get_CRMER_CD2_fast()`, `get_CRCL_CD_fast()`,
  `get_CRCL_CD_fast_earlyLate()`, `get_CRCL_CD_full()`, `get_go_contact_irf()`,
  `pullAnmTsepComputeSMOXval()`: build trial-type-pooled and shuffled ROI activity tensors, compute
  coding-direction projections and cross-validated SMO (Sensory / Motor-Choice / Outcome-Reward)
  scores, and fit go-cue vs. lick-contact impulse response functions.
- **Animal-level aggregation** — `get_allAnmParams()` / `get_allAnmParams2()` (verbose extended-phase
  vs. quiet standard-phase variants), `update_allAnmParams()`, `check_anms_ROIs()`, `check_trials()`,
  `full_get_SMO()`: assemble the master per-animal parameter dictionary (tracker, cue/contact times,
  shift-phase masks, error flags, ROI-keep masks).
- **Statistics** — `summary_stats()`, `pstar()`, `multi_test_correction()`,
  `prep_manual_contrasts()`, `standard_bootstrap()`, `vectorized_bootstrap()`,
  `simple_trial_trace_summary_stats()`, `bootstrap_summary_stats()`.
- **Plotting** — `clean_subplots()`, `display_clean_subplots()`, `plot_CD_projections()`,
  `enforce_equal_axes()`, `figure_grid()`, `remove_all_clipping()`, `set_dynamic_suptitle()`,
  `zoom_lookup()`, `add_plot_scaleBar()`, `generate_cmap()`, `convert2RGB()`, `export_colorbar()`:
  standardized formatting that keeps panels visually consistent (transparent backgrounds, Arial,
  editable-text PDF export).
- **Raster / movie panels** — `aligned_ROI_cluster_trace_raster_fig()`, `ROI_cluster_traceSets()`,
  `trial_alignments_byROI_pages()`, `image_and_trace_movie_preview()`, `image_and_trace_panels()`,
  `draw_ROI_borders()`, `add_image_scaleBar()`.

### Figure notebooks

| Notebook | Data it loads | Description |
|---|---|---|
| `Fig1.ipynb` | full `masterData` + `examples/fig1_behav_examples/` | Licking behavior. Session/trial count stats, example camera frames with the lick-angle construction, 3D lick projections and the LDA left/right boundary, per-animal lick-angle distributions, motor-error rate (MER) and directional-error rate (DER) across the sensorimotor shift. |
| `Fig2.ipynb` | `examples/fig2_indiv_examples/` | Reconstructs two example ROIs (one dendrite, one soma) from the trimmed exports and re-renders `image_and_trace_movie_preview()` frames — imaging field of view with ROI overlays alongside the simultaneous traces. |
| `Fig3.ipynb` | `traceAlignments/fig3/`, `traceAlignments/fig3_licking/`, `examples/fig3_indiv_examples/` | Example-ROI NMF footprints, trial-aligned activity rasters by cluster, and the matching lick rasters. |
| `Fig4.ipynb` | full `masterData` + `examples/fig4_indiv_examples/` | Modeled responses; go-cue vs. lick-contact impulse response functions per ROI; SMO (Sensory / Motor-Choice / Outcome-Reward) coding-direction preparation, cross-validated calculation, projections, selectivity and correlations. **Uses a Dask `LocalCluster`** for the heavy SMO passes. |
| `Fig5.ipynb` | full `masterData` + `examples/fig5_indiv_examples/` + `traceAlignments/fig5*/` | Core motor-error / coding-direction figure: CR–CA, CR–CL and AP–CA coding-direction comparisons, licking data, example rasters, and a supplement section on raw-trial correlations. |
| `Fig6.ipynb` | full `masterData` + `examples/fig6_indiv_examples/` | Coding direction across the shift: pre/post selectivity index, CD projections at finer time bins, a CR–CL scrolling panel, and pre-late / post-early / post-late window comparisons including percent-change summaries. |

#### Dask (Fig4 only)

Two cells in `Fig4.ipynb` spin up a local Dask cluster with lab-specific settings:

```python
dask.config.set(temporary_directory='/export/disk1/aaron/')   # also '/export/disk1/scratch/dask'
cluster = LocalCluster(n_workers=4, threads_per_worker=8, memory_limit="20GiB")
client  = Client(cluster)
```

Change `temporary_directory` to a scratch path that exists on your machine, and scale
`n_workers` / `threads_per_worker` / `memory_limit` to your hardware (the defaults assume ~80 GB of
usable RAM).

---

## Dependencies

Full pins are in [`requirements.txt`](requirements.txt) and [`Scheib2026_env.yml`](Scheib2026_env.yml).
Summary of what is actually imported:

| Package | Used by |
|---|---|
| `numpy`, `scipy`, `matplotlib` | everything |
| `seaborn`, `statsmodels`, `scikit-learn` | `js_manuscript_final.py`; `sklearn.mixture.GaussianMixture` also directly in Fig1 |
| `tqdm` | progress bars throughout |
| `tifffile` (+ `imagecodecs`) | Fig1/Fig2 example `.tiff` / `.tif` reads |
| `pandas`, `pillow` | transitive (seaborn, matplotlib) |
| `dask`, `distributed` | Fig4 only |
| `ipykernel`, `ipython`, `notebook`, `ipympl` | Jupyter runtime; `ipympl` backs the `%matplotlib widget` cell in Fig1 |

Optional, not needed for any committed notebook:

- `opencv-python` — only for `jsm.readFrames()` (`cv2` is not imported at module level).
- `joblib` — only for `jsm.run_parallel_bootstrap()` (`Parallel`/`delayed` are not imported at module
  level). `jsm.vectorized_bootstrap()` is the NumPy-only equivalent and is what the notebooks use.
- `ipynbname` — only if you enable the commented-out auto-`PARENT_DIR` block.

**NumPy must stay on the 1.x line.** `js_manuscript_final.py` installs a compatibility shim at import
time so that pickles written under NumPy 2.x can be read under NumPy 1.x:

```python
sys.modules.setdefault('numpy._core', numpy.core)
```

Installing NumPy ≥ 2 changes pickle module paths and defeats this shim.

---

## Known limitations

- **Two Fig1 cells still require lab-internal raw data.** The cells that set
  `behPath = '/home/range6-raid17/kerlinlab/raw/behavior'` enumerate recording days and derive
  `camIDs` from the raw session directory. The cell that follows falls back to the repo copies in
  `Data/examples/fig1_behav_examples/` when they exist, but it still needs `camIDs` from the earlier
  cell, so it cannot run standalone without setting `camIDs` (and `days`, `dayPath`) by hand.
- **`Data/examples/fig2_hist_examples/` and `fig2_zstack_examples/` are unused by the committed
  notebooks.** `Fig2.ipynb` covers only the movie-preview panels; the histology and z-stack panels
  are not in the repo.
- **`PDF_export_active` does not gate every write.** Most cells route through `maybe_pdf()`, but some
  still call `PdfPages(os.path.join(figSaveDir, figName))` directly and will attempt to write even
  when the flag is `False`. Create `figSaveDir` before running a full notebook top-to-bottom.
- **Absolute lab paths ship as the defaults.** `PARENT_DIR`, `DATA_DIR`, `figSaveDir`,
  `movieParams['previewPDFDir']` (Fig2) and the Dask `temporary_directory` (Fig4) all point at the
  original lab server and must be edited.
- **Save-side cells write into `Data/`.** The alignment/example cells in Fig3 and Fig5 are the
  original export cells; the load cells beside them are what you want for reproduction. Re-running an
  export cell requires the full upstream dataset and will overwrite files under `Data/`.
- **`js_manuscript_final.py` has two `clean_subplots()` and two `summary_stats()` definitions**
  (lines 2642/6598 and 941/6853). The later definition wins at import time.

---

## Citation

[Scheib et al., 2026 Distinct sensorimotor encoding in tuft dendrites and somata associated with action, correction, and learning. eLife](https://elifesciences.org/reviewed-preprints/111876)


====================================================================================================
## Raw data inventory

Due to the large size of the raw data files they are not included here in the repo but here is a summary of all files that are available upon special request

22 animals, 167 imaging sessions (somas / dendrites) over 167 recording days, 336 raw directories, **147,203 files**, **20.07 TB** (20,073,459,670,400 bytes).

Indexed 2026-08-14 17:30:19 from the raw paths recorded in `summaryInfo` (`imagingRawPath` / `behaviorRawPath` per session). 
Sessions are counted per compartment; where a directory is shared by two sessions (a behavior day recorded for both compartments) 
it is counted once in the animal and total rows, so the per-session rows can sum higher.

| Animal | Sessions | Imaging files | Imaging size | Behavior files | Behavior size | Total files | Total size |
|---|---|---|---|---|---|---|---|
| `B00002121749` | 5 | 1,347 | 70.31 GB | 2,707 | 242.23 GB | 4,054 | 312.54 GB |
| `B00002121774` | 9 | 2,934 | 705.04 GB | 5,339 | 467.68 GB | 8,273 | 1.17 TB |
| `B00002121777` | 10 | 3,039 | 826.72 GB | 6,489 | 580.75 GB | 9,528 | 1.41 TB |
| `B00002213772` | 11 | 3,384 | 698.98 GB | 6,795 | 608.72 GB | 10,179 | 1.31 TB |
| `B00002213773` | 3 | 1,015 | 259.06 GB | 2,039 | 182.8 GB | 3,054 | 441.86 GB |
| `B00002213784` | 9 | 2,586 | 532.05 GB | 4,914 | 435.38 GB | 7,500 | 967.44 GB |
| `B00002213785` | 9 | 2,336 | 505.5 GB | 4,419 | 390.52 GB | 6,755 | 896.02 GB |
| `B00002213889` | 7 | 1,913 | 325.4 GB | 3,843 | 343.92 GB | 5,756 | 669.32 GB |
| `B00002213908` | 9 | 2,753 | 517.28 GB | 5,531 | 495.47 GB | 8,284 | 1.01 TB |
| `B00002213909` | 5 | 1,569 | 364.17 GB | 3,149 | 282.14 GB | 4,718 | 646.31 GB |
| `B00002213920` | 8 | 2,298 | 418.65 GB | 4,614 | 413.07 GB | 6,912 | 831.72 GB |
| `B00002213921` | 6 | 1,697 | 357.19 GB | 3,413 | 305.43 GB | 5,110 | 662.61 GB |
| `B00002213932` | 8 | 2,534 | 598.81 GB | 5,088 | 455.87 GB | 7,622 | 1.05 TB |
| `B00002213943` | 7 | 1,840 | 406.54 GB | 3,103 | 279.47 GB | 4,943 | 686.01 GB |
| `B00002213944` | 2 | 571 | 216.52 GB | 1,139 | 102.72 GB | 1,710 | 319.24 GB |
| `B00002213970` | 10 | 2,805 | 944.15 GB | 5,589 | 503.51 GB | 8,394 | 1.45 TB |
| `B00002213982` | 8 | 2,488 | 530.44 GB | 4,958 | 446.82 GB | 7,446 | 977.26 GB |
| `B00002213985` | 7 | 2,229 | 465.53 GB | 4,445 | 400.6 GB | 6,674 | 866.14 GB |
| `B00002213997` | 9 | 2,845 | 750.6 GB | 5,681 | 511.99 GB | 8,526 | 1.26 TB |
| `B00002213998` | 7 | 1,759 | 684.3 GB | 3,509 | 316.11 GB | 5,268 | 1.0 TB |
| `B00002213999` | 9 | 2,833 | 606.59 GB | 5,653 | 509.47 GB | 8,486 | 1.12 TB |
| `B00002214001` | 9 | 2,710 | 546.22 GB | 5,301 | 468.73 GB | 8,011 | 1.01 TB |
| **All animals** | **167** | **49,485** | **11.33 TB** | **97,718** | **8.74 TB** | **147,203** | **20.07 TB** |


## License

Released under the [MIT License](LICENSE).
