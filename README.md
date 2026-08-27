# Scheib_2026

Code and data repository for the Scheib et al., 2026 manuscript in eLife from
the Kerlin Lab at the University of Minnesota, Twin Cities.

The repository contains the per-animal behavioral and two-photon calcium-imaging data used in the
study, the pre-computed trace alignments and example-ROI exports the figures are built from, and the
six Jupyter notebooks (plus one shared Python library) that generate every manuscript figure panel.

The notebooks are **self-contained**: they import only `js_manuscript_final.py` and standard
scientific-Python packages. No lab-internal packages are required.

**Companion raw-data repository:** the true raw and motion-registered two-photon imaging data behind
the two Fig2 example ROIs (one dendrite, one soma) are published separately at
**[kerlinlab/Scheib_2026_data](https://github.com/kerlinlab/Scheib_2026_data)** (Git LFS). Nothing in
this repository needs it — `Fig2.ipynb` runs off the trimmed exports in `Data/examples/` — but it is
there if you want the full movies rather than the frames the panels draw.

---

## Citation

[Scheib et al., 2026 Distinct sensorimotor encoding in tuft dendrites and somata associated with action, correction, and learning. eLife](https://elifesciences.org/reviewed-preprints/111876)

Raw data for the Fig2 examples is archived separately:
[kerlinlab/Scheib_2026_data](https://github.com/kerlinlab/Scheib_2026_data).

---

## Repository structure

```
Scheib_2026/
├── Data/                              # 6.64 GB, 469 files
│   ├── <animal ID>/                   # 22 animal folders (e.g. B00002213999) — 5.28 GB, 397 files
│   │   ├── behavior/                  # Bpod behavioral data, split as one '_meta' + one '_day_<N>' .pkl
│   │   ├── somas/                     # Per-session somatic ROI calcium data (.pkl, 20 of 22 animals)
│   │   └── dendrites/                 # Per-session dendritic ROI calcium data (.pkl, 21 of 22 animals)
│   ├── examples/                      # Per-figure example ROIs / images used by individual panels
│   │   ├── fig1_behav_examples/       # Side + bottom camera .tiff frames + rawBehavior.pkl for Fig1
│   │   ├── fig2_hist_examples/        # Histology crop + image data (example animal; see Known limitations)
│   │   ├── fig2_indiv_examples/       # Movie-preview frames + traces for two example ROIs
│   │   ├── fig2_zstack_examples/      # Structural z-stack .tif projections (example animal; see Known limitations)
│   │   ├── fig3_indiv_examples/       # Example ROI traces, NMF footprints, ROI manifest
│   │   ├── fig4_indiv_examples/       # Example ROI traces (go-cue aligned)
│   │   ├── fig5_indiv_examples/       # Example ROI traces (contact-time aligned)
│   │   └── fig6_indiv_examples/       # Example ROI traces (go-cue aligned)
│   └── traceAlignments/               # Pre-computed, split cluster-level trace alignments
│       ├── fig3/  fig3_licking/
│       └── fig5/  fig5_licking/
├── Figures/
│   ├── js_manuscript_final.py         # Shared analysis/plotting library (14,945 lines, 162 functions)
│   ├── Fig1.ipynb … Fig6.ipynb        # One notebook per manuscript figure
├── requirements.txt                   # pip dependency pins
├── Scheib2026_env.yml                 # conda environment definition
├── .gitignore
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

> The notebooks are saved against a kernel named `scheib2026_env` (display name `Scheib2026_env`,
> Python 3.8.20) — the exact name the `ipykernel install` line above creates. No notebook edits are
> needed; if the kernel is missing, just pick `Scheib2026_env` from the Jupyter kernel picker.

### 2. Point the notebooks at your copy of the data

**No paths need editing.** Every notebook's `init` cell (cell 1) derives them from the notebook's own
location:

```python
try:
    import ipynbname
    NOTEBOOK_PATH = ipynbname.path()
    PARENT_DIR = NOTEBOOK_PATH.parent.parent       # the repository root
    DATA_DIR = os.path.join(PARENT_DIR, "Data")
except:
    print("UNABLE TO DETERMINE CURRENT SCRIPT LOCAITON PLEASE PROVIDE PATHS MANUALLY")

FIG_EXPORT_DIR = os.path.join(PARENT_DIR, "figurePanels", "Fig1")
PDF_export_active = False
if PDF_export_active:
    os.makedirs(FIG_EXPORT_DIR, exist_ok=True)
```

- **`PARENT_DIR`** is the repository root: `ipynbname.path()` returns `<repo>/Figures/FigN.ipynb`, so
  `.parent.parent` walks up past `Figures/`. **`ipynbname` must be installed** (it is pinned in
  `requirements.txt` / `Scheib2026_env.yml`) — the `try/except` only prints a warning, but
  `PARENT_DIR` is then left undefined and the next line raises `NameError`. Commented-out
  `DATA_DIR` / `PARENT_DIR` lines are provided just below for setting them by hand in that case.
- **`DATA_DIR`** resolves to `<repo>/Data` as long as the notebooks stay in `Figures/`. Move a
  notebook elsewhere and you must set the paths manually.
- **`FIG_EXPORT_DIR`** (formerly `figSaveDir`) is `<repo>/figurePanels/FigN`, where PDF panels are
  written. You no longer need to create it yourself — the init cell creates it when
  `PDF_export_active` is `True`.
- **`PDF_export_active`** gates all PDF writing. It ships as `False`, so figures render inline and
  nothing is written to disk. Set it to `True` to export the panels. Every panel cell routes through
  the `maybe_pdf()` / `_NoPdf` helper defined in the same cell, except one Fig4 cell that uses an
  equivalent `contextlib.nullcontext()` guard.

The init cell also tries an `importlib.util.spec_from_file_location` load of `js_manuscript_final.py`
before falling back to `import js_manuscript_final as jsm`. The guard it tests
(`os.path.join(NOTEBOOK_PATH, "js_manuscript_final.py")`) joins onto the notebook *file* path, so it
never matches and the plain `import` is always what runs. That works because the library sits next to
the notebooks in `Figures/` — but keep them together, or put `Figures/` on `PYTHONPATH`.

### 3. Run a notebook

All figures are produced by running notebooks interactively, cell by cell. None take command-line
arguments. Start with the `init` cell, then the data-loading cell (`masterData_loading()` in Fig1,
Fig4, Fig5 and Fig6; the `examples/` or `traceAlignments/` load cell in Fig2 and Fig3), then whichever
panel section you want.

Fonts: panels are exported with `font.family='sans-serif'`, `font.sans-serif=['Arial']` and
`pdf.fonttype=42` (editable text in the PDF). **Arial must be installed locally** for the panels to
render and export as intended.

---

## Data/

### Per-animal raw data

22 subjects, folders named with barcode-style IDs (e.g. `B00002213999`), 167 imaging sessions in
total (92 dendrite + 75 soma). Every animal has `behavior/`; 21 have `dendrites/` and 20 have
`somas/`. The three animals missing a compartment are marked `'NOT RECORDED'` in
`masterData_loading()`: `B00002121749` (no dendrites), `B00002213944` and `B00002213998` (no somas).

| File | Contents |
|---|---|
| `behavior/<ID>_behavior_meta.pkl` | One per animal. All day-invariant behavior fields — `anmID`, `recorded_days`, `shift_idx`, `twoCams`, `mouth`, `port_locs` (pre/post-shift lickport positions), `cameraCal`, `shiftInfo` — plus a `_split_index` listing the day part files. |
| `behavior/<ID>_behavior_day_<N>.pkl` | One per behavior day. The day-indexed Bpod payload: `mats` (session info / state machine), `trajs` (lick trajectories), `raw`. Written in the same order as `recorded_days`. |
| `somas/<ID>_somas.pkl` | Small per-animal metadata file: `daysRelativeToShift`, `roiScores`. |
| `somas/<ID>_somas_session_N.pkl` | One per imaging session (N runs 0–4). Per-ROI deconvolved calcium "impulse trace" activity aligned to trial frames, somatic ROIs. |
| `dendrites/<ID>_dendrites*.pkl` | Same structure, dendritic ROIs. Session indices run 0–7. |

**Behavior is stored split, not monolithic.** A single `<ID>_behavior.pkl` was too large to commit, so
each animal's behavior dict is written as one small `_meta` file plus one file per behavior day
(`jsm.save_behavior_split()`), after `jsm.reduce_behavior_for_repo()` prunes the Bpod struct down to
the fields the analysis reads and downcasts float64 → float32. `jsm.load_behavior_split()` reassembles
a dict identical to the old monolithic one. `masterData_loading()` auto-detects the split set via the
`_meta` file (`jsm.find_behavior_split_extension()`) and falls back to a monolithic
`<ID>_behavior.pkl` if no `_meta` file is present, so both layouts work. The committed files are
uncompressed `.pkl`; the loaders also accept `.pkl.gz` / `.pkl.xz` / `.pkl.bz2`.

> `load_behavior_split(..., days=[...])` is for inspecting one day. A `days` filter makes `mats` /
> `trajs` / `raw` shorter than `recorded_days` and breaks the positional alignment `get_tracker()` and
> the rest of the analysis assume — leave it as `None` for reproduction.

`jsm.masterData_loading(DATA_DIR)` walks this tree and assembles the in-memory `masterData`
dictionary the notebooks operate on. The per-animal shift metadata (`shiftDay`, `shiftTrial`,
`shiftDir`) and the list of sessions to load are hard-coded inside that function, so the animal set is
fixed — you do not need a separate consolidated dataset file.

**Memory note:** `masterData_loading()` loads every session of every animal at once. Budget well
above the 5.3 GB on-disk size of the animal folders. Fig1, Fig4, Fig5 and Fig6 call it; Fig2 and Fig3
run entirely off `examples/` and `traceAlignments/`.

### `Data/examples/`

Small, pre-extracted slices of the raw data that let individual example panels be reproduced without
re-deriving them from the full dataset. The trace arrays are stored trimmed to only the ROIs/frames a
panel actually reads and are scattered back into full-shape arrays at load time. Three loading
patterns are used:

- `fig1_behav_examples/` — camera `.tiff` frames read with `tifffile`, plus a `rawBehavior.pkl` read
  with `jsm.unpickler()`. `Fig1.ipynb` cells 12–13 hard-code the example animal, day, trial indices
  and `camIDs`, so this panel needs nothing from the lab servers.
- `fig2_indiv_examples/` — a `_movie_preview_index.pkl` per example ROI, unpacked with
  `jsm.resolve()` and restored with `jsm.restore_trace()`, `jsm.restore_licks()` and
  `jsm.restore_image_frames()`.
- `fig3/4/5/6_indiv_examples/` — a `_meta` file, one `_traces_*` part file per example ROI, and a
  `_ROI_manifest.json` naming the ROIs, all read with `jsm.load_all_grouped_aligned_traces_split()`.
  The manifest is what lets these notebooks pick their example ROIs without touching the full
  `ROI_clustering` structure.

### `Data/traceAlignments/`

Cluster-level trace alignments for Fig3 and Fig5, each stored as a `_meta` file plus one part file per
cluster (`fig3/` 14 parts, `fig5/` 8, `fig3_licking/` and `fig5_licking/` 4 each). They are read back
with `jsm.load_ROI_clustering_reduced_split()` and `jsm.load_all_grouped_aligned_traces_split()`. Both
loaders take filter arguments (`ROI_type_keys`, `ROI_scores`, `clusterList`, `trialGroupings`,
`trace_types`, …) so a notebook only reads the parts it needs.

The matching *writers* (`save_ROI_clustering_reduced_split` /
`save_all_grouped_aligned_traces_split`) live in the lab-internal `jsm_figs` module and are **not**
shipped here — the committed notebooks only ever load. Regenerating these files requires the full
upstream dataset and the lab pipeline.

---

## Figures/

### `js_manuscript_final.py`

A single shared library (imported everywhere as `jsm`), 14,945 lines holding 162 functions. Its only
third-party imports are `numpy`, `scipy`, `matplotlib`, `seaborn`, `statsmodels`, `scikit-learn` and
`tqdm`; it also sets the manuscript-wide matplotlib rcParams at import time. Major functional groups:

- **Data loading / IO** — `masterData_loading()`, `pickler()`, `unpickler()`,
  `load_ROI_clustering_reduced_split()`, `load_all_grouped_aligned_traces_split()`,
  `resolve_align_trace_type()`: assemble `masterData` from `Data/`, and read back the split,
  optionally gzip/bz2/lzma-compressed alignment and example files.
- **Behavior export / reassembly** — `report_behavior_sizes()`, `reduce_behavior_for_repo()`,
  `save_behavior_split()`, `find_behavior_split_extension()`, `load_behavior_split()`: the
  list-aware pruning and per-day splitting that make the behavior data committable, and the loader
  `masterData_loading()` calls to put it back together.
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

| Notebook | Cells | Data it loads | Description |
|---|---|---|---|
| `Fig1.ipynb` | 38 | full `masterData` + `examples/fig1_behav_examples/` | Licking behavior. Session/trial count stats, example camera frames with the lick-angle construction, 3D lick projections and the LDA left/right boundary, per-animal lick-angle distributions, motor-error rate (MER) and directional-error rate (DER) across the sensorimotor shift. |
| `Fig2.ipynb` | 5 | `examples/fig2_indiv_examples/` | Reconstructs two example ROIs (one dendrite, one soma) from the trimmed exports and re-renders `image_and_trace_movie_preview()` frames — imaging field of view with ROI overlays alongside the simultaneous traces. Does not load `masterData`. The raw and registered source data for both example sessions (`B00002213784`, 230123 dendrites and 230124 somas) is in the companion [Scheib_2026_data](https://github.com/kerlinlab/Scheib_2026_data) repo. |
| `Fig3.ipynb` | 18 | `traceAlignments/fig3/`, `traceAlignments/fig3_licking/`, `examples/fig3_indiv_examples/` | Example-ROI NMF footprints, trial-aligned activity rasters by cluster, and the matching lick rasters. Does not load `masterData`. |
| `Fig4.ipynb` | 47 | full `masterData` + `examples/fig4_indiv_examples/` | Modeled responses; go-cue vs. lick-contact impulse response functions per ROI; SMO (Sensory / Motor-Choice / Outcome-Reward) coding-direction preparation, cross-validated calculation, projections, selectivity and correlations. **Uses a Dask `LocalCluster`** for the heavy SMO passes. |
| `Fig5.ipynb` | 54 | full `masterData` + `examples/fig5_indiv_examples/` + `traceAlignments/fig5*/` | Core motor-error / coding-direction figure: CR–CA, CR–CL and AP–CA coding-direction comparisons, licking data, example rasters, and a supplement section on raw-trial correlations. |
| `Fig6.ipynb` | 60 | full `masterData` + `examples/fig6_indiv_examples/` | Coding direction across the shift: pre/post selectivity index, CD projections at finer time bins, a CR–CL scrolling panel, and pre-late / post-early / post-late window comparisons including percent-change summaries. |

#### Dask (Fig4 only)

Two cells in `Fig4.ipynb` (cells 15 and 26) spin up a local Dask cluster:

```python
dask.config.set(temporary_directory=SCRATCH_DIR)   # SCRATCH_DIR = <repo>/scratch, set in the init cell
cluster = LocalCluster(n_workers=4, threads_per_worker=8, memory_limit="20GiB")
client  = Client(cluster)
```

`SCRATCH_DIR` is derived from `PARENT_DIR` in the init cell, so no path editing is needed —
`distributed` creates `<repo>/scratch/dask-worker-space/` itself on first use. Repoint it at a local
disk if the repository lives on a network filesystem. Scale `n_workers` / `threads_per_worker` /
`memory_limit` to your hardware (the defaults assume ~80 GB of usable RAM).

---

## Dependencies

Full pins are in [`requirements.txt`](requirements.txt) and [`Scheib2026_env.yml`](Scheib2026_env.yml).
Summary of what is actually imported:

| Package | Used by |
|---|---|
| `numpy`, `scipy`, `matplotlib` | everything |
| `seaborn`, `statsmodels`, `scikit-learn` | `js_manuscript_final.py`; `sklearn.mixture.GaussianMixture` also directly in Fig1 |
| `tqdm` | progress bars throughout |
| `tifffile` (+ `imagecodecs`) | Fig1's example camera `.tiff` reads. Every notebook imports it in the init cell, but only Fig1 calls it. |
| `pandas`, `pillow` | transitive (seaborn, matplotlib) |
| `dask`, `distributed` | Fig4 only |
| `ipykernel`, `ipython`, `notebook` | Jupyter runtime |
| `ipynbname` | every notebook's init cell, to derive `PARENT_DIR` from the notebook's own location |

Optional, not needed for any committed notebook (installed by `requirements.txt` /
`Scheib2026_env.yml` anyway, so nothing extra is required):

- `opencv-python` — only for `jsm.readFrames()` (`cv2` is not imported at module level).
- `joblib` — only for `jsm.run_parallel_bootstrap()` (`Parallel`/`delayed` are not imported at module
  level). `jsm.vectorized_bootstrap()` is the NumPy-only equivalent and is what the notebooks use.
- `ipympl` — no notebook uses `%matplotlib widget` any more; pinned only to keep the environment
  matching the reference one.

**NumPy must stay on the 1.x line.** `js_manuscript_final.py` installs a compatibility shim at import
time so that pickles written under NumPy 2.x can be read under NumPy 1.x:

```python
sys.modules.setdefault('numpy._core', numpy.core)
sys.modules.setdefault('numpy._core.multiarray', numpy.core.multiarray)
sys.modules.setdefault('numpy._core.numeric', numpy.core.numeric)
```

Installing NumPy ≥ 2 changes pickle module paths and defeats this shim.

---

## Known limitations

- **`Data/examples/fig2_hist_examples/` and `fig2_zstack_examples/` are unused by the committed
  notebooks.** `Fig2.ipynb` covers only the movie-preview panels; the histology and z-stack panels
  are not in the repo. The files are shipped anyway (367 MB) for anyone reproducing those panels by
  hand.
- **Running the notebooks writes two untracked directories into the repository root.**
  `figurePanels/` (only when `PDF_export_active = True`) and, for Fig4, `scratch/dask-worker-space/`.
  Neither is in `.gitignore`, so they will show up as untracked files in `git status`.
- **`Data/` is read-only from the notebooks' point of view.** The alignment and example files were
  written by lab-internal `jsm_figs` savers that are not part of this repo; the committed notebooks
  contain only the matching load calls, so running any notebook top-to-bottom will not overwrite
  anything under `Data/`. The one exception is `Fig1.ipynb` cell 13, which `os.mkdir`s the
  `examples/fig1_behav_examples/<anm>/<day>/run1` path it then reads from — harmless when the
  directory already exists.

---

====================================================================================================
## Raw data inventory

Due to the large size of the raw data files they are not included here in the repo but here is a summary of all files that are available upon special request.

Two sessions **are** published: the raw and motion-registered imaging data for the Fig2 dendrite and
soma examples (`B00002213784` 230123 and 230124) is in
**[kerlinlab/Scheib_2026_data](https://github.com/kerlinlab/Scheib_2026_data)**. For any other
session, contact newmanza@umn.edu.

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


## NOTES:
SLeD refers to the lab specific environment (Synapses Learning and Dendrites) however most of required
tools should be available in the single package file.

The preparation of this readme was performed with assistance from Claude claude-opus-5 however most
of the actual code was written by Zachary Newman, Jackson Scheib and Aaron Kerlin with minimal AI
assistance.

Please send any questions to newmanza@umn.edu.

## License

Released under the [MIT License](LICENSE).
