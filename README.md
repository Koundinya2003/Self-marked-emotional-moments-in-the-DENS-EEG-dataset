# Self marked emotional moments in the DENS EEG dataset

Author: Aditya Koundinya

My preregistered secondary analysis, written entirely in R, of OpenNeuro **ds003751 v1.0.6**, the *Dataset on
Emotion with Naturalistic Stimuli (DENS)* (Mishra, Asif, Tiwary and Srinivasan). Participants watched one-minute
film clips during 128-channel EEG and clicked whenever they felt an emotion.

**Research question.** Do the seconds before a participant marks an emotional moment show distinct EEG spectral
activity compared with unmarked periods of the same clip, and does that activity scale with the clip's rated
valence and arousal?

**Status: project complete (2026-10-19).** Download and audit, per-subject preprocessing, feature extraction,
preregistered mixed models, seven robustness analyses, five exploratory families, one post hoc analysis, figures and
a Quarto paper rendered to PDF and HTML (`paper/paper.pdf`, 9 pages). The whole pipeline reruns with
`Rscript code/run_all.R`.

This file is the single source of truth for the project. It absorbs the preregistration, the decision log, the data
audit, the self review and the final summary that previously lived in `docs/` and `PROGRESS.md`.

| Hypothesis | Test |
|---|---|
| H1 | Posterior alpha (8 to 13 Hz) lower in marked than in control epochs |
| H2 | Frontal alpha asymmetry (ln right minus ln left) in marked epochs increases with valence |
| H3a | Frontal theta (4 to 8 Hz) in marked epochs increases with arousal |
| H3b | Global beta (13 to 30 Hz) in marked epochs increases with arousal |

**Main findings:**

* 31 of 40 EEG files in ds003751 v1.0.6 are truncated; 211 of 432 clip presentations could be analysed
  (30 participants with features, 24 in H1, 21 in H2/H3).
* **None of the four preregistered hypotheses was supported.** Posterior alpha was lower before clicks
  (b = -0.121 ln units, 95% CI [-0.219, -0.023], uncorrected p = .015) but not after FDR correction (p = .062).
  Frontal asymmetry was unrelated to valence, and theta and beta were unrelated to arousal.
* The pre-click alpha decrease vanished in an earlier window and was left-lateralised over sensorimotor cortex
  (post hoc), consistent with preparing the right-hand click rather than with the felt emotion.

## Contents

1. [Important data caveat](#important-data-caveat)
2. [Repository map](#repository-map)
3. [How to reproduce](#how-to-reproduce)
4. [Data audit](#data-audit)
5. [Preregistration](#preregistration)
6. [Decision log](#decision-log)
7. [Analysis pipeline](#analysis-pipeline)
8. [Confirmatory results](#confirmatory-results)
9. [Robustness analyses](#robustness-analyses)
10. [Exploratory and post hoc results](#exploratory-and-post-hoc-results)
11. [Deviations from the preregistration](#deviations-from-the-preregistration)
12. [Self review](#self-review)
13. [Limitations](#limitations)
14. [Reference verification](#reference-verification)
15. [Licence and data](#licence-and-data)


## Important data caveat

In version 1.0.6, 31 of the 40 EEG binary files (`.fdt`) are truncated (verified against the git-annex checksums of
the official DataLad mirror). The analysis uses every clip whose full window lies inside the recorded samples; the
full audit is in [Data audit](#data-audit) below.

## Repository map

| Path | Content |
|---|---|
| `code/00_setup_environment.sh` | builds the R environment (micromamba, conda-forge), Quarto, LaTeX |
| `code/01_download_metadata.R` | file inventory, metadata download, v1.0.6 checksums from the DataLad mirror |
| `code/02_data_audit.R` | audit tables (`results/audit_*.csv`) |
| `code/03_preprocess.R`, `code/03_run_parallel.sh` | per subject pipeline: download, filter, bad channels, ICA, epochs, spectra, heart rate |
| `code/04_features.R` | tidy feature tables and QC summary |
| `code/05_confirmatory.R` | preregistered mixed models, FDR, diagnostics |
| `code/06_robustness_exploratory.R` | robustness analyses R1 to R7 and exploratory analyses E1 to E5 |
| `code/07_figures.R` | figures (PNG and PDF) |
| `code/08_render_paper.R` | renders `paper/paper.qmd` (IEEE format, IEEEtran class) |
| `code/09_checks.R` | consistency checks between results files and the paper |
| `code/run_all.R` | runs everything |
| `code/functions/` | I/O, events, montage and regions, preprocessing, epochs, ECG, models |
| `results/` | audit tables, QC summary, features, model results (small CSV); `results/subjects/` per subject outputs |
| `figures/` | publication figures; `figures/diagnostics/` model diagnostics |
| `paper/` | `paper.qmd`, rendered `paper.pdf` and `paper.html`, `references.bib` |
| `env/environment.yml`, `renv.lock` | pinned software versions |

## How to reproduce

```bash
# 1. environment (R 4.4 + packages from conda-forge, Quarto, LaTeX)
code/00_setup_environment.sh
export PATH=/opt/quarto/opt/quarto/bin:/opt/mm/envs/dens/bin:$PATH
# 2. everything (download, audit, preprocessing, features, statistics, figures, paper)
Rscript code/run_all.R
```

Preprocessing downloads one subject at a time (about 5 GB in total over the run, at most two raw files on disk at once)
and deletes each raw file after use. With `results/subjects/*.rds` present (committed), steps 04 to 09 run without
the raw data: `Rscript code/run_all.R --from 4`.

**Software.** R 4.4.3, lme4 2.0.6, lmerTest 3.2.1, fastICA 1.2.8, signal 1.8.1, R.matlab 3.8.0, Quarto 1.5.57,
pdf-engine xelatex with the `IEEEtran` class and the `TeX Gyre Termes` main font. Exact versions are pinned in
`env/environment.yml` and `renv.lock`.


## Data audit

Produced by `code/01_download_metadata.R` and `code/02_data_audit.R`; all tables referred to are in `results/`.
Only metadata, `.set` headers and (for one subject) one `.fdt` file were inspected at this stage; no EEG outcome
measure was computed.

### Sources and version

* Files were downloaded anonymously from the public OpenNeuro S3 bucket
  (`https://s3.amazonaws.com/openneuro.org/ds003751/`). The bucket holds 204 objects, 5.06 GB in total
  (`results/file_inventory.csv`).
* `openneuro.org` (and therefore its download API) was blocked by the network policy of the build machine.
* The DataLad mirror `github.com/OpenNeuroDatasets/ds003751` was reachable with `git`. Its default branch `main`
  carries tag `1.0.6` (CHANGES: "1.0.6 2023-07-08, How to cite?").
* Every annexed `.set`/`.fdt` file in 1.0.6 has exactly the size of the S3 object (80 of 80) and the MD5 in the
  git-annex key matches the downloaded file (`results/annex_md5.csv`; every download is MD5 verified). The S3 copies
  of `README`, `CHANGES` and `dataset_description.json` are older (1.0.2) and were replaced by the 1.0.6 copies.
  All other small files (participants, events, beh, json sidecars) are byte identical between S3 and 1.0.6.

### Participants: 38 or 40?

| Source | Count |
|---|---|
| DENS data paper (Mishra et al., bioRxiv 2021): recorded / retained after excluding 3 for excessive movement | 43 / 40 |
| Subject folders with EEG in ds003751 | **40** |
| Rows in `participants.tsv` | 55 (54 unique IDs; `sub-mit081` listed twice) |
| IDs in `participants.tsv` without any data folder | 14 |
| Subjects with a behavioural ratings file (`beh/*_beh.tsv`) | **34** |
| Subjects with EEG that can be analysed after the truncation problem | **32** with at least one recorded video that contains a click |

Resolution: the dataset contains EEG for 40 participants (the number the authors report after exclusions).
`participants.tsv` is a superset that also lists people without data. Six EEG subjects have no ratings file
(`sub-mit074, 079, 080, 104, 106, 107`), so only 34 have ratings; the DENS paper itself mentions that "30 out of 34
remaining participants" answered a post experiment question. The "38" figure in some secondary sources could not be
traced to a definitive origin and does not correspond to any count in the files; the definitive counts are the
ones above. Demographics of the 40 subjects with EEG: mean age 23.3 years (SD 1.4), 3 female, all right handed.

### Recording parameters (all 40 subjects identical)

| Property | Value |
|---|---|
| Format | EEGLAB `.set` (MATLAB v5 MAT file, read with `R.matlab`) plus `.fdt` (float32, channels x samples) |
| Sampling rate | 250 Hz |
| Channels | 132: E1 to E128 (EGI HydroCel GSN 128 scalp electrodes), E129 (vertex, recording reference Cz, all zeros), ECG, EMG, EMG_2 |
| Reference | Cz / E129 (json sidecar: "EEGReference": "Cz"; `.set` ref field: E129) |
| Channel coordinates | template HydroCel GSN 129 positions, identical for every subject (`results/channel_locations.csv`) |
| Power line | 50 Hz (json). Line noise is strong in the raw data (50 Hz peak 40 to 50 dB above neighbouring bins) |
| Software filters | none ("n/a"); data are raw with large DC offsets |
| Recording length (header) | 22.1 to 54.0 min (median 29.6) |


### Critical finding: truncated EEG binaries in v1.0.6

For 31 of 40 subjects the `.fdt` file is shorter than the header (`pnts x nbchan x 4` bytes) implies. This is the
published state of v1.0.6 (checksums match the annex keys), not a download error. EEGLAB writes samples interleaved
by channel, so a truncated file holds all channels for the first `floor(bytes / 528)` samples. Inspection of
`sub-mit007` showed plausible continuous EEG up to the cut, which falls mid sample (456 remainder bytes), consistent
with an interrupted upload.

* 9 subjects have complete files (`mit003, 004, 061, 076, 108, 117, 121, 123, b2017007`).
* 7 subjects have only 9 to 37 s of data (`mit054, 067, 072, 111, 113, 114, 116`) and contribute nothing.
* The remaining 24 have between 6.0 and 31.3 minutes of recorded data.
* Of 432 video presentations, **211 lie completely inside the recorded data**; of 744 clicks, **261** do.

Consequence: the analysis uses every video whose full 60 s window lies inside the recorded samples. Videos are
lost from the end of the session, so later trials (and therefore some stimuli more than others, depending on the
randomised order) are under represented. This is a missing-at-random-by-position mechanism unrelated to the EEG
signal itself, but it reduces power and is a central limitation (decision D07).

### Event structure and what marks emotional moments

`*_events.tsv` equals the EEGLAB `EEG.event` structure for all 40 subjects. **The `onset` column is in samples at
250 Hz, not in seconds** (verified: click minus video onset, divided by 250, reproduces `MouseClick` in `beh.tsv`;
median of the per video maximum absolute difference 0.008 s over the 248 videos with clicks in both sources).

Trial sequence: `trno` (trial number) -> `fixL` (fixation) -> `stm` (video start; label `<stimulus>_<trial>`) ->
`clic` (zero or more mouse clicks during the video) -> `vlnc`, `arsl`, `dmns`, `lkng`, `fmrt`, `relv` (rating
screens) -> `cate` (emotion category per click, only when the video had clicks). Other codes: `base`/`bend` (rest
baseline), `IBEG`/`IEND` (impedance checks), `puse`/`boundary` (pauses removed from the recording; never inside a
video), `rsrt` (logged "restart" marker; 18 cases fall inside a video window with no data discontinuity), `quiz`,
`qdon`, `stop`, `net`.

* Videos last 60.5 s (from `stm` to `vlnc`; range 58.5 to 60.7 s). One video (`sub-mit108`, `24_11`) has no rating
  screen because the recording ended; it is excluded.
* Each session had 9 to 11 videos (two neutral, `neutral_1`, `neutral_2`, plus emotional clips).
* **Emotional moments** are the `clic` events: participants were instructed to click "only when they felt any
  discernible emotion" during the video (Mishra et al., 2021). After the ratings, each click was assigned an emotion
  category and an intensity ("not at all" to "very much", `FivePointScale`). All 744 clicks occur inside videos
  (703 in emotional videos, 41 in neutral videos); 0 to 10 per video, 6 to 54 per subject.
* **Timing of the felt moment**: 30 of 34 participants reported clicking 2 to 3 s after becoming aware of the
  emotion (Mishra et al., 2021). The marked epoch is therefore placed before the click (see
  [Preregistration](#preregistration)).
* Click times in `beh.tsv` match the EEG events to within 0.017 s except for one video (`sub-mit072`, `8_9`) with a
  constant 2.04 s offset; the EEG events are used throughout.

### Ratings

Per video ratings (in `beh.tsv`, matched to EEG trials by stimulus name): valence, arousal, dominance (continuous
1 to 9), liking, familiarity (discrete 1 to 5), relevance (continuous 1 to 5), emotion category and intensity per
click. Across all rated videos: mean valence 4.01, mean arousal 5.49. In the analysable marked epochs valence and
arousal are nearly uncorrelated (r = 0.03 across videos), which helps separate H2 from H3.

### Feasibility for the research question

Applying the preregistered epoch rules to the recorded videos, before any EEG feature was computed:

| Quantity | Count |
|---|---|
| Recorded videos containing at least one click | 124 |
| Marked epochs (pre-click windows meeting the rules) | 237 in 123 videos of 32 subjects |
| Control epochs from those same videos | 1025 |
| Marked epochs in videos with ratings (for H2, H3) | 173 in 95 videos of 26 subjects |

The event structure therefore supports the original question, with reduced power. No redesign was needed (D09).

### Peripheral channels

ECG (one lead) and two EMG channels (placement not documented) are present. The ECG is used for the exploratory
heart rate analysis and as a reference for cardiac ICA components; EMG is not analysed because its placement is
unknown.


## Preregistration

**Status: locked.** The preregistration was committed on its own with the message "Preregistration, locked before
analysis" (commit `ac8815b`) before any EEG outcome measure (band power, asymmetry or heart rate in any epoch) was
computed. Every later change is recorded as a deviation in the [Decision log](#decision-log) and in the paper.

What had been looked at before this point: dataset metadata, events and ratings files, the `.set` headers, and a
QC-only preprocessing pilot of two subjects (`sub-mit003`, `sub-mit039`) that produced bad channel lists, ICA
component metrics and channel-median spectra of whole video periods, without epochs, conditions or ratings
(decisions D11 to D13). Epoch counts per condition were computed from event times only.

### Study information

**Title.** Pre-click EEG signatures of self-marked emotional moments during naturalistic film viewing: a secondary
analysis of the DENS dataset (OpenNeuro ds003751 v1.0.6).

**Authors.** Analysis carried out autonomously by an AI research agent (Claude Code) for the repository owner.
Dataset authors: Mishra, Asif, Tiwary and Srinivasan.

**Description.** Participants watched about 60 s emotional and neutral film clips while 128-channel EEG and ECG were
recorded, and clicked a mouse whenever they felt a discernible emotion. They then rated each clip on valence,
arousal, dominance, liking, familiarity and relevance. We ask whether the seconds that precede a click (the
self-marked emotional moment) differ in EEG spectral power from unmarked periods of the same clip, and whether EEG
activity in marked moments scales with the clip's rated valence and arousal.

**Hypotheses (confirmatory family, four tests).**

* **H1 (alpha and engagement).** Posterior alpha power (8 to 13 Hz) is lower in marked epochs than in control epochs
  of the same videos. Directional prediction: negative condition effect (marked minus control) on ln alpha power.
* **H2 (frontal asymmetry and valence).** In marked epochs, frontal alpha asymmetry (FAA = ln alpha at the right
  frontal cluster minus ln alpha at the left frontal cluster) is positively related to the video's valence rating
  (higher valence, relatively more left frontal activity, hence larger FAA).
* **H3a (theta and arousal).** In marked epochs, frontal theta power (4 to 8 Hz) is positively related to the
  video's arousal rating.
* **H3b (beta and arousal).** In marked epochs, global beta power (13 to 30 Hz) is positively related to the video's
  arousal rating.

All four are tested two sided; the direction is part of the prediction and is reported.

### Design plan

* **Study type:** secondary analysis of existing observational data (within-subject design).
* **Blinding:** not applicable; analysis code is fixed before outcomes are computed.
* **Design:** within-subject comparison of epoch types (marked versus control) inside the same video presentation, and
  trial-level association of EEG measures with self-report ratings.
* **Randomisation:** video order was randomised by the original authors; nothing is randomised here except the ICA
  initialisation (fixed seed 1).

### Sampling plan

* **Existing data:** yes. Data were collected and published by others (OpenNeuro ds003751 v1.0.6). The analyst had not
  analysed the EEG of this dataset before; prior knowledge is limited to what is listed above.
* **Sample size:** all 40 participants with EEG. Because 31 of 40 `.fdt` files are truncated in v1.0.6, only videos
  whose complete window lies inside the recorded samples are used. Expected before artefact rejection: 32 subjects,
  123 videos, 237 marked and 1025 control epochs for H1; 26 subjects, 95 videos and 173 marked epochs for H2 and H3
  (6 subjects have no ratings file).
* **Sample size rationale:** fixed by data availability, so no a priori power analysis is meaningful. No sensitivity
  simulation is planned; every estimate is reported with its 95% confidence interval and non-significant results are
  interpreted through those intervals.
* **Stopping rule:** not applicable.


### Variables

**Manipulated variables.** None (observational). The factor of interest, **condition**, is defined from events. All
times relative to video onset (`stm` event). Video length = time to the first `vlnc` event (about 60.5 s), capped
at 60 s.

* **Usable span** of a video: from 2.0 s after onset to 0.5 s before its end.
* **Marked epoch:** [click - 4.5 s, click - 0.5 s] (4 s). Kept only if the window lies fully in the usable span and
  no other click occurred in the 5 s before this click. Rationale: participants reported clicking 2 to 3 s after
  becoming aware of an emotion, so the window is centred 2.5 s before the click and ends before movement execution.
* **Control epoch:** 4 s windows of the same video, non-overlapping, outside every exclusion zone
  [click - 7 s, click + 3 s], packed consecutively from the start of each allowed stretch of the usable span.
* **Marked early epoch (robustness only):** [click - 6.5 s, click - 2.5 s], no other click in the preceding 7 s.
* For H1 only videos with at least one retained marked epoch contribute control epochs. Clicks in neutral videos are
  treated like any other click (the participant marked a felt emotion).

**Measured variables.**

* EEG: 128 scalp channels plus the reference E129, 250 Hz.
* ECG (one lead).
* Ratings per video: valence, arousal (1 to 9). Other ratings are used only in exploratory analyses.

**Preprocessing (identical for every subject; code in `code/functions/`).**

1. Read `.set`/`.fdt` in R; use the recorded samples. Crop from 10 s before the first usable video to 5 s after the
   last. Remove channel means.
2. Band pass: zero-phase windowed-sinc FIR (Hamming), order 826, passband 1 to 40 Hz (-6 dB at 0.5 and 45 Hz),
   applied separately to each segment between `boundary` events. No notch filter (the low pass attenuates 50 Hz by
   69 dB). No downsampling (250 Hz fits in memory).
3. Bad channels (E1 to E128; E129 is the flat reference): flat (median 2 s window SD < 0.5 uV), deviant (robust z
   of log SD > 3), or poorly correlated (median over 2 s windows of the median correlation with the 6 nearest
   neighbours < 0.4).
4. Average reference over good channels including E129 (which recovers the vertex signal).
5. ICA: FastICA (logcosh, parallel), PCA reduced to min(70, rank) dimensions, trained on every second sample of all
   1 s windows of the cropped data whose maximum absolute amplitude is below 500 uV. Components are removed if:
   * eye: |r| >= 0.7 with the VEOG derivation (mean E8, E25 minus mean E126, E127) or the HEOG derivation (E128 minus
     E125), or a z score >= 3 across components with |r| >= 0.3;
   * muscle: ratio of mean power 20 to 40 Hz to mean power 2 to 20 Hz >= 1 and either >= 50% of the squared mixing
     weights on outer ring channels or >= 40% on a single channel;
   * cardiac: |r| >= 0.5 with the band passed ECG.

   All metrics are computed on the ICA training samples.
6. Bad channels interpolated with spherical splines (Perrin et al., 1989; m = 4, 50 Legendre terms); final average
   reference over all 129 channels.
7. Epochs extracted as above. **Epoch rejection:** epoch rejected if the maximum absolute amplitude in any feature
   channel (all channels except the 38 outer ring channels E1, E8, E14, E17, E21, E25, E32, E38, E43, E44, E48, E49,
   E56, E57, E63, E64, E68, E69, E73, E74, E81, E82, E88, E89, E94, E95, E99, E100, E107, E113, E114, E119, E120,
   E121, E125, E126, E127, E128) exceeds 100 uV.

**Indices.**

* **Spectra:** Welch, 2 s Hann windows, 50% overlap (3 windows per epoch), 0.5 Hz resolution, uV^2/Hz.
* **Band power** per channel: PSD summed over bins times 0.5 Hz; theta [4, 8) Hz, alpha [8, 13) Hz, beta [13, 30]
  Hz (also delta [1, 4) and total [1, 40] for exploratory use).
* **Region power** = mean of channel band power over the region's channels, then natural log.
* **Regions** (anchor 10-20 equivalents plus all non outer ring channels within 3.0 cm, template coordinates):
  * frontal (F3 E24, Fz E11, F4 E124): E3 E4 E11 E16 E19 E20 E23 E24 E27 E28 E117 E118 E123 E124
  * central (C3 E36, Cz E129, C4 E104): E7 E29 E30 E31 E35 E36 E41 E55 E80 E103 E104 E105 E106 E110 E111 E129
  * posterior (P3 E52, Pz E62, P4 E92, O1 E70, Oz E75, O2 E83): E42 E47 E51 E52 E60 E62 E65 E67 E70 E72 E75 E77
    E83 E85 E90 E92 E93 E97 E98
  * left lateral (F7 E33, T7 E45, P7 E58): E33 E45 E50 E58 E65
  * right lateral (F8 E122, T8 E108, P8 E96): E90 E96 E101 E108 E122
  * FAA left (F3 E24): E19 E20 E23 E24 E27 E28; FAA right (F4 E124): E3 E4 E117 E118 E123 E124
  * global: all 91 feature channels (non outer ring E1 to E129).
* **FAA** = ln(mean alpha, right cluster) - ln(mean alpha, left cluster).
* **Heart rate (exploratory):** R peaks from the ECG (5 to 20 Hz FIR, polarity by larger extreme, threshold 40% of
  the 99.5th percentile, 300 ms refractory); RR intervals outside 0.33 to 1.5 s or deviating more than 30% from the
  running median of 11 intervals are discarded; epoch HR = 60 / mean valid RR with midpoint in the epoch (at least 2).


### Analysis plan

**Statistical models (R, lme4 + lmerTest, REML, Satterthwaite df).**

* **H1:** `ln_alpha_posterior ~ condition + (1 + condition | subject) + (1 | stimulus) + (1 | subject:stimulus)`;
  epoch level; condition coded control = 0, marked = 1; marked and control epochs only.
* **H2:** `FAA ~ valence_z + (1 + valence_z | subject) + (1 | stimulus) + (1 | subject:stimulus)`; marked epochs with
  ratings.
* **H3a:** `ln_theta_frontal ~ arousal_z + (1 + arousal_z | subject) + (1 | stimulus) + (1 | subject:stimulus)`.
* **H3b:** `ln_beta_global ~ arousal_z + (same random structure)`.
* `valence_z`, `arousal_z`: ratings standardised across the analysed videos. The `subject:stimulus` intercept
  represents the video presentation (each stimulus is shown once per subject), so trial level predictors are not
  pseudo replicated across epochs. Random intercepts for subject and stimulus are always kept.
* **Convergence rule:** fit the model above; if it fails to converge or is singular, drop the random slope
  correlation, then the random slope, then any random intercept other than subject and stimulus whose variance is
  estimated as zero. The first model that converges without a singular fit is reported (if subject and stimulus
  intercepts alone are singular, that model is reported anyway).
* **Effect sizes:** H1: d = b / sqrt(sum of random intercept variances + residual variance) (Westfall, Kenny and
  Judd, 2014). H2, H3: standardised slope from the same model refitted with the outcome z scored (SD of outcome per
  SD of rating). 95% CIs from Satterthwaite t (Wald) for b; profile CIs reported if they differ materially.

**Inference criteria.** Benjamini-Hochberg FDR across the four confirmatory p values, q = 0.05. A hypothesis is
supported if its adjusted p < 0.05 and the effect has the predicted sign. A result with adjusted p >= 0.05 is
reported as not supported, interpreted through its CI.

**Data exclusion.**

* Subject excluded: no recorded usable video; more than 25 bad channels (about 20%); pipeline failure (logged).
* Epoch excluded: amplitude criterion (above); video without a rating screen (`sub-mit108`, `24_11`).
* H1 needs at least one retained marked and one retained control epoch in the same video for that video to
  contribute. H2/H3 need ratings.
* No outlier removal in the confirmatory models (see robustness).

**Missing data.** Missing ratings (6 subjects without beh file) exclude those subjects from H2/H3 only. Missing data
due to truncation are not imputed.

**Assumption checks.** Residual QQ plot, residuals versus fitted, QQ plots of random intercepts
(`performance::check_model` style plots saved to `figures/diagnostics/`). Log transformation is used for power
because band power is right skewed.

**Robustness analyses (for each confirmatory test).**

* R1 alternative filter: passband 0.5 to 40 Hz (same bad channels, same ICA weights and component selection).
* R2 stricter rejection: 75 uV absolute threshold.
* R3 excluding outlier subjects: subjects whose mean outcome (per condition for H1) lies more than 3 MAD from the
  median of subjects, or whose random effect conditional mode exceeds 3 SD.
* R4 early marked window [-6.5, -2.5] s (motor preparation control; H1 only, and H2/H3 on early windows).
* R5 complete-file subjects only (the 9 untruncated recordings).
* R6 intercept-only random structure (subject, stimulus, subject:stimulus).

Robustness results are descriptive; they are not FDR corrected and do not change the confirmatory verdict.

**R7 was added after the preregistration** as a consequence of deviation D20: it re-runs the confirmatory analyses
with the original whole-epoch rejection rule (100 uV) and no channel repair, so that the deviation can be judged
against the preregistered pipeline. It is a consequence of a deviation, not a preregistered robustness analysis.

**Exploratory analyses (clearly labelled, BH-FDR within each family).**

* E1: condition effect for every band (theta, alpha, beta) x region (frontal, central, posterior, left, right) and
  per-channel topographies of the marked minus control difference.
* E2: valence and arousal associations for every band x region in marked epochs; and with the trial-level difference
  (marked minus mean control of the same video).
* E3: heart rate: marked versus control epoch HR (H1 model structure) and the HR time course from -10 to +10 s around
  clicks.
* E4: click intensity (FivePointScale) and posterior alpha in marked epochs.
* E5: emotional versus neutral video epochs (control epochs only) for posterior alpha.

**Software.** R 4.4.3, packages pinned in `env/environment.yml` and `renv.lock`.


## Decision log

Every analytic and engineering decision, with its reason. Entries are numbered and dated. Decisions after the
preregistration commit that change a preregistered element are additionally flagged as **DEVIATION** and listed in
the paper (see [Deviations](#deviations-from-the-preregistration)).

| # | Date | Phase | Decision | Reason |
|---|------|-------|----------|--------|
| D01 | 2026-10-19 | Setup | R 4.4.3 and all R packages installed from conda-forge with micromamba (`env/environment.yml`), not from CRAN. | CRAN, Posit Package Manager and r-universe were blocked by the network policy (HTTP 403 on CONNECT). conda-forge was reachable and provides binary builds of the same CRAN releases. Ubuntu apt R (4.3) lacked `R.matlab` and `signal`. |
| D02 | 2026-10-19 | Setup | `eegUtils` is not used; filtering, re-referencing, bad channel detection, interpolation, ICA and topographic plotting are implemented directly in R with `signal`, `fastICA` and `interp`/`akima`. | `eegUtils` is not on conda-forge and CRAN/GitHub source downloads were blocked. Implementing the steps directly keeps the whole pipeline in R and every step is visible in `code/functions/`. |
| D03 | 2026-10-19 | Setup | Package versions are pinned by `env/environment.yml` plus an `renv.lock` snapshot generated from the conda library (`renv::snapshot(type = "all")`). | The task asks for renv. Because packages came from conda-forge, `renv::restore()` on a machine with CRAN access reproduces the same versions, and the conda file reproduces the exact binaries. |
| D04 | 2026-10-19 | Setup | Quarto 1.5.57 from its GitHub release; LaTeX from the Ubuntu archive. | Only these routes were reachable. |
| D05 | 2026-10-19 | Audit | Data are downloaded anonymously over HTTPS from the public OpenNeuro S3 bucket (`https://s3.amazonaws.com/openneuro.org/ds003751/`). | First source in the prescribed order, and it worked. openneuro.org itself was blocked by the network policy. |
| D06 | 2026-10-19 | Audit | Version 1.0.6 is identified through the DataLad mirror (git tag 1.0.6 on branch `main`); S3 subject files are identical to it (sizes and MD5). The 1.0.6 top-level documents replace the older S3 copies. Each EEG download is MD5 verified against the git-annex key. | The S3 bucket only keeps the latest object of each file and its top-level documents are from 1.0.2; the mirror proves which version the binaries belong to. |
| D07 | 2026-10-19 | Audit | 31 of 40 `.fdt` files are truncated in 1.0.6. The analysis uses the recorded part of every file and keeps every video whose complete window lies inside the recorded samples (211 of 432 videos, 261 of 744 clicks, 32 subjects with a usable clicked video). | No complete copy exists in any reachable source (the truncated files are the published, checksummed state). Using the recorded part is the only way to use most participants; the truncation depends only on position in the session, not on the EEG. The subset is stated prominently in the paper. |
| D08 | 2026-10-19 | Audit | Event onsets are treated as samples at 250 Hz; click timing is taken from the EEG `clic` events, not from `beh.tsv`. | Verified against `MouseClick` (section 5 of the audit); EEG events exist for all subjects, `beh.tsv` only for 34, and one video in beh has a 2 s offset. |
| D09 | 2026-10-19 | Design | The original research question is kept (no redesign). | 237 marked and 1025 matched control epochs are available before artefact rejection; 173 marked epochs have ratings. |
| D10 | 2026-10-19 | Design | Marked epoch = [click - 4.5 s, click - 0.5 s]; control epochs = non-overlapping 4 s windows of the same video outside [click - 7 s, click + 3 s]; first 2 s and last 0.5 s of each video excluded; a click is used only if no other click occurred in the preceding 5 s. | Participants reported clicking 2 to 3 s after awareness (Mishra et al., 2021), so the window is centred on the felt moment; ending 0.5 s before the click limits movement execution activity; the exclusion zone keeps controls free of click related motor and visual activity; the onset buffer avoids the video onset response. An earlier window [-6.5, -2.5] s is a preregistered robustness check against motor preparation. |

| D11 | 2026-10-19 | Preprocessing pilot | Before the preregistration, the preprocessing was piloted on two subjects (`sub-mit003`, complete; `sub-mit039`, truncated) in QC-only mode (no epochs, no band power per condition, no ratings touched). | To verify the R implementation and set QC thresholds on data quality grounds only, without access to any outcome. |
| D12 | 2026-10-19 | Preprocessing pilot | Bad channel neighbour criterion changed from "correlation with the neighbour mean" to "median pairwise correlation with the 6 nearest neighbours"; ICA component metrics are computed on the artefact free ICA training samples only. | Pilot: one channel at 2900 uV SD made 7 neighbours look bad, and gross artefact windows outside the training set made 38 of 40 components correlate with the EOG. |
| D13 | 2026-10-19 | Preprocessing pilot | ICA: FastICA with 70 components (PCA reduced); muscle rule = 20-40 Hz / 2-20 Hz power ratio >= 1 AND focal (>= 40% of squared weights on one channel) or peripheral (>= 50% on the outer ring) topography; eye rule = \|r\| >= 0.7 with VEOG/HEOG derivations or z >= 3 across components with \|r\| >= 0.3; cardiac rule = \|r\| >= 0.5 with ECG. | A simulation confirmed the ICA code recovers and removes a known source (r = 0.9999). In the pilot, a ratio threshold of 0.5 with 40 components flagged mixed components and removing them raised alpha power; 70 components plus a focality requirement removed clearly artefactual components (alpha retained at 0.95 and 0.91 of the pre-cleaning value in the two pilot subjects). Thresholds were fixed on QC grounds only. |
| D14 | 2026-10-19 | Preprocessing | No notch filter. Band pass FIR (Hamming, order 826, -6 dB at 0.5 and 45 Hz, passband 1 to 40 Hz) applied separately to each segment between boundary events. | The FIR attenuates 50 Hz by 69 dB, which removes the strong line noise (40 to 50 dB) below the feature range; a notch would be redundant. Segment wise filtering avoids transients at recording discontinuities. |
| D15 | 2026-10-19 | Design | Regions = anchor channels (10-20 equivalents on the HydroCel net) plus all non outer ring channels within 3.0 cm; "left" and "right" are lateral fronto-temporo-parietal regions (F7/T7/P7 and F8/T8/P8 anchors); FAA clusters are F3 (E24) and F4 (E124) plus neighbours within 3.0 cm (6 mirror symmetric channels each). | A geometric rule is reproducible; 3.0 cm gives 5 to 19 channels per region; frontal, central and posterior already cover the midline, so the left/right regions are defined laterally to be distinct from them. |
| D16 | 2026-10-19 | Preprocessing | Bug found and fixed before any outcome was produced: `montage.R` still had a 2.5 cm region radius while the preregistration specifies 3.0 cm. The run was stopped, outputs deleted, the radius set to 3.0 cm and a hard check added that stops the pipeline unless the regions equal the preregistered channel lists. Not a deviation (the preregistered value is used). | Region definitions must match the preregistration exactly. |
| D17 | 2026-10-19 | Paper | References verified against Crossref or DataCite records (see [Reference verification](#reference-verification)); MNE-Python not cited because its registry record could not confirm the author list. | The task forbids unverifiable citations. |
| D18 | 2026-10-19 | Preprocessing | All subjects run in two parallel worker processes (`code/03_run_parallel.sh`), each processing one subject at a time and deleting the raw file immediately after it is loaded. | 4 cores and 15 GB RAM allow two subjects at once (about 3 to 5 GB each); disk use stays below one raw file per worker. |
| D19 | 2026-10-19 | Preprocessing (bug fix) | If every periocular channel of one side of an EOG derivation was bad, the derivation was constant, its correlations were NA, and NA propagated into the component selection (all 70 components removed in 3 early subjects). Fixed: bad but non-flat periocular channels are used for the derivation when no good one exists; NA correlations are set to 0; the pipeline stops if any classification is NA. All subjects reprocessed. | Bug, not a change of rule. |

| D20 | 2026-10-19 | Preprocessing | **DEVIATION.** Epoch handling changed from "reject the epoch if any of the 91 scalp channels exceeds 100 uV" to epoch-level channel repair: scalp channels above 100 uV in an epoch are flagged; with 1 to 9 flagged channels (10%), all channels above 100 uV are spline interpolated within that epoch and the epoch is re-referenced and kept if no scalp channel then exceeds 100 uV; with more than 9 flagged channels the epoch is rejected. The preregistered rule is kept as an additional robustness analysis (R7), and the 75 uV robustness analysis (R2) uses the same repair logic at 75 uV. | QC only, before any statistics: the preregistered rule kept 0 of 138 marked and control epochs in the first three subjects reprocessed (sub-mit003, 007, 023), including a complete, low-noise recording, because sporadic single-channel transients (typical electrode "pops" of saline nets) exceed 100 uV in almost every 4 s epoch of 91 channels while the median channel SD is about 5.5 uV. Epoch-wise repair of a limited number of channels follows established practice (FASTER, autoreject) and keeps the 100 uV criterion. |
| D21 | 2026-10-19 | Preprocessing | **DEVIATION.** "Do no harm" safeguard for ICA: a component flagged by the eye, muscle or cardiac rule is removed only if removing it alone does not increase the median (over clean 1 s video windows) of the median scalp channel SD. | QC: in sub-mit007 and sub-mit023 some components flagged by the adaptive eye criterion were mixtures; removing them raised typical amplitudes and beta power (beta power after cleaning 1.7 and 2.8 times the pre-cleaning value). With the safeguard the cleaning ratios returned to about 1 (alpha 0.93 to 1.00). The rule uses EEG amplitude only and is applied identically to every subject. |
| D22 | 2026-10-19 | Preprocessing | **DEVIATION.** Spherical spline interpolation regularised with lambda = 1e-5 on the diagonal of G (the MNE-Python default). | Without regularisation the interpolation weights summed to 6 to 11 in absolute value (noise gain 1.3 to 2.3), so interpolated channels amplified neighbour noise and repaired epochs failed the amplitude check; with lambda = 1e-5 the absolute weights sum to 2 to 3 (gain 0.4 to 0.7). |
| D23 | 2026-10-19 | Statistics | Post hoc analysis PH1 added (labelled as not preregistered): alpha and beta lateralisation ln(C4 cluster) - ln(C3 cluster), and C3/C4 alpha separately, marked vs control, in the main and early windows. | The confirmatory H1 effect disappeared in the preregistered early window (R4) and was largest over the left lateral region (E1), a pattern expected if the pre-click window captures preparation of the right hand mouse click (all participants right handed). Testing this alternative explanation is necessary for interpreting H1; it is reported separately from confirmatory and preregistered exploratory results. |
| D24 | 2026-10-19 | Statistics | Effect size bug fixed before reporting: `westfall_d` returned NA when the residual row of `VarCorr` was compared with `==` (NA comparison); `%in%` is used instead. | Bug. |
| D25 | 2026-10-19 | Self review | Added conservative sanity analyses (participant-level paired t test for H1; clip-level models for H2/H3) to `code/09_checks.R` and the paper, reported as self-review checks, not as confirmatory tests. | Selected random structures were simplified by the preregistered rule (H1 lost the random slope; H2/H3a lost the clip presentation intercept), which can make Satterthwaite df anticonservative or pseudo-replicate clip ratings. |


## Analysis pipeline

1. **Download and audit** (`code/01_download_metadata.R`, `code/02_data_audit.R`): file inventory, MD5 verification
   against the DataLad annex keys, participant/event/rating audit tables.
2. **Preprocessing** (`code/03_preprocess.R`, run by `code/03_run_parallel.sh`): per subject, crop to the usable
   span, remove channel means, band pass filter, detect bad channels, average reference, FastICA with rule-based
   component selection, spherical spline interpolation of bad channels, final average reference, extract marked and
   control epochs, epoch-level channel repair and rejection, Welch spectra per epoch per region, ECG heart rate.
   Subjects whose outputs already exist in `results/subjects/` are skipped, so the step is resumable.
3. **Features** (`code/04_features.R`): tidy per-epoch and per-region tables (`results/features_epoch.csv`,
   `results/features_long.csv`) and the QC summary (`results/qc_summary.csv`).
4. **Confirmatory models** (`code/05_confirmatory.R`): the four preregistered mixed models with the preregistered
   convergence rule, Westfall d effect sizes, Benjamini-Hochberg FDR across the four p values, assumption-check
   plots saved to `figures/diagnostics/`.
5. **Robustness and exploratory** (`code/06_robustness_exploratory.R`): R1 to R7 and E1 to E5, plus the post hoc
   motor lateralisation analysis (D23).
6. **Figures** (`code/07_figures.R`): all figures in PNG and PDF.
7. **Paper** (`code/08_render_paper.R`): renders `paper/paper.qmd` with inline R that reads the results CSV files, so
   every number in the paper is computed, not typed.
8. **Checks** (`code/09_checks.R`): consistency checks between the results files and the rendered paper, and the
   conservative sanity analyses of D25.

## Confirmatory results

From `results/confirmatory.csv`. Benjamini-Hochberg FDR is applied across the four confirmatory p values at
q = 0.05. A hypothesis is supported only if its adjusted p < 0.05 **and** the effect has the predicted sign.

| Hypothesis | Estimate | 95% CI | Uncorrected p | FDR p | Effect size | N obs / subj | Supported |
|---|---|---|---|---|---|---|---|
| H1 posterior alpha, marked vs control | b = -0.121 | [-0.219, -0.023] | .015 | .062 | d = -0.16 | 451 / 24 | **No** |
| H2 FAA ~ valence | b = -0.015 | [-0.123, 0.094] | .789 | .789 | std. slope = -0.03 | 70 / 21 | **No** |
| H3a frontal theta ~ arousal | b = -0.063 | [-0.190, 0.064] | .324 | .433 | std. slope = -0.11 | 70 / 21 | **No** |
| H3b global beta ~ arousal | b = 0.061 | [-0.067, 0.190] | .323 | .433 | std. slope = 0.16 | 70 / 21 | **No** |

**No confirmatory hypothesis was supported.** H1 is the only one with an effect in the predicted direction: posterior
alpha was lower in marked than in control epochs, but the effect did not survive FDR correction across the family of
four tests. H2, H3a and H3b are null; their confidence intervals are centred near zero.

**Random structures selected by the preregistered convergence rule** (`results/model_selection.csv`): H1 ended at
`(1 | subject) + (1 | stim) + (1 | subject:stim)`; H2 and H3a ended at the minimal structure
`(1 | subject) + (1 | stim)`; H3b kept the full structure including the random slope.

### Self-review sanity checks

From `results/self_review_checks.csv`. These are conservative re-analyses, reported as checks and **not** as
confirmatory tests.

| Check | Estimate | 95% CI | p | N |
|---|---|---|---|---|
| H1 paired t test on participant means | -0.155 | [-0.290, -0.019] | .027 | 24 participants |
| H2 clip-level model (marked epochs averaged per clip) | -0.029 | [-0.124, 0.066] | .537 | 44 clips |
| H3a clip-level model | -0.082 | [-0.197, 0.033] | .157 | 44 clips |
| H3b clip-level model | 0.055 | [-0.060, 0.169] | .339 | 44 clips |

The participant-level test for H1 reaches the same qualitative conclusion as the preregistered model (effect in the
predicted direction, not significant after BH correction: with this p value the smallest of four p values would
adjust to about .11). The clip-level models for H2 and H3a are still null, so the null results are not an artefact of
the epoch-level pseudo-replication.


## Robustness analyses

From `results/robustness.csv`. Descriptive only: they are not FDR corrected and do not change the confirmatory
verdict.

| Analysis | H1 b (p) | H2 b (p) | H3a b (p) | H3b b (p) |
|---|---|---|---|---|
| R0 confirmatory | -0.121 (.015) | -0.015 (.789) | -0.063 (.324) | 0.061 (.323) |
| R1 alternative filter (0.5 to 40 Hz) | -0.168 (.015) | -0.030 (.406) | 0.005 (.912) | 0.036 (.519) |
| R2 stricter rejection (75 uV) | -0.152 (.006) | -0.061 (.314) | -0.004 (.951) | 0.033 (.566) |
| R3 excluding outlier subjects | -0.066 (.216) | -0.015 (.789) | -0.063 (.324) | 0.027 (.620) |
| R4 early marked window [-6.5, -2.5] s | -0.024 (.625) | 0.003 (.966) | -0.007 (.889) | 0.046 (.443) |
| R5 complete-file subjects only | -0.275 (.142) | -0.051 (.509) | -0.148 (.089) | 0.023 (.787) |
| R6 intercept-only random structure | -0.121 (.015) | -0.015 (.789) | -0.063 (.324) | 0.058 (.308) |
| R7 preregistered rejection rule (added after D20) | -0.092 (.165) | -0.115 (.361) | -0.013 (.875) | 0.060 (.431) |

Reading of the robustness table:

* The H1 effect survives the alternative filter, the stricter 75 uV threshold and the intercept-only structure, but
  **not** the early window (R4), not the outlier exclusion (R3), and only with a very wide interval in the nine
  complete-file subjects (R5, 7 participants, df 6). It is weaker but same-signed under the preregistered rejection
  rule (R7, b = -0.092, p = .165) and that estimate lies inside the confirmatory CI, so the D20 deviation did not
  create the pattern.
* R4 is the decisive robustness result: shifting the marked window 2 s earlier removes the effect entirely.
* H2, H3a and H3b are null in every single robustness analysis.

## Exploratory and post hoc results

From `results/exploratory.csv` (BH-FDR within each family) and `results/posthoc_motor.csv`.

**E1, condition effect by band and region** (24 participants, 451 epochs, 53 videos). The only surviving condition
effect within the family is the alpha decrease over the **left lateral** region (b = -0.148, 95% CI [-0.239, -0.057],
uncorrected p = .002, FDR p = .023). The right lateral region shows a similar but weaker effect (b = -0.128, p = .028,
FDR p = .139). Frontal, central, posterior and all theta and beta regions are null.

**E2, rating associations.** No valence or arousal association with band power survives FDR correction in any band
and region, either in marked epochs (E2a) or in the trial-level marked-minus-control difference (E2b). This is the
exploratory counterpart of the null H2 and H3 results.

**E3, heart rate.** Marked versus control epoch heart rate: b = -1.000 bpm, 95% CI [-2.053, 0.053], p = .063, N = 448,
30 participants. The HR change in the 0 to 5 s after the click versus the -10 to -6 s baseline: b = 0.016, p = .969,
N = 29. No heart rate difference between marked and control epochs.

**E4, click intensity.** Click intensity (FivePointScale) and posterior alpha in marked epochs: b = 0.054, 95% CI
[-0.028, 0.136], p = .195, N = 70. Null.

**E5, neutral versus emotional videos.** In control epochs only, neutral clips have **higher** posterior alpha than
emotional clips: b = 0.168, 95% CI [0.028, 0.308], uncorrected p = .025, N = 979, 30 participants, 128 videos.

**PH1, post hoc motor lateralisation (not preregistered, uncorrected).** Label ln(C4 cluster) - ln(C3 cluster), and
C3/C4 alpha separately, marked versus control, in the main and early windows.

| Measure | Main window b (p) | Early window b (p) |
|---|---|---|
| alpha lateralisation (ln C4 - ln C3) | 0.113 (.021) | 0.043 (.406) |
| beta lateralisation (ln C4 - ln C3) | 0.102 (.007) | 0.058 (.125) |
| ln alpha C3 cluster | -0.116 (.025) | -0.076 (.139) |
| ln alpha C4 cluster | -0.002 (.971) | -0.033 (.551) |

The alpha decrease is confined to the **left** C3 cluster and is accompanied by a leftward beta lateralisation, and
both disappear in the early window. This is the pattern expected if the pre-click window captures preparation of the
right-hand mouse click (all participants are right handed) rather than the felt emotion. This analysis was chosen
after seeing the confirmatory and robustness results, so it is reported as an interpretation aid and is never used to
claim support for a hypothesis.


## Deviations from the preregistration

Three deviations (D20, D21, D22) were made after the preregistration commit. All three were decided on **quality
control grounds before any statistical model was fitted** and before any condition contrast was computed; none was
made after seeing a hypothesis result. They are analyst degrees of freedom and are reported in the paper's Appendix A.

| # | Preregistered rule | What was done | Why | Consequence |
|---|---|---|---|---|
| D20 | Reject the whole epoch if any of the 91 scalp channels exceeds 100 uV | Epoch-level channel repair: 1 to 9 flagged channels are spline interpolated within the epoch, the epoch is re-referenced and kept if no scalp channel then exceeds 100 uV; more than 9 flagged channels reject the epoch | The preregistered rule kept 0 of 138 marked and control epochs in the first three subjects reprocessed, including a complete low-noise recording: single-channel electrode "pops" exceed 100 uV in almost every 4 s epoch of 91 channels while the median channel SD is about 5.5 uV. Epoch-wise repair of a limited number of channels follows established practice (FASTER, autoreject) and keeps the 100 uV criterion | The preregistered rule is reported as robustness analysis R7. Its H1 estimate (-0.092) has the same sign and lies inside the confirmatory CI, so the deviation did not create the H1 pattern |
| D21 | Remove every ICA component flagged by the eye, muscle or cardiac rule | "Do no harm" safeguard: a flagged component is removed only if removing it alone does not increase the median (over clean 1 s video windows) of the median scalp channel SD | In sub-mit007 and sub-mit023 some components flagged by the adaptive eye criterion were mixtures; removing them raised typical amplitudes and beta power (1.7 and 2.8 times the pre-cleaning value). With the safeguard the cleaning ratios returned to about 1 (alpha 0.93 to 1.00) | The rule uses EEG amplitude only and is applied identically to every subject. It may leave some artefact in the data, which mainly affects beta power (H3b) |
| D22 | Spherical spline interpolation without regularisation | lambda = 1e-5 on the diagonal of G (the MNE-Python default) | Without regularisation the interpolation weights summed to 6 to 11 in absolute value (noise gain 1.3 to 2.3), so interpolated channels amplified neighbour noise and repaired epochs failed the amplitude check; with lambda = 1e-5 the absolute weights sum to 2 to 3 (gain 0.4 to 0.7) | Makes the D20 repair numerically stable |

The preregistered analysis code was also used unchanged wherever possible: the region radius bug (D16), the NA EOG
correlation bug (D19) and the Westfall d bug (D24) were all fixed **before** any outcome was produced, and all
subjects were reprocessed after D19.


## Self review

A skeptical reviewer pass. Automated parts are in `code/09_checks.R` (all checks pass); conservative re-analyses are
in `results/self_review_checks.csv`.

### 1. Data leakage and circular analysis

* **Epoch selection independent of outcomes?** Yes. Epochs are defined from event times only (`code/functions/epochs.R`);
  rejection depends only on amplitude; the H1 analysis set depends only on rejection flags and condition labels.
  Ratings are not used to select H1 data.
* **Regions chosen after looking at effects?** No. Region channel lists were written into the preregistration and the
  pipeline stops if they differ (`check_regions`; check 5).
* **QC-driven deviations (D20 to D22) made with knowledge of results?** They were made after looking at amplitude
  distributions, rejection counts and cleaning ratios of the first processed participants, before any model was
  fitted. They are nevertheless analyst degrees of freedom and are labelled as deviations; the preregistered rule is
  reported (R7) and gives a similar H1 estimate.
* **Post hoc analysis.** The motor lateralisation analysis was chosen after seeing the confirmatory and robustness
  results. It is labelled as post hoc, uncorrected, and used only to interpret H1.
* **Rating standardisation.** Ratings are z scored across the analysed clips (one value per clip presentation, not per
  epoch), so clips with many epochs do not dominate the scaling.

### 2. Degrees of freedom and model structure

* H1: the preregistered rule ended at the intercept structure (random slope singular). Epoch-level Satterthwaite df
  (412) may be anticonservative. The participant-level paired t test (24 participants) gives the same conclusion.
* H2, H3a: the rule ended at the minimal structure, which drops the clip presentation intercept, so several marked
  epochs of one clip share a rating (pseudo-replication). Clip-level models are still null.
* H3b kept the full structure (df 12.8, appropriately conservative).
* Singular fits were reported, not hidden (`results/model_selection.csv`).

### 3. Mislabelled conditions or channels

* Marked epochs end exactly 0.5 s before their click; no control epoch intersects any click exclusion zone; marked
  and control epochs never overlap; no epoch starts in the first 2 s (checks 3).
* Event onsets were verified to be in samples (click minus onset matched `MouseClick` within 0.017 s in 247 of 248
  clips), and sample indices are converted with the 1-based EEGLAB convention minus the crop offset.
* Hemispheres: in the template coordinates +Y is left (E24 = F3 has Y = +5.1). FAA = ln(right cluster around E124) -
  ln(left cluster around E24), as preregistered. In topographies x = radius x sin(theta), which places E24 on the left.
* The reference channel E129 is all zeros in the raw files and is recovered by average referencing; it was excluded
  from bad channel detection, which would otherwise always flag it.

### 4. Numbers in the paper versus results files

* Every result number in `paper/paper.qmd` is produced by inline R from `results/*.csv` or the saved models; the four
  confirmatory estimates were string-matched in the rendered HTML (check 6), and BH adjusted p values were recomputed
  (check 1). Saved model objects reproduce the table estimates and N (check 2).
* Literature values typed in the text (43 recorded and 3 excluded participants; 30 of 34 participants clicking 2 to 3 s
  after awareness) are quoted from the DENS paper and cited.


### 5. Problems found during the review and fixed

| Problem | Fix |
|---|---|
| Region radius 2.5 cm in code versus 3.0 cm in the preregistration | Fixed before any output (D16); hard check added |
| NA EOG correlations made the pipeline remove all ICA components in 3 participants | Bug fixed, all participants reprocessed (D19) |
| Westfall d returned NA (NA comparison in `VarCorr`) | Fixed (D24) |
| QC figure showed 0% rejection for participants without marked epochs | Shown as missing |
| Heart rate text reported the overall minimum (after the click) as the pre-click dip | Text now reports the pre-click minimum and post-click maximum separately |
| Unicode micro sign rendered as `<U+00B5>` in LaTeX tables | Replaced with "uV" inside tables |
| Abstract over 250 words | Shortened to 250 |

### 6. Confidence in the findings

* **High confidence** in the data audit (truncation, sample-based onsets, participant counts): directly verified
  against checksums and timing.
* **High confidence** that there is no strong valence-FAA or arousal-theta/beta relationship detectable in this
  subset; **moderate confidence** in the null itself (small samples: 70 marked epochs, 21 participants).
* **Moderate confidence** that the H1 pattern reflects motor preparation: the post hoc lateralisation analysis is
  consistent with it but was chosen after seeing the data.
* **Moderate confidence** in artefact handling: R-based ICA rules without ICLabel, plus the three deviations D20 to
  D22. The preregistered rejection rule gives a similar but less precise H1 estimate.

### 7. What a human should check first

1. **Deviations D20 to D22**: epoch-level channel repair, the ICA do-no-harm safeguard and spline regularisation.
2. **ICA component selection quality**: inspect `figures/fig_qc_subjects.png` and, ideally, compare with ICLabel in
   EEGLAB for a few participants.
3. **The truncation finding**: worth reporting to the dataset authors / OpenNeuro, since complete files would double
   the usable data.
4. **The motor-confound interpretation** of H1 and the post hoc analysis (`results/posthoc_motor.csv`).
5. That no confirmatory claim is made beyond what `results/confirmatory.csv` shows (all four hypotheses not
   supported).

## Limitations

* **Motor confound built into the design:** every marked epoch precedes a click, no control epoch does. This is the
  main interpretive limitation and is discussed in the paper.
* **Heavy data loss:** truncated files in the published dataset, and about half of the epochs rejected even with
  channel repair. Only 21 participants contribute to H2/H3 and 24 to H1.
* **Automated ICA rules are not validated against expert labels:** the do-no-harm rule protects against harmful
  removals but cannot guarantee that all artefacts were removed. Beta-band results are the most vulnerable.
* **H1 residuals are heavy tailed:** conclusions were unchanged under outlier exclusion and in the participant-level
  test.
* **Missing at random by session position:** the truncation removes later trials, so some stimuli are under
  represented depending on the randomised order.
* **No power analysis:** the sample size is fixed by data availability, so the null results are statements about what
  this subset can detect, not proof of no effect.


## Reference verification

Every reference in `paper/references.bib` was checked on 2026-10-19 against its registry record: Crossref
(`https://api.crossref.org/works/<DOI>`, fetched through a web scraping service because Crossref and doi.org were
blocked from the build machine) or DataCite (`https://api.datacite.org/dois/<DOI>`) for the OpenNeuro dataset.
Title, author list and order, venue, year, volume, issue, pages and DOI were compared with the bib entry. **All 34
references matched.**

| Key | DOI | Registry | Result |
|---|---|---|---|
| klimesch2012 | 10.1016/j.tics.2012.10.007 | Crossref | matches |
| jensen2010 | 10.3389/fnhum.2010.00186 | Crossref | matches (article 186 from DOI suffix) |
| schubring2019 | 10.1111/psyp.13386 | Crossref | matches |
| decesarei2011 | 10.1111/j.1469-8986.2011.01204.x | Crossref | matches |
| davidson1992 | 10.1016/0278-2626(92)90065-T | Crossref | matches |
| harmonjones2010 | 10.1016/j.biopsycho.2009.08.010 | Crossref | matches |
| reznik2018 | 10.1111/psyp.12965 | Crossref | matches |
| smith2017 | 10.1016/j.ijpsycho.2016.11.005 | Crossref | matches |
| aftanas2001 | 10.1016/S0304-3940(01)01703-7 | Crossref | matches |
| ray1985 | 10.1126/science.3992243 | Crossref | matches |
| knyazev2007 | 10.1016/j.neubiorev.2006.10.004 | Crossref | matches |
| sonkusare2019 | 10.1016/j.tics.2019.05.004 | Crossref | matches |
| koelstra2012 | 10.1109/T-AFFC.2011.15 | Crossref | matches |
| perrin1989 | 10.1016/0013-4694(89)90180-6 | Crossref | matches |
| hyvarinen2000 | 10.1016/S0893-6080(00)00026-5 | Crossref | matches |
| piontonachini2019 | 10.1016/j.neuroimage.2019.05.026 | Crossref | matches |
| bates2015 | 10.18637/jss.v067.i01 | Crossref | matches |
| kuznetsova2017 | 10.18637/jss.v082.i13 | Crossref | matches |
| westfall2014 | 10.1037/xge0000014 | Crossref | matches |
| welch1967 | 10.1109/TAU.1967.1161901 | Crossref | matches |
| benjamini1995 | 10.1111/j.2517-6161.1995.tb02031.x | Crossref | matches |
| barr2013 | 10.1016/j.jml.2012.11.001 | Crossref | matches |
| poldrack2006 | 10.1016/j.tics.2005.12.004 | Crossref | matches |
| pfurtscheller1999 | 10.1016/S1388-2457(99)00141-8 | Crossref | matches |
| kreibig2010 | 10.1016/j.biopsycho.2010.03.010 | Crossref | matches |
| muthukumaraswamy2013 | 10.3389/fnhum.2013.00138 | Crossref | matches |
| saarimaki2021 | 10.3389/fnhum.2021.675068 | Crossref | matches |
| dmochowski2012 | 10.3389/fnhum.2012.00112 | Crossref | matches |
| delorme2004 | 10.1016/j.jneumeth.2003.10.009 | Crossref | matches |
| mishra2021dens | 10.1101/2021.08.04.455041 | Crossref (bioRxiv) | matches; content also read from the PDF |
| mishra2021theta | 10.1101/2021.08.05.455187 | Crossref (bioRxiv) | matches |
| asif2023 | 10.1109/ACCESS.2023.3266804 | Crossref search | matches (IEEE Access 11, 39913-39925) |
| dens2023openneuro | 10.18112/openneuro.ds003751.v1.0.6 | DataCite | matches (OpenNeuro, 2023) |

Left out: Gramfort et al. (2013, MNE-Python, 10.3389/fnins.2013.00267). Its Crossref record lists only the first
author, so the full author list could not be checked from a registry; the adaptive eye criterion is described in
words instead of being attributed by citation.

## Licence and data

Code: MIT. Data: ds003751 is distributed by OpenNeuro under CC0; raw data are never stored in this repository.

