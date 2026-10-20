# Self marked emotional moments in the DENS EEG dataset

A preregistered secondary analysis of the public **Dataset on Emotion with Naturalistic Stimuli (DENS)**
([OpenNeuro ds003751, version 1.0.6](https://openneuro.org/datasets/ds003751/versions/1.0.6); Mishra, Asif, Tiwary
and Srinivasan). Participants watched one minute film clips during 128 channel EEG and clicked a mouse whenever they
felt an emotion. The project asks whether the seconds before such a click carry a distinct EEG signature, and whether
that activity follows the clip's rated valence and arousal. The whole analysis is written in R.

**Status.** Analysis complete; none of the four preregistered hypotheses was supported. Paper:
[`paper/paper.pdf`](paper/paper.pdf) (9 pages, IEEE format), also as [`paper/paper.html`](paper/paper.html).

## Research question and hypotheses

Do the seconds before a participant marks an emotional moment show distinct EEG spectral activity compared with
unmarked periods of the same clip, and does that activity scale with the clip's rated valence and arousal?

<table>
<tr><th></th><th>Preregistered hypothesis</th><th>Prediction</th></tr>
<tr><td>H1</td><td>Posterior alpha (8 to 13 Hz) is lower in marked than in control epochs of the same clips</td><td>negative condition effect</td></tr>
<tr><td>H2</td><td>In marked epochs, frontal alpha asymmetry (ln right minus ln left) increases with valence</td><td>positive slope</td></tr>
<tr><td>H3a</td><td>In marked epochs, frontal theta (4 to 8 Hz) increases with arousal</td><td>positive slope</td></tr>
<tr><td>H3b</td><td>In marked epochs, global beta (13 to 30 Hz) increases with arousal</td><td>positive slope</td></tr>
</table>

All four are tested two sided with Benjamini Hochberg FDR correction across the family (q = 0.05). A hypothesis
counts as supported only if its adjusted p is below .05 and the effect has the predicted sign.

**Marked and control epochs.** Participants in the DENS study reported clicking 2 to 3 s after becoming aware of an
emotion, so the marked epoch is the 4 s window from 4.5 s to 0.5 s before a click (kept only if no other click fell in
the preceding 5 s). Control epochs are nonoverlapping 4 s windows of the same clip outside an exclusion zone from 7 s
before to 3 s after every click. The first 2 s and last 0.5 s of each clip are never used. A second marked window,
6.5 s to 2.5 s before the click, serves as a robustness check against motor preparation.

## Data

40 participants with EEG (mean age 23.3 years, 3 women, all right handed), 9 to 11 clips per session (two neutral),
128 scalp channels on an EGI HydroCel net plus the vertex reference, ECG and two EMG channels, recorded at 250 Hz
without software filtering. After each clip, participants rated valence, arousal and dominance (1 to 9), liking and
familiarity (1 to 5) and relevance, and assigned an emotion category and intensity to every click.

Problems found in the published files:

* **31 of the 40 EEG binaries (`.fdt`) are truncated** in version 1.0.6. The file sizes and MD5 sums match the
  git annex keys of the official DataLad mirror, so this is the published state, not a download error. Nine
  recordings are complete, seven hold only 9 to 37 s of data, and the rest stop between 6 and 31 minutes in.
* Event onsets in `*_events.tsv` are in samples at 250 Hz, not in seconds (verified against the click times in
  `beh.tsv`, which agree to within 0.017 s except for one clip with a constant 2.04 s offset).
* `participants.tsv` lists 55 rows (54 unique IDs, one duplicated) for 40 participants with EEG.
* Six participants with EEG have no ratings file, so only 34 have ratings.
* One clip (participant mit108, clip `24_11`) has no rating screen because the recording ended.

Because of the truncation, the analysis keeps every clip whose full 60 s window lies inside the recorded samples:
211 of 432 clip presentations and 261 of 744 clicks. Clips are lost from the end of each session, so the missing data
depend on position in the session rather than on the EEG. After preprocessing and epoch rejection, 30 participants
contribute features, 24 enter H1 and 21 enter H2 and H3.

## Workflow

<table>
<tr><td><code>code/01_download_metadata.R</code></td><td>file inventory, metadata download, MD5 verification against the DataLad mirror</td></tr>
<tr><td><code>code/02_data_audit.R</code></td><td>participant, event, click and rating audit tables (<code>results/audit_*.csv</code>)</td></tr>
<tr><td><code>code/03_preprocess.R</code>, <code>code/03_run_parallel.sh</code></td><td>per participant preprocessing: download, filter, bad channels, ICA, epochs, spectra, heart rate; resumable</td></tr>
<tr><td><code>code/04_features.R</code></td><td>tidy feature tables and QC summary</td></tr>
<tr><td><code>code/05_confirmatory.R</code></td><td>preregistered mixed models, FDR correction, model diagnostics</td></tr>
<tr><td><code>code/06_robustness_exploratory.R</code></td><td>robustness analyses R1 to R7, exploratory analyses E1 to E5, post hoc analysis PH1</td></tr>
<tr><td><code>code/07_figures.R</code></td><td>all figures as PNG and PDF</td></tr>
<tr><td><code>code/08_render_paper.R</code></td><td>renders <code>paper/paper.qmd</code>; every number in the paper is computed from the results files</td></tr>
<tr><td><code>code/09_checks.R</code></td><td>consistency checks between results files and the paper, plus conservative sanity analyses</td></tr>
<tr><td><code>code/run_all.R</code></td><td>runs steps 1 to 9 in order</td></tr>
<tr><td><code>code/functions/</code></td><td>I/O, events, montage and regions, preprocessing, epochs, ECG, models</td></tr>
</table>

**Preprocessing** (identical for every participant). Zero phase FIR band pass 1 to 40 Hz applied per segment between
boundary events (no notch; the low pass attenuates the 50 Hz line noise by 69 dB). Bad channels detected as flat,
deviant (robust z of log SD above 3) or poorly correlated with their six nearest neighbours. Average reference.
FastICA with up to 70 components; components removed by eye, muscle and cardiac rules (correlation with EOG
derivations, high frequency power ratio with focal or peripheral topography, correlation with the ECG), subject to a
safeguard that a component is only removed if doing so does not raise typical channel amplitudes. Bad channels
interpolated with regularised spherical splines. Epochs with 1 to 9 scalp channels above 100 µV have those channels
repaired within the epoch; epochs still above 100 µV, or with more than 9 such channels, are rejected.

**Features.** Welch spectra (2 s Hann windows, 50% overlap, 0.5 Hz resolution); band power per channel; region power
as the natural log of the mean over the region's channels. Regions are anchor electrodes (10 20 system equivalents on
the HydroCel net) plus all non peripheral channels within 3.0 cm, fixed in the preregistration and checked by the
pipeline. Frontal alpha asymmetry uses six channel clusters around F3 (E24) and F4 (E124).

**Models** (lme4, lmerTest, REML, Satterthwaite df):

* H1: `ln_alpha_posterior ~ condition + (1 + condition | subject) + (1 | stimulus) + (1 | subject:stimulus)`
* H2: `FAA ~ valence_z + (1 + valence_z | subject) + (1 | stimulus) + (1 | subject:stimulus)`
* H3a and H3b: the same structure with `arousal_z` predicting frontal theta and global beta

When a model is singular or fails to converge, the preregistered fallback drops the slope correlation, then the
slope, then any zero variance intercept other than subject and stimulus. The structures finally used are in
`results/model_selection.csv`.

## Results

<table>
<tr><th>Hypothesis</th><th>Estimate [95% CI]</th><th>p</th><th>FDR p</th><th>Effect size</th><th>Observations / participants</th><th>Supported</th></tr>
<tr><td>H1 posterior alpha, marked vs control</td><td>−0.121 [−0.219, −0.023]</td><td>.015</td><td>.062</td><td>d = −0.16</td><td>451 / 24</td><td>no</td></tr>
<tr><td>H2 frontal asymmetry ~ valence</td><td>−0.015 [−0.123, 0.094]</td><td>.789</td><td>.789</td><td>std. slope −0.03</td><td>70 / 21</td><td>no</td></tr>
<tr><td>H3a frontal theta ~ arousal</td><td>−0.063 [−0.190, 0.064]</td><td>.324</td><td>.433</td><td>std. slope −0.11</td><td>70 / 21</td><td>no</td></tr>
<tr><td>H3b global beta ~ arousal</td><td>0.061 [−0.067, 0.190]</td><td>.323</td><td>.433</td><td>std. slope 0.16</td><td>70 / 21</td><td>no</td></tr>
</table>

Source: `results/confirmatory.csv`. H1 is the only effect in the predicted direction; it does not survive FDR
correction. H2, H3a and H3b have intervals centred near zero.

**Conservative checks** (`results/self_review_checks.csv`, not confirmatory). A participant level paired t test for
H1 gives −0.155 [−0.290, −0.019], p = .027 (24 participants), which would also not survive correction. Clip level
models for H2, H3a and H3b (44 clips) remain null, so the null results are not an artefact of treating epochs from
the same clip as independent.

**Robustness** (`results/robustness.csv`, descriptive, not corrected):

<table>
<tr><th>Analysis</th><th>H1 b (p)</th><th>H2 b (p)</th><th>H3a b (p)</th><th>H3b b (p)</th></tr>
<tr><td>Confirmatory</td><td>−0.121 (.015)</td><td>−0.015 (.789)</td><td>−0.063 (.324)</td><td>0.061 (.323)</td></tr>
<tr><td>R1 filter 0.5 to 40 Hz</td><td>−0.168 (.015)</td><td>−0.030 (.406)</td><td>0.005 (.912)</td><td>0.036 (.519)</td></tr>
<tr><td>R2 stricter rejection, 75 µV</td><td>−0.152 (.006)</td><td>−0.061 (.314)</td><td>−0.004 (.951)</td><td>0.033 (.566)</td></tr>
<tr><td>R3 outlier participants excluded</td><td>−0.066 (.216)</td><td>−0.015 (.789)</td><td>−0.063 (.324)</td><td>0.027 (.620)</td></tr>
<tr><td>R4 earlier marked window</td><td>−0.024 (.625)</td><td>0.003 (.966)</td><td>−0.007 (.889)</td><td>0.046 (.443)</td></tr>
<tr><td>R5 complete recordings only</td><td>−0.275 (.142)</td><td>−0.051 (.509)</td><td>−0.148 (.089)</td><td>0.023 (.787)</td></tr>
<tr><td>R6 intercepts only</td><td>−0.121 (.015)</td><td>−0.015 (.789)</td><td>−0.063 (.324)</td><td>0.058 (.308)</td></tr>
<tr><td>R7 preregistered rejection rule</td><td>−0.092 (.165)</td><td>−0.115 (.361)</td><td>−0.013 (.875)</td><td>0.060 (.431)</td></tr>
</table>

The H1 effect holds under a different filter, a stricter threshold and simpler random effects, but disappears when
the marked window is moved 2 s earlier (R4) and when outlying participants are excluded (R3). H2, H3a and H3b are null
in every analysis.

**Exploratory and post hoc** (`results/exploratory.csv`, `results/posthoc_motor.csv`):

* E1, band by region: the only condition effect surviving FDR within the family is lower alpha over the left lateral
  region (b = −0.148, FDR p = .023).
* E2, ratings: no valence or arousal association with band power survives correction in any band or region.
* E3, heart rate: marked minus control −1.0 bpm [−2.1, 0.1], p = .063; no change after the click.
* E4, click intensity and posterior alpha: b = 0.054, p = .195.
* E5, neutral clips show higher posterior alpha than emotional clips in control epochs (b = 0.168, uncorrected
  p = .025).
* PH1 (chosen after seeing the results, uncorrected): the alpha decrease is confined to the left central (C3) cluster
  and comes with leftward beta lateralisation; both vanish in the earlier window. This matches preparation of the
  right hand mouse click more than the felt emotion.

## Limitations

* **Motor confound in the design.** Every marked epoch precedes a click and no control epoch does. The earlier window
  and the lateralisation pattern suggest the H1 effect reflects movement preparation.
* **Heavy data loss.** Truncated files in the published dataset, plus about half of the epochs rejected even with
  channel repair, leave 24 participants for H1 and 21 for H2 and H3.
* **Automated ICA rules** were not validated against expert labels or ICLabel. Residual artefact would mostly affect
  beta (H3b).
* **Missing data by session position.** Truncation removes later clips, so some stimuli are under represented
  depending on the randomised order.
* **No power analysis.** The sample is fixed by data availability; the null results describe what this subset can
  detect, not the absence of an effect.
* **H1 residuals are heavy tailed**; conclusions were unchanged under outlier exclusion and at the participant level.

## Deviations from the preregistration

All three were decided on quality control grounds before any statistical model was fitted, and are listed in
Appendix A of the paper.

1. **Epoch rejection (D20).** The preregistered rule (reject an epoch if any of 91 scalp channels exceeds 100 µV) kept
   0 of 138 epochs in the first three participants processed, because single channel electrode pops exceed 100 µV in
   almost every epoch. It was replaced by within epoch repair of up to 9 channels. The original rule is reported as
   R7 and gives a same signed H1 estimate inside the confirmatory interval.
2. **ICA safeguard (D21).** A flagged component is removed only if removing it does not raise typical channel
   amplitudes, after some flagged components proved to be mixtures whose removal inflated beta power.
3. **Spline regularisation (D22).** A small ridge term (λ = 0.00001) stabilises spherical spline interpolation, which
   otherwise amplified neighbour noise.

Three coding errors (region radius, missing EOG correlations, an effect size function returning NA) were found and
fixed before any outcome was produced; all participants were reprocessed afterwards.

The full preregistration text, decision log (D01 to D25) and parameter details are in the paper and its appendix.
The preregistration is recorded as committed before any EEG outcome was computed. The published history of this
repository is a single commit, so that earlier commit is not visible here and the timing cannot be checked from the
repository alone.

## Reproducing the analysis

Requirements: R 4.4.3 with lme4, lmerTest, fastICA, signal and R.matlab (versions pinned in `renv.lock`), Quarto
1.5.57 and XeLaTeX with the IEEEtran class for the paper.

```bash
bash code/00_setup_environment.sh
Rscript code/run_all.R
```

Preprocessing downloads one participant at a time from the public OpenNeuro S3 bucket (about 5 GB in total, at most
two raw files on disk at once) and deletes each raw file after use. Because the per participant outputs in
`results/subjects/` are committed, steps 4 to 9 can be rerun without the raw data; `code/run_all.R` documents the
option for starting at a later step.

**Known gaps.** `code/00_setup_environment.sh` builds the environment from `env/environment.yml`, which is not
included in this repository; until it is added, use `renv.lock` (for example with `renv::restore()`). The licence
below is stated but no `LICENSE` file is included yet.

## Repository map

<table>
<tr><td><code>code/</code></td><td>numbered pipeline scripts and <code>functions/</code></td></tr>
<tr><td><code>results/</code></td><td>audit tables, QC summary, feature tables, model results (CSV), fitted models (<code>models/</code>), per participant outputs (<code>subjects/</code>)</td></tr>
<tr><td><code>figures/</code></td><td>paper figures as PNG and PDF; model diagnostics in <code>diagnostics/</code></td></tr>
<tr><td><code>paper/</code></td><td>Quarto source, rendered PDF and HTML, bibliography and citation style</td></tr>
<tr><td><code>renv.lock</code></td><td>pinned R package versions</td></tr>
</table>

## Data, licence and citation

Data: Mishra, S., Asif, M., Tiwary, U. S. and Srinivasan, N. (2023). *Dataset on Emotion with Naturalistic Stimuli
(DENS)* (version 1.0.6). OpenNeuro. doi:10.18112/openneuro.ds003751.v1.0.6. Distributed under CC0; no raw data are
stored in this repository.

Dataset paper: Mishra et al. (2021), bioRxiv, doi:10.1101/2021.08.04.455041.

Code: MIT.
