# fel-lamp-calibration

Absolute radiometric calibration of an Ocean Optics USB4000 spectrometer
against a NIST-traceable FEL standard lamp (Eppley Laboratory, S.O.
52435, lamp EN-66; certificate dated 15 July 1992, 250-2400 nm,
calibrated at 50 cm, 7.90 A DC).

## Problem

The most intuitive model-selection criterion — lowest error on the
calibration data itself — picks the wrong interpolation model for this
lamp's spectral irradiance certificate. Only leave-one-out
cross-validation reveals it: a Gaussian process fit that looks best
in-sample (0.59% RMS) generalizes six times worse than the runner-up
(3.51% vs. 1.40% LOOCV). This repository documents and automates that
check across six candidate models.

## Models

| # | Model | Parameters |
|---|-------|------------|
| 1 | Black body (pure Planck function) | A, T |
| 2 | Black body x linear emissivity | A, T, alpha |
| 3 | Black body x tabulated tungsten emissivity | A, T |
| 4 | Model 3 x linear residual correction | A, T, beta |
| 5 | Two-stage log-linear (Huang, Cebula & Hilsenrath, 1998) | A, T, h2-h6 |
| 6 | Gaussian process regression (log-linear space) | 3 hyperparameters |

Each model is fit by weighted nonlinear least squares (weight = 1/F^2,
i.e. constant relative uncertainty) to the 31-point certificate, then
scored by in-sample RMS, AIC, AICc, BIC, and leave-one-out
cross-validation (LOOCV). LOOCV is the criterion that determines the
selected model (Model 5); the others are reported as supporting
evidence, along with a nested F-test for the two nested pairs (1 vs. 2,
3 vs. 4).

## Data sources

- Eppley Laboratory, Inc. *Certificate of Calibration of a Standard of
  Spectral Irradiance*, S.O. 52435, Lamp Serial No. EN-66, 1992.
- Tungsten emissivity (Models 3 and 4): Table III of Branstetter, J.R.,
  *Radiant Heat Transfer Between Nongray Parallel Plates of Tungsten*,
  NASA TN D-1088, 1961 — a public-domain compilation of De Vos (1954,
  Physica 20, 669) corrected for scattered light by Larrabee (1959, J.
  Opt. Soc. Am. 49, 619).
- Functional form of Model 5: Huang, L.K., Cebula, R.P., Hilsenrath, E.,
  *New procedure for interpolating NIST FEL lamp irradiances*,
  Metrologia 35, 381-386, 1998.
- Model-selection criteria: Burnham, K.P., Anderson, D.R., *Model
  Selection and Multimodel Inference*, 2nd ed., Springer, 2002.
- Model 6: Rasmussen, C.E., Williams, C.K.I., *Gaussian Processes for
  Machine Learning*, MIT Press, 2006.

## Requirements

MATLAB with the Optimization Toolbox (`lsqnonlin`) and the Statistics
and Machine Learning Toolbox (`fitrgp`, Model 6 only). Figure export
uses `print -RGBImage`, available from R2020a.

## Usage

```matlab
ABSOLUTE_CALIBRATION
```

No external input files are required — the certificate and the
tungsten emissivity table are hardcoded in the script (Section 2).

## Outputs

| File | Description |
|------|-------------|
| `lampfile_EN66_Model5.lmp` | 250-2400 nm, 1 nm step, for SpectraSuite's Absolute Irradiance wizard |
| `model_comparison.csv` | Per-model metrics (k, RMS, AIC, AICc, LOOCV) |
| `model_comparison_color.png` / `.pdf` | Publication figure, six models and residuals |
| `model_comparison_bw.png` | Same figure, true grayscale (pixel-level RGB-to-luminance conversion) |
| `six_models_comparison*.png` | Comparison figures saved during development |

## Legacy

`legacy/` contains the six earlier iterations of this analysis, in
order (`Calibration01.m` through `Calibration06.m`): a plot-only look
at the certificate (01); a single black-body-times-linear-emissivity
fit via the Curve Fitting Toolbox (02); a first three-model comparison
using illustrative, non-tabulated emissivity values (03); a four-model
comparison after switching to the real NASA TN D-1088 emissivity data
(04-05, near-identical, differing only in text encoding); and a
five-model comparison adding the Huang log-linear model (06). The final
script added Model 6 (Gaussian process regression), leave-one-out
cross-validation as the deciding criterion, and the publication-figure
and lamp-file output.
