# fel-lamp-calibration

Absolute radiometric calibration of an Ocean Optics USB4000 spectrometer
against a NIST-traceable FEL standard lamp (Eppley Laboratory, S.O.
52435, lamp EN-66; certificate dated 15 July 1992, 250-2400 nm,
calibrated at 50 cm, 7.90 A DC).

## Problem

The most intuitive model-selection criterion — lowest error on the
calibration data itself — picks the wrong interpolation model for this
lamp's spectral irradiance certificate. Only leave-one-out
cross-validation reveals it: a Gaussian process fit with in-sample RMS
close to the best model's (0.63% vs. 0.57%) generalizes about 3.5 times
worse under LOOCV (3.51% vs. 1.01%), with every parameter re-estimated
at each fold. This repository documents and automates that
check across six candidate models.

## Models

| # | Model | Parameters |
|---|-------|------------|
| 1 | Gray body (Planck function with free amplitude) | A, T |
| 2 | Gray body x linear emissivity | A, T, alpha |
| 3 | Gray body x tabulated tungsten emissivity | A, T |
| 4 | Model 3 x linear residual correction | A, T, beta |
| 5 | Two-stage log-linear, adapted from Huang, Cebula & Hilsenrath (1998) | A, T, h2-h6 |
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

```matlab
ABSOLUTE_IRRADIANCE_VALIDATION
```

Unlike `ABSOLUTE_CALIBRATION`, this script reads external inputs: for
each working distance, a SpectraSuite `.IrradCal` calibration file and
a raw High-Speed-Acquisition counts export (100 acquisitions) from the
same instrument. The `RUNS` table in Section 1 of the script points at
the files in `example_data/` for 50, 57 and 60 cm (the 57 cm run reuses
the 50 cm calibration), with the integration time of each acquisition;
edit it to point at your own files. Validation statistics are computed
between 400 and 800 nm (`STAT_BAND`); saturated pixels are flagged and
excluded, and the region below 300 nm is shaded in the figure.

## Outputs

| File | Description |
|------|-------------|
| `lampfile_EN66_Model5.lmp` | 250-2400 nm, 1 nm step, for SpectraSuite's Absolute Irradiance wizard |
| `model_comparison.csv` | Per-model metrics (k, RMS, AIC, AICc, ΔAICc relative to Model 5, LOOCV) |
| `model_comparison_color.png` / `.pdf` | Publication figure, six models and residuals |
| `model_comparison_bw.png` | Same figure, true grayscale (pixel-level RGB-to-luminance conversion) |
| `six_models_comparison*.png` | Comparison figures saved during development |
| `irradiance_results_all_distances.csv` | Per-pixel table for 50, 57 and 60 cm: wavelength, mean irradiance, both confidence bands, and saturation/out-of-range/negative-signal flags (`ABSOLUTE_IRRADIANCE_VALIDATION`) |
| `certificate_comparison_all_distances.csv` | Ratio of measured to certificate irradiance at the certificate wavelengths, per distance, with the expected (50/d)^2 (`ABSOLUTE_IRRADIANCE_VALIDATION`) |
| `distance_validation_color.png` / `.pdf`, `distance_validation_bw.png` | Measured-to-certificate ratio vs. wavelength at the three distances (`ABSOLUTE_IRRADIANCE_VALIDATION`) |

## Distance and the calibration

A SpectraSuite absolute-irradiance calibration converts detector counts
into the spectral irradiance at the fiber, wherever the fiber is. When
the calibration is built at a distance d other than the certificate's
50 cm, the certified irradiance must be scaled to that distance by
(50/d)^2, either through the wizard's Reference/Corrected Distance
fields or by pre-scaling the lamp file (both give the same result); the
coefficients written to the .IrradCal file then do not depend on d.

The example calibrations confirm this. Between 400 and 800 nm, outside
the band saturated at 50 cm, the coefficients of the 60 cm file agree
with those of the 50 cm file to within 1-5% in 50 nm bands (median 4%),
against the factor 1.44 expected if the scaling were missing. Measuring
the lamp at d with such a calibration therefore returns the certificate times
(50/d)^2: the ratio `ABSOLUTE_IRRADIANCE_VALIDATION.m` reports at 57 and
60 cm (0.770 and 0.698, against 0.769 and 0.694 predicted) is the
expected inverse-square attenuation, not a software fault. Do not
rescale measurements of other sources by (d/50)^2.

The 4% difference between the two independently built calibrations is
the reproducibility of the setup (alignment and distance). It is the
dominant term of the irradiance-scale uncertainty, about 3% at 1 sigma.

Known limitation: the 50 cm example scan, taken at the USB4000's minimum
integration time of 3.8 ms, saturates between about 533 and 680 nm, and
the coefficients of `example_calibration_50cm.IrradCal` are biased high
by 7-9% in that band. Prefer `example_calibration_60cm.IrradCal`, and do
not use the 50 cm file between 533 and 680 nm.

Correction notice: v1.2.1 described this behaviour as a failure of
SpectraSuite's Distance Correction to persist into the .IrradCal file
and recommended multiplying measured irradiance by (d/50)^2. That
interpretation was wrong and is withdrawn in v1.2.2.

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
