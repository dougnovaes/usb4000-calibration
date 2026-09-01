# wavelength-calibration

Wavelength-axis calibration of an Ocean Optics USB4000 spectrometer
(3648-pixel linear CCD) from Hg and Ne emission-lamp lines, fitting the
cubic polynomial lambda(pixel) = a0 + a1*p + a2*p^2 + a3*p^3 that the
instrument's own firmware accepts.

## Problem

A spectrometer's wavelength-calibration coefficients can appear
correctly written in the instrument's own configuration panel and still
not be the ones actually active — a coefficient can silently retain its
factory value after a write that looked successful. Standard validation
(comparing against a pre-extracted peak table) lacks the resolution to
catch this. This repository extracts the coefficients actually encoded
in a fresh acquisition's own wavelength axis and compares them directly
against what was meant to be written, before reporting any other
diagnostic.

## Method

- Peak positions are located to sub-pixel precision (a parabola through
  the three central points of each peak, not a whole-profile Gaussian,
  which is pulled by the instrument's asymmetric coma tail).
- A prominence filter rejects noise ripples that a height threshold
  alone would misidentify as peaks.
- Unresolved line blends are corrected to their effective wavelength
  rather than the catalogue value of the dominant component.
- The fit is weighted by each line's measured 90% bandwidth, which
  doubles as a blend detector (a broadened line gets less weight) and
  as a positional error bar.
- Hg lines are identified first, from a small, unambiguous, high-SNR
  set. A cubic fit to a robust subset of them — outliers removed by
  leave-one-out residuals, not raw residuals, which a high-leverage
  point can hide from — is then extrapolated to identify the entire
  pattern of Ne lines at once, never by nearest-line matching, which
  always finds a plausible match and so cannot expose a wrong
  hypothesis about the lamp's contents.
- After the coefficients are written to the instrument, a fresh
  acquisition is used to (a) recover the polynomial actually encoded in
  its own wavelength axis, independent of any assumption about what
  should be there, and (b) re-locate its peaks with the same sub-pixel
  method used for fitting, so the validation is not limited by the
  ~0.2 nm integer-pixel quantization of the instrument's own reported
  peak positions.

Full methodological detail and the specific failure cases that
motivated each step are documented as numbered notes in
`WAVELENGTH_CALIBRATION.m`'s header.

## Data flow

| Stage | Inputs | Outputs |
|-------|--------|---------|
| Fit | `hg_spectrum.txt`, `ne_spectrum.txt`, `hg_bandwidths.txt`, `ne_bandwidths.txt` | `calibration_coefficients.txt`, `calibration_residuals.csv` |
| Validation (post-write) | `validation_hg_spectrum.txt`, `validation_ne_spectrum.txt` | `validation_residuals.csv` |

`hg_spectrum.txt`/`ne_spectrum.txt` are two-column SpectraSuite exports
(wavelength under the factory calibration, counts), acquired under the
factory calibration. The validation spectra are the same format,
acquired after writing the new coefficients, under whatever calibration
is actually active on the instrument.

`Factory coefficients - BACKUP (NOTE 5).txt` is not read by the script;
it is a manual backup of the original factory calibration, recorded
before it was overwritten, kept for provenance. `obsolete_data/` holds
raw SpectraSuite exports and a superseded validation format not used by
the current script; see that folder for what each file was for.

## Reference result

Bandwidth-weighted fit, robust Hg anchor, 32 lines (254-744 nm):

```
a0 =  1.7635567179E+02
a1 =  2.1606649029E-01
a2 = -3.6985346136E-06
a3 = -5.5381923002E-10
weighted RMS = 0.0339 nm   (factory calibration: 1.0419 nm, ~31x worse)
```

## Requirements

Base MATLAB. No additional toolboxes.

## Usage

```matlab
WAVELENGTH_CALIBRATION
```

## Legacy

`legacy/` contains four earlier iterations, in order:

- `CALIBRACAO_COMPRIMENTO_ONDA_Ar.m` — the earliest version, using
  hardcoded peak positions and misidentifying the second lamp as argon.
- `CALIBRACAO_COMPRIMENTO_ONDA_Ne.m` — corrected to neon, with the full
  sub-pixel peak-detection pipeline.
- `CALIBRACAO_COMPRIMENTO_DE_ONDA.m` — adds bandwidth-based weighting
  and blend classification; the most complete Portuguese-language
  version.
- `WAVELENGTH_CALIBRATION_legacy.m` — the English translation, before
  the robust anchor (leave-one-out outlier rejection) and the
  post-write validation stage described above were added.
