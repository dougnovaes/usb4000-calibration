# usb4000-calibration

Two independent calibration pipelines for the same instrument: an Ocean
Optics USB4000 spectrometer, at the Institute of Physics, University of
Sao Paulo. `wavelength-calibration/` states its end application
directly in its own header: spectral characterization of a cold
dielectric-barrier-discharge (DBD) plasma. `fel-lamp-calibration/`
calibrates the same instrument's absolute radiometric response and
does not state a specific downstream application of its own.

- **[wavelength-calibration/](wavelength-calibration/)** — calibrates the
  pixel-to-wavelength axis of the spectrometer from Hg and Ne emission
  lamp lines.
- **[fel-lamp-calibration/](fel-lamp-calibration/)** — calibrates the
  absolute spectral-irradiance response of the same spectrometer against
  a NIST-traceable FEL standard lamp.

The two pipelines solve unrelated calibration problems (wavelength axis
vs. radiometric response) and are kept in one repository because they
calibrate the same physical instrument — in the wavelength-calibration
script's own words, the two are otherwise "unrelated ... beyond
sharing a project." They also share a leave-one-out cross-validation
discipline — used for model selection in one case and for robust
outlier rejection in the other. Each subfolder is self-contained, with
its own README, inputs, outputs, and a `legacy/` folder documenting the
earlier iterations that led to the current script.

## License

MIT — see [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff).
