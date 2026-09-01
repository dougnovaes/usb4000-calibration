%% WAVELENGTH_CALIBRATION_Calibrated
%  Wavelength-axis calibration of an Ocean Optics USB4000 spectrometer from
%  Hg and Ne emission lines — fit, robust anchor, bandwidth-weighted least
%  squares, and validation of the WRITTEN coefficients against a fresh,
%  sub-pixel-accurate acquisition.
%
%  ------------------------------------------------------------------------
%  WHAT CHANGED IN THIS VERSION
%  ------------------------------------------------------------------------
%  Sections 1-10 (the fit itself) are unchanged from the previous version
%  and reproduce the same numbers — this file's news is entirely in
%  Section 11. Three validation rounds against peak tables (not raw
%  spectra) all showed a ~0.5 nm residual that did not shrink with an
%  air-conditioned room or a longer settling time, which ruled out simple
%  thermal drift. Two things were found once raw validation spectra
%  finally became available, and both are now built into Section 11
%  permanently rather than living in a one-off analysis (NOTE 13, NOTE 14
%  below):
%
%    1. The "Wavelength" column SpectraSuite reports for a peak is the
%       calibration evaluated at an INTEGER pixel — not a sub-pixel value.
%       Comparing it directly against the catalogue compares a ~0.2 nm
%       quantization step against a fit good to ~0.03 nm. Section 11 now
%       re-detects peaks from the raw validation spectrum with the exact
%       same sub-pixel routine used for fitting (NOTE 3), instead of
%       trusting a pre-computed peak table.
%
%    2. Most of the remaining ~0.5 nm was not drift at all: fitting a
%       cubic directly to a validation spectrum file's OWN wavelength
%       column (which is fab_active(pixel), whatever is actually loaded
%       on the instrument) recovers the coefficients really in effect,
%       independent of any assumption about what should be there. Doing
%       this for both the Hg and the Ne validation spectra gave IDENTICAL
%       coefficients to 8+ significant figures — strong evidence they
%       reflect the instrument, not noise — and a0, a1, a2 matched the
%       adopted fit closely while a3 matched the OLD FACTORY value instead:
%       the "3rd Coefficient" field had not actually been overwritten.
%       Section 11 now performs this extraction automatically (NOTE 14)
%       and flags any coefficient that looks like it was not written,
%       before it ever gets a chance to look like drift.
%
%    3. NOTE 15 (new): a validation "spectrum" whose own wavelength column
%       is not a smooth function of pixel — a mislabeled file, a table
%       exported under some other mode, or a file emptied/corrupted in
%       transit — produces an axis self-fit residual of tens of nm or
%       worse (or NaN outright), and everything downstream of that
%       (active-coefficient values, the Hg/Ne cross-check, the
%       percent-off flags) becomes meaningless noise rather than a
%       diagnosis. Section 11 now checks the axis self-fit residual
%       BEFORE reporting any coefficient comparison for that file, and
%       skips the comparison (reporting NaN) rather than printing
%       thousand-percent "differences" that describe the file, not the
%       instrument.
%
%  ------------------------------------------------------------------------
%  CONTEXT
%  ------------------------------------------------------------------------
%  Stage A of an absolute-calibration procedure for a USB4000, used as
%  training before the same work on Jobin-Yvon THR1000 and McPherson 207
%  monochromators. End application: spectral characterisation of a cold
%  dielectric-barrier-discharge (DBD) plasma. (The FEL lamp's own,
%  separate radiometric calibration — six interpolation models, LOOCV,
%  lamp-file generation — lives in ABSOLUTE_CALIBRATION.m and is
%  unrelated to this file beyond sharing a project.)
%
%  The detector does not measure wavelength — it only records which of its
%  3648 pixels received light, and how much. The pixel-to-wavelength
%  relation is a property of that particular optical assembly and must be
%  measured by pointing the instrument at sources of known wavelength.
%  Because the grating equation relates wavelength to diffraction ANGLE
%  rather than directly to position on a flat detector, the relation
%  carries residual curvature — hence a cubic polynomial.
%
%  ------------------------------------------------------------------------
%  INPUTS
%  ------------------------------------------------------------------------
%  hg_spectrum.txt, ne_spectrum.txt   (required, for the fit)
%      Two tab-separated columns exported from SpectraSuite under the
%      FACTORY calibration: wavelength [nm], counts. One row per pixel.
%
%  hg_bandwidths.txt, ne_bandwidths.txt   (optional, strongly recommended)
%      pixel, 90% bandwidth [nm], from the Peak Info panel. NaN for
%      unstable readings. Falls back to uniform weights if absent (NOTE 9).
%
%  validation_hg_spectrum.txt, validation_ne_spectrum.txt   (optional)
%      SAME two-column format as the fitting spectra, but acquired AFTER
%      writing the fitted coefficients, under whatever calibration is
%      CURRENTLY active on the instrument (not necessarily the adopted
%      one — that is exactly what Section 11 checks). This is the
%      preferred validation input; when present it drives both NOTE 13
%      and NOTE 14 automatically. Must be a REAL post-write acquisition —
%      see NOTE 15 for what happens when it is not.
%
%  validation_hg.txt, validation_ne.txt   (optional, legacy fallback)
%      Pre-extracted peak tables (pixel, wavelength, bandwidth, centre
%      wavelength), used only if the raw validation spectra above are not
%      found. Section 11 prints an explicit resolution warning in this
%      mode — see NOTE 13.
%
%  ------------------------------------------------------------------------
%  OUTPUTS
%  ------------------------------------------------------------------------
%  calibration_coefficients.txt, calibration_residuals.csv  — the fit
%  validation_residuals.csv                                 — Section 11
%  a diagnostic figure (6 panels, +1 more when validation data is present)
%
%  ------------------------------------------------------------------------
%  REFERENCE RESULT (for verification)
%  ------------------------------------------------------------------------
%    Bandwidth-weighted fit, robust Hg anchor, 32 lines (254-744 nm):
%        a0 =  1.7635567179E+02
%        a1 =  2.1606649029E-01
%        a2 = -3.6985346136E-06
%        a3 = -5.5381923002E-10
%        weighted RMS ................ 0.0339 nm
%        factory calibration .......... 1.0419 nm   (~31x improvement)
%        Ne pattern hit rate .......... 25 of 27 peaks (93%)
%        anchor: 253.6521 nm excluded (leverage 0.94); anchor offset
%                scatter 0.183 -> 0.066 nm
%
%    Section 11, raw validation spectra, coefficients as WRITTEN before
%    the a3 issue below was found:
%        active-coefficient extraction: a0,a1,a2 matched the adopted fit;
%        a3 matched FACTORY (-5.9958e-10) instead of adopted
%        (-5.5382e-10) — the field had not been overwritten.
%        sub-pixel residual using the file's own (uncorrected) axis: 0.38 nm
%        sub-pixel residual after evaluating at the ADOPTED coefficients: 0.05-0.08 nm
%        (the two sub-pixel methods compared during development — parabola
%        in wavelength-space vs. in pixel-index-space, see NOTE 3 — agree
%        to within this range, which is itself the practical floor of this
%        validation: below about 0.1 nm, the choice of sub-pixel method
%        matters as much as anything else being measured.)
%
%    STATUS AT TIME OF WRITING: the adopted coefficients above (all four)
%    have been re-entered into the instrument and confirmed directly on
%    the instrument's own Wavelength Calibration panel (all four terms
%    match to the panel's displayed precision). Re-run Section 11 with a
%    fresh, REAL pair of post-write validation spectra to confirm at
%    sub-pixel resolution; if the residual is within the 0.05-0.1 nm
%    range quoted above, the calibration is closed.
%
%  ------------------------------------------------------------------------
%  METHODOLOGICAL NOTES
%  ------------------------------------------------------------------------
%  NOTE 1 (vacuum vs. air), NOTE 2 (query the full NIST ASD, not the
%  curated table), NOTE 4 (SpectraSuite's peak threshold is global), NOTE
%  5 (back up factory coefficients), NOTE 6 (saturation), NOTE 8 (blended
%  lines need an effective wavelength, not the dominant component's
%  catalogue value), and NOTE 10 (identify Hg first, extrapolate to Ne as
%  a whole pattern, never peak-by-peak) are unchanged from the previous
%  version and are not repeated in full here; see the inline code
%  comments at each step for the specifics.
%
%  NOTE 3 — Sub-pixel position: parabola on the core, not a Gaussian
%    The peak pixel alone limits precision to about ±0.1 nm. A parabola
%    through the three central points, fitted in PIXEL-INDEX space (the
%    physical detector coordinate — this is what Sections 3-4 use, and
%    what Section 11 now also uses, for consistency), gives ~1/10 pixel.
%    It is not pulled by the asymmetric red-side coma tail the way a
%    whole-profile Gaussian would be. During development, fitting the
%    same parabola directly in WAVELENGTH space (using a possibly-wrong
%    active calibration as the x-axis) gave a validation RMS of 0.046 nm
%    against 0.083 nm for the index-space method on the same data — a
%    reminder that once real drift and coefficient errors are removed,
%    the remaining ~0.05-0.1 nm is sub-pixel METHOD noise, not signal.
%
%  NOTE 7 — Prominence filter (indispensable)
%    A height threshold alone lets noise ripples on a strong line's flank
%    pass as peaks. PROMINENCE (height above the higher of the two
%    flanking saddles) separates them: a ripple showed 1.8% of its own
%    height; a real peak, close to 100%. Threshold: 30%.
%
%  NOTE 9 — Why the 90% bandwidth is READ, not recomputed
%    A version that computed peak width from the raw counts correlated
%    +0.12 with what SpectraSuite reports and was useless as a blend
%    detector — it only sees the peak core, exactly where a blend does
%    not show. The instrument-reported value is read from a file; if
%    absent, weights fall back to uniform rather than a misleading
%    substitute.
%
%  NOTE 11 — The anchor needs a ROBUST subset, built from LEAVE-ONE-OUT
%            residuals, not raw ones
%    Stage 1 (Hg) and stage 2 (Ne) have opposite requirements: the anchor
%    needs reliable POSITION (one bad line poisons the extrapolation),
%    the final fit wants COVERAGE (bandwidth weighting already handles
%    poor lines). Raw-residual screening fails on exactly the point that
%    matters: the 253.65 nm line has leverage 0.94, so the fit bends to
%    pass through it and its raw residual is the third-smallest of eight
%    — while an innocent line gets blamed instead. The leave-one-out
%    residual r_LOO,i = r_i/(1-h_i) undoes this: 253.65 nm jumps to 10.7x
%    the median absolute deviation and is correctly isolated. Shared by
%    the anchor here (Section 5) and, previously, by the validation drift
%    model — now mostly superseded by NOTE 14's more direct diagnosis,
%    but kept as a fallback (Section 11c) for when the residual isn't
%    fully explained by a coefficient mismatch.
%
%  NOTE 12 — Reading a validation residual: drift vs. a wiring error vs. a
%            coefficient that silently stayed at its old value
%    All three produce a "the written calibration disagrees with the
%    catalogue" symptom, and they are easy to conflate. A wiring/entry
%    error (wrong coefficient typed in) tends to reproduce almost exactly
%    at the low-pixel end (where all cubics agree, having similar a0) and
%    diverge increasingly toward the other end — which is also what mild
%    thermal drift of the LEADING coefficients would look like. The
%    decisive discriminator turned out to be NOTE 14, not the shape of the
%    residual: extracting the coefficients actually driving the file's own
%    axis and comparing them one at a time to what was supposed to be
%    written settles the question directly, in a way that no amount of
%    residual-shape analysis can.
%
%  NOTE 13 — Validate at the SAME sub-pixel precision as the fit, or the
%            comparison is not fair
%    Comparing the fit (sub-pixel, ~0.03-0.05 nm) against a validation
%    read at integer-pixel resolution (~0.2 nm quantization) makes any
%    fit look roughly ten times worse than it is, for no physical reason.
%    Confirmed directly: repeating the SAME acquisition's peak table
%    twice gave IDENTICAL "Wavelength" values at 26-29 of ~31-33 lines
%    across three separate rounds, and the handful that did change moved
%    by almost exactly one pixel of local dispersion (0.99x, measured) —
%    the unmistakable signature of a value that only updates when the
%    true peak crosses an integer-pixel boundary, not a continuously
%    varying measurement. Section 11 (a) below re-detects peaks from the
%    raw validation spectrum with the fit's own sub-pixel routine
%    whenever that spectrum is available, and only falls back to the
%    coarser peak-table comparison — with an explicit warning — when it
%    is not.
%
%  NOTE 14 — Extract what's actually active before diagnosing anything else
%    A spectrum file's wavelength column is fab_active(pixel) for
%    whatever fab_active is currently loaded on the instrument — it does
%    not have to be the fit you think you wrote. Fitting a cubic to that
%    column against the row index recovers fab_active's coefficients
%    directly, with no assumption about what SHOULD be there. Applied
%    here: a0, a1, a2 matched the adopted fit to 4+ significant figures;
%    a3 matched the OLD FACTORY value (-5.9958e-10) far more closely than
%    the adopted one (-5.5382e-10) — an 8% relative gap, reproduced
%    identically from the independent Hg and Ne validation spectra. The
%    "3rd Coefficient" field on the instrument had silently kept its old
%    value. This check costs one polynomial fit and should be the FIRST
%    thing Section 11 does with a raw validation spectrum, before any
%    line-by-line residual or drift analysis — it answers "is this even
%    the calibration I think it is?" before asking "how good is it?".
%
%  NOTE 15 — A validation file that fails its own axis self-fit is not a
%            coefficient mismatch, it's not real data
%    NOTE 14's extraction only means something if the file's wavelength
%    column IS a cubic-in-pixel to begin with (i.e. it really did come
%    off this spectrometer's own axis). A genuine post-write acquisition
%    self-fits to well under 0.1 nm, matching the in-sample fit RMS. A
%    file that does not — self-fit residual of tens of nm, or NaN because
%    the column is degenerate — is not encoding any fab_active worth
%    comparing to anything; running report_coef_check on it manufactures
%    "percent off ADOPTED" numbers in the thousands or billions that
%    describe the file's unrelatedness to a real spectrum, not an
%    instrument problem. MAX_AXIS_SELFFIT_RESID (Section 1) gates this:
%    above threshold, the coefficient comparison for that file is skipped
%    and a_active is reported as NaN instead of being printed as if it
%    meant something.
%
%  ------------------------------------------------------------------------
%  REFERENCES
%  ------------------------------------------------------------------------
%    Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024).
%    NIST Atomic Spectra Database (ver. 5.12). physics.nist.gov/asd
%    DOI: 10.18434/T4W30F
%    Ocean Optics. SpectraSuite Installation and Operation Manual.
%
%  Author: Douglas Oliveira Novaes                             Licence: MIT
%  ========================================================================

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURATION
%  ========================================================================
HG_SPECTRUM  = 'hg_spectrum.txt';
NE_SPECTRUM  = 'ne_spectrum.txt';
HG_BANDWIDTH = 'hg_bandwidths.txt';
NE_BANDWIDTH = 'ne_bandwidths.txt';

HG_VAL_SPECTRUM = 'validation_hg_spectrum.txt';   % preferred (NOTE 13/14)
NE_VAL_SPECTRUM = 'validation_ne_spectrum.txt';
HG_VALIDATION   = 'validation_hg.txt';            % legacy fallback
NE_VALIDATION   = 'validation_ne.txt';

DEGREE       = 3;      N_PARAB    = 3;      FRAC_THRESH = 0.02;
PROM_MIN     = 0.30;   SEP_MIN    = 5;      SHIFT_MAX   = 1.0;
TOL_HG       = 2.0;    TOL_NE     = 0.30;   TOL_VALID   = 1.50;
RATE_MIN     = 0.60;   CLEAN_RATIO= 1.3;    PIX_TOL     = 3;
MAD_REJECT   = 4;      DRIFT_DEGREE = 2;    N_PIX       = 3648;
COEF_REL_TOL = 0.02;   % NOTE 14: flag a coefficient as "maybe not written"
                        % if it differs from the adopted fit by more than
                        % this fraction of the adopted value's magnitude.
                        % The real a3 bug measured 8.3% -- 2% catches it
                        % with ample margin while ignoring ordinary rounding.
MAX_AXIS_SELFFIT_RESID = 1.0;  % nm -- NOTE 15: a validation file's own
                        % wavelength column must self-fit a cubic to
                        % below this to be trusted as a real post-write
                        % acquisition. Sections 1-10 operate at
                        % 0.03-0.05 nm RMS, so 1.0 nm is already a very
                        % generous margin above any real instrument axis
                        % and comfortably rejects a mislabeled/corrupt file.

FACTORY_COEF = [178.085770, 2.1517764E-1, -3.3423826E-6, -5.9958050E-10];

%% ========================================================================
%  2. REFERENCE LINE TABLES (NIST ASD, air wavelengths)
%  ========================================================================
HG_LINES = [ ...
    253.6521, 253.6521; ...
    313.1716, 313.1555; ...
    365.1205, 365.0158; ...
    404.6565, 404.6565; ...
    435.8335, 435.8335; ...
    546.0750, 546.0750; ...
    576.9610, 576.9610; ...
    579.0670, 579.0670];

NE_LINES = [585.2488; 588.1895; 594.4834; 597.5534; 603.0000; 607.4338; ...
            609.6163; 614.3062; 616.3594; 621.7281; 626.6495; 630.4789; ...
            633.4428; 638.2991; 640.2248; 650.6528; 653.2882; 659.8953; ...
            667.8276; 671.7043; 692.9467; 703.2413; 717.3938; 724.5167; ...
            743.8899];
ALL_LINES = [HG_LINES(:,1); NE_LINES];

%% ========================================================================
%  3. READ SPECTRA
%  ========================================================================
[pix_hg, cnt_hg] = read_spectrum(HG_SPECTRUM, FACTORY_COEF);
[pix_ne, cnt_ne] = read_spectrum(NE_SPECTRUM, FACTORY_COEF);
fprintf('Spectra loaded: Hg (%d points), Ne (%d points)\n', numel(cnt_hg), numel(cnt_ne));

%% ========================================================================
%  4. PEAK DETECTION AND SUB-PIXEL POSITION
%  ========================================================================
[p_hg, rej_hg] = find_peaks(pix_hg, cnt_hg, FRAC_THRESH, PROM_MIN, SEP_MIN, N_PARAB, SHIFT_MAX);
[p_ne, rej_ne] = find_peaks(pix_ne, cnt_ne, FRAC_THRESH, PROM_MIN, SEP_MIN, N_PARAB, SHIFT_MAX);
if any(rej_hg), fprintf('NOTE: %d Hg peak(s) rejected.\n', sum(rej_hg)); end
if any(rej_ne), fprintf('NOTE: %d Ne peak(s) rejected.\n', sum(rej_ne)); end
p_hg = p_hg(~rej_hg);  p_ne = p_ne(~rej_ne);
fprintf('Usable peaks: Hg = %d, Ne = %d\n', numel(p_hg), numel(p_ne));

%% ========================================================================
%  5. TWO-STAGE IDENTIFICATION (NOTE 10) WITH A ROBUST ANCHOR (NOTE 11)
%  ========================================================================
lam_fac = polyval(fliplr(FACTORY_COEF), p_hg);
[id_hg, ok_hg] = match(lam_fac, HG_LINES(:,1), TOL_HG);
p_hg = p_hg(ok_hg);  lam_hg = id_hg(ok_hg);
fprintf('\nHg lines identified: %d\n', numel(p_hg));
assert(numel(p_hg) >= DEGREE+2, 'Too few Hg lines (%d) to anchor the fit.', numel(p_hg));

offset = polyval(fliplr(FACTORY_COEF), p_hg) - lam_hg;
fprintf('Hg factory-NIST offset: %.3f to %.3f nm (scatter %.3f nm)\n', min(offset), max(offset), std(offset));

anchor = robust_subset(p_hg, lam_hg, DEGREE, MAD_REJECT);
if ~all(anchor)
    fprintf('Excluded from the ANCHOR (still used in the final fit):');
    fprintf(' %.4f', lam_hg(~anchor)); fprintf(' nm\n');
    fprintf('Anchor offset scatter: %.3f nm (was %.3f with all lines)\n', std(offset(anchor)), std(offset));
end

[c_hg, ~, mu_hg] = polyfit(p_hg(anchor), lam_hg(anchor), DEGREE);
lam_extrap = polyval(c_hg, (p_ne - mu_hg(1))/mu_hg(2));
[id_ne, ok_ne] = match(lam_extrap, NE_LINES, TOL_NE);
p_ne = p_ne(ok_ne);  lam_ne = id_ne(ok_ne);
rate = numel(p_ne)/numel(ok_ne);
fprintf('Ne lines identified: %d of %d peaks (%.0f%%, NO free parameter)\n', numel(p_ne), numel(ok_ne), 100*rate);
if rate < RATE_MIN
    warning(['Only %.0f%% of the Ne peaks landed on catalogued lines. Likely ' ...
        'causes: (a) a false peak (NOTE 7); (b) an insufficiently robust ' ...
        'anchor (NOTE 11); (c) wrong lamp species (NOTE 10).'], 100*rate);
end

%% ========================================================================
%  6. BANDWIDTHS, CLASSIFICATION AND WEIGHTS (NOTE 9)
%  ========================================================================
pixel  = [p_hg(:);  p_ne(:)];
target = [lam_hg(:); lam_ne(:)];
is_Hg  = [true(numel(p_hg),1); false(numel(p_ne),1)];
[pixel, ord] = sort(pixel);  target = target(ord);  is_Hg = is_Hg(ord);

bw = [read_bandwidths(HG_BANDWIDTH, p_hg, PIX_TOL); read_bandwidths(NE_BANDWIDTH, p_ne, PIX_TOL)];
bw = bw(ord);

have_bw = any(isfinite(bw));
if have_bw
    keep_bw = isfinite(bw);
    bw_min  = min(bw(keep_bw));
    ratio   = bw / bw_min;
    weight  = (bw_min ./ bw).^2;
    weight(~keep_bw) = 0;
    if any(~keep_bw)
        fprintf('NOTE: %d line(s) have no usable bandwidth and were given zero weight.\n', sum(~keep_bw));
    end
else
    warning('No bandwidth file found. Falling back to UNIFORM weights (NOTE 9).');
    keep_bw = true(size(bw));  bw_min = NaN;  ratio = nan(size(bw));  weight = ones(size(bw));
end

pixel = pixel(keep_bw);  target = target(keep_bw);  is_Hg = is_Hg(keep_bw);
bw = bw(keep_bw);  ratio = ratio(keep_bw);  weight = weight(keep_bw);
n = numel(pixel);
clean = have_bw & (ratio <= CLEAN_RATIO);

if have_bw
    fprintf('\n=== 90%% BANDWIDTH CLASSIFICATION ===\n');
    fprintf('Smallest bandwidth = %.3f nm  ->  instrumental FWHM ~ %.2f nm\n', bw_min, bw_min/0.39);
    fprintf('%10s %9s %7s %8s  %s\n','lambda','bandwidth','ratio','weight','class');
    for i = 1:n
        fprintf('%10.4f %9.2f %6.1fx %8.3f  %s\n', target(i), bw(i), ratio(i), weight(i), ...
                bandwidth_class(ratio(i), CLEAN_RATIO));
    end
    fprintf('Clean lines: %d of %d\n', sum(clean), n);
end

%% ========================================================================
%  7. FIT — THREE SCHEMES, FOR COMPARISON
%  ========================================================================
[a_wt, r_wt, rms_wt] = fit_poly(pixel, target, weight,      DEGREE);
[a_un, r_un, rms_un] = fit_poly(pixel, target, ones(n,1),   DEGREE);
if sum(clean) >= DEGREE+2
    [a_cl, r_cl, rms_cl] = fit_poly(pixel(clean), target(clean), ones(sum(clean),1), DEGREE);
else
    a_cl = nan(1,DEGREE+1); rms_cl = NaN;
end
rms_fac = rms(polyval(fliplr(FACTORY_COEF), pixel) - target);

fprintf('\n=== COMPARISON OF FITTING SCHEMES ===\n');
fprintf('%-28s %4s %11s %12s\n','scheme','n','RMS (nm)','range (nm)');
fprintf('%-28s %4d %11.4f %5.0f-%.0f\n','bandwidth-weighted (NOTE 9)', n, rms_wt, min(target), max(target));
fprintf('%-28s %4d %11.4f %5.0f-%.0f\n','uniform', n, rms_un, min(target), max(target));
if sum(clean) >= DEGREE+2
    fprintf('%-28s %4d %11.4f %5.0f-%.0f\n','clean lines only', sum(clean), rms_cl, min(target(clean)), max(target(clean)));
    fprintf('\nWeighted fit evaluated on clean lines only: %.4f nm\n', rms(r_wt(clean)));
end
fprintf('%-28s %4d %11.4f\n','factory calibration', n, rms_fac);

a = a_wt;  residual = r_wt;    % ADOPTED coefficients -- used throughout Section 11
fprintf('\n=== ADOPTED COEFFICIENTS (weighted) ===\n');
fprintf('Intercept        (a0) = %.6f\n',  a(1));
fprintf('1st Coefficient  (a1) = %.10E\n', a(2));
fprintf('2nd Coefficient  (a2) = %.10E\n', a(3));
fprintf('3rd Coefficient  (a3) = %.10E\n', a(4));

%% ========================================================================
%  8. CROSS-VALIDATION AND LEVERAGE
%  ========================================================================
loo = zeros(n,1);
for i = 1:n
    m = true(n,1); m(i) = false;
    ai = fit_poly(pixel(m), target(m), weight(m), DEGREE);
    loo(i) = polyval(fliplr(ai), pixel(i)) - target(i);
end
mu = mean(pixel); sd = std(pixel);
X  = vandermonde((pixel-mu)/sd, DEGREE);
% Alavancagem do ajuste PONDERADO (hat matrix de WLS), não a de OLS: H deve
% refletir o mesmo ajuste que 'a'/'residual' representam (fit_poly com
% 'weight'), senão o diagnóstico de alavancagem da NOTE 11 descreve um
% ajuste diferente do que foi de fato adotado. O traço continua = DEGREE+1,
% então "expected mean = (DEGREE+1)/n" abaixo permanece válido sem alteração.
Wd = diag(weight(:));
H  = diag(X*pinv(X'*Wd*X)*X'*Wd);

fprintf('\n=== CROSS-VALIDATION ===\n');
fprintf('LOOCV RMS = %.4f nm (in-sample: %.4f nm), ratio = %.2f\n', rms(loo), rms(residual), rms(loo)/rms(residual));
fprintf('Leverage: max = %.3f (lambda = %.1f nm), expected mean = %.3f\n', max(H), target(H==max(H)), (DEGREE+1)/n);
if rms(loo)/rms(residual) > 1.5 && max(H) > 3*(DEGREE+1)/n
    fprintf(['  -> high ratio EXPLAINED BY LEVERAGE, not overfitting (NOTE 11). ' ...
             'Add reference lines in the sparsely covered region to reduce it.\n']);
end

%% ========================================================================
%  9. EXPORT (fit results)
%  ========================================================================
source = repmat("Ne", n, 1); source(is_Hg) = "Hg";
cls = strings(n,1);
for i = 1:n, cls(i) = string(bandwidth_class(ratio(i), CLEAN_RATIO)); end
in_anchor = false(n,1);
for i = 1:n
    if is_Hg(i)
        j = find(abs(p_hg-pixel(i))<1e-6,1);
        if ~isempty(j), in_anchor(i) = anchor(j); end
    end
end
T = table(pixel, target, polyval(fliplr(a),pixel), residual, bw, ratio, weight, source, cls, in_anchor, ...
    'VariableNames', {'pixel','lambda_NIST_nm','predicted_nm','residual_nm', ...
                      'bandwidth90_nm','bandwidth_ratio','weight','source','class','in_anchor'});
writetable(T,'calibration_residuals.csv');

fid = fopen('calibration_coefficients.txt','w');
fprintf(fid,'# Wavelength calibration - Ocean Optics USB4000\n# Generated %s\n#\n', datestr(now,'yyyy-mm-dd HH:MM'));
fprintf(fid,'[WEIGHTED]    n=%d  RMS=%.4f nm  range=%.1f-%.1f nm\n', n, rms_wt, min(target), max(target));
fprintf(fid,'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n\n', a_wt);
fprintf(fid,'[FACTORY_BACKUP]  RMS=%.4f nm\n', rms_fac);
fprintf(fid,'Intercept\t%.6f\n1st\t%.7E\n2nd\t%.7E\n3rd\t%.7E\n', FACTORY_COEF);
fclose(fid);
fprintf('\nFit results written: calibration_residuals.csv, calibration_coefficients.txt\n');

%% ========================================================================
%  10. FIGURE (fit diagnostics) — panels 1-6; 7 added in Section 11
%  ========================================================================
have_raw_val   = isfile(HG_VAL_SPECTRUM) || isfile(NE_VAL_SPECTRUM);
have_table_val = isfile(HG_VALIDATION)   || isfile(NE_VALIDATION);
have_validation = have_raw_val || have_table_val;
n_rows = 2 + have_validation;
figure('Units','normalized','OuterPosition',[0.02 0.05 0.96 0.90]);
cHg=[0 0.45 0.74]; cNe=[0.85 0.33 0.10]; grey=[0.6 0.6 0.6];

% pp_full cobre todo o detector -- usado só para a curva de FÁBRICA, que é
% um polinômio conhecido, não um ajuste extrapolado. pp_data cobre apenas a
% faixa de pixel com linhas de referência (254-744 nm); pp_lo/pp_hi cobrem
% o restante do detector, onde o polinômio ADOTADO é extrapolação pura --
% cúbicas são notoriamente instáveis fora do intervalo ajustado, então essa
% faixa é traçada pontilhada em vez de contínua, para não sugerir cobertura
% validada onde não há nenhum marcador Hg/Ne.
pp_full = linspace(0, N_PIX, 500);
pp_data = linspace(min(pixel), max(pixel), 400);
pp_lo   = linspace(0, min(pixel), 60);
pp_hi   = linspace(max(pixel), N_PIX, 60);

subplot(n_rows,3,1); hold on; grid on; box on;
h_fit = plot(pp_data, polyval(fliplr(a),pp_data), 'k-', 'LineWidth',1.4);
plot(pp_lo, polyval(fliplr(a),pp_lo), 'k:', 'LineWidth',1.0, 'HandleVisibility','off');
plot(pp_hi, polyval(fliplr(a),pp_hi), 'k:', 'LineWidth',1.0, 'HandleVisibility','off');
h_fac = plot(pp_full, polyval(fliplr(FACTORY_COEF),pp_full), '--','Color',grey,'LineWidth',1.1);
h_hg  = plot(pixel(is_Hg), target(is_Hg), 'o','MarkerFaceColor',cHg,'MarkerEdgeColor','k','MarkerSize',6);
h_ne  = plot(pixel(~is_Hg),target(~is_Hg),'s','MarkerFaceColor',cNe,'MarkerEdgeColor','k','MarkerSize',6);
xlabel('Pixel'); ylabel('Wavelength (nm)');
title('Calibration curve \lambda(pixel)  (dotted = extrapolated beyond reference lines)');
lg1 = legend([h_fit,h_fac,h_hg,h_ne], {'fit','fact','Hg','Ne'}, 'Location','southeast');
lg1.ItemTokenSize = [15, 18];
lg1.Position(3) = lg1.Position(3) * 1.0;
lg = findobj(gcf, 'Type', 'Legend');
disp(lg(end).String)

subplot(n_rows,3,2); hold on; grid on; box on;
yline(0,'k-');
sz = 20 + 120*weight/max(weight);
scatter(pixel(is_Hg),  residual(is_Hg),  sz(is_Hg),  cHg,'filled','MarkerEdgeColor','k');
scatter(pixel(~is_Hg), residual(~is_Hg), sz(~is_Hg), cNe,'filled','MarkerEdgeColor','k');
yline( rms_wt,':k'); yline(-rms_wt,':k');
xlabel('Pixel'); ylabel('Residual (nm)'); title(sprintf('Fit residuals (weighted RMS = %.4f nm)', rms_wt));

subplot(n_rows,3,3); hold on; grid on; box on;
if have_bw
    plot(pixel, ratio, 'ko-','MarkerFaceColor','k','MarkerSize',4);
    yline(CLEAN_RATIO,'--','clean','Color',[0 .6 0]); yline(2.5,'--','blend','Color',[.8 .5 0]);
    yline(6.0,'--','severe','Color',[.8 0 0]);
    set(gca,'YScale','log'); xlabel('Pixel'); ylabel('bandwidth / min');
    title('Blend diagnostic from the 90% bandwidth');
else
    text(0.5,0.5,{'no bandwidth file','(uniform weights)'},'HorizontalAlignment','center','Units','normalized'); axis off;
end

subplot(n_rows,3,4); hold on; grid on; box on;
plot(pix_hg, max(cnt_hg,1), '-','Color',cHg,'LineWidth',0.6);
for i = find(is_Hg)', xline(pixel(i),':','Color','k'); end
set(gca,'YScale','log'); xlabel('Pixel'); ylabel('Counts'); title('Hg spectrum and lines used');

subplot(n_rows,3,5); hold on; grid on; box on;
plot(pix_ne, max(cnt_ne,1), '-','Color',cNe,'LineWidth',0.6);
for i = find(~is_Hg)', xline(pixel(i),':','Color','k'); end
set(gca,'YScale','log'); xlabel('Pixel'); ylabel('Counts'); title('Ne spectrum and lines used');

subplot(n_rows,3,6); hold on; grid on; box on;
plot(pp_data, polyval(polyder(fliplr(a)), pp_data), 'k-','LineWidth',1.4);
plot(pp_lo, polyval(polyder(fliplr(a)), pp_lo), 'k:','LineWidth',1.0);
plot(pp_hi, polyval(polyder(fliplr(a)), pp_hi), 'k:','LineWidth',1.0);
xlabel('Pixel'); ylabel('d\lambda/dpixel (nm/pixel)');
title('Dispersion (derivative of the polynomial; dotted = extrapolated)');

% lg = findobj(gcf, 'Type', 'Legend');
% lg(end).String   % a legenda do subplot 1 deve imprimir {'fit','factory','Hg','Ne'}

%% ========================================================================
%  11. VALIDATION AGAINST A POST-WRITE ACQUISITION
%  ========================================================================
if ~have_validation
    fprintf(['\nNo validation files found (%s, %s, or the legacy peak-table ' ...
        'pair). Skipping Section 11.\nRe-acquire the Hg and Ne lamps after ' ...
        'writing the coefficients to validate them.\n'], HG_VAL_SPECTRUM, NE_VAL_SPECTRUM);

elseif have_raw_val
    % ---- preferred path: raw spectra, sub-pixel + active-coefficient check
    fprintf('\n\n=== SECTION 11: VALIDATION (raw spectra, sub-pixel, NOTE 13/14/15) ===\n');
    % val_desc documenta O QUE o RMS abaixo está de fato respondendo -- aqui,
    % "o ajuste é bom?", pois meas_v é avaliado nos coeficientes ADOTADOS,
    % não nos que porventura estejam gravados no instrumento agora.
    val_desc = 'ADOPTED coefficients (fit quality, independent of what is currently written)';

    pixel_v = []; meas_idx_v = []; is_Hg_v = false(0,1); cnt_v_all = {}; pixv_all = {}; labels={};
    active_coefs = [];

    if isfile(HG_VAL_SPECTRUM)
        [pv, cv] = read_spectrum_novalidate(HG_VAL_SPECTRUM);
        [a_act_hg, resid_axis_hg] = extract_active_coefficients(HG_VAL_SPECTRUM, DEGREE);
        fprintf('Active coefficients from %s (axis self-fit residual %.4f nm):\n', HG_VAL_SPECTRUM, resid_axis_hg);
        % NOTE 15: only trust this file's coefficient comparison if its own
        % wavelength column actually self-fits a cubic well -- otherwise
        % the "active coefficients" and every %-off-ADOPTED flag below
        % describe the file's unrelatedness to a real spectrum, not the
        % instrument.
        if ~isfinite(resid_axis_hg) || resid_axis_hg > MAX_AXIS_SELFFIT_RESID
            warning(['%s: axis self-fit residual (%.2f nm) far exceeds a real ' ...
                'post-write acquisition (expect < %.1f nm) -- this file is ' ...
                'probably not a genuine spectrum export. Coefficient comparison ' ...
                'skipped.'], HG_VAL_SPECTRUM, resid_axis_hg, MAX_AXIS_SELFFIT_RESID);
            a_act_hg = nan(1,4);
        else
            report_coef_check(a_act_hg, a, COEF_REL_TOL);
        end
        active_coefs = [active_coefs; a_act_hg];
        [pk, rj] = find_peaks(pv, cv, FRAC_THRESH, PROM_MIN, SEP_MIN, N_PARAB, SHIFT_MAX);
        pk = pk(~rj);
        pixel_v = [pixel_v; pk(:)]; is_Hg_v = [is_Hg_v; true(numel(pk),1)]; %#ok<AGROW>
    end
    if isfile(NE_VAL_SPECTRUM)
        [pv, cv] = read_spectrum_novalidate(NE_VAL_SPECTRUM);
        [a_act_ne, resid_axis_ne] = extract_active_coefficients(NE_VAL_SPECTRUM, DEGREE);
        fprintf('Active coefficients from %s (axis self-fit residual %.4f nm):\n', NE_VAL_SPECTRUM, resid_axis_ne);
        if ~isfinite(resid_axis_ne) || resid_axis_ne > MAX_AXIS_SELFFIT_RESID
            warning(['%s: axis self-fit residual (%.2f nm) far exceeds a real ' ...
                'post-write acquisition (expect < %.1f nm) -- this file is ' ...
                'probably not a genuine spectrum export. Coefficient comparison ' ...
                'skipped.'], NE_VAL_SPECTRUM, resid_axis_ne, MAX_AXIS_SELFFIT_RESID);
            a_act_ne = nan(1,4);
        else
            report_coef_check(a_act_ne, a, COEF_REL_TOL);
        end
        active_coefs = [active_coefs; a_act_ne];
        [pk, rj] = find_peaks(pv, cv, FRAC_THRESH, PROM_MIN, SEP_MIN, N_PARAB, SHIFT_MAX);
        pk = pk(~rj);
        pixel_v = [pixel_v; pk(:)]; is_Hg_v = [is_Hg_v; false(numel(pk),1)]; %#ok<AGROW>
    end

    if size(active_coefs,1) == 2
        if any(isnan(active_coefs(:)))
            fprintf(['Cross-check, Hg vs. Ne active coefficients: skipped -- at least ' ...
                'one file failed the axis self-fit check above.\n']);
        else
            d = max(abs(active_coefs(1,:) - active_coefs(2,:)) ./ max(abs(active_coefs(2,:)),eps));
            fprintf('Cross-check, Hg vs. Ne active coefficients: max relative difference %.2e', d);
            if d < 1e-4, fprintf('  (agree -- both read off the same instrument state)\n');
            else, fprintf('  (DISAGREE -- check the two files were taken back to back)\n'); end
        end
    end

    % Sub-pixel positions were found in PIXEL-INDEX space (same as the fit,
    % NOTE 3); evaluate them at the ADOPTED coefficients regardless of what
    % is actually active on the instrument -- this is what makes the
    % comparison fair even when NOTE 14 above found a mismatch.
    meas_v = polyval(fliplr(a), pixel_v);
    [id_v, ok_v] = match(meas_v, ALL_LINES, TOL_VALID);
    n_unmatched = sum(~ok_v);
    if n_unmatched > 0
        fprintf('%d point(s) matched no catalogued line within %.2f nm and were dropped.\n', n_unmatched, TOL_VALID);
    end
    pixel_v = pixel_v(ok_v); meas_v = meas_v(ok_v); is_Hg_v = is_Hg_v(ok_v); id_v = id_v(ok_v);
    nv = numel(pixel_v);
    res_v = meas_v - id_v;

    fprintf('\nSub-pixel validation RMS (evaluated at ADOPTED coefficients): %.4f nm  (n=%d)\n', rms(res_v), nv);
    fprintf('min/max residual: %+.4f / %+.4f nm\n', min(res_v), max(res_v));
    fprintf('For comparison, the in-sample fit RMS was %.4f nm.\n', rms_wt);

    bw_v = nan(nv,1);  ctr_v = nan(nv,1);  have_ctr = false;   % not available from raw spectra alone

else
    % ---- legacy fallback: pre-extracted peak tables (NOTE 13 caveat) ----
    fprintf('\n\n=== SECTION 11: VALIDATION (peak tables -- LOW RESOLUTION, see NOTE 13) ===\n');
    % val_desc aqui é o OPOSTO do ramo acima: meas_v vem da coluna
    % "Wavelength" que o próprio instrumento reportou para cada pico, sob
    % QUALQUER calibração esteja ativa agora -- não necessariamente a
    % adotada. Este RMS responde "o que está gravado bate com o catálogo?",
    % não "meu ajuste é bom?".
    val_desc = 'ACTIVE/as-written calibration reported by the instrument (NOT necessarily the adopted fit)';

    [pv_hg, wv_hg, bwv_hg, ctrv_hg] = read_validation(HG_VALIDATION);
    [pv_ne, wv_ne, bwv_ne, ctrv_ne] = read_validation(NE_VALIDATION);
    pixel_v  = [pv_hg(:); pv_ne(:)];
    meas_v0  = [wv_hg(:); wv_ne(:)];
    bw_v     = [bwv_hg(:); bwv_ne(:)];
    ctr_v    = [ctrv_hg(:); ctrv_ne(:)];
    is_Hg_v  = [true(numel(pv_hg),1); false(numel(pv_ne),1)];

    [id_v, ok_v] = match(meas_v0, ALL_LINES, TOL_VALID);
    n_unmatched = sum(~ok_v);
    if n_unmatched > 0
        fprintf('%d point(s) matched no catalogued line within %.2f nm and were dropped.\n', n_unmatched, TOL_VALID);
    end
    pixel_v = pixel_v(ok_v); meas_v = meas_v0(ok_v); bw_v = bw_v(ok_v);
    ctr_v = ctr_v(ok_v); is_Hg_v = is_Hg_v(ok_v); id_v = id_v(ok_v);
    nv = numel(pixel_v);
    res_v = meas_v - id_v;
    have_ctr = any(isfinite(ctr_v));

    fprintf('Raw (quantized) validation RMS: %.4f nm  (n=%d)\n', rms(res_v), nv);
    fprintf('min/max residual: %+.4f / %+.4f nm\n', min(res_v), max(res_v));

    if have_ctr
        gap = abs(meas_v - ctr_v);
        suspect_ctr = isfinite(gap) & (gap > 0.30);
        if any(suspect_ctr)
            fprintf('%d point(s) show |Wavelength - Center Wavelength| > 0.30 nm:\n', sum(suspect_ctr));
            for i = find(suspect_ctr)', fprintf('   pixel %d: gap = %.3f nm\n', pixel_v(i), gap(i)); end
        end
    end

    % Fallback smooth-drift model (NOTE 11's machinery reused, NOTE 12) --
    % only meaningful once a coefficient-level error has been ruled out,
    % which the peak-table path cannot do (no raw spectrum to extract from).
    drift_ok = robust_subset(pixel_v, res_v, DRIFT_DEGREE, MAD_REJECT);
    [c_drift_s, ~, mu_drift] = polyfit(pixel_v(drift_ok), res_v(drift_ok), DRIFT_DEGREE);
    c_drift = scaled_to_raw(c_drift_s, mu_drift, DRIFT_DEGREE);
    rms_after_drift = rms(res_v(drift_ok) - polyval(c_drift, pixel_v(drift_ok)));
    fprintf('\nSmooth drift model (degree %d): explains %d of %d points, RMS after removal %.4f nm\n', ...
            DRIFT_DEGREE, sum(drift_ok), nv, rms_after_drift);
    if any(~drift_ok)
        fprintf('Points NOT explained by the trend (pixel):'); fprintf(' %d', pixel_v(~drift_ok)); fprintf('\n');
    end
end

if have_validation
    fonte_v = repmat("Ne", numel(pixel_v), 1); fonte_v(is_Hg_v) = "Hg";
    Tv = table(pixel_v, id_v, meas_v, res_v, fonte_v, 'VariableNames', ...
        {'pixel','lambda_NIST_nm','measured_nm','residual_nm','source'});
    writetable(Tv,'validation_residuals.csv');
    fprintf('\nValidation results written: validation_residuals.csv\n');

    % Same handle-based legend pattern as subplot 1 (Section 10): the
    % zero-line is excluded via HandleVisibility so it can never contend
    % for a legend slot, and only the two real series are passed by handle.
    subplot(n_rows,3,7); hold on; grid on; box on; yline(0,'k-','HandleVisibility','off');
    h1 = plot(pixel_v(is_Hg_v),  res_v(is_Hg_v),  'o', 'Color',cHg,'MarkerFaceColor',cHg,'MarkerSize',5);
    h2 = plot(pixel_v(~is_Hg_v), res_v(~is_Hg_v), 's', 'Color',cNe,'MarkerFaceColor',cNe,'MarkerSize',5);
    legend([h1 h2], {'Hg','Ne'}, 'Location','best');
    xlabel('Pixel'); ylabel('Residual (nm)');
    title(sprintf('Validation vs. catalogue (RMS %.4f nm)  --  %s', rms(res_v), val_desc), ...
          'FontSize', 9);
    fprintf('\nValidation mode: %s\n', val_desc);
end

fprintf(['\nREMINDER: Sections 1-10 compute a fit; whether that fit is the one ' ...
    'currently written\nto the instrument, coefficient by coefficient, is what ' ...
    'Section 11 checks (NOTE 14/15). Keep\nthe factory backup.\n']);

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [pixel, counts] = read_spectrum(filename, factory_coef)
% Reads a two-column SpectraSuite export and verifies its wavelength axis
% matches factory_coef evaluated at integer pixel (used for the FITTING
% spectra, which are acquired under the known factory calibration).
    M = readmatrix(filename);
    counts = M(:,2);  pixel = (0:numel(counts)-1)';
    dev = max(abs(M(:,1) - polyval(fliplr(factory_coef), pixel)));
    if dev > 0.01
        warning(['The wavelength axis of %s does not match the factory ' ...
            'calibration given (max deviation %.4f nm).'], filename, dev);
    end
end

function [pixel, counts] = read_spectrum_novalidate(filename)
% Same as read_spectrum, but for a VALIDATION spectrum whose active
% calibration is exactly what we do not want to assume (NOTE 14) -- so no
% consistency check against any particular coefficient set is performed
% here. The row index is still the pixel by construction of the export.
    M = readmatrix(filename);
    counts = M(:,2);  pixel = (0:numel(counts)-1)';
end

function [a_active, max_axis_resid] = extract_active_coefficients(filename, degree)
% NOTE 14/15. Recovers the polynomial actually encoding a spectrum file's
% own wavelength column, by fitting it directly against the row index --
% with no assumption about which coefficients SHOULD be there. Returns
% ascending-order coefficients [a0 a1 a2 a3], matching the convention used
% for the ADOPTED coefficients elsewhere in this script, and the max
% residual of that self-fit -- the caller uses this residual (NOTE 15) to
% decide whether a_active is trustworthy before comparing it to anything.
    M = readmatrix(filename);
    wl = M(:,1);  idx = (0:numel(wl)-1)';
    mu = mean(idx); sd = std(idx);
    c_s = polyfit((idx-mu)/sd, wl, degree);
    coef_raw = scaled_to_raw(c_s, [mu,sd], degree);   % descending powers
    a_active = fliplr(coef_raw);                       % -> [a0 a1 a2 a3]
    max_axis_resid = max(abs(polyval(coef_raw, idx) - wl));
end

function report_coef_check(a_active, a_adopted, rel_tol)
% Prints a_active alongside the adopted fit, coefficient by coefficient,
% and flags any that differ by more than rel_tol of the adopted value's
% magnitude -- the check that caught a3 silently holding its factory
% value (NOTE 14). Callers must gate this on the axis self-fit residual
% first (NOTE 15) -- this function does not know whether a_active came
% from a genuine spectrum.
    names = {'a0 (Intercept)','a1 (1st Coeff.)','a2 (2nd Coeff.)','a3 (3rd Coeff.)'};
    for i = 1:4
        rel = abs(a_active(i)-a_adopted(i)) / max(abs(a_adopted(i)), eps);
        flag = '';
        if rel > rel_tol
            flag = sprintf('  <-- %.0f%% off ADOPTED: likely NOT written correctly', 100*rel);
        end
        fprintf('  %-16s active=%.6E  adopted=%.6E%s\n', names{i}, a_active(i), a_adopted(i), flag);
    end
end

function bw = read_bandwidths(filename, peak_pos, pix_tol)
    bw = nan(numel(peak_pos),1);
    if ~isfile(filename), fprintf('Bandwidth file not found: %s\n', filename); return; end
    M = readmatrix(filename, 'CommentStyle', '#');
    M = M(all(~isnan(M(:,1)),2), :);
    for i = 1:numel(peak_pos)
        [d, j] = min(abs(M(:,1) - peak_pos(i)));
        if d <= pix_tol, bw(i) = M(j,2); end
    end
    fprintf('%s: %d of %d peaks matched to a bandwidth entry\n', filename, sum(isfinite(bw)), numel(peak_pos));
end

function [pixel, wavelength, bandwidth, center] = read_validation(filename)
% Legacy peak-table reader (NOTE 13 fallback path only).
    if ~isfile(filename)
        pixel = []; wavelength = []; bandwidth = []; center = []; return
    end
    M = readmatrix(filename, 'CommentStyle', '#');
    M = M(all(~isnan(M(:,1:2)),2), :);
    pixel = M(:,1);  wavelength = M(:,2);  bandwidth = M(:,3);
    if size(M,2) >= 4, center = M(:,4); else, center = nan(size(pixel)); end
    fprintf('%s: %d points read\n', filename, numel(pixel));
end

function [pos, rejected] = find_peaks(pixel, counts, frac_thresh, prom_min, sep_min, n_parab, shift_max)
% Sub-pixel peak detection in PIXEL-INDEX space (NOTE 3), shared by the
% fitting spectra (Sections 3-4) and, now, the raw validation spectra
% (Section 11) -- the same method on both sides is what makes the
% validation comparison fair.
    thresh = frac_thresh * max(counts);  N = numel(counts);  h = floor(n_parab/2);
    cand = find(counts(2:end-1) > thresh & counts(2:end-1) >= counts(1:end-2) & ...
                counts(2:end-1) >= counts(3:end)) + 1;
    keep = false(size(cand));
    for k = 1:numel(cand)
        i = cand(k);  pk = counts(i);
        j = i-1; minL = pk; while j >= 1 && counts(j) <= pk, minL = min(minL,counts(j)); j = j-1; end
        j = i+1; minR = pk; while j <= N && counts(j) <= pk, minR = min(minR,counts(j)); j = j+1; end
        keep(k) = (pk - max(minL,minR)) >= prom_min*pk;
    end
    cand = cand(keep);
    cand = sort(cand);  keep = true(size(cand));
    for i = 1:numel(cand)-1
        if ~keep(i), continue; end
        j = i+1;
        while j <= numel(cand) && cand(j)-cand(i) < sep_min
            if counts(cand(j)) > counts(cand(i)), keep(i) = false; else, keep(j) = false; end
            j = j+1;
        end
    end
    cand = cand(keep);
    pos = zeros(numel(cand),1);  rejected = false(numel(cand),1);
    for k = 1:numel(cand)
        i = cand(k);
        neigh = counts(max(1,i-3):min(N,i+3));
        if sum(abs(neigh - counts(i)) < 1e-6) > 1, rejected(k) = true; end
        if i-h < 1 || i+h > N, pos(k) = pixel(i); rejected(k) = true; continue; end
        c = polyfit(pixel(i-h:i+h), counts(i-h:i+h), 2);
        if c(1) >= 0
            pos(k) = pixel(i); rejected(k) = true;
        else
            pos(k) = -c(2)/(2*c(1));
            if abs(pos(k)-pixel(i)) > shift_max, rejected(k) = true; end
        end
    end
end

function [matched, ok] = match(lambda_meas, table, tol)
    matched = nan(size(lambda_meas));  ok = false(size(lambda_meas));
    for i = 1:numel(lambda_meas)
        [d,j] = min(abs(table - lambda_meas(i)));
        if d <= tol, matched(i) = table(j); ok(i) = true; end
    end
    idx = find(ok);  [u,~,g] = unique(matched(ok));
    for k = 1:numel(u)
        grp = idx(g==k);
        if numel(grp) > 1
            [~,best] = min(abs(lambda_meas(grp)-u(k)));
            drop = grp; drop(best) = [];
            ok(drop) = false;  matched(drop) = NaN;
        end
    end
end

function [a, r, rms_w] = fit_poly(pixel, target, weight, degree)
    mu = mean(pixel);  sd = std(pixel);
    X  = vandermonde((pixel-mu)/sd, degree);
    W  = diag(weight(:));
    c  = (X'*W*X) \ (X'*W*target(:));
    r  = X*c - target(:);
    rms_w = sqrt(sum(weight(:).*r.^2)/sum(weight(:)));
    coef_raw = scaled_to_raw(c, [mu,sd], degree);
    a = fliplr(coef_raw);
    assert(max(abs(polyval(coef_raw,pixel) - X*c)) < 1e-6, 'Scaled-to-raw coefficient conversion is inconsistent.');
end

function X = vandermonde(s, degree)
    X = zeros(numel(s), degree+1);
    for k = 0:degree, X(:,k+1) = s(:).^(degree-k); end
end

function coef_raw = scaled_to_raw(coef_scaled, mu, degree)
% Converts polynomial coefficients from a centred/scaled predictor
% s=(x-mu(1))/mu(2) back to coefficients in raw x (descending powers,
% ready for polyval against raw x directly). Shared by fit_poly, the
% robust anchor/drift model, and extract_active_coefficients -- a single
% implementation, so a scale-mismatch bug (found and fixed in an earlier
% version of this script, where a similar conversion was duplicated and
% one copy forgot it) cannot recur in one place without recurring in all.
    coef_raw = zeros(1, degree+1);
    for k = 1:degree+1
        g = degree-k+1;  term = 1;
        for j = 1:g, term = conv(term, [1/mu(2), -mu(1)/mu(2)]); end
        coef_raw(end-g:end) = coef_raw(end-g:end) + coef_scaled(k)*term;
    end
end

function keep = robust_subset(x, y, degree, k)
% Leave-one-out-robust outlier rejection (NOTE 11). Shared by the Hg
% anchor and the legacy drift-model fallback in Section 11.
    keep = true(numel(x),1);
    for it = 1:numel(x)
        if sum(keep) <= degree+2, break; end
        [c,~,mu] = polyfit(x(keep), y(keep), degree);
        s = (x - mu(1))/mu(2);
        r = polyval(c, s) - y;
        Xall = vandermonde(s, degree);  Xk = Xall(keep,:);
        h = min(max(sum((Xall/(Xk'*Xk)).*Xall, 2),0),0.999);
        rl = r ./ (1 - h);
        m  = median(abs(rl(keep) - median(rl(keep))));
        if m <= 0, break; end
        bad = keep & (abs(rl) > k*m);
        if ~any(bad), break; end
        [~,j] = max(abs(rl).*bad);
        keep(j) = false;
    end
end

function c = bandwidth_class(ratio, clean_ratio)
    if     isnan(ratio),         c = 'n/a';
    elseif ratio <= clean_ratio, c = 'clean';
    elseif ratio <= 2.5,         c = 'broadened';
    elseif ratio <= 6.0,         c = 'blend';
    else,                        c = 'SEVERE BLEND';
    end
end