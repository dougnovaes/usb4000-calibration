%% WAVELENGTH_CALIBRATION
% Wavelength-axis calibration of an Ocean Optics USB4000 spectrometer from
% Hg and Ne emission lines, using sub-pixel peak localisation and a
% weighted fit driven by the measured 90% bandwidth of each line.
%
% =========================================================================
% SUMMARY
% =========================================================================
%   Fits  lambda(pixel) = a0 + a1*p + a2*p^2 + a3*p^3  — the polynomial the
%   USB4000 firmware accepts — to catalogued emission lines.
%
%   What distinguishes this from a naive fit:
%     1. peak positions from a parabola through the top of each peak, not
%        from the peak pixel: ~10x better positional precision (NOTE 3);
%     2. a prominence filter that rejects noise ripples on the flank of
%        strong lines (NOTE 7);
%     3. effective-wavelength correction for unresolved blends (NOTE 8);
%     4. weighted least squares driven by the instrument-reported 90%
%        bandwidth of each line, which doubles as an error bar and as a
%        blend detector (NOTE 9).
%
% =========================================================================
% CONTEXT
% =========================================================================
%   Stage A of an absolute-calibration procedure for a USB4000, used as
%   training before the same work on Jobin-Yvon THR1000 and McPherson 207
%   monochromators. End application: spectral characterisation of a cold
%   dielectric-barrier-discharge (DBD) plasma.
%
%   The detector does not measure wavelength — it only records which of its
%   3648 pixels received light, and how much. The pixel-to-wavelength
%   relation is a property of that particular optical assembly and must be
%   measured by pointing the instrument at sources of known wavelength.
%   Because the grating equation relates wavelength to diffraction ANGLE
%   rather than directly to position on a flat detector, the relation
%   carries residual curvature — hence a cubic polynomial.
%
% =========================================================================
% INPUTS
% =========================================================================
%   hg_spectrum.txt, ne_spectrum.txt   (required)
%       Two tab-separated columns exported from SpectraSuite:
%       column 1 = wavelength under the FACTORY calibration [nm]
%       column 2 = counts (with Electric Dark Correction enabled)
%       One row per pixel, 3648 rows.
%
%   hg_bandwidths.txt, ne_bandwidths.txt   (optional but recommended)
%       Two columns: detector pixel index, and the "90% Bandwidth" [nm]
%       reported by the SpectraSuite Show Peak Info panel. Use NaN for
%       peaks whose reported value is unstable. Lines beginning with '#'
%       are comments.
%
%       If these files are absent the script falls back to UNIFORM weights
%       and says so. It does NOT attempt to recompute the bandwidth — see
%       NOTE 9 for why.
%
% =========================================================================
% OUTPUTS
% =========================================================================
%   calibration_coefficients.txt  — coefficients for all fitting schemes,
%                                   plus the factory backup
%   calibration_residuals.csv     — per-line table
%   a six-panel diagnostic figure
%
% =========================================================================
% REFERENCE RESULT (for verification)
% =========================================================================
%   Bandwidth-weighted fit, 31 lines (254-744 nm):
%       a0 =  1.7635567179E+02
%       a1 =  2.1606649029E-01
%       a2 = -3.6985346142E-06
%       a3 = -5.5381922989E-10
%       weighted RMS ............... 0.0339 nm
%       RMS over clean lines only ... 0.0185 nm
%       factory calibration ......... 1.0419 nm   (~31x improvement)
%
%   Clean lines only, 11 lines (405-725 nm):
%       a0 =  1.7701734777E+02
%       a1 =  2.1485460168E-01
%       a2 = -3.0222362698E-06
%       a3 = -6.7253730018E-10
%       RMS ........................ 0.0133 nm
%
% =========================================================================
% METHODOLOGICAL NOTES
% =========================================================================
% Each note records a mistake made during development and the defence
% adopted. They are here because they are the kind of trap that does not
% appear in any manual and cost several iterations to diagnose.
%
% -------------------------------------------------------------------------
% NOTE 1 — Vacuum vs. air wavelengths
% -------------------------------------------------------------------------
%   Reference tables often list VACUUM wavelengths above 200 nm without
%   flagging it. The USB4000 measures in air. When querying the NIST ASD,
%   select the output option that already returns air wavelengths, which
%   removes the manual conversion (factor ~1 - 2.8e-4 in the visible) and
%   the error that comes with it.
%
% -------------------------------------------------------------------------
% NOTE 2 — Do not filter by intensity when querying the catalogue
% -------------------------------------------------------------------------
%   The curated "Strong Lines" table of the NIST Handbook is a small subset
%   and leaves most detected peaks without a reliable match. Request the
%   full list (leave "Relative intensity minimum" blank) and filter later
%   using your own criterion.
%
% -------------------------------------------------------------------------
% NOTE 3 — Sub-pixel position: parabola on the core, not a Gaussian on the
%          whole profile
% -------------------------------------------------------------------------
%   Using the peak pixel limits precision to about +-0.1 nm (half a pixel),
%   which is comparable to the spacing between catalogued lines — the
%   identification is then ambiguous by construction. A parabola through
%   the three central points gives the position to about 1/10 pixel
%   (~0.02 nm).
%
%   Why a parabola on the core and NOT a Gaussian over the whole profile:
%   USB4000 peaks have an asymmetric red-side tail (coma, typical of
%   compact Czerny-Turner designs). A Gaussian fitted to the whole profile
%   is pulled by that tail; a parabola through the three central points
%   sees only the core, where the asymmetry is still negligible.
%
%   Measured on this dataset: the centroid sits systematically ~0.38 pixel
%   (~0.08 nm) redward of the parabola vertex. That is a systematic bias —
%   it does not average out with more measurements.
%
%   Net effect: the scatter of the (factory - NIST) offset across the Hg
%   lines dropped from ~0.20 nm (integer pixel) to 0.057 nm (sub-pixel).
%
% -------------------------------------------------------------------------
% NOTE 4 — The SpectraSuite peak threshold is global
% -------------------------------------------------------------------------
%   The software uses a single horizontal threshold for the entire graph.
%   With peaks of very different heights, a threshold suited to the tallest
%   peak makes the shorter ones report an absurd 90% bandwidth (hundreds of
%   nm) even when they are real and narrow. When using the GUI, zoom into a
%   local window and re-position the threshold there. This script does not
%   depend on that: it detects peaks directly from the raw spectrum.
%
% -------------------------------------------------------------------------
% NOTE 5 — Back up the factory coefficients
% -------------------------------------------------------------------------
%   Writing new coefficients OVERWRITES the factory ones. Record them first
%   (see FACTORY_COEF).
%
% -------------------------------------------------------------------------
% NOTE 6 — Saturation
% -------------------------------------------------------------------------
%   A clipped peak has a flat top and the three-point parabola returns a
%   meaningless vertex. The script detects consecutive identical values at
%   the top and rejects the peak. If the strongest line saturates, lower
%   the integration time and re-acquire: losing the brightest line is
%   expensive, particularly in the UV.
%
% -------------------------------------------------------------------------
% NOTE 7 — Prominence filter (indispensable)
% -------------------------------------------------------------------------
%   A height threshold alone is NOT enough. In an earlier version, a noise
%   ripple on the rising flank of a strong line (height 1700 counts, above
%   the 1284-count threshold) was taken as a peak and matched to that same
%   line, corrupting the whole fit — the Ne identification rate collapsed
%   from 18/19 to 4/25.
%
%   The defence is PROMINENCE: how far a peak rises above the higher of the
%   two saddles flanking it. For that ripple the prominence was 31 counts,
%   or 1.8% of its own height; for a real peak the prominence/height ratio
%   is close to 100%. A threshold of 30% separates the two cases with a
%   wide margin and still preserves partially resolved doublets (Ne
%   638.3/640.2 nm has a ratio of 0.97).
%
% -------------------------------------------------------------------------
% NOTE 8 — Unresolved blends: use the EFFECTIVE wavelength
% -------------------------------------------------------------------------
%   When two or more catalogued lines fall inside the instrumental width
%   (~0.95 nm FWHM here), the detector sees a single peak, and the position
%   of that peak is NEITHER the dominant component NOR the intensity-
%   weighted centroid — it is the maximum of the summed profiles.
%
%   Two Hg lines in this set are blends:
%     313.1555 (int. 3000) + 313.1844 (4000)             -> eff. 313.1716
%     365.0158 (9000) + 365.4842 (3000) + 366.2887 (500) -> eff. 365.1205
%
%   The correction on the 365 nm line is +0.105 nm — larger than the RMS of
%   the entire fit. The effective values were obtained by summing the
%   component profiles at the measured instrumental width and locating the
%   maximum of the sum.
%
% -------------------------------------------------------------------------
% NOTE 9 — Why the 90% bandwidth is READ, not recomputed
% -------------------------------------------------------------------------
%   The bandwidth serves two purposes at once:
%
%   (a) BLEND DETECTOR, without consulting a catalogue. The smallest
%       bandwidth in this set is 0.37 nm; that is the bare instrumental
%       response of an isolated line. Anything above it is broadened by
%       something:
%           ratio = bandwidth / min(bandwidth)
%           <=1.3x  clean        1.3-2.5x  broadened
%           2.5-6x  blend        >6x       severe blend
%       Applied here, this criterion independently recovered the two Hg
%       blends found via the catalogue in NOTE 8 (313 and 365 come out at
%       2.2-2.3x) and also flagged the 253.65 nm line (2.3x).
%
%   (b) WEIGHT in the fit, as w = (min_bandwidth / bandwidth)^2, i.e. a
%       positional uncertainty proportional to the width. This replaces the
%       binary include/exclude decision for suspect lines with a continuous
%       grading, which is statistically more defensible and preserves
%       spectral coverage.
%
%   IMPORTANT — why it is read from a file rather than computed here.
%   An earlier version of this script computed the width at 90% of peak
%   height directly from the spectrum. That quantity is NOT what
%   SpectraSuite reports and is useless as a blend detector: it measures
%   only the peak core, which is precisely where a blend does not show. The
%   correlation between the two quantities across this dataset is +0.12,
%   i.e. none, and the disagreement is extreme in the worst cases — the
%   616.36 nm line, which SpectraSuite flags as a severe blend (9.83 nm),
%   came out as the cleanest line of all (1.0x) under the computed metric.
%   Attempts to reproduce the SpectraSuite definition from the raw spectrum
%   (contiguous width above the detection threshold) also failed. The
%   SpectraSuite manual does not define the quantity; it only states that
%   the threshold is set "at the level needed to isolate the desired
%   peaks", which implies the metric is threshold-relative.
%
%   The honest resolution is to treat the instrument-reported value as the
%   authoritative input and read it. If the file is missing, the script
%   uses uniform weights rather than a misleading substitute.
%
%   Concrete case: the 253.65 nm line is the only anchor below 313 nm, but
%   its offset is an outlier (+1.56 nm against 0.93-1.15 nm for the rest),
%   most likely from self-absorption — it is the Hg resonance line, and in
%   a low-pressure lamp the vapour reabsorbs the core of the profile,
%   distorting it. With uniform weights it degrades the RMS from 0.027 to
%   0.044 nm; with bandwidth weighting it receives 19% of the maximum
%   weight, still anchors the UV end, and the weighted RMS is 0.034 nm.
%
% -------------------------------------------------------------------------
% NOTE 10 — Two-stage identification, never peak by peak
% -------------------------------------------------------------------------
%   The Hg lines are identified first (few, strong, well separated,
%   unambiguous). A cubic is fitted to the Hg lines ALONE and extrapolated
%   into the Ne region; the ENTIRE pattern of Ne peaks must then land on
%   catalogued lines at once.
%
%   It is this global check — many peaks simultaneously, not one at a time
%   — that gives confidence in the identification. Nearest-neighbour
%   matching peak by peak ALWAYS finds some line and therefore hides
%   errors: in the 580-730 nm range there is a catalogued argon line every
%   0.78 nm on average.
%
%   This is exactly how the second lamp was found to be NEON and not argon:
%   under the argon hypothesis no smooth model explained more than 7 of 13
%   peaks; under neon, 18 of 19 matched within 0.25 nm with no free
%   parameter at all. The negative test was equally decisive — the
%   strongest Ar I lines (696.5 / 706.7 / 750.4 / 763.5 / 811.5 nm) are
%   ABSENT from the spectrum.
%
% =========================================================================
% REFERENCES
% =========================================================================
%   Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024).
%   NIST Atomic Spectra Database (ver. 5.12). National Institute of
%   Standards and Technology, Gaithersburg, MD.
%   https://physics.nist.gov/asd     DOI: 10.18434/T4W30F
%
%   Ocean Optics. SpectraSuite Installation and Operation Manual.
%
% =========================================================================
% STATUS
% =========================================================================
%   The coefficients computed here have NOT been written to the instrument.
%   Before writing: confirm from the lamp label that it is an Hg-Ne lamp,
%   and preserve the factory backup.
%
% Author: (fill in)                                        Licence: MIT 
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURATION
%  ========================================================================
HG_SPECTRUM  = 'hg_spectrum.txt';
NE_SPECTRUM  = 'ne_spectrum.txt';
HG_BANDWIDTH = 'hg_bandwidths.txt';   % optional; see NOTE 9
NE_BANDWIDTH = 'ne_bandwidths.txt';

DEGREE     = 3;      % polynomial degree (3 = what the firmware accepts)
N_PARAB    = 3;      % points in the top parabola (NOTE 3)
FRAC_THRESH= 0.02;   % height threshold, as a fraction of the tallest peak
PROM_MIN   = 0.30;   % minimum prominence / height (NOTE 7)
SEP_MIN    = 5;      % minimum separation between peaks, in pixels
SHIFT_MAX  = 1.0;    % max vertex-to-maximum distance [pixel] before reject
TOL_HG     = 2.0;    % matching tolerance [nm] for Hg lines
TOL_NE     = 0.30;   % matching tolerance [nm] for Ne lines (NOTE 10)
RATE_MIN   = 0.60;   % minimum Ne pattern hit rate before warning
CLEAN_RATIO= 1.3;    % bandwidth ratio below which a line counts as clean
PIX_TOL    = 3;      % max pixel distance when matching a bandwidth entry
N_PIX      = 3648;
MAD_REJECT = 4;      % anchor outlier rejection, in MAD units (NOTE 11)

% Factory coefficients — BACKUP (NOTE 5), ascending order [a0 a1 a2 a3]
FACTORY_COEF = [178.085770, 2.1517764E-1, -3.3423826E-6, -5.9958050E-10];

%% ========================================================================
%  2. REFERENCE LINE TABLES (NIST ASD, air wavelengths)
%  ========================================================================
% Column 1: wavelength used in the fit (blend-corrected, NOTE 8)
% Column 2: catalogue wavelength of the dominant component, for the record
HG_LINES = [ ...
    253.6521, 253.6521; ...  % resonance line; self-absorption (NOTE 9)
    313.1716, 313.1555; ...  % BLEND 313.1555 + 313.1844
    365.1205, 365.0158; ...  % BLEND 365.0158 + 365.4842 + 366.2887
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

%% ========================================================================
%  3. READ SPECTRA
%  ========================================================================
[pix_hg, cnt_hg] = read_spectrum(HG_SPECTRUM, FACTORY_COEF);
[pix_ne, cnt_ne] = read_spectrum(NE_SPECTRUM, FACTORY_COEF);
fprintf('Spectra loaded: Hg (%d points), Ne (%d points)\n', ...
        numel(cnt_hg), numel(cnt_ne));

%% ========================================================================
%  4. PEAK DETECTION AND SUB-PIXEL POSITION
%  ========================================================================
[p_hg, rej_hg] = find_peaks(pix_hg, cnt_hg, FRAC_THRESH, PROM_MIN, ...
                            SEP_MIN, N_PARAB, SHIFT_MAX);
[p_ne, rej_ne] = find_peaks(pix_ne, cnt_ne, FRAC_THRESH, PROM_MIN, ...
                            SEP_MIN, N_PARAB, SHIFT_MAX);
if any(rej_hg), fprintf('NOTE: %d Hg peak(s) rejected.\n', sum(rej_hg)); end
if any(rej_ne), fprintf('NOTE: %d Ne peak(s) rejected.\n', sum(rej_ne)); end
p_hg = p_hg(~rej_hg);  p_ne = p_ne(~rej_ne);
fprintf('Usable peaks: Hg = %d, Ne = %d\n', numel(p_hg), numel(p_ne));

%% ========================================================================
%  5. TWO-STAGE IDENTIFICATION (NOTE 10)
%  ========================================================================
% --- Stage 1: Hg -------------------------------------------------------
lam_fac = polyval(fliplr(FACTORY_COEF), p_hg);
[id_hg, ok_hg] = match(lam_fac, HG_LINES(:,1), TOL_HG);
p_hg = p_hg(ok_hg);  lam_hg = id_hg(ok_hg);
fprintf('\nHg lines identified: %d\n', numel(p_hg));
assert(numel(p_hg) >= DEGREE+2, ...
    'Too few Hg lines (%d) to anchor the fit.', numel(p_hg));

% Consistency check: the (factory - NIST) offset must vary SMOOTHLY with
% pixel, because the factory calibration and the true curve are both smooth
% polynomials. A point off the trend betrays a wrong identification.
offset = polyval(fliplr(FACTORY_COEF), p_hg) - lam_hg;
% --- Stage 1b: robust anchor subset (NOTE 11) --------------------------
% The anchor and the final fit have opposite requirements: the anchor needs
% only lines of RELIABLE POSITION, because one bad line poisons the
% extrapolation that identifies the Ne pattern; the final fit wants maximum
% spectral COVERAGE, and the weighting already handles poor lines. A line
% may therefore be excluded here and still enter the final fit.
anchor = robust_anchor(p_hg, lam_hg, DEGREE, MAD_REJECT);
if ~all(anchor)
    fprintf('Excluded from the ANCHOR (still used in the final fit):');
    fprintf(' %.4f', lam_hg(~anchor));
    fprintf(' nm\n');
    fprintf('Anchor offset scatter: %.3f nm (was %.3f with all lines)\n', ...
            std(offset(anchor)), std(offset));
end

% --- Stage 2: Ne, by extrapolating the robust Hg anchor ----------------
[c_hg, ~, mu_hg] = polyfit(p_hg(anchor), lam_hg(anchor), DEGREE);
lam_extrap = polyval(c_hg, (p_ne - mu_hg(1))/mu_hg(2));
[id_ne, ok_ne] = match(lam_extrap, NE_LINES, TOL_NE);
p_ne = p_ne(ok_ne);  lam_ne = id_ne(ok_ne);
rate = numel(p_ne)/numel(ok_ne);
fprintf('Ne lines identified: %d of %d peaks (%.0f%%, NO free parameter)\n', ...
        numel(p_ne), numel(ok_ne), 100*rate);
if rate < RATE_MIN
    warning(['Only %.0f%% of the Ne peaks landed on catalogued lines. ' ...
        'When the identification is correct this rate exceeds 90%%. ' ...
        'Likely causes: (a) a false peak corrupting the Hg cubic ' ...
        '(NOTE 7); (b) the lamp is not Ne (NOTE 10). Do NOT write these ' ...
        'coefficients without investigating.'], 100*rate);
end

%% ========================================================================
%  6. BANDWIDTHS, CLASSIFICATION AND WEIGHTS (NOTE 9)
%  ========================================================================
pixel = [p_hg(:);  p_ne(:)];
target= [lam_hg(:); lam_ne(:)];
is_Hg = [true(numel(p_hg),1); false(numel(p_ne),1)];
[pixel, ord] = sort(pixel);  target = target(ord);  is_Hg = is_Hg(ord);

bw = [read_bandwidths(HG_BANDWIDTH, p_hg, PIX_TOL); ...
      read_bandwidths(NE_BANDWIDTH, p_ne, PIX_TOL)];
bw = bw(ord);

have_bw = any(isfinite(bw));
if have_bw
    keep   = isfinite(bw);
    bw_min = min(bw(keep));
    ratio  = bw / bw_min;
    weight = (bw_min ./ bw).^2;
    weight(~keep) = 0;                 % unusable bandwidth -> excluded
    if any(~keep)
        fprintf(['NOTE: %d line(s) have no usable bandwidth and were ' ...
                 'given zero weight.\n'], sum(~keep));
    end
else
    warning(['No bandwidth file found. Falling back to UNIFORM weights. ' ...
        'The script deliberately does not recompute the bandwidth from ' ...
        'the spectrum — see NOTE 9.']);
    keep = true(size(bw));  bw_min = NaN;  ratio = nan(size(bw));
    weight = ones(size(bw));
end

pixel = pixel(keep);  target = target(keep);  is_Hg = is_Hg(keep);
bw = bw(keep);  ratio = ratio(keep);  weight = weight(keep);
n = numel(pixel);
clean = have_bw & (ratio <= CLEAN_RATIO);

if have_bw
    fprintf('\n=== 90%% BANDWIDTH CLASSIFICATION ===\n');
    fprintf('Smallest bandwidth = %.3f nm  ->  instrumental FWHM ~ %.2f nm\n', ...
            bw_min, bw_min/0.39);
    fprintf('%10s %9s %7s %8s  %s\n','lambda','bandwidth','ratio','weight','class');
    for i = 1:n
        fprintf('%10.4f %9.2f %6.1fx %8.3f  %s\n', target(i), bw(i), ...
                ratio(i), weight(i), classify(ratio(i), CLEAN_RATIO));
    end
    fprintf('Clean lines: %d of %d\n', sum(clean), n);
end

%% ========================================================================
%  7. FIT — THREE SCHEMES, FOR COMPARISON
%  ========================================================================
% NUMERICAL CONDITIONING: fitting in raw pixel gives an ill-conditioned
% Vandermonde matrix (cond ~ 1e11 at degree 3 with pixels up to 2800). The
% predictor is centred and scaled internally (cond ~ 9) and the
% coefficients are converted back to raw pixel at the end, because that is
% the form the instrument needs. The conversion is checked by assertion.
[a_wt,   r_wt,   rms_wt  ] = fit_poly(pixel, target, weight,      DEGREE);
[a_un,   r_un,   rms_un  ] = fit_poly(pixel, target, ones(n,1),   DEGREE);
if sum(clean) >= DEGREE+2
    [a_cl, r_cl, rms_cl] = fit_poly(pixel(clean), target(clean), ...
                                    ones(sum(clean),1), DEGREE);
else
    a_cl = nan(1,DEGREE+1); rms_cl = NaN;
end
rms_fac = rms(polyval(fliplr(FACTORY_COEF), pixel) - target);

fprintf('\n=== COMPARISON OF FITTING SCHEMES ===\n');
fprintf('%-28s %4s %11s %12s\n','scheme','n','RMS (nm)','range (nm)');
fprintf('%-28s %4d %11.4f %5.0f-%.0f\n','bandwidth-weighted (NOTE 9)', ...
        n, rms_wt, min(target), max(target));
fprintf('%-28s %4d %11.4f %5.0f-%.0f\n','uniform', n, rms_un, ...
        min(target), max(target));
if sum(clean) >= DEGREE+2
    fprintf('%-28s %4d %11.4f %5.0f-%.0f\n','clean lines only', sum(clean), ...
            rms_cl, min(target(clean)), max(target(clean)));
    fprintf('\nWeighted fit evaluated on clean lines only: %.4f nm\n', ...
            rms(r_wt(clean)));
end
fprintf('%-28s %4d %11.4f\n','factory calibration', n, rms_fac);

% Recommended choice: the weighted fit (full coverage, each point counting
% for what it is worth). Change here to adopt another scheme.
a = a_wt;  residual = r_wt;
fprintf('\n=== ADOPTED COEFFICIENTS (weighted) ===\n');
fprintf('Intercept        (a0) = %.6f\n',  a(1));
fprintf('1st Coefficient  (a1) = %.10E\n', a(2));
fprintf('2nd Coefficient  (a2) = %.10E\n', a(3));
fprintf('3rd Coefficient  (a3) = %.10E\n', a(4));

%% ========================================================================
%  8. CROSS-VALIDATION AND LEVERAGE
%  ========================================================================
% The RMS above is measured on the same points used to fit. Leave-one-out
% cross-validation removes one point, refits and predicts the held-out
% point. CAUTION when interpreting: a high LOOCV/in-sample ratio is not
% always overfitting — it can be LEVERAGE, when few isolated points anchor
% one end of the range. Compare the ratio against the maximum leverage
% before concluding.
loo = zeros(n,1);
for i = 1:n
    m = true(n,1); m(i) = false;
    ai = fit_poly(pixel(m), target(m), weight(m), DEGREE);
    loo(i) = polyval(fliplr(ai), pixel(i)) - target(i);
end
mu = mean(pixel);  sd = std(pixel);
X  = vandermonde((pixel-mu)/sd, DEGREE);
H  = diag(X*pinv(X'*X)*X');

fprintf('\n=== CROSS-VALIDATION ===\n');
fprintf('LOOCV RMS = %.4f nm (in-sample: %.4f nm), ratio = %.2f\n', ...
        rms(loo), rms(residual), rms(loo)/rms(residual));
fprintf('Leverage: max = %.3f (lambda = %.1f nm), expected mean = %.3f\n', ...
        max(H), target(H==max(H)), (DEGREE+1)/n);
if rms(loo)/rms(residual) > 1.5 && max(H) > 3*(DEGREE+1)/n
    fprintf(['  -> high ratio EXPLAINED BY LEVERAGE: the sampling is\n' ...
             '     unbalanced (few points anchoring one end of the range).\n' ...
             '     This is not overfitting. To reduce it, add reference\n' ...
             '     lines in the sparsely covered region.\n']);
end

%% ========================================================================
%  9. FIGURES
%  ========================================================================
figure('Units','normalized','OuterPosition',[0.03 0.08 0.94 0.84]);
cHg=[0 0.45 0.74]; cNe=[0.85 0.33 0.10]; grey=[0.6 0.6 0.6];
pp = linspace(0, N_PIX, 500);

subplot(2,3,1); hold on; grid on; box on;
plot(pp, polyval(fliplr(a),pp), 'k-','LineWidth',1.4);
plot(pp, polyval(fliplr(FACTORY_COEF),pp), '--','Color',grey,'LineWidth',1.1);
plot(pixel(is_Hg), target(is_Hg), 'o','MarkerFaceColor',cHg,'MarkerEdgeColor','k','MarkerSize',6);
plot(pixel(~is_Hg),target(~is_Hg),'s','MarkerFaceColor',cNe,'MarkerEdgeColor','k','MarkerSize',6);
xlabel('Pixel'); ylabel('Wavelength (nm)');
title('Calibration curve \lambda(pixel)');
legend({'fit','factory','Hg','Ne'},'Location','southeast');

subplot(2,3,2); hold on; grid on; box on;
yline(0,'k-');
sz = 20 + 120*weight/max(weight);
scatter(pixel(is_Hg),  residual(is_Hg),  sz(is_Hg),  cHg,'filled','MarkerEdgeColor','k');
scatter(pixel(~is_Hg), residual(~is_Hg), sz(~is_Hg), cNe,'filled','MarkerEdgeColor','k');
yline( rms_wt,':k'); yline(-rms_wt,':k');
xlabel('Pixel'); ylabel('Residual (nm)');
title(sprintf('Residuals (weighted RMS = %.4f nm); area \\propto weight', rms_wt));

subplot(2,3,3); hold on; grid on; box on;
if have_bw
    plot(pixel, ratio, 'ko-','MarkerFaceColor','k','MarkerSize',4);
    yline(CLEAN_RATIO,'--','clean','Color',[0 .6 0]);
    yline(2.5,'--','blend','Color',[.8 .5 0]);
    yline(6.0,'--','severe','Color',[.8 0 0]);
    set(gca,'YScale','log');
    xlabel('Pixel'); ylabel('bandwidth / min');
    title('Blend diagnostic from the 90% bandwidth');
else
    text(0.5,0.5,{'no bandwidth file','(uniform weights)'}, ...
        'HorizontalAlignment','center','Units','normalized');
    axis off; title('Blend diagnostic unavailable');
end

subplot(2,3,4); hold on; grid on; box on;
plot(pix_hg, max(cnt_hg,1), '-','Color',cHg,'LineWidth',0.6);
for i = find(is_Hg)', xline(pixel(i),':','Color','k'); end
set(gca,'YScale','log'); xlabel('Pixel'); ylabel('Counts');
title('Hg spectrum and lines used');

subplot(2,3,5); hold on; grid on; box on;
plot(pix_ne, max(cnt_ne,1), '-','Color',cNe,'LineWidth',0.6);
for i = find(~is_Hg)', xline(pixel(i),':','Color','k'); end
set(gca,'YScale','log'); xlabel('Pixel'); ylabel('Counts');
title('Ne spectrum and lines used');

subplot(2,3,6); grid on; box on;
plot(pp, polyval(polyder(fliplr(a)), pp), 'k-','LineWidth',1.4);
xlabel('Pixel'); ylabel('d\lambda/dpixel (nm/pixel)');
title('Dispersion (derivative of the polynomial)');

%% ========================================================================
%  10. EXPORT
%  ========================================================================
source = repmat("Ne", n, 1); source(is_Hg) = "Hg";
cls = strings(n,1);
for i = 1:n, cls(i) = string(classify(ratio(i), CLEAN_RATIO)); end
T = table(pixel, target, polyval(fliplr(a),pixel), residual, bw, ratio, ...
          weight, source, cls, 'VariableNames', ...
    {'pixel','lambda_NIST_nm','predicted_nm','residual_nm', ...
     'bandwidth90_nm','bandwidth_ratio','weight','source','class'});
writetable(T,'calibration_residuals.csv');

fid = fopen('calibration_coefficients.txt','w');
fprintf(fid,'# Wavelength calibration - Ocean Optics USB4000\n');
fprintf(fid,'# Generated %s\n#\n', datestr(now,'yyyy-mm-dd HH:MM'));
fprintf(fid,'# Sources: Hg and Ne lamps; reference wavelengths from the\n');
fprintf(fid,'# NIST ASD v5.12 (DOI 10.18434/T4W30F), air wavelengths.\n');
fprintf(fid,'# Peak positions: %d-point sub-pixel parabola.\n', N_PARAB);
fprintf(fid,'# Weights: w = (min_bandwidth / bandwidth_90pct)^2\n#\n');
fprintf(fid,'[WEIGHTED]    n=%d  RMS=%.4f nm  range=%.1f-%.1f nm\n', ...
        n, rms_wt, min(target), max(target));
fprintf(fid,'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n\n', a_wt);
if sum(clean) >= DEGREE+2
fprintf(fid,'[CLEAN_ONLY]  n=%d  RMS=%.4f nm  range=%.1f-%.1f nm\n', ...
        sum(clean), rms_cl, min(target(clean)), max(target(clean)));
fprintf(fid,'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n\n', a_cl);
end
fprintf(fid,'[UNIFORM]     n=%d  RMS=%.4f nm\n', n, rms_un);
fprintf(fid,'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n\n', a_un);
fprintf(fid,'[FACTORY_BACKUP]  RMS=%.4f nm\n', rms_fac);
fprintf(fid,'Intercept\t%.6f\n1st\t%.7E\n2nd\t%.7E\n3rd\t%.7E\n', FACTORY_COEF);
fclose(fid);

fprintf('\nFiles written: calibration_residuals.csv, calibration_coefficients.txt\n');
fprintf(['\nREMINDER: these coefficients have NOT been written to the ' ...
    'instrument.\nBefore writing, confirm the lamp is Hg-Ne and preserve ' ...
    'the factory backup.\n']);

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [pixel, counts] = read_spectrum(filename, factory_coef)
% Reads the two-column SpectraSuite export. The wavelength axis in the file
% is the FACTORY calibration evaluated at integer pixel, so the row index
% is the pixel itself. This is verified and a warning is issued otherwise.
    M = readmatrix(filename);
    counts = M(:,2);
    pixel  = (0:numel(counts)-1)';
    dev = max(abs(M(:,1) - polyval(fliplr(factory_coef), pixel)));
    if dev > 0.01
        warning(['The wavelength axis of %s does not match the factory ' ...
            'calibration given (max deviation %.4f nm). Check that ' ...
            'FACTORY_COEF belongs to the instrument that produced this ' ...
            'spectrum.'], filename, dev);
    end
end

function bw = read_bandwidths(filename, peak_pos, pix_tol)
% Reads a two-column file (pixel, 90% bandwidth in nm) and returns the
% bandwidth for each detected peak, matched by pixel proximity. Returns NaN
% for peaks with no entry, and NaN if the file does not exist.
    bw = nan(numel(peak_pos),1);
    if ~isfile(filename)
        fprintf('Bandwidth file not found: %s\n', filename);
        return
    end
    M = readmatrix(filename, 'NumHeaderLines', 0, 'CommentStyle', '#');
    M = M(all(~isnan(M(:,1)),2), :);          % drop malformed rows
    for i = 1:numel(peak_pos)
        [d, j] = min(abs(M(:,1) - peak_pos(i)));
        if d <= pix_tol, bw(i) = M(j,2); end
    end
    fprintf('%s: %d of %d peaks matched to a bandwidth entry\n', ...
            filename, sum(isfinite(bw)), numel(peak_pos));
end

function [pos, rejected] = find_peaks(pixel, counts, frac_thresh, ...
                                      prom_min, sep_min, n_parab, shift_max)
% Detects peaks and returns their sub-pixel position (vertex of the
% parabola through the n_parab central points, NOTE 3). Rejects low-
% prominence peaks (NOTE 7), saturated peaks (NOTE 6), and peaks whose
% vertex falls too far from the maximum.
    thresh = frac_thresh * max(counts);
    N = numel(counts);
    h = floor(n_parab/2);

    cand = find(counts(2:end-1) > thresh & ...
                counts(2:end-1) >= counts(1:end-2) & ...
                counts(2:end-1) >= counts(3:end)) + 1;

    % --- PROMINENCE filter (NOTE 7) ---
    keep = false(size(cand));
    for k = 1:numel(cand)
        i = cand(k);  pk = counts(i);
        j = i-1; minL = pk;
        while j >= 1 && counts(j) <= pk, minL = min(minL,counts(j)); j = j-1; end
        j = i+1; minR = pk;
        while j <= N && counts(j) <= pk, minR = min(minR,counts(j)); j = j+1; end
        keep(k) = (pk - max(minL,minR)) >= prom_min*pk;
    end
    cand = cand(keep);

    % --- suppress neighbours, keeping the tallest of each group ---
    cand = sort(cand);  keep = true(size(cand));
    for i = 1:numel(cand)-1
        if ~keep(i), continue; end
        j = i+1;
        while j <= numel(cand) && cand(j)-cand(i) < sep_min
            if counts(cand(j)) > counts(cand(i)), keep(i) = false;
            else,                                 keep(j) = false; end
            j = j+1;
        end
    end
    cand = cand(keep);

    pos = zeros(numel(cand),1);  rejected = false(numel(cand),1);
    for k = 1:numel(cand)
        i = cand(k);
        % saturation: two or more identical values at the top (NOTE 6)
        neigh = counts(max(1,i-3):min(N,i+3));
        if sum(abs(neigh - counts(i)) < 1e-6) > 1, rejected(k) = true; end
        if i-h < 1 || i+h > N
            pos(k) = pixel(i); rejected(k) = true; continue;
        end
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
% Matches each measured wavelength to the nearest catalogued line within
% tol, preventing two measurements from claiming the same line.
    matched = nan(size(lambda_meas));
    ok = false(size(lambda_meas));
    for i = 1:numel(lambda_meas)
        [d,j] = min(abs(table - lambda_meas(i)));
        if d <= tol, matched(i) = table(j); ok(i) = true; end
    end
    idx = find(ok);
    [u,~,g] = unique(matched(ok));
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
% Weighted least squares with a centred/scaled predictor, converting the
% coefficients back to raw pixel (verified by assertion).
    mu = mean(pixel);  sd = std(pixel);
    X  = vandermonde((pixel-mu)/sd, degree);
    W  = diag(weight(:));
    c  = (X'*W*X) \ (X'*W*target(:));     % descending powers
    r  = X*c - target(:);
    rms_w = sqrt(sum(weight(:).*r.^2)/sum(weight(:)));

    coef_raw = zeros(1, degree+1);
    for k = 1:degree+1
        g = degree-k+1;  term = 1;
        for j = 1:g, term = conv(term, [1/sd, -mu/sd]); end
        coef_raw(end-g:end) = coef_raw(end-g:end) + c(k)*term;
    end
    a = fliplr(coef_raw);                 % [a0 a1 a2 a3]
    assert(max(abs(polyval(coef_raw,pixel) - X*c)) < 1e-6, ...
           'Scaled-to-raw coefficient conversion is inconsistent.');
end

function X = vandermonde(s, degree)
% Vandermonde matrix in descending powers.
    X = zeros(numel(s), degree+1);
    for k = 0:degree, X(:,k+1) = s(:).^(degree-k); end
end

function c = classify(ratio, clean_ratio)
% Blend class from the bandwidth ratio (NOTE 9).
    if     isnan(ratio),         c = 'n/a';
    elseif ratio <= clean_ratio, c = 'clean';
    elseif ratio <= 2.5,         c = 'broadened';
    elseif ratio <= 6.0,         c = 'blend';
    else,                        c = 'SEVERE BLEND';
    end
end

function keep = robust_anchor(pixel, target, degree, k)
% Iteratively drops the line with the largest LEAVE-ONE-OUT residual, as
% long as it exceeds k times the median absolute deviation of those
% residuals. Stops when no line qualifies or only degree+2 lines remain.
%
% Why the LOO residual and not the raw one (NOTE 11). A bad line sitting at
% a high-leverage position — typically at the edge of the sampled range —
% bends the fit until it passes through itself, so its RAW residual becomes
% small and rejection never catches it; the blame is pushed onto interior
% points instead. Measured on this dataset: the 253.65 nm line has leverage
% 0.940 and a raw residual of only +0.028 nm (1.5x MAD, third smallest of
% eight), while the innocent 313.17 nm line shows 3.9x MAD and would be
% removed in its place. Dividing by (1 - h) undoes the self-fitting: the
% amplification is 16.7x for that line, which then stands out at 10.7x MAD.
%
% The MAD is used instead of the standard deviation because a single gross
% outlier inflates the standard deviation enough to hide itself.
    keep = true(numel(pixel),1);
    for it = 1:numel(pixel)
        if sum(keep) <= degree+2, break; end
        [c,~,mu] = polyfit(pixel(keep), target(keep), degree);
        s = (pixel - mu(1))/mu(2);
        r = polyval(c, s) - target;

        % leverage of the kept points, evaluated for every point
        X  = zeros(numel(s), degree+1);
        for q = 0:degree, X(:,q+1) = s(:).^(degree-q); end
        Xk = X(keep,:);
        h  = sum((X/(Xk'*Xk)).*X, 2);            % h_ii for all points
        h  = min(max(h,0), 0.999);               % guard against 1/(1-h)

        rl = r ./ (1 - h);                        % leave-one-out residual
        m  = median(abs(rl(keep) - median(rl(keep))));
        if m <= 0, break; end
        bad = keep & (abs(rl) > k*m);
        if ~any(bad), break; end
        [~,j] = max(abs(rl).*bad);
        keep(j) = false;
    end
end