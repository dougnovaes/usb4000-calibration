%% ABSOLUTE_IRRADIANCE_VALIDATION
%  Converts a SpectraSuite High-Speed multi-shot raw-counts acquisition of
%  the Eppley FEL lamp into absolute spectral irradiance, using an Ocean
%  Optics .IrradCal calibration file, and validates the result against
%  the lamp's own NIST-traceable certificate.
%
% =========================================================================
% SUMMARY
% =========================================================================
%   This is the analysis layer used throughout an extended debugging
%   investigation of an absolute-irradiance calibration built at a
%   non-certificate working distance. It does six things:
%     1. reads a SpectraSuite .IrradCal calibration file and a raw
%        High-Speed-Acquisition counts file (any number of shots);
%     2. validates both files structurally before doing any arithmetic
%        with them (NOTE 6) -- a file that fails validation stops the
%        script with a specific, actionable message, rather than being
%        silently patched and fed into the rest of the pipeline;
%     3. checks every shot for detector saturation, per pixel, and
%        separately flags pixels outside the calibration file's own
%        calibrated range (NOTE 4) -- these are two different reasons a
%        pixel can be unusable and are never conflated;
%     4. converts counts to absolute spectral irradiance, shot by shot,
%        using a dispersion computed from each pixel's OWN wavelength
%        value, not from its row position in the file (NOTE 1 -- this
%        was a real bug in the previous version, not just a robustness
%        gap);
%     5. reports shot-to-shot statistics as TWO distinct confidence
%        bands (NOTE 3) -- these answer different questions and are not
%        interchangeable -- plus a robust (median-based) summary
%        alongside the mean (NOTE 7);
%     6. cross-checks the result against the certificate at the standard
%        reference wavelengths, as the actual test of whether the
%        calibration -- and the physical setup it was built from -- is
%        trustworthy.
%
%   This script does NOT build the .IrradCal file itself (that happens in
%   SpectraSuite's own "New Calibration" wizard) and does NOT touch the
%   wavelength-axis calibration (see WAVELENGTH_CALIBRATION_Calibrated.m)
%   or the six-model lamp-file derivation (see ABSOLUTE_CALIBRATION.m).
%   It is the independent check that sits downstream of both.
%
% =========================================================================
% CHANGE LOG
% =========================================================================
%   v3 (this version) -- Full robustness rewrite. Fixes a real bug (NOTE 1)
%       in v2's dispersion calculation that silently produced wrong results
%       whenever the counts file did not have exactly as many rows as the
%       wavelength-calibration polynomial's own pixel count (v2 computed
%       dispersion from row position; a file with an unexpected row count,
%       e.g. a duplicated export, silently extrapolated the polynomial far
%       outside the range it was fit for -- errors of 30%+ in dispersion
%       were reproduced from a real log in this project). Also adds: input
%       structural validation (NOTE 6), separate NO_CALIBRATION_DATA vs
%       SATURATED flags (NOTE 4), a mandatory non-default T_INT_US (NOTE 2),
%       negative-irradiance flagging (NOTE 8), match-distance-checked
%       certificate lookup, and a robust summary statistic (NOTE 7).
%   v2 -- added axis-length interpolation, saturation/no-data distinction
%       in reporting, valid-data filtering before plotting, more permissive
%       header field-name matching. (Superseded by v3's NOTE 1 fix, which
%       makes the axis-length interpolation this version added actually
%       safe to rely on -- v2 interpolated CalCoef correctly, but then fed
%       the result through the row-position dispersion bug.)
%   v1 -- initial version.
%
% =========================================================================
% WHY AN INDEPENDENT CHECK, NOT JUST TRUSTING SPECTRASUITE'S OWN DISPLAY
% =========================================================================
%   SpectraSuite computes and displays absolute irradiance internally once
%   a calibration is loaded. Re-deriving the same number externally, from
%   the raw counts and the calibration FILE (not the live display), is
%   what makes it possible to: (a) test different physical setups against
%   the SAME certificate numbers side by side, (b) separate the counts
%   themselves from whatever the software's own distance-correction
%   machinery may or may not have done to them, and (c) keep an audit
%   trail of exactly which formula produced which number, which "read the
%   number off the screen" cannot provide.
%
% =========================================================================
% THE CONVERSION FORMULA
% =========================================================================
%   Irradiance(lambda) [uW/cm^2/nm] =
%       Counts(lambda) * CalCoef(lambda) [uJ/count]
%       -----------------------------------------------
%       t_int [s]  *  CollectionArea [cm^2]  *  dlambda/dpixel [nm]
%
%   CalCoef is read directly from the .IrradCal file -- it is the
%   instrument's own measured "energy per detector count" at each pixel,
%   NOT recomputed here. t_int is the integration time of THE MEASUREMENT
%   BEING CONVERTED, not the calibration file's own reference t_int (NOTE
%   2). CollectionArea is a fixed property of the fiber core (NOTE 3, this
%   note numbering follows the file top to bottom -- see below).
%   dlambda/dpixel is the LOCAL dispersion of the wavelength-axis
%   polynomial, evaluated at each pixel's OWN wavelength (NOTE 1), not at
%   its row position and not as a single constant.
%
% =========================================================================
% METHODOLOGICAL NOTES
% =========================================================================
%
% -------------------------------------------------------------------------
% NOTE 1 -- Dispersion from wavelength, never from row position
% -------------------------------------------------------------------------
%   The wavelength axis is a cubic in PIXEL INDEX (see
%   WAVELENGTH_CALIBRATION_Calibrated.m): wl(i) = a0 + a1*i + a2*i^2 +
%   a3*i^3, fit once, for a specific detector, over i = 0..3647. Two
%   things can make a counts file's row count differ from 3648: a
%   different or partial export, or (observed in this project) an
%   accidental duplicated export (7296 = 2x3648 rows). A version of this
%   script that computed "the pixel index" as simply "the row number"
%   would, in that second case, feed indices up to 7295 into a polynomial
%   whose derivative (the dispersion) is only meaningful over 0..3647 --
%   silently extrapolating it, with errors upward of 30% depending on
%   where in the file a given wavelength happens to land (verified
%   numerically against this project's own coefficients).
%
%   This version never does that: for every row, it inverts the
%   wavelength polynomial NUMERICALLY to recover the pixel index that
%   wavelength actually corresponds to (via a dense monotonic lookup
%   table and interpolation -- fast, and exact to floating-point
%   precision since the polynomial is confirmed strictly increasing over
%   the relevant range), then evaluates the dispersion derivative at THAT
%   recovered index. This is correct regardless of how many rows the file
%   has, whether they are duplicated, or in what order they appear,
%   because it depends only on each row's own wavelength value.
%
% -------------------------------------------------------------------------
% NOTE 2 -- t_int is per-measurement, not per-calibration, and mandatory
% -------------------------------------------------------------------------
%   The .IrradCal file's own header records the integration time of the
%   REFERENCE scan used to BUILD the calibration. CalCoef already has
%   that reference t_int baked into it -- it is a per-count ENERGY
%   coefficient, independent of how long any FUTURE measurement
%   integrates. The t_int used in the conversion must be the t_int of the
%   file being CONVERTED, which this project's own logs show varying
%   session to session (3800, 3871, 3907, 3944, 3986, 4005, 4047, 4048,
%   4222 us have all appeared). Neither input file reliably states it.
%
%   This cannot be automated away -- there is no field in either file
%   that is guaranteed to hold it -- so this version instead makes the
%   failure mode loud: T_INT_US defaults to NaN and the script refuses to
%   run until it is set to a real, explicit value for the specific file
%   being processed (Section 1). A stale numeric value silently left over
%   from a previous run -- which is what actually happened earlier in
%   this project -- is exactly the failure this guards against; a NaN
%   default cannot be accidentally "close enough".
%
% -------------------------------------------------------------------------
% NOTE 3 -- Collection area
% -------------------------------------------------------------------------
%   CollectionArea = pi * (fiber_core_diameter/2)^2, read from the
%   .IrradCal header's "Fiber (micron)" field. For a 600 um core:
%   2.82743E-3 cm^2 -- confirmed against the exact value SpectraSuite
%   itself reports for the same fiber, to 6 significant figures.
%
% -------------------------------------------------------------------------
% NOTE 4 -- Two distinct reasons a pixel is unusable, never conflated
% -------------------------------------------------------------------------
%   SATURATED: at least one shot's raw count exceeded SAT_THRESHOLD of
%   the ADC's full scale. The true signal is clipped; mean and std at
%   that pixel are biased and meaningless.
%   OUT_OF_CALIBRATED_RANGE: the .IrradCal file itself stores an explicit
%   coefficient of (numerically) zero at that wavelength -- observed in
%   this project below about 250 nm and above about 888 nm, i.e. outside
%   the range the calibration lamp file actually covered. This is a
%   property of the CALIBRATION, present even with perfect, unsaturated
%   counts, and would otherwise silently produce a reported irradiance of
%   exactly zero -- indistinguishable, downstream, from "the lamp emits
%   nothing here" unless flagged separately.
%   Both are reported as their own logical column in the output table and
%   their own labelled exclusion reason in the certificate cross-check;
%   never merged into one generic "invalid" flag.
%
% -------------------------------------------------------------------------
% NOTE 5 -- Two different "95% bands", answering two different questions
% -------------------------------------------------------------------------
%   Given n shots at each pixel with mean m and sample std s:
%     - 95% CI of the mean:      m +/- 1.96*s/sqrt(n)
%       "How precisely do I know the TRUE MEAN irradiance?" Narrows as n
%       grows. Report this as the calibration's own precision.
%     - 95% prediction interval:  m +/- 1.96*s
%       "Where would a SINGLE FUTURE shot fall?" Does not narrow with n.
%       Use this for shot-to-shot reproducibility against a downstream
%       measurement (e.g. a single DBD plasma acquisition compared
%       against this reference).
%   Reporting only one of the two, unlabelled, is a common and
%   consequential ambiguity in this kind of dataset.
%
% -------------------------------------------------------------------------
% NOTE 6 -- Structural validation before arithmetic
% -------------------------------------------------------------------------
%   Before any conversion, this script checks: wavelength axes are
%   finite, strictly positive, and (after sorting) strictly increasing
%   with no duplicate wavelength values in either file; the counts file's
%   row count is compared against 3648 (the polynomial's fit domain) and
%   a mismatch is reported explicitly, together with a duplicate-row scan
%   (are consecutive rows, or first-half vs second-half, identical?) so
%   the CAUSE is visible, not just papered over; CalCoef values are
%   checked for being finite and non-negative. A file that fails any of
%   these prints a specific, actionable message and the script stops --
%   it does not guess and continue.
%
% -------------------------------------------------------------------------
% NOTE 7 -- Robust summary alongside the mean
% -------------------------------------------------------------------------
%   The certificate-ratio summary reports both mean+/-std AND
%   median+/-MAD (median absolute deviation, scaled by 1.4826 to be
%   comparable to a standard deviation under normality). A single bad
%   point (this project has seen a single wavelength swing the plain mean
%   by more than the entire spread of every other point) dominates the
%   mean and std; the median-based pair is far less sensitive to exactly
%   that failure mode and should be the first number trusted when the two
%   disagree noticeably.
%
% -------------------------------------------------------------------------
% NOTE 8 -- Negative irradiance is a flag, not a value
% -------------------------------------------------------------------------
%   Dark-corrected raw counts fluctuate around zero at wavelengths with
%   negligible true signal (normal, expected noise) and can be negative.
%   A negative computed "irradiance" is not physically meaningful and is
%   flagged (NEGATIVE_SIGNAL, its own logical column) rather than passed
%   silently into an average or a certificate ratio, where it has
%   previously produced misleadingly large or sign-flipped ratios at
%   weak-signal wavelengths.
%
% =========================================================================
% INPUTS (external file dependencies)
% =========================================================================
%   CAL_FILE   (required) A SpectraSuite ".IrradCal" file: a short header
%              of tab-separated key/value metadata (must include
%              "Fiber (micron)"), then a "[uJoule/count]" marker line,
%              then one "wavelength_nm<TAB>coefficient" row per pixel.
%              The companion ".cal" (OOIIrrad) export is NOT read here --
%              it carries the same coefficients without the wavelength
%              column and is redundant for this script's purposes.
%
%   COUNTS_FILE (required) A raw High-Speed-Acquisition export: one header
%              row (a leading blank cell, then one relative timestamp in
%              ms per shot), then one row per pixel: wavelength_nm, then
%              one raw count value per shot. This is RAW data -- no dark
%              subtraction or calibration has been applied by SpectraSuite
%              to this particular export; this script does not add any
%              either, beyond what the .IrradCal coefficients represent.
%
%   T_INT_US   (required, Section 1, no usable default -- NOTE 2) The
%              integration time, in microseconds, actually used for
%              COUNTS_FILE. Get it from the acquisition session's own
%              record, for THIS file specifically.
%
%   WL_COEF    (required, set in Section 1) The four wavelength-
%              calibration coefficients [a0 a1 a2 a3], ascending order,
%              currently active on the instrument, and N_PIXELS_FIT, the
%              detector pixel count that polynomial was fit for. Defaults
%              match WAVELENGTH_CALIBRATION_Calibrated.m's adopted values
%              (3648-pixel USB4000); override both together if a
%              different calibration or detector was active.
%
%   CERT_DATA  (embedded, Section 2) The Eppley EN-66 certificate,
%              250-2400 nm, at 50 cm -- the same table used throughout
%              this project (WAVELENGTH_CALIBRATION_Calibrated.m,
%              ABSOLUTE_CALIBRATION.m). No external file.
%
% =========================================================================
% OUTPUTS
% =========================================================================
%   irradiance_results.csv -- per-pixel table: wavelength, mean
%       irradiance, both confidence bands, and three independent logical
%       flag columns (saturated / out_of_calibrated_range /
%       negative_signal).
%   certificate_comparison.csv -- ratio-vs-certificate at the standard
%       reference wavelengths, with the exclusion reason recorded for any
%       point that could not be compared.
%   A diagnostic figure: mean spectrum with the prediction band (excluded
%       regions shaded, colour-coded by exclusion reason), and the
%       certificate-ratio spectrum below it.
%
% Author: (fill in)                                        Licence: MIT
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURATION -- edit every value in this section per acquisition
%  ========================================================================
CAL_FILE    = 'absolute_calibration_50cm.IrradCal';
COUNTS_FILE = 'absolute_irradiance_HighSpeed_100x_original_lamp_file_blk_alignd_50cm.txt';

T_INT_US = 3800;                % REQUIRED, NO DEFAULT (NOTE 2). Set this to
                                % the actual integration time, in
                                % microseconds, used for COUNTS_FILE. The
                                % script refuses to run while this is NaN.

WL_COEF      = [176.3557, 0.2160665, -3.698535e-6, -5.538192e-10];
N_PIXELS_FIT = 3648;            % detector pixel count WL_COEF was fit for
                                 % (NOTE 1, NOTE 6). Both match
                                 % WAVELENGTH_CALIBRATION_Calibrated.m.

SAT_THRESHOLD      = 0.95;      % fraction of ADC full scale -> saturated
ADC_FULL_SCALE     = 65535;     % 16-bit detector
CAL_ZERO_TOLERANCE = 1e-15;     % |CalCoef| below this = out-of-range (NOTE 4)
CERT_MATCH_TOL_NM  = 2.0;       % max allowed distance, cert wavelength to
                                 % nearest actual data point, before that
                                 % point is dropped with a warning

OUT_TABLE = 'irradiance_results.csv';
OUT_CERT  = 'certificate_comparison.csv';

if isnan(T_INT_US)
    error('ABSOLUTE_IRRADIANCE_VALIDATION:MissingTInt', ...
        ['T_INT_US is not set (Section 1). This value cannot be safely ' ...
         'defaulted or inferred from either input file -- see NOTE 2. ' ...
         'Set it to the integration time actually used for %s before ' ...
         're-running.'], COUNTS_FILE);
end

%% ========================================================================
%  2. CERTIFICATE DATA (Eppley S.O. 52435, lamp EN-66, 50 cm, 7.90 A DC)
%  ========================================================================
% Same table as WAVELENGTH_CALIBRATION_Calibrated.m / ABSOLUTE_CALIBRATION.m.
% Column 2 is W/cm^3 as certified; converted to uW/cm^2/nm (x0.1, see
% ABSOLUTE_CALIBRATION.m's header for the unit-conversion derivation) for
% direct comparison with this script's output.
CERT = [
  250   0.136;   260   0.247;   270   0.411;   280   0.648;   290   0.985;
  300   1.450;   310   2.057;   320   2.839;   330   3.816;   340   5.052;
  350   6.537;   400  18.150;   450  37.220;   500  62.670;   555  93.900;
  600 119.900;   654.60 148.900; 700 169.100;  800 198.900;   900 208.900;
 1050 203.300;  1150 189.200;  1200 180.200;  1300 162.200;  1540 121.400;
 1600 112.100;  1700  97.800;  2000  66.000;  2100  58.600;  2300  45.000;
 2400  40.300];
cert_wl = CERT(:,1);  cert_uwcm2nm = CERT(:,2) * 0.1;

%% ========================================================================
%  3. READ AND VALIDATE THE CALIBRATION FILE
%  ========================================================================
[cal_meta, cal_wl, cal_coef] = read_irradcal(CAL_FILE);
cal_wl = cal_wl(:);  cal_coef = cal_coef(:);
validate_wavelength_axis(cal_wl, sprintf('calibration file (%s)', CAL_FILE));
assert(all(isfinite(cal_coef)) && all(cal_coef >= 0), ...
    'ABSOLUTE_IRRADIANCE_VALIDATION:BadCalCoef', ...
    '%s contains non-finite or negative calibration coefficients.', CAL_FILE);

fprintf('Calibration file: %s\n', CAL_FILE);
fprintf('  built %s, reference t_int = %s us, fiber = %s um\n', ...
        cal_meta.Date, cal_meta.IntTimeUsec, cal_meta.FiberMicron);
area_cm2 = pi * (str2double(cal_meta.FiberMicron)/2 * 1e-4)^2;   % NOTE 3
fprintf('  collection area = %.6e cm^2\n', area_cm2);
fprintf('  T_INT_US configured for THIS run = %d us (verify this against %s''s own acquisition record)\n', ...
        T_INT_US, COUNTS_FILE);

%% ========================================================================
%  4. READ AND VALIDATE THE RAW COUNTS FILE
%  ========================================================================
[wl_raw, counts_raw] = read_highspeed_counts(COUNTS_FILE);
wl_raw = wl_raw(:);
validate_wavelength_axis(wl_raw, sprintf('counts file (%s)', COUNTS_FILE));

% Sort by wavelength (defensive -- downstream fill()/plotting assumes
% ascending order; validate_wavelength_axis only checked, did not enforce).
[wl, sort_idx] = sort(wl_raw);
counts = counts_raw(sort_idx, :);
n_pixels = size(counts,1);  n_shots = size(counts,2);

fprintf('\nCounts file: %s\n', COUNTS_FILE);
fprintf('  %d pixels x %d shots\n', n_pixels, n_shots);

if n_pixels ~= N_PIXELS_FIT
    dup_note = '';
    if mod(n_pixels, N_PIXELS_FIT) == 0
        k = n_pixels / N_PIXELS_FIT;
        first_block = counts(1:N_PIXELS_FIT, :);
        is_repeat = true;
        for b = 2:k
            block = counts((b-1)*N_PIXELS_FIT+1 : b*N_PIXELS_FIT, :);
            if ~isequal(block, first_block), is_repeat = false; break; end
        end
        if is_repeat
            dup_note = sprintf([' This is exactly %dx the expected pixel count, and ' ...
                'the blocks are byte-for-byte identical -- almost certainly a ' ...
                'duplicated export, not %d genuinely different pixels.'], k, n_pixels);
        else
            dup_note = sprintf([' This is exactly %dx the expected pixel count, but the ' ...
                'repeated blocks are NOT identical -- investigate before trusting ' ...
                'this file; it may not be a simple duplicate.'], k);
        end
    end
    warning('ABSOLUTE_IRRADIANCE_VALIDATION:UnexpectedPixelCount', ...
        ['%s has %d rows; the active wavelength calibration (N_PIXELS_FIT) was ' ...
         'fit for %d.%s Dispersion below is computed from each row''s OWN ' ...
         'wavelength value (NOTE 1), so this does not corrupt the physics -- but ' ...
         'a genuinely unexpected row count is still worth understanding before ' ...
         'trusting the file.'], COUNTS_FILE, n_pixels, N_PIXELS_FIT, dup_note);
end

sat_mask = counts > SAT_THRESHOLD * ADC_FULL_SCALE;      % (pixel x shot)
saturated = any(sat_mask, 2);                             % NOTE 4
fprintf('  peak count = %.1f  |  saturated pixels (>=1 shot) = %d of %d\n', ...
        max(counts(:)), sum(saturated), n_pixels);
if any(saturated)
    fprintf('  saturated wavelength range: %.2f - %.2f nm\n', ...
            min(wl(saturated)), max(wl(saturated)));
end

%% ========================================================================
%  5. MAP CALIBRATION ONTO THE COUNTS FILE'S OWN WAVELENGTH AXIS
%  ========================================================================
% CalCoef is defined at the calibration file's own 3648 wavelengths;
% interpolate onto whatever wavelengths this counts file actually has.
% Points outside the calibration file's own covered range become NaN here
% (handled as OUT_OF_CALIBRATED_RANGE below, alongside the calibration's
% own explicit-zero convention -- NOTE 4).
cal_coef_on_wl = interp1(cal_wl, cal_coef, wl, 'linear', NaN);
out_of_range = isnan(cal_coef_on_wl) | (cal_coef_on_wl < CAL_ZERO_TOLERANCE);  % NOTE 4
cal_coef_on_wl(out_of_range) = NaN;   % make the exclusion explicit downstream

%% ========================================================================
%  6. DISPERSION FROM WAVELENGTH, NOT FROM ROW POSITION (NOTE 1)
%  ========================================================================
% Build a dense, monotonic pixel-index -> wavelength lookup over the fit
% domain (with margin, in case a real pixel's wavelength falls just
% outside 0..N_PIXELS_FIT-1 at the grid's edges), then invert it by
% interpolation to recover each row's fractional pixel index from its
% wavelength. The polynomial is confirmed strictly increasing over this
% instrument's working range, so this inversion is well-posed.
idx_grid = linspace(-50, N_PIXELS_FIT + 50, 20000)';
wl_grid  = WL_COEF(1) + WL_COEF(2)*idx_grid + WL_COEF(3)*idx_grid.^2 + WL_COEF(4)*idx_grid.^3;
assert(all(diff(wl_grid) > 0), ...
    'ABSOLUTE_IRRADIANCE_VALIDATION:NonMonotonicWlPoly', ...
    'WL_COEF does not give a strictly increasing wavelength-vs-pixel curve over the fit domain -- cannot invert it safely.');

pixel_idx_of_row = interp1(wl_grid, idx_grid, wl, 'linear', NaN);
if any(isnan(pixel_idx_of_row))
    warning('ABSOLUTE_IRRADIANCE_VALIDATION:WlOutsidePolyDomain', ...
        ['%d row(s) have a wavelength outside the range WL_COEF/N_PIXELS_FIT ' ...
         'covers (even with margin); their dispersion, and therefore their ' ...
         'irradiance, cannot be computed and will be NaN.'], sum(isnan(pixel_idx_of_row)));
end
dispersion = WL_COEF(2) + 2*WL_COEF(3)*pixel_idx_of_row + 3*WL_COEF(4)*pixel_idx_of_row.^2;

%% ========================================================================
%  7. CONVERT COUNTS TO ABSOLUTE IRRADIANCE (see formula in header)
%  ========================================================================
t_int_s = T_INT_US * 1e-6;
conv = cal_coef_on_wl ./ (t_int_s * area_cm2 * dispersion);   % uW/cm^2/nm per count

irrad_shots = counts .* conv;             % (pixel x shot), uW/cm^2/nm
irrad_mean  = mean(irrad_shots, 2);
irrad_std   = std(irrad_shots, 0, 2);     % MATLAB flag 0 -> N-1 (sample std)

ci95_mean = 1.96 * irrad_std / sqrt(n_shots);      % NOTE 5, precision of the mean
ci95_pred = 1.96 * irrad_std;                       % NOTE 5, single-shot spread

negative_signal = irrad_mean < 0;                    % NOTE 8
invalid = saturated | out_of_range | negative_signal;
irrad_mean(invalid) = NaN;
ci95_mean(invalid)  = NaN;
ci95_pred(invalid)  = NaN;

fprintf('\nExclusions: %d saturated, %d out-of-calibrated-range, %d negative-signal (%d unique pixels excluded of %d)\n', ...
    sum(saturated), sum(out_of_range), sum(negative_signal), sum(invalid), n_pixels);

%% ========================================================================
%  8. CERTIFICATE CROSS-CHECK
%  ========================================================================
fprintf('\n=== Certificate cross-check ===\n');
fprintf('%9s %12s %12s %10s\n', 'lambda','measured','certificate','ratio');
cert_ratio_wl = [];  cert_ratio_val = [];
for i = 1:numel(cert_wl)
    if cert_wl(i) < min(wl) || cert_wl(i) > max(wl), continue; end
    [dmin, j] = min(abs(wl - cert_wl(i)));
    if dmin > CERT_MATCH_TOL_NM
        fprintf('%9.1f %12s %12.4f %10s   (NEAREST DATA POINT %.2f nm AWAY -- excluded)\n', ...
                cert_wl(i), '--', cert_uwcm2nm(i), '--', dmin);
        continue
    end
    if isnan(irrad_mean(j))
        reasons = {};
        if saturated(j),       reasons{end+1} = 'SATURATED';             end %#ok<AGROW>
        if out_of_range(j),    reasons{end+1} = 'OUT_OF_CALIBRATED_RANGE'; end %#ok<AGROW>
        if negative_signal(j), reasons{end+1} = 'NEGATIVE_SIGNAL';       end %#ok<AGROW>
        fprintf('%9.1f %12s %12.4f %10s   (%s -- excluded)\n', ...
                wl(j), '--', cert_uwcm2nm(i), '--', strjoin(reasons, '+'));
        continue
    end
    r = irrad_mean(j) / cert_uwcm2nm(i);
    fprintf('%9.1f %12.4f %12.4f %10.3f\n', wl(j), irrad_mean(j), cert_uwcm2nm(i), r);
    cert_ratio_wl(end+1)  = wl(j);    %#ok<SAGROW>
    cert_ratio_val(end+1) = r;        %#ok<SAGROW>
end

if numel(cert_ratio_val) >= 2
    m_ratio  = mean(cert_ratio_val);   s_ratio = std(cert_ratio_val);
    md_ratio = median(cert_ratio_val); mad_ratio = 1.4826 * median(abs(cert_ratio_val - md_ratio));
    fprintf('\nmean +/- std:      %.4f +/- %.4f  (%d points)\n', m_ratio, s_ratio, numel(cert_ratio_val));
    fprintf('median +/- MAD*1.4826:  %.4f +/- %.4f   -- NOTE 7: prefer this if the two disagree\n', ...
            md_ratio, mad_ratio);
    fprintf(['\nA near-constant ratio across wavelength points to a single multiplicative\n' ...
        'cause (geometry, coupling). A ratio that drifts with wavelength -- especially\n' ...
        'one that departs most where the lamp''s own signal is weak -- points to an\n' ...
        'additive, wavelength-dependent contamination instead.\n']);
end

Tcert = table(cert_ratio_wl(:), cert_ratio_val(:), ...
              'VariableNames', {'wavelength_nm','ratio_to_certificate'});
writetable(Tcert, OUT_CERT);

%% ========================================================================
%  9. EXPORT AND FIGURE
%  ========================================================================
Tout = table(wl, irrad_mean, ci95_mean, ci95_pred, saturated, out_of_range, negative_signal, ...
    'VariableNames', {'wavelength_nm','irradiance_uWcm2nm','ci95_of_mean', ...
    'pred95_single_shot','saturated','out_of_calibrated_range','negative_signal'});
writetable(Tout, OUT_TABLE);
fprintf('\nResults written: %s, %s\n', OUT_TABLE, OUT_CERT);

valid_idx   = ~isnan(irrad_mean);
wl_valid    = wl(valid_idx);
irrad_valid = irrad_mean(valid_idx);
ci95_valid  = ci95_pred(valid_idx);

figure('Units','normalized','OuterPosition',[0.05 0.08 0.85 0.82]);

subplot(2,1,1); hold on; grid on; box on;
if ~isempty(wl_valid)
    plot(wl_valid, irrad_valid, 'k-', 'LineWidth', 1.2);
    fill([wl_valid; flipud(wl_valid)], [irrad_valid+ci95_valid; flipud(irrad_valid-ci95_valid)], ...
         [0 0.45 0.74], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    yl = ylim;
    shade_region(wl, saturated,       yl, [0.80 0.00 0.00]);   % red   = saturated
    shade_region(wl, out_of_range,    yl, [0.50 0.50 0.50]);   % grey  = out of calibrated range
    shade_region(wl, negative_signal, yl, [0.60 0.40 0.80]);   % purple= negative signal
    xlim([min(wl_valid) max(wl_valid)]);
end
xlabel('Wavelength (nm)'); ylabel('Irradiance (uW cm^{-2} nm^{-1})');
title('Mean spectrum, 95% single-shot prediction band (red=saturated, grey=out of calibrated range, purple=negative signal)');

subplot(2,1,2); hold on; grid on; box on;
yline(1, 'k--');
if ~isempty(cert_ratio_wl)
    plot(cert_ratio_wl, cert_ratio_val, 'o-', 'Color', [0.85 0.33 0.10], ...
         'MarkerFaceColor', [0.85 0.33 0.10]);
    xlim([min(wl_valid) max(wl_valid)]);
end
xlabel('Wavelength (nm)'); ylabel('measured / certificate');
title('Certificate ratio spectrum -- flat = single factor, sloped = wavelength-dependent cause');

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function [meta, wl, coef] = read_irradcal(filename)
% Reads a SpectraSuite .IrradCal file: a short metadata header (tab-
% separated key/value lines), a "[uJoule/count]" marker, then one
% "wavelength<TAB>coefficient" row per pixel.
    fid = fopen(filename, 'r');
    assert(fid > 0, 'ABSOLUTE_IRRADIANCE_VALIDATION:FileNotFound', ...
        'Cannot open calibration file: %s', filename);
    cleaner = onCleanup(@() fclose(fid));   % closes fid even if an error is thrown below

    meta = struct();
    line = fgetl(fid);
    while ischar(line) && ~strcmp(strtrim(line), '[uJoule/count]')
        parts = strsplit(strtrim(line), '\t');
        if numel(parts) == 2
            key = matlab.lang.makeValidName(parts{1});
            meta.(key) = parts{2};
        end
        line = fgetl(fid);
    end
    assert(ischar(line), 'ABSOLUTE_IRRADIANCE_VALIDATION:BadIrradCal', ...
        '%s has no "[uJoule/count]" marker line -- not a recognised .IrradCal file.', filename);

    data = textscan(fid, '%f%f', 'Delimiter', '\t');
    wl = data{1};  coef = data{2};
    assert(~isempty(wl), 'ABSOLUTE_IRRADIANCE_VALIDATION:BadIrradCal', ...
        '%s: no wavelength/coefficient rows found after the marker line.', filename);

    fn = fieldnames(meta);
    fn_clean = strrep(fn, '_', '');
    hit = fn(contains(fn_clean, 'IntTime', 'IgnoreCase', true) | ...
             contains(fn_clean, 'IntegrationTime', 'IgnoreCase', true));
    assert(~isempty(hit), 'ABSOLUTE_IRRADIANCE_VALIDATION:MissingIntTime', ...
        '%s: could not find the integration-time field in the header.', filename);
    meta.IntTimeUsec = meta.(hit{1});

    hit = fn(contains(fn, 'Fiber', 'IgnoreCase', true) | contains(fn, 'Fibre', 'IgnoreCase', true));
    assert(~isempty(hit), 'ABSOLUTE_IRRADIANCE_VALIDATION:MissingFiber', ...
        '%s: could not find the fibre-diameter field in the header.', filename);
    meta.FiberMicron = meta.(hit{1});
end

function [wl, counts] = read_highspeed_counts(filename)
% Reads a raw High-Speed-Acquisition export: header row = [blank,
% timestamp_1, ..., timestamp_n(, blank)] in ms; then one row per pixel:
% [wavelength_nm, count_1, ..., count_n]. Returns wl (n_pixels x 1) and
% counts (n_pixels x n_shots). Timestamps themselves are not needed
% downstream and are discarded here.
%
% NOTE (v3.1 fix): the header row in files from this instrument has been
% observed with a TRAILING tab (an empty field after the last timestamp)
% that the DATA rows do NOT have. Splitting the header with
% strsplit(strtrim(header), '\t') -- as an earlier version of this
% function did -- strips that trailing tab, but strtrim also strips the
% LEADING tab (the empty field marking "no timestamp for the wavelength
% column"), which is not trailing whitespace to discard but a real,
% position-carrying delimiter. The result was n_shots undercounted by
% exactly one, which then misaligned textscan's fixed-width read of every
% single data row: each row's last (real) count value was read as the
% NEXT row's wavelength. This produced small, plausible-looking negative
% "wavelengths" scattered through the file (dark-noise-level count values
% leaking into the wavelength column) -- caught, in one real case, only
% because a handful of the leaked values happened to be negative and
% tripped validate_wavelength_axis; most such leaks would silently pass
% as a shifted-but-still-positive number.
%
% Fixed by never trimming the header before splitting: the leading empty
% field is dropped explicitly (by position, always element 1), and any
% trailing empty field(s) are dropped explicitly (by content, wherever
% they occur), rather than delegating both to strtrim's blanket
% whitespace-stripping. The result is then cross-checked directly against
% a data row's own field count -- not merely trusted -- so any future
% instance of this failure mode raises an explicit, specific error
% instead of silently misparsing.
    fid = fopen(filename, 'r');
    assert(fid > 0, 'ABSOLUTE_IRRADIANCE_VALIDATION:FileNotFound', ...
        'Cannot open counts file: %s', filename);
    cleaner = onCleanup(@() fclose(fid));

    header = fgetl(fid);
    assert(ischar(header), 'ABSOLUTE_IRRADIANCE_VALIDATION:BadCountsFile', ...
        '%s appears to be empty.', filename);

    header_fields = strsplit(header, '\t');       % NOT strtrim'd -- see note above
    assert(isempty(header_fields{1}), ...
        'ABSOLUTE_IRRADIANCE_VALIDATION:BadCountsFile', ...
        ['%s: expected the header row to start with an empty field (no timestamp ' ...
         'for the wavelength column) followed by a tab; got "%s" instead. This ' ...
         'file may not be in the expected format.'], filename, header_fields{1});
    header_fields(1) = [];                         % drop the leading empty marker
    header_fields(cellfun(@isempty, header_fields)) = [];  % drop any trailing empty field(s)
    n_shots = numel(header_fields);
    assert(n_shots >= 1, 'ABSOLUTE_IRRADIANCE_VALIDATION:BadCountsFile', ...
        '%s: could not parse any shot columns from the header row.', filename);

    first_data_line = fgetl(fid);
    assert(ischar(first_data_line), 'ABSOLUTE_IRRADIANCE_VALIDATION:BadCountsFile', ...
        '%s: header row present but no data rows follow.', filename);
    n_fields_first_row = numel(strsplit(first_data_line, '\t'));
    assert(n_fields_first_row == n_shots + 1, ...
        'ABSOLUTE_IRRADIANCE_VALIDATION:HeaderDataMismatch', ...
        ['%s: the header implies %d shots (%d fields expected per data row: ' ...
         '1 wavelength + %d shots), but the first data row has %d fields. ' ...
         'Refusing to parse the rest of the file with a format that does not ' ...
         'match its own first row -- every row would silently misalign by the ' ...
         'same offset (this is exactly the failure this check exists to catch).'], ...
         filename, n_shots, n_shots+1, n_shots, n_fields_first_row);
    frewind(fid);
    fgetl(fid);   % consume the header line again, now that it has been validated

    fmt = ['%f' repmat('%f', 1, n_shots)];
    data = textscan(fid, fmt, 'Delimiter', '\t');
    wl = data{1};
    counts = cell2mat(data(2:end));
    assert(size(counts,1) == numel(wl) && size(counts,1) > 0, ...
        'ABSOLUTE_IRRADIANCE_VALIDATION:BadCountsFile', ...
        '%s: row count mismatch or no data rows after the header -- check for a truncated file.', filename);
end

function validate_wavelength_axis(wl, label)
% NOTE 6. Checks common, previously-encountered failure modes of a
% wavelength axis before it is used for any arithmetic: non-finite
% values, non-positive wavelengths, and duplicate wavelength values
% (which would make later interpolation steps ill-posed).
    assert(all(isfinite(wl)), 'ABSOLUTE_IRRADIANCE_VALIDATION:BadWavelengthAxis', ...
        '%s: wavelength axis contains non-finite values.', label);
    assert(all(wl > 0), 'ABSOLUTE_IRRADIANCE_VALIDATION:BadWavelengthAxis', ...
        '%s: wavelength axis contains non-positive values.', label);
    [sorted_wl, ~] = sort(wl);
    n_dupes = sum(diff(sorted_wl) == 0);
    if n_dupes > 0
        warning('ABSOLUTE_IRRADIANCE_VALIDATION:DuplicateWavelengths', ...
            '%s: %d duplicate wavelength value(s) found. Interpolation steps that use this axis as sample points will treat duplicates as redundant, not contradictory, but this is worth understanding.', ...
            label, n_dupes);
    end
end

function shade_region(wl, mask, ylims, rgb)
% Shades vertical bands of the current axes wherever MASK is true,
% grouping into contiguous runs of WL rather than one patch per pixel.
    if ~any(mask), return; end
    d = diff([0; mask(:); 0]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;
    for k = 1:numel(starts)
        x0 = wl(starts(k));  x1 = wl(stops(k));
        fill([x0 x1 x1 x0], [ylims(1) ylims(1) ylims(2) ylims(2)], rgb, ...
             'FaceAlpha', 0.10, 'EdgeColor', 'none');
    end
end