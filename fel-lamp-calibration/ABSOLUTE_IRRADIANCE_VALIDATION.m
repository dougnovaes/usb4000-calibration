%% ABSOLUTE_IRRADIANCE_VALIDATION -- multi-distance field validation
%  Converts SpectraSuite High-Speed (100-shot) raw-count acquisitions of the
%  Eppley FEL lamp EN-66 into absolute spectral irradiance with SpectraSuite
%  .IrradCal calibration files, and compares the result with the lamp
%  certificate (S.O. 52435, 50 cm) at 50, 57 and 60 cm.
%
% =========================================================================
% REVISION OF 2026-10-30 (coauthor review)
% =========================================================================
%   R1. read_highspeed_counts rewritten. The previous parser (a) took the
%       timestamp header line for a data row (its first token, 0.0, parses
%       as a number), creating a bogus pixel at 0 nm; (b) padded short rows
%       with zeros and then silently truncated the sorted table to 3648
%       rows, dropping real pixels at the top of the range. In addition,
%       example_HighSpeed100x_50cm.txt as published in release v1.2.1
%       carried 11 rows (repaired in v1.2.2)
%       split across two lines: a CR+LF was inserted at every 262144-byte
%       boundary of the file, a file-transfer/synchronisation artefact, not
%       instrument behaviour. The old parser read each fragment as a row of
%       its own (one of them as a fake pixel at 516.39 nm). The new parser
%       reads the header explicitly, rejoins split rows -- reporting each
%       one -- and stops on any row whose field count is still wrong.
%       Certificate-matched results are unaffected (no split row falls on a
%       certificate wavelength), but the per-pixel table was.
%   R2. Statistics are reported both over all matched certificate points
%       and over STAT_BAND (400-800 nm). Below about 300 nm the lamp signal
%       is two to three orders of magnitude under its visible peak and
%       close to the detector noise floor; the ratio there is not a
%       meaningful test of the calibration.
%   R3. Each run carries its distance; the inverse-square expectation
%       (50/d)^2 is printed next to the measured median and drawn in the
%       figure, and the 250-300 nm region is shaded.
%   R4. New Section 10: direct comparison of the two calibration files
%       (coefficient ratio 60 cm / 50 cm by band) and conversion of the
%       60 cm acquisition with the 50 cm calibration.
%   R5. validate_wavelength_axis also rejects non-positive wavelengths.
%
% =========================================================================
% WHICH CALIBRATION IS USED AT EACH DISTANCE
% =========================================================================
%   50 cm and 60 cm each have their own .IrradCal and their own 100-shot
%   acquisition. 57 cm has its own acquisition but no .IrradCal of its own;
%   it is converted with the 50 cm calibration. A calibration coefficient
%   (energy per detector count) describes the instrument and fibre, not the
%   lamp distance, so a valid calibration converts counts into the spectral
%   irradiance at the measurement plane. At a distance d the lamp delivers
%   the certificate value times (50/d)^2 -- the inverse-square law, as the
%   certificate's own instructions state (not to be used below ~45 cm).
%
%   The 50 cm acquisitions saturate the detector between about 533 and
%   680 nm even at the USB4000's minimum integration time (3.8 ms); those
%   pixels are flagged and excluded below.
% =========================================================================
clc;
DATA_DIR = 'example_data';
if ~isfolder(DATA_DIR)
    error('Data directory not found: %s. Run this script from its own directory.', DATA_DIR);
end

RUNS = struct( ...
    'label',       {'50 cm (certificate)', '57 cm', '60 cm'}, ...
    'dist_cm',     {50, 57, 60}, ...
    'cal_file',    {fullfile(DATA_DIR, 'example_calibration_50cm.IrradCal'), ...
                     fullfile(DATA_DIR, 'example_calibration_50cm.IrradCal'), ...   % reused, see header
                     fullfile(DATA_DIR, 'example_calibration_60cm.IrradCal')}, ...
    'counts_file', {fullfile(DATA_DIR, 'example_HighSpeed100x_50cm.txt'), ...
                     fullfile(DATA_DIR, 'example_HighSpeed100x_57cm.txt'), ...
                     fullfile(DATA_DIR, 'example_HighSpeed100x_60cm.txt')}, ...
    't_int_us',    {3800, 4049, 4555});  % confirmed from each acquisition's header

WL_COEF      = [176.3557, 0.2160665, -3.698535e-6, -5.538192e-10];
N_PIXELS_FIT = 3648;
SAT_THRESHOLD = 0.95;  ADC_FULL_SCALE = 65535;
CAL_ZERO_TOLERANCE = 1e-15;  CERT_MATCH_TOL_NM = 2.0;
STAT_BAND   = [400 800];   % R2: band used for the validation statistics (nm)
LOW_SN_BAND = [250 300];   % R3: shaded in the figure (nm)

OUT_TABLE_COMBINED = 'irradiance_results_all_distances.csv';
OUT_CERT_COMBINED  = 'certificate_comparison_all_distances.csv';
OUT_FIG_BASE       = 'distance_validation';
FIG_DPI = 300;

%% ========================================================================
%  2. CERTIFICATE DATA (Eppley S.O. 52435, lamp EN-66, 50 cm, 7.90 A DC)
%  ========================================================================
%  Column 2 in W/cm^3 as certified; x0.1 converts to uW cm^-2 nm^-1.
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
%  3-8. ONE CALL PER DISTANCE
%  ========================================================================
results_all = cell(numel(RUNS),1);
for k = 1:numel(RUNS)
    fprintf('\n=====================================================\n  %s\n=====================================================\n', RUNS(k).label);
    results_all{k} = process_one_distance(RUNS(k), WL_COEF, N_PIXELS_FIT, ...
        SAT_THRESHOLD, ADC_FULL_SCALE, CAL_ZERO_TOLERANCE, CERT_MATCH_TOL_NM, ...
        cert_wl, cert_uwcm2nm, STAT_BAND);
end

%% ========================================================================
%  9. COMBINED EXPORT AND COMPARISON FIGURE (Figure 3 of the manuscript)
%  ========================================================================
Tcert_all = table();
Tirr_all  = table();
for k = 1:numel(RUNS)
    r = results_all{k};
    if ~isempty(r.Tcert)
        Tcert_all = [Tcert_all; addvars(r.Tcert, repmat(string(RUNS(k).label), height(r.Tcert), 1), 'NewVariableNames', 'distance')]; %#ok<AGROW>
    end
    if ~isempty(r.Tirr)
        Tirr_all  = [Tirr_all;  addvars(r.Tirr,  repmat(string(RUNS(k).label), height(r.Tirr), 1),  'NewVariableNames', 'distance')]; %#ok<AGROW>
    end
end
writetable(Tcert_all, OUT_CERT_COMBINED);
writetable(Tirr_all,  OUT_TABLE_COMBINED);
fprintf('\nArquivos combinados salvos: %s, %s\n', OUT_TABLE_COMBINED, OUT_CERT_COMBINED);

STYLE = struct('color', {[0.00 0.45 0.70], [0.90 0.60 0.00], [0.00 0.60 0.50]}, ...
               'marker', {'o', 's', '^'}, 'line', {'-', '--', ':'});

fig = figure('Units', 'normalized', 'OuterPosition', [0.08 0.15 0.8 0.65], 'Color', 'w');
ax = axes('Parent', fig);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
set(ax, 'FontName', 'Helvetica', 'FontSize', 13, 'LineWidth', 1.0);

plotted_count = 0;
for k = 1:numel(RUNS)
    r = results_all{k};
    if ~isempty(r.Tcert) && height(r.Tcert) > 0
        plot(ax, r.Tcert.wavelength_nm, r.Tcert.ratio_to_certificate, ...
             'Color', STYLE(k).color, 'LineStyle', STYLE(k).line, 'Marker', STYLE(k).marker, ...
             'MarkerFaceColor', STYLE(k).color, 'MarkerSize', 8, 'LineWidth', 1.3, ...
             'DisplayName', RUNS(k).label);
        % R3: inverse-square expectation for this distance
        yline(ax, (50/RUNS(k).dist_cm)^2, ':', 'Color', STYLE(k).color, 'LineWidth', 1.0, ...
              'HandleVisibility', 'off');
        plotted_count = plotted_count + 1;
    end
end

xlabel(ax, 'Wavelength (nm)', 'FontSize', 13);
ylabel(ax, 'measured / certificate (50 cm)', 'FontSize', 13);
title(ax, 'Field validation: certificate ratio vs. working distance', 'FontSize', 13);

% R3: shade the low signal-to-noise region, behind the data
yl = ylim(ax);
hs = patch(ax, [LOW_SN_BAND(1) LOW_SN_BAND(2) LOW_SN_BAND(2) LOW_SN_BAND(1)], ...
           [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
           'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
uistack(hs, 'bottom');
ylim(ax, yl);

if plotted_count > 0
    legend(ax, 'Location', 'best', 'Box', 'off', 'FontSize', 13);
else
    warning('ABSOLUTE_IRRADIANCE_VALIDATION:EmptyPlot', ...
        'Nenhum dado válido foi encontrado para traçar o gráfico.');
end

exportgraphics(fig, [OUT_FIG_BASE '_color.png'], 'Resolution', FIG_DPI);
exportgraphics(fig, [OUT_FIG_BASE '_color.pdf'], 'ContentType', 'vector');
img = print(fig, '-RGBImage', ['-r' num2str(FIG_DPI)]);
gray = uint8(0.299*double(img(:,:,1)) + 0.587*double(img(:,:,2)) + 0.114*double(img(:,:,3)));
imwrite(gray, [OUT_FIG_BASE '_bw.png']);
fprintf('Figuras salvas com sucesso: %s_color.png/.pdf, %s_bw.png\n', OUT_FIG_BASE, OUT_FIG_BASE);

%% ========================================================================
%  10. CONSISTENCY BETWEEN CALIBRATION FILES (R4)
%  ========================================================================
%  A calibration coefficient (uJ/count) does not depend on the lamp
%  distance. Files built at 50 and 60 cm therefore agree (ratio ~1) if the
%  60 cm reference irradiance was the certificate scaled by (50/60)^2, and
%  differ by (60/50)^2 = 1.44 if it was the unscaled certificate.
[~, wl50, c50] = read_irradcal(RUNS(1).cal_file);
[~, wl60, c60] = read_irradcal(RUNS(3).cal_file);
wl50 = wl50(:);  c50 = c50(:);
c60_on_50 = interp1(wl60(:), c60(:), wl50, 'linear', NaN);
ok_coef = (c50 > CAL_ZERO_TOLERANCE) & (c60_on_50 > CAL_ZERO_TOLERANCE);
ratio_coef = c60_on_50 ./ c50;
BANDS = [250 350; 350 450; 450 546; 546 668; 668 800; 800 889];
fprintf('\n=== Calibration-file consistency: coefficient ratio 60 cm / 50 cm ===\n');
fprintf('    (1.00 if the 60 cm reference was scaled by (50/60)^2; %.2f if it was not)\n', (60/50)^2);
for b = 1:size(BANDS,1)
    s = ok_coef & wl50 >= BANDS(b,1) & wl50 < BANDS(b,2);
    if any(s)
        fprintf('  %4.0f-%4.0f nm: median %.4f   (16th-84th percentile %.3f-%.3f)\n', ...
            BANDS(b,1), BANDS(b,2), median(ratio_coef(s)), prctile(ratio_coef(s),16), prctile(ratio_coef(s),84));
    end
end

cross = RUNS(3);
cross.label    = '60 cm counts converted with the 50 cm calibration';
cross.cal_file = RUNS(1).cal_file;
fprintf('\n=====================================================\n  %s\n=====================================================\n', cross.label);
process_one_distance(cross, WL_COEF, N_PIXELS_FIT, SAT_THRESHOLD, ADC_FULL_SCALE, ...
    CAL_ZERO_TOLERANCE, CERT_MATCH_TOL_NM, cert_wl, cert_uwcm2nm, STAT_BAND);

%% ========================================================================
%  LOCAL FUNCTION -- one distance
%  ========================================================================
function r = process_one_distance(cfg, WL_COEF, N_PIXELS_FIT, SAT_THRESHOLD, ...
        ADC_FULL_SCALE, CAL_ZERO_TOLERANCE, CERT_MATCH_TOL_NM, cert_wl, cert_uwcm2nm, STAT_BAND)
    if isnan(cfg.t_int_us)
        error('ABSOLUTE_IRRADIANCE_VALIDATION:MissingTInt', ...
            'T_INT_US not defined for %s. Fill in the RUNS structure.', cfg.label);
    end

    [cal_meta, cal_wl, cal_coef] = read_irradcal(cfg.cal_file);
    cal_wl = cal_wl(:); cal_coef = cal_coef(:);
    validate_wavelength_axis(cal_wl, sprintf('arquivo de calibração (%s)', cfg.cal_file));
    assert(all(isfinite(cal_coef)) && all(cal_coef >= 0), ...
        'ABSOLUTE_IRRADIANCE_VALIDATION:BadCalCoef', ...
        '%s contains non-finite or negative calibration coefficients.', cfg.cal_file);

    fiber_um = cal_meta.FiberMicron;
    if ischar(fiber_um) || isstring(fiber_um)
        fiber_um = str2double(fiber_um);
    else
        fiber_um = double(fiber_um);
    end
    area_cm2 = pi * ((fiber_um / 2) * 1e-4)^2;

    [wl_raw, counts_raw, n_repaired] = read_highspeed_counts(cfg.counts_file);   % R1
    [wl, sort_idx] = sort(wl_raw(:));
    counts = counts_raw(sort_idx, :);
    validate_wavelength_axis(wl, sprintf('arquivo de contagens (%s)', cfg.counts_file));
    n_pixels = size(counts, 1);
    n_shots  = size(counts, 2);
    if n_pixels ~= N_PIXELS_FIT
        warning('ABSOLUTE_IRRADIANCE_VALIDATION:UnexpectedPixelCount', ...
            '%s (%s): %d linhas de dados, esperado %d.', cfg.label, cfg.counts_file, n_pixels, N_PIXELS_FIT);
    end
    if n_repaired > 0
        fprintf('%s: %d linha(s) partida(s) em duas foram reunidas (ver avisos acima).\n', cfg.label, n_repaired);
    end

    sat_mask  = counts > SAT_THRESHOLD * ADC_FULL_SCALE;
    saturated = any(sat_mask, 2);

    cal_coef_on_wl = interp1(cal_wl, cal_coef, wl, 'linear', NaN);
    out_of_range   = isnan(cal_coef_on_wl) | (cal_coef_on_wl < CAL_ZERO_TOLERANCE);
    cal_coef_on_wl(out_of_range) = NaN;

    idx_grid = linspace(-50, N_PIXELS_FIT + 50, 20000)';
    wl_grid  = WL_COEF(1) + WL_COEF(2)*idx_grid + WL_COEF(3)*idx_grid.^2 + WL_COEF(4)*idx_grid.^3;
    pixel_idx_of_row = interp1(wl_grid, idx_grid, wl, 'linear', NaN);
    dispersion = WL_COEF(2) + 2*WL_COEF(3)*pixel_idx_of_row + 3*WL_COEF(4)*pixel_idx_of_row.^2;

    t_int_s = cfg.t_int_us * 1e-6;
    conv = cal_coef_on_wl ./ (t_int_s * area_cm2 * dispersion);
    irrad_shots = counts .* conv;
    irrad_mean  = mean(irrad_shots, 2);
    irrad_std   = std(irrad_shots, 0, 2);
    ci95_mean = 1.96 * irrad_std / sqrt(n_shots);
    ci95_pred = 1.96 * irrad_std;
    negative_signal = irrad_mean < 0;
    invalid = saturated | out_of_range | negative_signal;
    irrad_mean(invalid) = NaN;
    ci95_mean(invalid)  = NaN;
    ci95_pred(invalid)  = NaN;

    fprintf('%s: %d saturados, %d fora da faixa, %d negativos (de %d pixels)\n', ...
        cfg.label, sum(saturated), sum(out_of_range), sum(negative_signal), n_pixels);
    if any(saturated)
        fprintf('  faixa saturada: %.2f-%.2f nm\n', min(wl(saturated)), max(wl(saturated)));
    end

    cert_ratio_wl  = [];   % wavelength of the matched PIXEL (used for plotting)
    cert_ratio_cwl = [];   % CERTIFICATE wavelength (used for band selection, R2)
    cert_ratio_val = [];
    for i = 1:numel(cert_wl)
        if cert_wl(i) < min(wl) || cert_wl(i) > max(wl), continue; end
        [dmin, j] = min(abs(wl - cert_wl(i)));
        if dmin > CERT_MATCH_TOL_NM || isnan(irrad_mean(j)), continue; end
        cert_ratio_wl(end+1)  = wl(j); %#ok<AGROW>
        cert_ratio_cwl(end+1) = cert_wl(i); %#ok<AGROW>
        cert_ratio_val(end+1) = irrad_mean(j) / cert_uwcm2nm(i); %#ok<AGROW>
    end

    r.Tcert = table(cert_ratio_wl(:), cert_ratio_cwl(:), cert_ratio_val(:), ...
        'VariableNames', {'wavelength_nm','certificate_wavelength_nm','ratio_to_certificate'});
    r.Tirr  = table(wl, irrad_mean, ci95_mean, ci95_pred, saturated, out_of_range, negative_signal, ...
        'VariableNames', {'wavelength_nm','irradiance_uWcm2nm','ci95_of_mean', ...
        'pred95_single_shot','saturated','out_of_calibrated_range','negative_signal'});

    expected = (50/cfg.dist_cm)^2;   % R3
    r.expected = expected;
    if numel(cert_ratio_val) >= 2
        r.median_ratio = median(cert_ratio_val);
        r.mad_ratio    = 1.4826 * median(abs(cert_ratio_val - r.median_ratio));
        fprintf('  todos os pontos:   mediana %.4f +/- %.4f (n=%d)\n', ...
            r.median_ratio, r.mad_ratio, numel(cert_ratio_val));
    end
    % R2: select by CERTIFICATE wavelength. Selecting by the matched pixel
    % wavelength dropped the 400 nm point (nearest pixel 399.95 nm).
    in_band = cert_ratio_cwl >= STAT_BAND(1) & cert_ratio_cwl <= STAT_BAND(2);
    if sum(in_band) >= 2
        vb = cert_ratio_val(in_band);
        r.median_band = median(vb);
        r.mad_band    = 1.4826 * median(abs(vb - r.median_band));
        r.n_band      = numel(vb);
        fprintf('  %d-%d nm:       mediana %.4f +/- %.4f (n=%d) | esperado (50/d)^2 = %.4f | desvio %+.1f%%\n', ...
            STAT_BAND(1), STAT_BAND(2), r.median_band, r.mad_band, r.n_band, expected, ...
            100*(r.median_band/expected - 1));
    else
        fprintf('  %d-%d nm: pontos insuficientes (n=%d)\n', STAT_BAND(1), STAT_BAND(2), sum(in_band));
    end
end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================
function [meta, wl, coef] = read_irradcal(filename)
    %READ_IRRADCAL  Read an Ocean Optics .IrradCal calibration file
    fid = fopen(filename, 'r');
    if fid < 0
        error('Cannot open file: %s', filename);
    end
    meta = struct();
    while true
        line = fgetl(fid);
        if ischar(line) && ~isempty(strtrim(line))
            parts = strsplit(line, char(9));
            if numel(parts) >= 2
                key = strtrim(parts{1});
                val = strtrim(parts{2});
                if strcmpi(key, 'Fiber (micron)')
                    field_name = 'FiberMicron';
                elseif strcmpi(key, 'Int. Time (usec)')
                    field_name = 'IntTime';
                elseif strcmpi(key, 'Lamp SN')
                    field_name = 'LampSN';
                else
                    field_name = regexprep(key, '[\s\(\)\.]', '');
                end
                num_val = str2double(val);
                if ~isnan(num_val)
                    meta.(field_name) = num_val;
                else
                    meta.(field_name) = val;
                end
            end
        else
            break;
        end
    end
    frewind(fid);
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        if contains(line, '['), break; end
    end
    wl   = [];
    coef = [];
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);
        if isempty(line) || startsWith(line, '%') || startsWith(line, '[')
            continue;
        end
        parts = strsplit(line);
        if numel(parts) >= 2
            w = str2double(parts{1});
            c = str2double(parts{2});
            if ~(isnan(w) || isnan(c))
                wl(end+1)   = w; %#ok<AGROW>
                coef(end+1) = c; %#ok<AGROW>
            end
        end
    end
    fclose(fid);
end

function [wl, counts, n_repaired] = read_highspeed_counts(filename)
    %READ_HIGHSPEED_COUNTS  Read a SpectraSuite High-Speed Acquisition export (R1)
    %  Line 1: one relative timestamp (ms) per shot, optionally preceded by
    %  an empty field. Following lines: wavelength_nm, then one count per
    %  shot. A row split across two lines (see header, R1) is rejoined and
    %  reported; any row whose field count is still wrong stops the script.
    txt = fileread(filename);
    txt = strrep(txt, sprintf('\r'), '');
    lines = strsplit(txt, sprintf('\n'), 'CollapseDelimiters', false);
    while ~isempty(lines) && isempty(strtrim(lines{end}))
        lines(end) = [];
    end
    assert(numel(lines) >= 2, 'ABSOLUTE_IRRADIANCE_VALIDATION:BadCountsFile', ...
        '%s: no data rows.', filename);

    n_shots = numel(sscanf(lines{1}, '%f'));   % header = timestamps only
    assert(n_shots >= 1, 'ABSOLUTE_IRRADIANCE_VALIDATION:BadCountsFile', ...
        '%s: could not read the timestamp header.', filename);

    n_lines = numel(lines);
    data = nan(n_lines - 1, n_shots + 1);
    n_rows = 0;  n_repaired = 0;  i = 2;
    while i <= n_lines
        vals = sscanf(lines{i}, '%f');
        if numel(vals) < n_shots + 1 && i < n_lines
            joined = sscanf([lines{i} lines{i+1}], '%f');
            if numel(joined) == n_shots + 1
                warning('ABSOLUTE_IRRADIANCE_VALIDATION:SplitRow', ...
                    '%s: row split across lines %d-%d rejoined (%.2f nm).', filename, i, i+1, joined(1));
                vals = joined;  n_repaired = n_repaired + 1;  i = i + 1;
            end
        end
        if numel(vals) ~= n_shots + 1
            error('ABSOLUTE_IRRADIANCE_VALIDATION:BadCountsFile', ...
                '%s, line %d: %d numeric fields, expected %d (1 wavelength + %d shots).', ...
                filename, i, numel(vals), n_shots + 1, n_shots);
        end
        n_rows = n_rows + 1;
        data(n_rows, :) = vals(:)';
        i = i + 1;
    end
    data   = data(1:n_rows, :);
    wl     = data(:, 1);
    counts = data(:, 2:end);
end

function validate_wavelength_axis(wl, label)
    %VALIDATE_WAVELENGTH_AXIS  Finite, positive, strictly increasing (R5)
    if any(~isfinite(wl))
        error('Non-finite wavelengths found in %s', label);
    end
    if any(wl <= 0)
        error('Non-positive wavelengths found in %s', label);
    end
    if ~all(diff(wl) > 0)
        error('Wavelength axis not strictly increasing in %s', label);
    end
end
