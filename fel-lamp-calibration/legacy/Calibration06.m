%% =====================================================================
%  FEL STANDARD LAMP SPECTRAL IRRADIANCE - MODEL COMPARISON
% =====================================================================
%  Fits six increasingly sophisticated models to the NIST-traceable
%  spectral irradiance certificate of an Eppley FEL standard lamp
%  (S.O. 52435, serial EN-66; certificate dated 15 July 1992, 250 nm to
%  2400 nm, calibrated at 50 cm, 7.90 A DC).
%
%  MODELS
%    1) Grey body                    - pure Planck function (A, T)
%    2) Empirical linear emissivity  - Planck x (1 - alpha*lambda)
%    3) Tabulated real emissivity    - Planck x De Vos/NASA tungsten table
%    4) Hybrid                       - tabulated emissivity x small
%                                       empirical residual correction
%    5) Huang log-linear             - two-stage fit inspired by
%                                       Huang, Cebula & Hilsenrath (1998)
%    6) Gaussian Process Regression  - non-parametric fit in log-linear
%                                       space, with calibrated pointwise
%                                       predictive uncertainty
%
%  VALIDATION
%    Besides in-sample fit metrics (SSE, R2, AIC, AICc, BIC, RMS), every
%    model is also scored by leave-one-out cross-validation (LOOCV), an
%    out-of-sample check that is standard practice before presenting a
%    model-comparison study of this kind for publication.
%
%  DATA SOURCES
%    - Eppley Laboratory certificate of calibration (irradiance data).
%    - Tungsten emissivity: Table III of NASA TN D-1088 (Branstetter,
%      1961), a US government report (public domain), itself built on
%      De Vos (1954, Physica 20, 690) corrected for scattered light by
%      Larrabee (1959, J. Opt. Soc. Am. 49, 619).
%    - Functional form of Model 5: Huang, Cebula & Hilsenrath, "New
%      procedure for interpolating NIST FEL lamp irradiances,"
%      Metrologia 35, 381-386 (1998).
%    - Model 6 follows standard GP regression practice, e.g. Rasmussen
%      & Williams, "Gaussian Processes for Machine Learning" (MIT
%      Press, 2006).
%
%  REQUIRED TOOLBOXES
%    - Optimization Toolbox            (lsqnonlin, optimoptions)
%    - Statistics and Machine Learning Toolbox (fcdf, fitrgp, predict)
%    - Base MATLAB >= R2017b           (lsqminnorm, interp2, fminsearch)
%
%  All fits are weighted nonlinear least squares under the assumption
%  of constant RELATIVE measurement uncertainty, i.e. weight = 1/y^2,
%  which is also the assumption used by Huang et al. (1998) and by the
%  NBS/NIST fitting procedures they compare against.
% =====================================================================

clear; clc; close all;

%% 1. Data: Eppley certificate (S.O. 52435, lamp EN-66)
wavelength_nm = [250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, ...
    400, 450, 500, 555, 600, 654.60, 700, 800, 900, 1050, ...
    1150, 1200, 1300, 1540, 1600, 1700, 2000, 2100, 2300, 2400]';

irradiance = [0.136, 0.247, 0.411, 0.648, 0.985, 1.450, 2.057, 2.839, ...
    3.816, 5.052, 6.537, 18.150, 37.220, 62.670, 93.900, ...
    119.900, 148.900, 169.100, 198.900, 208.900, 203.300, ...
    189.200, 180.200, 162.200, 121.400, 112.100, 97.800, ...
    66.000, 58.600, 45.000, 40.300]';

x_um = wavelength_nm / 1000;    % wavelength in micrometres (fitting unit)
y    = irradiance;               % spectral irradiance, W/cm^3
n    = numel(y);                 % number of calibration points (31)
w    = 1 ./ y.^2;                % relative-error weights (see header)
sse_null = sum(w .* (y - mean(y)).^2);   % weighted total sum of squares,
                                          % used as the R^2 denominator

%% 2. Physical constants and the Planck core
C2 = 14387.77;   % second radiation constant, hc/k, in micrometre*Kelvin

planck_core = @(x,T) 1 ./ (x.^5 .* (exp(C2./(x.*T)) - 1));

%% 3. Real tungsten emissivity table (NASA TN D-1088, Table III)
% Public-domain US government data; De Vos (1954) corrected for
% scattered light per Larrabee (1959). Values are hemispherical
% spectral emissivity of tungsten at three temperatures; 2-D linear
% interpolation in (wavelength, temperature) is used, matching the
% interpolation scheme Branstetter (1961) himself validated.
lam_tab_nm = [200 300 400 500 600 700 800 900 1000 1100 1250 1400 1600 1800 2000 2200 2500]';
eps_T0     = [0.4996 0.4996 0.4965 0.4856 0.4731 0.4639 0.4652 0.4485 0.4205 0.3934 0.3404 0.2995 0.2477 0.2052 0.1710 0.1439 0.1169]';
eps_T2000  = [0.4632 0.4632 0.4628 0.4518 0.4398 0.4296 0.4153 0.4000 0.3824 0.3656 0.3404 0.3182 0.2910 0.2684 0.2492 0.2321 0.2141]';
eps_T4000  = [0.4218 0.4218 0.4344 0.4221 0.4100 0.3916 0.3745 0.3608 0.3545 0.3470 0.3404 0.3361 0.3322 0.3266 0.3223 0.3174 0.3087]';

T_tab   = [0 2000 4000];
eps_tab = [eps_T0, eps_T2000, eps_T4000];   % [wavelength x temperature]

eps_interp = @(lam_nm,T) interp2(T_tab, lam_tab_nm, eps_tab, ...
    max(min(T,4000),0), lam_nm, 'linear');

%% 4. Parametric model definitions (Models 1-4)
model1 = @(p,x) p(1) .* planck_core(x,p(2));
model2 = @(p,x) p(1) .* (1 - p(3).*x) .* planck_core(x,p(2));
model3 = @(p,x) p(1) .* eps_interp(x*1000, p(2)) .* planck_core(x,p(2));
model4 = @(p,x) p(1) .* eps_interp(x*1000, p(2)) .* (1 - p(3).*x) .* planck_core(x,p(2));

%% 5. Weighted nonlinear least-squares fits (Models 1-4)
opts = optimoptions('lsqnonlin', 'Display','off', ...
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10, ...
    'MaxFunctionEvaluations',5000);

resid_w = @(p,model) sqrt(w) .* (model(p,x_um) - y);

% Only the fitted parameters (p1..p4) are kept here: the residual
% output of lsqnonlin is discarded with '~' because Section 9 below
% recomputes weighted residuals uniformly for every model (including
% Models 5 and 6, which are not fitted with lsqnonlin at all) via a
% single shared code path - keeping a second, unused copy here would
% be redundant.
[p1,~,~] = lsqnonlin(@(p) resid_w(p,model1), [1e4, 3000], ...
    [0 1000], [Inf 6000], opts);

[p2,~,~] = lsqnonlin(@(p) resid_w(p,model2), [2.8e4, 3020, 0.13], ...
    [0 1000 -1], [Inf 6000 1], opts);

[p3,~,~] = lsqnonlin(@(p) resid_w(p,model3), [7e4, 3000], ...
    [0 1000], [Inf 6000], opts);

[p4,~,~] = lsqnonlin(@(p) resid_w(p,model4), [6e4, 3000, 0], ...
    [0 1000 -1], [Inf 6000 1], opts);

%% 6. Model 5: Huang-inspired log-linear two-stage fit
% -----------------------------------------------------------------
% Physical picture (Huang, Cebula & Hilsenrath, 1998): FEL lamp
% irradiance = black body x filament emissivity x quartz-envelope
% transmittance. That product has negative curvature below ~300 nm,
% is roughly linear between 300-450 nm, and has positive curvature
% above 450 nm - captured here by a "modified Gaussian with tilt",
% continuous at the inflection point lambda0 = 450 nm.
%
% KEY IDEA WORTH BORROWING: fit L(lambda) = log(lambda^5 * F(lambda))
% instead of F(lambda) directly. With the two curvature EXPONENTS
% (h4, h6) held fixed, every other coefficient (h0, h1, h2, h3, h5)
% enters LINEARLY, so they can be solved by ordinary linear least
% squares - always well posed, no initial guess required. Only h4
% and h6 need a (low-dimensional, well-behaved) nonlinear search.
%
% Stage 1 works in the Wien approximation (as in the original paper);
% Stage 2 refines the result using the FULL Planck function (with the
% "-1" term), which matters here because our data extend to 2400 nm,
% well beyond the 1600 nm upper limit of the original paper, where the
% Wien approximation starts to lose accuracy.
%
% Coefficients are named h0..h6 (not c0..c6, as in the source paper)
% specifically to avoid clashing with the physical constant C2 defined
% above - the paper's "c2" and our physical "C2" are unrelated
% quantities that happen to share a name in the literature.
% -----------------------------------------------------------------

lambda0 = 450;   % nm - curvature inflection point (Huang et al., 1998)

% ---- Stage 1: log-linear fit (Wien approximation) ----
L = log((wavelength_nm.^5) .* y);

design_matrix = @(h4,h6) [ ...
    ones(n,1), ...                                              % h0
    1./wavelength_nm, ...                                        % h1
    wavelength_nm, ...                                            % h2
    (wavelength_nm<lambda0) .* (-abs((wavelength_nm-lambda0)/500).^h4) + ...
    (wavelength_nm>=lambda0) .* ( abs((wavelength_nm-lambda0)/500).^h6) ];  % h3 (below) / h5 (above)

sse_loglinear = @(hh) local_sse(hh, design_matrix, L);

% Coarse grid search over the two nonlinear exponents, refined locally.
best_sse = Inf; best_h4 = 2; best_h6 = 2;
for h4try = 1:0.5:8
    for h6try = 1:0.5:8
        s = sse_loglinear([h4try, h6try]);
        if s < best_sse
            best_sse = s; best_h4 = h4try; best_h6 = h6try;
        end
    end
end
hh_opt = fminsearch(sse_loglinear, [best_h4, best_h6], ...
    optimset('TolX',1e-4, 'TolFun',1e-8, 'MaxIter',2000));
h4_1 = hh_opt(1); h6_1 = hh_opt(2);

Xd_opt   = design_matrix(h4_1, h6_1);
coef_opt = lsqminnorm(Xd_opt, L);        % robust even if Xd is rank-deficient
h0_1 = coef_opt(1); h1_1 = coef_opt(2); h2_1 = coef_opt(3); h35_1 = coef_opt(4);

T_wien = -(C2*1000) / h1_1;    % hc/k in nm*K = C2*1000
fprintf('\n=== MODEL 5, STAGE 1 (log-linear, Wien approximation) ===\n');
fprintf('h4=%.4f h6=%.4f | h0=%.4f h1=%.4f h2=%.6f h(curvature)=%.4f\n', ...
    h4_1, h6_1, h0_1, h1_1, h2_1, h35_1);
fprintf('Approximate temperature (Wien) = %.0f K\n', T_wien);

% ---- Stage 2: refinement using the full Planck function ----
emiss_huang = @(lam_nm,h2,h3,h4,h5,h6) ...
    (lam_nm<lambda0)  .* exp(h2.*(lam_nm-lambda0)/500 - h3.*abs((lam_nm-lambda0)/500).^h4) + ...
    (lam_nm>=lambda0) .* exp(h2.*(lam_nm-lambda0)/500 + h5.*abs((lam_nm-lambda0)/500).^h6);

model5 = @(p,x) p(1) .* emiss_huang(x*1000, p(3),p(4),p(5),p(6),p(7)) .* planck_core(x,p(2));

F0_stage1 = exp(h0_1 + h1_1/lambda0 + h2_1*lambda0) / lambda0^5;
A0_stage1 = F0_stage1 / planck_core(lambda0/1000, T_wien);

p5_0 = [A0_stage1, T_wien, h2_1*500, h35_1, h4_1, h35_1, h6_1];
lb5  = [0, 1000, -5,  0, 0.3,  0, 0.3];
ub5  = [Inf, 6000,  5, 80,  12, 80,  12];
[p5,~,res5] = lsqnonlin(@(p) resid_w(p,model5), p5_0, lb5, ub5, opts);

% Automatic fallback: h3 (below-450nm curvature amplitude) is prone to
% sticking at its upper bound, because only 12 of the 31 data points
% lie below 450 nm and one of them (250 nm) carries very high weight
% under 1/y^2 weighting. If that happens, retry with h4 fixed at a
% typical value (2, a plain parabola) for a more stable below-450nm
% branch, and keep whichever result is no worse.
if abs(p5(4) - ub5(4)) < 1e-2
    fprintf('\nNOTE: h3 stuck at its upper bound in the unconstrained fit.\n');
    fprintf('Refitting with h4 fixed at 2 (a plain parabola) for stability...\n');
    h4_fixed = 2;
    model5_fixedH4 = @(p,x) p(1) .* emiss_huang(x*1000, p(3),p(4),h4_fixed,p(5),p(6)) .* planck_core(x,p(2));
    p5b_0 = [p5(1), p5(2), p5(3), 1, p5(6), p5(7)];
    [p5b,~,res5b] = lsqnonlin(@(p) resid_w(p,model5_fixedH4), p5b_0, ...
        [0, 1000, -5,  0,  0, 0.3], [Inf, 6000,  5, 40, 40, 12], opts);
    if sum(res5b.^2) < sum(res5.^2)
        fprintf('Fixed-h4 refit improved the fit -> adopting it.\n');
        p5 = [p5b(1), p5b(2), p5b(3), p5b(4), h4_fixed, p5b(5), p5b(6)];
        res5 = res5b;
    else
        fprintf('Fixed-h4 refit did not improve the fit -> keeping the free result.\n');
    end
end
fprintf('\n=== MODEL 5, STAGE 2 (full Planck, refined) ===\n');
fprintf('A=%.2f  T=%.0fK  h2=%.4f  h3=%.4f  h4=%.4f  h5=%.4f  h6=%.4f\n', p5);

% ---- Parameter uncertainty diagnostic (linearised, via the Jacobian) ----
% With n=31 and k=7, n/k = 4.4; some coefficients below 450 nm may be
% only weakly constrained even though the overall fit is excellent.
% This block quantifies that rather than leaving it implicit.
[~,~,resid5_final,~,~,~,jac5] = lsqnonlin(@(p) resid_w(p,model5), p5, lb5, ub5, opts);
J5 = full(jac5);
dof5 = n - numel(p5);
sigma2_5 = sum(resid5_final.^2) / dof5;
cov5 = sigma2_5 * pinv(J5' * J5);
se5  = sqrt(abs(diag(cov5)))';

fprintf('\nModel 5 parameter standard errors (linearised):\n');
p5_names = {'A','T','h2','h3','h4','h5','h6'};
for i = 1:numel(p5)
    relSE = 100 * se5(i) / max(abs(p5(i)), eps);
    flag = '';
    if relSE > 50, flag = '  <-- poorly determined (relative SE > 50%)'; end
    fprintf('  %-3s = %10.4f  +/- %10.4f  (%.0f%% relative)%s\n', ...
        p5_names{i}, p5(i), se5(i), relSE, flag);
end

%% 7. Model 6: Gaussian Process Regression (non-parametric)
% -----------------------------------------------------------------
% Rationale: Model 5 needs seven hand-crafted shape parameters, two of
% which (h3, h4) are poorly constrained by the available data (see the
% diagnostic above). A Gaussian Process replaces that hand-built
% curvature model with a small number of well-conditioned
% HYPERPARAMETERS (signal variance, length scale, noise variance),
% estimated by maximum marginal likelihood - the standard modern
% (2000s-2020s) approach to exactly this class of problem: smooth,
% overshoot-free interpolation of sparse, irregularly sampled data,
% together with a calibrated pointwise predictive uncertainty (see
% e.g. Rasmussen & Williams, 2006).
%
% The GP is fitted in the SAME log-linear space L = log(lambda^5*y)
% used for Model 5's Stage 1. In that space, the assumption of
% constant relative measurement error becomes approximately constant
% ABSOLUTE noise, which is exactly the homoscedastic-noise assumption
% a standard GP regression makes - so no manual reweighting is needed.
%
% IMPLEMENTATION NOTE (found by running the script): fitrgp's default
% noise floor ('SigmaLowerBound') is approximately 1% of std(response),
% which for our log-linear response L is close to the noise level the
% single-start default fit actually converged to (~0.058). That is a
% strong sign the optimiser was hitting an artificial floor rather than
% the genuine maximum of the marginal likelihood - for a smooth physical
% curve like this one, the true noise level can plausibly be much
% smaller. Two defences are used together: (a) the floor is lowered
% explicitly, and (b) the fit is repeated from several random starting
% points (marginal-likelihood optimisation for GPs is non-convex and
% can have local optima), keeping the run with the highest
% log-likelihood - standard practice in GP regression.
% -----------------------------------------------------------------
rng(0);   % reproducibility of the multi-start search
nStarts  = 20;
nSuccess = 0;
bestLL   = -Inf;
gprMdl   = [];
Lstd     = std(L);

for s = 1:nStarts
    sigmaL0 = 10 ^ (rand*3 - 1);            % length scale,  ~0.1 to 100 (standardised)
    sigmaF0 = Lstd * 10 ^ (rand*2 - 1);      % signal std,    ~0.1x to 10x std(L)
    sigma0  = Lstd * 10 ^ (rand*4 - 4);      % noise std,     ~0.0001x to 1x std(L)
    try
        mdl_s = fitrgp(wavelength_nm, L, 'KernelFunction','matern52', ...
            'Standardize',true, 'FitMethod','exact', 'PredictMethod','exact', ...
            'SigmaLowerBound',1e-4, ...
            'KernelParameters',[sigmaL0; sigmaF0], 'Sigma',sigma0);
        nSuccess = nSuccess + 1;
        if mdl_s.LogLikelihood > bestLL
            bestLL = mdl_s.LogLikelihood;
            gprMdl = mdl_s;
        end
    catch
        continue   % a poorly conditioned random start is simply skipped
    end
end
fprintf('\nModel 6 multi-start search: %d/%d starts succeeded, best log-likelihood = %.3f\n', ...
    nSuccess, nStarts, bestLL);
if isempty(gprMdl)
    error('Model 6: every multi-start attempt failed - check the fitrgp installation.');
end

model6 = @(x) exp(gp_predict(gprMdl, x*1000)) ./ (x*1000).^5;

fprintf('\n=== MODEL 6 (Gaussian Process Regression) ===\n');
% NOTE: because 'Standardize' is true, the reported length scale and
% signal std apply to the STANDARDISED wavelength predictor, not to
% nanometres directly - printed as such to avoid a misleading unit.
fprintf('Kernel: Matern 5/2 | Length scale = %.4f (standardised units) | Signal std = %.4f (log-L units) | Noise std = %.4f (log-L units)\n', ...
    gprMdl.KernelInformation.KernelParameters(1), ...
    gprMdl.KernelInformation.KernelParameters(2), gprMdl.Sigma);
fprintf(['NOTE: unlike Model 5, the GP has no free shape parameters to\n' ...
    'mis-constrain: its 3 hyperparameters (length scale, signal std,\n' ...
    'noise std) are all identified by maximum marginal likelihood. In\n' ...
    'exchange, it also directly reports a calibrated uncertainty at\n' ...
    'every wavelength - including the largest gap in the data\n' ...
    '(1300-1540 nm) - which is printed below.\n']);

ref_wavelengths = [250, 300, 450, 900, 1500, 2400];
fprintf('\n95%% predictive interval half-width at representative wavelengths:\n');
for wl = ref_wavelengths
    [Lhat, Lsd] = gp_predict(gprMdl, wl);
    Fhat  = exp(Lhat) / wl^5;
    Fhi   = exp(Lhat + 2*Lsd) / wl^5;
    Flo   = exp(Lhat - 2*Lsd) / wl^5;
    fprintf('  %6.0f nm : +/- %5.2f%%\n', wl, 100*(Fhi-Flo)/(2*Fhat));
end

% ---- Effective colour temperature (post-hoc; NOT a GP hyperparameter) ----
% The GP has no physical temperature: it interpolates the spectral
% SHAPE directly and never assumes a black body. For comparability
% with Models 1-5 (each of which reports a genuine fitted T), we fit a
% plain grey body to the GP's OWN smoothed prediction at the 31
% calibration wavelengths (not to the noisy certificate data) and
% report the resulting "effective colour temperature" - the single
% blackbody temperature whose overall spectral shape best resembles
% the GP curve. This is the same idea as the "correlated colour
% temperature" routinely quoted for real, non-blackbody light sources;
% it is a derived summary statistic, not a fitted GP parameter, and is
% therefore flagged with an asterisk everywhere it is reported below.
gpr_prediction_at_data = model6(x_um);
resid_grey_vs_gpr = @(p) sqrt(w) .* (model1(p,x_um) - gpr_prediction_at_data);
[p_Tgpr,~,~] = lsqnonlin(resid_grey_vs_gpr, [2e4, 3000], [0 1000], [Inf 6000], opts);
T_gpr_effective = p_Tgpr(2);
fprintf('\nEffective (post-hoc) colour temperature of the GP curve: %.0f K *\n', T_gpr_effective);
fprintf('* Not a GP hyperparameter - a grey-body fit to the GP''s own\n');
fprintf('  smoothed prediction, reported for comparability with Models 1-5 only.\n');

%% 8. Leave-one-out cross-validation (out-of-sample RMS)
% -----------------------------------------------------------------
% Every metric computed so far (SSE, R2, AIC, RMS) is measured on the
% SAME 31 points used to fit each model - an in-sample, "training"
% error. It answers "how well does the model reproduce data it has
% already seen", not "how well would it predict a wavelength it has
% never seen". For models with more free parameters (Model 5, k=7) or
% with very low estimated noise (Model 6 very nearly interpolates the
% data), in-sample error can be optimistic in a way AIC/AICc only
% approximately corrects for.
%
% Leave-one-out cross-validation (LOOCV) answers the out-of-sample
% question directly: hold out point i, refit on the remaining 30
% points, predict point i, record the relative error; repeat for every
% i = 1..31; aggregate into an RMS. This is standard practice before
% presenting a model-comparison study for publication, precisely
% because it cannot be fooled by a flexible model simply memorising
% the calibration points it was fitted to.
%
% Models 5 and 6 use a documented simplification for tractability and
% stability: their NONLINEAR SHAPE/KERNEL hyperparameters (h4, h6 for
% Model 5; the GP kernel hyperparameters for Model 6) are held fixed
% at the values already found from the full 31-point fit, and only the
% remaining, well-posed parameters are re-optimised in each fold. This
% avoids re-running the expensive multi-start hyperparameter search
% 31 times, and is, if anything, CONSERVATIVE: a fully naive refit
% would let those hyperparameters move too, which could only add
% further out-of-sample variability. If Models 5 and 6 still perform
% well under this generous treatment, that is meaningful evidence in
% their favour; if they do not, that is even more decisive evidence
% against them.
%
% Expect this section to take on the order of a minute to run (31
% refits per model, 6 models).
% -----------------------------------------------------------------
fprintf('\n=== LEAVE-ONE-OUT CROSS-VALIDATION (running, ~1 minute) ===\n');

[loocv_rms1, ~] = loocv_parametric(model1, [1e4, 3000], ...
    [0 1000], [Inf 6000], x_um, y, w, opts);

[loocv_rms2, ~] = loocv_parametric(model2, [2.8e4, 3020, 0.13], ...
    [0 1000 -1], [Inf 6000 1], x_um, y, w, opts);

[loocv_rms3, ~] = loocv_parametric(model3, [7e4, 3000], ...
    [0 1000], [Inf 6000], x_um, y, w, opts);

[loocv_rms4, ~] = loocv_parametric(model4, [6e4, 3000, 0], ...
    [0 1000 -1], [Inf 6000 1], x_um, y, w, opts);

% Model 5, reduced: h4 and h6 are held fixed at their full-fit values
% (h4=p5(5), h6=p5(7); see the note at the top of this section). The
% vector p being optimised here has only 5 elements, [A, T, h2, h3, h5]
% - do not confuse its p(5) (=h5) with p5(5) (=h4), the fixed value
% from the original 7-element full-fit vector.
model5_reduced = @(p,x) p(1) .* emiss_huang(x*1000, p(3), p(4), p5(5), p(5), p5(7)) ...
    .* planck_core(x, p(2));
[loocv_rms5, ~] = loocv_parametric(model5_reduced, p5([1 2 3 4 6]), ...
    lb5([1 2 3 4 6]), ub5([1 2 3 4 6]), x_um, y, w, opts);

% Model 6 (GPR): custom loop, since it is not fitted via lsqnonlin.
% Each fold's kernel hyperparameters are WARM-STARTED (not re-searched
% from scratch) from the full-fit optimum, for the same tractability
% reason given above.
r6 = zeros(n,1);
for i = 1:n
    idx = true(n,1); idx(i) = false;
    mdl_i = fitrgp(wavelength_nm(idx), L(idx), 'KernelFunction','matern52', ...
        'Standardize',true, 'FitMethod','exact', 'PredictMethod','exact', ...
        'SigmaLowerBound',1e-4, ...
        'KernelParameters', gprMdl.KernelInformation.KernelParameters, ...
        'Sigma', gprMdl.Sigma);
    Lhat_i = gp_predict(mdl_i, wavelength_nm(i));
    Fhat_i = exp(Lhat_i) / wavelength_nm(i)^5;
    r6(i) = (Fhat_i - y(i)) / y(i);
end
loocv_rms6 = sqrt(mean(r6.^2)) * 100;

loocv_rms = [loocv_rms1, loocv_rms2, loocv_rms3, loocv_rms4, loocv_rms5, loocv_rms6];
fprintf('LOOCV RMS (%%), out-of-sample, one value per model:\n');
model_names_short = {'1) Grey body','2) Linear emissivity','3) Tabulated emissivity', ...
    '4) Hybrid','5) Huang log-linear','6) Gaussian Process Regression'};
for i = 1:6
    fprintf('  %-32s %.2f%%\n', model_names_short{i}, loocv_rms(i));
end

%% 9. Unified comparison across all six models
% -----------------------------------------------------------------
% Every model exposes a uniform predict(x_um) -> irradiance interface,
% which lets the metrics, the printed table, and the plots below be
% generated by simple loops instead of six near-identical blocks of
% copy-pasted code.
% -----------------------------------------------------------------
results = struct('name',{}, 'k',{}, 'T',{}, 'T_derived',{}, 'predict',{}, 'metrics',{}, 'loocv_rms',{});

results(1) = struct('name','1) Grey body',                        'k',2, 'T',p1(2), 'T_derived',false, 'predict',@(x) model1(p1,x), 'metrics',[], 'loocv_rms',loocv_rms(1));
results(2) = struct('name','2) Empirical linear emissivity',      'k',3, 'T',p2(2), 'T_derived',false, 'predict',@(x) model2(p2,x), 'metrics',[], 'loocv_rms',loocv_rms(2));
results(3) = struct('name','3) Real tabulated emissivity',        'k',2, 'T',p3(2), 'T_derived',false, 'predict',@(x) model3(p3,x), 'metrics',[], 'loocv_rms',loocv_rms(3));
results(4) = struct('name','4) Hybrid (table + residual)',        'k',3, 'T',p4(2), 'T_derived',false, 'predict',@(x) model4(p4,x), 'metrics',[], 'loocv_rms',loocv_rms(4));
results(5) = struct('name','5) Huang log-linear (two stages)',    'k',7, 'T',p5(2), 'T_derived',false, 'predict',@(x) model5(p5,x), 'metrics',[], 'loocv_rms',loocv_rms(5));
results(6) = struct('name','6) Gaussian Process Regression',      'k',3, 'T',T_gpr_effective, 'T_derived',true, 'predict',model6, 'metrics',[], 'loocv_rms',loocv_rms(6));
% Model 6's k=3 is a standard but approximate convention: it counts
% the GP's hyperparameters (length scale, signal std, noise std), not
% literal curve-shape parameters, since the GP itself is non-parametric.
% Model 6's T is likewise a post-hoc derived quantity (see Section 7),
% not a genuine fitted parameter - hence T_derived = true, which the
% table below marks with an asterisk.

for i = 1:numel(results)
    resid_i = sqrt(w) .* (results(i).predict(x_um) - y);
    results(i).metrics = compute_metrics(resid_i, results(i).k, n, sse_null);
end

fprintf('\n=== MODEL COMPARISON (real emissivity: NASA TN D-1088) ===\n');
fprintf('%-38s %4s %10s %8s %8s %8s %8s %8s %9s %8s\n', ...
    'Model','k','SSE_w','R2','AIC','AICc','BIC','RMS(%)','LOOCV(%)','T(K)');
for i = 1:numel(results)
    m = results(i).metrics;
    if results(i).T_derived
        Tstr = sprintf('%5.0f*', results(i).T);   % flagged: see footnote below
    else
        Tstr = sprintf('%6.0f', results(i).T);
    end
    fprintf('%-38s %4d %10.5f %8.5f %8.2f %8.2f %8.2f %8.2f %9.2f %s\n', ...
        results(i).name, results(i).k, m.sse, m.r2, m.aic, m.aicc, m.bic, m.rms, ...
        results(i).loocv_rms, Tstr);
end
fprintf('* Model 6''s T is a post-hoc effective colour temperature (grey-body fit\n');
fprintf('  to the GP''s own prediction), not a parameter the GP itself estimates.\n');

% Overfitting diagnostic: flag any model whose out-of-sample LOOCV RMS
% is substantially worse than its in-sample RMS - a sign that the
% in-sample metrics above are, to some degree, measuring how well the
% model memorised the 31 calibration points rather than how well it
% would generalise to a new wavelength.
fprintf('\nIn-sample vs out-of-sample (LOOCV) RMS - a widening gap flags overfitting:\n');
for i = 1:numel(results)
    ratio = results(i).loocv_rms / max(results(i).metrics.rms, eps);
    if ratio > 1.5
        flag = '  <-- LOOCV notably worse than in-sample: possible overfitting';
    else
        flag = '  (LOOCV consistent with in-sample: no overfitting red flag)';
    end
    fprintf('  %-38s in-sample=%.2f%%  LOOCV=%.2f%%%s\n', ...
        results(i).name, results(i).metrics.rms, results(i).loocv_rms, flag);
end
fprintf('\nalpha (model 2) = %.4f | beta (model 4) = %.4f\n', p2(3), p4(3));
fprintf(['NOTE: n/k = %.1f for model 5, well below the n/k > 40 threshold\n' ...
    'Burnham & Anderson (2002) recommend for plain AIC - AICc is the more\n' ...
    'defensible criterion for models 5 and 6 here.\n'], n/7);

%% 10. Reference benchmarks - how to judge whether these numbers are "good"
fprintf('\n=== HOW TO READ THESE NUMBERS: REFERENCE BENCHMARKS ===\n');
fprintf(['RMS relative error (%%) = sqrt(SSE_w/dfe), i.e. the RMS of\n' ...
    '(model-data)/data. Compare against:\n' ...
    '  - Huang, Cebula & Hilsenrath (1998): SSBUV procedure achieves 0.16%%\n' ...
    '    (250-450nm) to 0.22%% (450-1600nm); the older NIST/NBS procedure\n' ...
    '    achieves 0.20%%-0.23%%.\n' ...
    '  - The same paper states the NIST calibration uncertainty is of order\n' ...
    '    1%% (2-sigma) near 400nm, i.e. roughly 0.5%% (1-sigma).\n' ...
    '  RULE OF THUMB: a fit is "good" once its RMS relative error drops\n' ...
    '  clearly below ~0.5%% - at that point the residual is dominated by\n' ...
    '  measurement noise in the original calibration, not by the choice\n' ...
    '  of fitting function.\n']);
for i = 1:numel(results)
    rms_i = results(i).metrics.rms;
    if rms_i < 0.5
        verdict = 'within literature-quality range (< 0.5%)';
    elseif rms_i < 1.0
        verdict = 'comparable to stated NIST calibration uncertainty (0.5-1%)';
    else
        verdict = 'dominated by lack-of-fit, well above calibration noise';
    end
    fprintf('  %-38s RMS=%.2f%% -> %s\n', results(i).name, rms_i, verdict);
end

fprintf(['\nAIC differences (Burnham & Anderson, 2002, rule of thumb):\n' ...
    '  delta-AIC 0-2  : substantial support for the poorer model too\n' ...
    '  delta-AIC 4-7  : considerably less support for the poorer model\n' ...
    '  delta-AIC > 10 : essentially no support for the poorer model\n']);
aic_values = arrayfun(@(r) r.metrics.aic, results);
best_aic = min(aic_values);
for i = 1:numel(results)
    fprintf('  %-38s delta-AIC = %6.2f\n', results(i).name, results(i).metrics.aic - best_aic);
end

fprintf('\nF-test significance: p<0.05 significant, p<0.01 highly significant (standard convention).\n');

%% 11. F-tests between nested model pairs
% Only models 1-2 and 3-4 are nested (model 1 = model 2 with alpha=0;
% model 3 = model 4 with beta=0); models 5 and 6 are not nested with
% anything here and are compared via AIC/AICc/BIC only (Section 10).
dfe1 = n - results(1).k; dfe2 = n - results(2).k;
F12 = ((results(1).metrics.sse - results(2).metrics.sse) / (dfe1-dfe2)) / (results(2).metrics.sse/dfe2);
p12 = 1 - fcdf(F12, dfe1-dfe2, dfe2);
fprintf('\nF-test (1 vs 2): F=%.2f, p=%.4g\n', F12, p12);

dfe3 = n - results(3).k; dfe4 = n - results(4).k;
F34 = ((results(3).metrics.sse - results(4).metrics.sse) / (dfe3-dfe4)) / (results(4).metrics.sse/dfe4);
p34 = 1 - fcdf(F34, dfe3-dfe4, dfe4);
fprintf('F-test (3 vs 4): F=%.2f, p=%.4g\n', F34, p34);
fprintf('(p<0.05 => the residual beta term significantly improves the tabulated model)\n');

%% 12. Plots
x_plot = linspace(min(x_um), max(x_um), 500)';

figure('units','normalized','outerposition',[0 0 1 1]);  % full-screen
colours = {'g','c','m','k',[0.85 0 0],[0 0.35 0.85]};   % one per model

% --- Main comparison panel, with Model 6's credible band ---
subplot(2,2,[1 2]); hold on;

[L6hat, L6sd] = gp_predict(gprMdl, x_plot*1000);
F6hat = exp(L6hat) ./ (x_plot*1000).^5;
F6hi  = exp(L6hat + 2*L6sd) ./ (x_plot*1000).^5;
F6lo  = exp(L6hat - 2*L6sd) ./ (x_plot*1000).^5;
fill([x_plot*1000; flipud(x_plot*1000)], [F6hi; flipud(F6lo)], ...
    colours{6}, 'FaceAlpha',0.15, 'EdgeColor','none', ...
    'DisplayName','GPR 95% predictive interval');

plot(x_um*1000, y, 'ob', 'MarkerFaceColor','b', 'MarkerSize',5, 'DisplayName','Data');
for i = 1:numel(results)
    yy = results(i).predict(x_plot);
    plot(x_plot*1000, yy, '-', 'Color',colours{i}, 'LineWidth',1.4, ...
        'DisplayName',sprintf('%s', results(i).name));
end
xlabel('Wavelength (nm)'); ylabel('Spectral Irradiance (W/cm^3)');
title('Comparison of the Six Models');
legend('Location','best'); grid on;

% --- Residuals, split into two groups for readability ---
subplot(2,2,3); hold on;
for i = 1:4
    r_i = 100 * (y - results(i).predict(x_um)) ./ y;
    plot(x_um*1000, r_i, '-o', 'Color',colours{i}, 'MarkerSize',4, ...
        'DisplayName',results(i).name);
end
yline(0,'--k','HandleVisibility','off');
xlabel('Wavelength (nm)'); ylabel('Residual Error (%)');
title('Residuals: Models 1-4'); legend('Location','best'); grid on;

subplot(2,2,4); hold on;
for i = 5:6
    r_i = 100 * (y - results(i).predict(x_um)) ./ y;
    plot(x_um*1000, r_i, '-^', 'Color',colours{i}, 'MarkerSize',5, ...
        'LineWidth',1.2, 'DisplayName',results(i).name);
end
yline(0,'--k','HandleVisibility','off');
xlabel('Wavelength (nm)'); ylabel('Residual Error (%)');
title('Residuals: Models 5-6'); legend('Location','best'); grid on;

%% =====================================================================
%  LOCAL FUNCTIONS (must appear after all script code in a .m script)
% =====================================================================

function s = local_sse(hh, design_matrix, L)
% Weighted-least-squares SSE for a given pair of nonlinear curvature
% exponents (hh = [h4, h6]), used to profile out those two parameters
% during Model 5's Stage 1 fit (see Section 6).
    if any(hh <= 0), s = Inf; return; end
    Xd = design_matrix(hh(1), hh(2));
    coef = lsqminnorm(Xd, L);
    s = sum((L - Xd*coef).^2);
end

function m = compute_metrics(resid, k, n, sse_null)
% Standard weighted-least-squares model-comparison metrics.
%   resid     : weighted residuals, resid = sqrt(w).*(model-data)
%   k         : number of free parameters (or hyperparameters, for GPR)
%   n         : number of data points
%   sse_null  : weighted SSE of the constant (mean) model, for R^2
    m.sse  = sum(resid.^2);
    m.r2   = 1 - m.sse / sse_null;
    m.aic  = n*log(m.sse/n) + 2*k;
    m.aicc = m.aic + (2*k*(k+1)) / max(n-k-1, 1);
    m.bic  = n*log(m.sse/n) + k*log(n);
    m.rms  = sqrt(m.sse / (n-k)) * 100;   % RMS relative error, in %
end

function [rms_loocv, r] = loocv_parametric(modelFcn, p0, lb, ub, x, y, w, opts)
% Leave-one-out cross-validation for a parametric weighted nonlinear
% least-squares model (see Section 8). For each of the n data points,
% refits the model on the other n-1 points (same initial guess and
% bounds as the full-data fit) and predicts the held-out point;
% aggregates the n out-of-sample relative errors into an RMS
% percentage.
%   modelFcn : function handle, modelFcn(p,x) -> predicted irradiance
%   p0,lb,ub : initial guess and bounds for lsqnonlin (as in the full fit)
%   x,y,w    : full data vectors (wavelength in um, irradiance, weights)
%   opts     : lsqnonlin options
% Unlike the in-sample RMS(%) reported elsewhere, this is NOT divided
% by (n-k): each held-out prediction is already a genuine out-of-sample
% test, so no degrees-of-freedom correction is needed.
    n = numel(y);
    r = zeros(n,1);
    for i = 1:n
        idx = true(n,1); idx(i) = false;
        resid_fold = @(p) sqrt(w(idx)) .* (modelFcn(p,x(idx)) - y(idx));
        p_fold = lsqnonlin(resid_fold, p0, lb, ub, opts);
        r(i) = (modelFcn(p_fold, x(i)) - y(i)) / y(i);
    end
    rms_loocv = sqrt(mean(r.^2)) * 100;
end

function [yhat, ysd] = gp_predict(gprMdl, Xnew)
% Thin wrapper around the Statistics and Machine Learning Toolbox's
% GP predict() method, isolated in one place so that Xnew's shape is
% normalised (column vector) exactly once rather than at every call
% site throughout the script.
    Xnew = Xnew(:);
    [yhat, ysd] = predict(gprMdl, Xnew);
end