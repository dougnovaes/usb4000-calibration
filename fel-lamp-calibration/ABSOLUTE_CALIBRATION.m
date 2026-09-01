%% ABSOLUTE_CALIBRATION — FEL lamp spectral-irradiance interpolation: full
%  six-model comparison, cross-validation, model selection, and lamp
%  file generation
%
%  (Renamed from Calibration07 to avoid confusion with the wavelength-
%  axis calibration work, which lives in WAVELENGTH_CALIBRATION_Calibrated.m
%  and is an unrelated calibration of the same instrument.)
%
% =========================================================================
% WHAT CHANGED IN THIS VERSION
% =========================================================================
%   All changes below are robustness/correctness hardening found by a
%   line-by-line review; none of them is expected to change the selected
%   model (5, log-linear Huang) or the overall conclusions, and the
%   REFERENCE RESULTS table further down still describes the same run
%   that produced the currently-adopted lamp file. Re-run and diff
%   against that table if exact reproducibility matters to you.
%
%   1. T upper bound for Models 3-4 capped at the tungsten emissivity
%      table's actual range (4000 K, T_MAX_TABLE), instead of 6000 K.
%      eps_table() silently returns NaN outside its table (by design,
%      via interp2(...,'linear',NaN)), so any evaluation of the model at
%      T>4000 K during the optimizer's search -- not just at the final
%      answer -- was capable of poisoning lsqnonlin's residual with NaN
%      without any warning. The final answers here (T~2991 K, ~3024 K)
%      were always safely inside range, so this had not yet bitten in
%      practice, but the same 6000 K bound was also used, unguarded, in
%      all 31 per-fold LOOCV refits (fit_tab/fit_hyb), where a bad
%      initial guess is more likely. eps_table() itself now also clamps
%      T to [0, T_MAX_TABLE] as a second, independent line of defense.
%
%   2. fit_bb/fit_lin/fit_tab/fit_hyb (the four LOOCV per-fold refit
%      functions) previously each re-derived the Planck function from
%      scratch with the radiation constant C2 hardcoded as a literal
%      (14387.77), duplicated in four places, instead of using the
%      single `planck` handle (backed by the one C2 defined in Section 1)
%      that the main fits already use. They now receive `planck` as a
%      parameter, so there is exactly one definition of the Planck
%      function and one place C2 could ever need updating.
%
%   3. The GPR multi-start loop (Section 7) called fitrgp() with
%      identical default initial hyperparameters on every one of its
%      N_GPR_STARTS iterations. fitrgp's default initialization is
%      deterministic given the data, so this was, in effect, a single
%      fit run 20 times -- "20/20 starts succeeded" was literally true
%      but did not represent genuine multi-start robustness. Each
%      iteration now perturbs the initial length scale, signal std, and
%      noise std by an independent random lognormal factor (still fully
%      reproducible: seeded once via rng(1) before the loop), and the
%      printed success count now reflects how many of the N_GPR_STARTS
%      attempts actually converged, rather than assuming all of them did.
%
%   4. jacobian_fd (used only for Model 5's reported standard errors, not
%      for the fit itself) switched from forward finite differences,
%      O(step) truncation error, to centered differences, O(step^2) --
%      doubles the number of residual evaluations per parameter, which
%      is free here (31-point residual, not an expensive simulation).
%
%   5. Removed the unused `w` argument from fit_bb/fit_lin/fit_tab/
%      fit_hyb. The certificate's 1/F^2 weighting was never applied
%      through this parameter -- it was already baked into the relative
%      residual (model-F)/F used inside each function, since summing the
%      square of that residual is algebraically identical to a 1/F^2-
%      weighted sum of squares. `w` was accepted and silently ignored;
%      MATLAB's code analyzer correctly flagged this as an unused input
%      argument in four places. The (now dead) `w = 1./F.^2` computation
%      in Section 2 is removed along with it, replaced by a comment
%      explaining where that weighting actually lives.
%
%   6. add_result's `ffun` argument (the model's own curve-evaluation
%      handle) was accepted, stored nowhere, and never reused -- Section
%      10 separately re-derived every model's formula a second time to
%      plot it, a duplication that could drift out of sync with the
%      fitting formula. `ffun` is now stored in the `results` struct and
%      Section 10 calls it directly, so each model's formula exists in
%      exactly one place. This also resolves the "unused input argument"
%      flag MATLAB raised on `ffun`.
%
%   7. Model 6 (GPR)'s AIC/AICc used k=3 (its declared hyperparameter
%      count), but a near-interpolating Gaussian process -- exactly the
%      behavior this file's own "WHY LOOCV" note describes -- has a much
%      higher EFFECTIVE degree of freedom than that. A printed caveat
%      now accompanies the metrics table so the AIC/AICc columns are not
%      read as a fair comparison for GPR; LOOCV, not this table, is what
%      actually decides the winner (unchanged).
%
%   8. loocv_model5's `idxset = [1 2 3 4 6]` (which of Model 5's 7
%      parameters get refit per LOOCV fold) is now documented inline --
%      it was a bare magic index before.
%
%   9. The 0.1 unit-conversion factor in the lamp-file export (W/cm^3 ->
%      µW·cm⁻²·nm⁻¹) is now derived in a comment instead of left as a
%      bare literal.
%
%  10. fopen() for the lamp file now checks for a -1 (failure) return
%      before writing, instead of writing to it unconditionally.
%
%  11. loocv_gpr's per-fold fallback (a fold whose fixed-hyperparameter
%      refit fails, silently re-fit with free hyperparameters instead)
%      now emits an explicit warning and a summary count when triggered,
%      instead of silently changing methodology for that one fold.
%
%   Unchanged from the previous version, kept for history: Sections 1-9
%   and 11's *numerical methodology* is the same; Section 10 (figures)
%   already carried the six-series color/style/grayscale rework done in
%   the prior revision (unique color+linestyle+marker per series, a true
%   pixel-level grayscale export, models 1-4 and 5-6 split into separate
%   panels).
%
% =========================================================================
% SUMMARY
% =========================================================================
%   Certificate (31 points, 250-2400 nm) -> six candidate models of
%   increasing physical sophistication -> weighted nonlinear least
%   squares fit of each -> leave-one-out cross-validation of each (the
%   step that actually decides the winner) -> AIC/AICc/BIC and a nested
%   F-test as supporting evidence -> Model 5 selected -> lamp file for
%   SpectraSuite's Absolute Irradiance wizard.
%
% =========================================================================
% THE SIX MODELS
% =========================================================================
%   1. Black body                          F = A * B(lambda,T)
%   2. Black body x linear emissivity      F = A*(1-alpha*lambda_um)*B
%   3. Black body x tabulated emissivity   F = A*eps(lambda,T)*B
%   4. Model 3 x linear residual           F = A*eps(lambda,T)*(1-beta*lambda_um)*B
%   5. Two-stage log-linear (Huang 1998)   see below
%   6. Gaussian process regression         non-parametric, log-linear space
%
%   B(lambda,T) = 1 / ( lambda_um^5 * (exp(C2/(lambda_um*T)) - 1) ), the
%   Planck function core, with lambda_um = lambda/1000 (nm -> um).
%
%   UNIT NOTE, verified against this study's own established results
%   before writing this script: the linear correction terms
%   (1-alpha*lambda) in Models 2 and 4 use lambda in MICROMETRES, not
%   nanometres. Nanometres makes the bracket wildly negative over
%   250-2400 nm; micrometres reproduces the reference alpha=0.1283,
%   beta=-0.1258 exactly.
%
% =========================================================================
% TABLE III — WHERE THE REAL TUNGSTEN EMISSIVITY DATA CAME FROM
% =========================================================================
%   Models 3 and 4 need real, measured tungsten emissivity. A file filed
%   under a "De Vos 1954" name in this project's source material is
%   actually NASA TN D-1088 (Branstetter, 1961), which compiles De Vos
%   (1954) and Larrabee (1959) into its own Table III -- hemispherical
%   spectral emissivity of tungsten at 0, 2000, 4000 K, 0.2-20 um. The
%   T3_* arrays below were transcribed directly from that table. Fitting
%   Model 3 against them reproduces the previously established T=2991 K,
%   RMS=4.34% to four significant figures -- the cross-check that
%   confirms the transcription, not just its plausibility. NOTE: this
%   table's temperature range (0-4000 K) is now enforced as a hard bound
%   on T for Models 3-4, see item 1 above.
%
% =========================================================================
% MODEL 5: A KNOWN INSTABILITY
% =========================================================================
%   Only 12 of 31 certificate points lie below lambda0=450 nm, and one
%   (250 nm) carries a very large 1/F^2 weight. h3 controls curvature
%   below lambda0 and is weakly constrained: a fresh Stage-2 solve can
%   land on h3 pinned at its own bound, with h2, h4 correspondingly
%   different from another run's values, at an equally good weighted RMS.
%   This matches the standard errors already on record for this model:
%   h2, h3, h5, h6 are all "poorly determined" (64-327% relative SE),
%   while A and T -- the physically interpretable parameters -- are not
%   (6% and 1%). PRACTICAL CONSEQUENCE: this script's lamp file (Section
%   11) is generated from the pinned REFERENCE parameters, not from
%   whatever a fresh Stage 2 happens to converge to on a given run --
%   the predicted CURVE is far better constrained than any individual
%   h-parameter.
%
% =========================================================================
% WHY LOOCV, NOT JUST RMS, DECIDES THE WINNER
% =========================================================================
%   In-sample RMS always improves with more parameters. Model 6 (GPR) is
%   the demonstration: its in-sample RMS is the second-best of six, which
%   is why it was, for a time, the leading candidate -- until LOOCV
%   exposed near-interpolation (a very low fitted noise level) that
%   generalises poorly (RMS 0.6% -> LOOCV 3.5%, ~6x worse). Model 5's
%   LOOCV barely moves relative to its in-sample value despite having the
%   most nominal parameters of any candidate.
%
%   Model 5's curvature exponents (h4,h6) and Model 6's kernel
%   hyperparameters are held fixed at their full-data values during
%   LOOCV, refitting only the better-posed remaining parameters each
%   fold -- a documented, conservative simplification that can only add
%   variability to these two models, never hide an overfitting problem.
%
% =========================================================================
% REFERENCE RESULTS (for verification)
% =========================================================================
%   Model                    k   RMS(%)  LOOCV(%)     AIC     AICc   T(K)
%   1  Black body             2    6.40     6.77   -168.53  -168.10  3069
%   2  Linear emissivity      3    2.74     3.05   -220.19  -219.30  3022
%   3  Tabulated real         2    4.34     4.47   -192.56  -192.13  2991
%   4  Hybrid                 3    2.07     2.29   -237.68  -236.79  3024
%   5  Log-linear (Huang)     7  0.54-0.57  0.78-1.40  ~-315   ~-311  3023
%   6  GPR                    3  0.59-0.63    3.51    ~-312   ~-311  3069*
%   (Model 5's exact RMS/LOOCV vary run to run per the h3 instability
%   above; Model 6's "T" is a post-hoc black-body fit to its own smoothed
%   prediction, not a hyperparameter. Model 6's AIC/AICc use k=3 but are
%   not directly comparable to the parametric models -- see item 7 above.)
%
%   F-test, Model 1 vs 2 (nested, alpha=0): F=130.08, p=4.9e-12
%   F-test, Model 3 vs 4 (nested, beta=0):  F=100.02, p=9.6e-11
%
%   Selected model: 5, on every run to date. Lamp file: 250-2400 nm,
%   1 nm step, 2151 rows; F(900nm) ~ 210.4 W/cm^3 (certificate: 208.900).
%
% =========================================================================
% REFERENCES
% =========================================================================
%   The Eppley Laboratory, Inc. Certificate of Calibration of a Standard
%   of Spectral Irradiance, S.O. 52435, Lamp Serial No. EN-66. 1992.
%   De Vos, J.C. Physica, 20:669-714, 1954.
%   Larrabee, R.D. J. Opt. Soc. Am., 49(6):619-625, 1959.
%   Branstetter, J.R. NASA Technical Note D-1088, 1961. (Table III.)
%   Huang, L.K., Cebula, R.P., Hilsenrath, E. Metrologia, 35:381-386, 1998.
%   Burnham, K.P., Anderson, D.R. Model Selection and Multimodel
%     Inference, 2nd ed. Springer, 2002.
%   Rasmussen, C.E., Williams, C.K.I. Gaussian Processes for Machine
%     Learning. MIT Press, 2006.
%   Okabe, M., Ito, K. Color Universal Design (CUD). 2008. (figure palette)
%
% =========================================================================
% OUTPUT
% =========================================================================
%   lampfile_EN66_Model5.lmp        — as before
%   model_comparison.csv            — per-model metrics table
%   model_comparison_color.png/.pdf — publication figure, color
%   model_comparison_bw.png         — the same figure, true grayscale
%
% Requires the Statistics and Machine Learning Toolbox (fitrgp, Model 6
% only) and the Optimization Toolbox (lsqnonlin). Figure export uses
% `print -RGBImage`, available from R2020a.
%
% Author: Douglas Oliveira Novaes                              Licence: MIT
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURATION
%  ========================================================================
OUTPUT_LAMPFILE = 'lampfile_EN66_Model5.lmp';
OUTPUT_TABLE    = 'model_comparison.csv';
OUTPUT_FIG_BASE = 'model_comparison';
GRID_STEP_NM    = 1.0;  GRID_MIN_NM = 250;  GRID_MAX_NM = 2400;
LAMBDA0         = 450.0;
C2              = 14387.77;
N_GPR_STARTS    = 20;
GPR_SIGMA_MIN   = 1e-4;
FIG_DPI         = 300;

REF5 = struct('A',26896.62, 'T',3022.9479, 'h2',-0.1157, 'h3',80.0000, ...
              'h4',7.0723, 'h5',0.0119, 'h6',1.8880);

%% ========================================================================
%  2. DATA
%  ========================================================================
CERT = [
  250   0.136;   260   0.247;   270   0.411;   280   0.648;   290   0.985;
  300   1.450;   310   2.057;   320   2.839;   330   3.816;   340   5.052;
  350   6.537;   400  18.150;   450  37.220;   500  62.670;   555  93.900;
  600 119.900;   654.60 148.900; 700 169.100;  800 198.900;   900 208.900;
 1050 203.300;  1150 189.200;  1200 180.200;  1300 162.200;  1540 121.400;
 1600 112.100;  1700  97.800;  2000  66.000;  2100  58.600;  2300  45.000;
 2400  40.300];
wl = CERT(:,1);  F = CERT(:,2);  wl_um = wl/1000;  n = numel(wl);
% NOTE: the certificate's 1/F^2 weighting is NOT a separate array here --
% it is already baked into the relative residual (model-F)/F used by every
% r1..r5 below: sum(((model-F)/F).^2) = sum((model-F).^2 / F.^2), which IS
% a 1/F^2-weighted sum of squares. A `w = 1./F.^2` array used to be defined
% here and threaded into fit_bb/fit_lin/fit_tab/fit_hyb, but those
% functions never referenced it -- it was dead weight (pun intended);
% removed along with the unused parameter (see changelog item 5).

T3_LAM_UM = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.25 1.4 1.6 1.8 2.0 ...
             2.2 2.5 3.0 3.5 4.0 5.0 6.0 8.0 10.0 12.0 16.0 20.0]';
T3_TEMPS  = [0 2000 4000];
T3_EPS = [ ...
 0.4996 0.4632 0.4218;  0.4996 0.4632 0.4218;  0.4965 0.4628 0.4344;
 0.4856 0.4518 0.4221;  0.4731 0.4398 0.4100;  0.4639 0.4296 0.3916;
 0.4652 0.4153 0.3745;  0.4485 0.4000 0.3608;  0.4205 0.3824 0.3545;
 0.3934 0.3656 0.3470;  0.3404 0.3404 0.3404;  0.2995 0.3182 0.3361;
 0.2477 0.2910 0.3322;  0.2052 0.2684 0.3266;  0.1710 0.2492 0.3223;
 0.1439 0.2321 0.3174;  0.1169 0.2141 0.3087;  0.0882 0.1890 0.2850;
 0.0710 0.1730 0.2500;  0.0610 0.1630 0.2355;  0.0550 0.1470 0.2160;
 0.0475 0.1360 0.2025;  0.0425 0.1200 0.1790;  0.0365 0.1085 0.1615;
 0.0350 0.1000 0.1460;  0.0300 0.0875 0.1320;  0.0270 0.0770 0.1190];

% Upper validity limit of the tungsten emissivity table above. Models 3-4
% (and their LOOCV refits) must never let T's optimizer search wander past
% this, or eps_table() silently returns NaN (changelog item 1).
T_MAX_TABLE = max(T3_TEMPS);

fprintf('Certificate: %d points, %.0f-%.0f nm\n', n, min(wl), max(wl));

%% ========================================================================
%  3. FORWARD MODELS
%  ========================================================================
planck = @(x_um,T) 1.0 ./ (x_um.^5 .* (exp(C2./(x_um*T)) - 1));
% eps_table clamps T to [0, T_MAX_TABLE] as a second, independent line of
% defense on top of the optimizer bounds in Section 4/8 (changelog item 1):
% even if a bound were ever loosened by mistake, this clamp still prevents
% a silent NaN from an out-of-table T reaching the optimizer.
eps_table = @(lam_um,T) interp2(T3_TEMPS, T3_LAM_UM, T3_EPS, min(max(T,0),T_MAX_TABLE), lam_um, 'linear', NaN);
emiss5 = @(lam,h2,h3,h4,h5,h6) ...
    (lam <  LAMBDA0) .* exp(h2.*(lam-LAMBDA0)/500 - h3*abs((lam-LAMBDA0)/500).^h4) + ...
    (lam >= LAMBDA0) .* exp(h2.*(lam-LAMBDA0)/500 + h5*abs((lam-LAMBDA0)/500).^h6);
model5fun = @(lam,A,T,h2,h3,h4,h5,h6) A.*emiss5(lam,h2,h3,h4,h5,h6).*planck(lam/1000,T);
fitopts = optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',5000,'FunctionTolerance',1e-12);

%% ========================================================================
%  4. FIT MODELS 1-4
%  ========================================================================
results = struct([]);

r1 = @(p) (p(1).*planck(wl_um,p(2)) - F)./F;
p1 = lsqnonlin(r1, [50, 3000], [0,1000], [1e6,6000], fitopts);
results = add_result(results, 'Black body', 2, p1, r1, @(lam,p) p(1).*planck(lam/1000,p(2)), n);

r2 = @(p) (p(1).*(1-p(3).*wl_um).*planck(wl_um,p(2)) - F)./F;
p2 = lsqnonlin(r2, [50, 3000, 0.1], [0,1000,-2], [1e6,6000,2], fitopts);
results = add_result(results, 'Linear emissivity', 3, p2, r2, ...
    @(lam,p) p(1).*(1-p(3).*(lam/1000)).*planck(lam/1000,p(2)), n);

% T bounded to T_MAX_TABLE (4000 K), not 6000 K -- eps_table has no data
% above this, see changelog item 1.
r3 = @(p) (p(1).*eps_table(wl_um,p(2)).*planck(wl_um,p(2)) - F)./F;
p3 = lsqnonlin(r3, [50000, 3000], [0,1000], [1e7,T_MAX_TABLE], fitopts);
results = add_result(results, 'Tabulated real', 2, p3, r3, ...
    @(lam,p) p(1).*eps_table(lam/1000,p(2)).*planck(lam/1000,p(2)), n);

r4 = @(p) (p(1).*eps_table(wl_um,p(2)).*(1-p(3).*wl_um).*planck(wl_um,p(2)) - F)./F;
p4 = lsqnonlin(r4, [50000, 3000, -0.1], [0,1000,-2], [1e7,T_MAX_TABLE,2], fitopts);
results = add_result(results, 'Hybrid', 3, p4, r4, ...
    @(lam,p) p(1).*eps_table(lam/1000,p(2)).*(1-p(3).*(lam/1000)).*planck(lam/1000,p(2)), n);

fprintf('\nModels 1-4 fitted. Quick check against REFERENCE RESULTS:\n');
fprintf('  M1: T=%.0f (ref 3069)   M2: T=%.0f a=%.4f (ref 3022, 0.1283)\n', p1(2), p2(2), p2(3));
fprintf('  M3: T=%.0f (ref 2991)   M4: T=%.0f b=%.4f (ref 3024, -0.1258)\n', p3(2), p4(2), p4(3));

%% ========================================================================
%  5. F-TESTS FOR NESTED MODELS
%  ========================================================================
f_test(results(1), results(2), n, 'Model 1 vs 2 (alpha=0)');
f_test(results(3), results(4), n, 'Model 3 vs 4 (beta=0)');

%% ========================================================================
%  6. MODEL 5 — TWO-STAGE LOG-LINEAR FIT (Huang et al. 1998)
%  ========================================================================
L = log(wl.^5 .* F);
d = (wl - LAMBDA0)/500;  below = wl < LAMBDA0;
design5 = @(h4,h6) [ones(n,1), 1./wl, wl, below.*(-abs(d).^h4) + (~below).*(abs(d).^h6)];
sse5_1 = @(p) local_sse_stage1(p, design5, L);

best = struct('sse',Inf);
for h4g = 0.5:0.5:10
    for h6g = 0.5:0.25:4
        s = sse5_1([h4g,h6g]);
        if s < best.sse, best = struct('sse',s,'h4',h4g,'h6',h6g); end
    end
end
p1_5 = fminsearch(sse5_1, [best.h4, best.h6], optimset('Display','off'));
h4_1 = p1_5(1); h6_1 = p1_5(2);
X1 = design5(h4_1, h6_1);  coef1 = lsqminnorm(X1, L);

p0_5 = [exp(coef1(1)), REF5.T, coef1(3)*500, 50, h4_1, 0.05, h6_1];
lb5  = [0,    1000, -5,   0, 0.1, -5, 0.1];
ub5  = [1e6,  6000,  5, 200,  15,  5,  15];
r5 = @(p) (model5fun(wl,p(1),p(2),p(3),p(4),p(5),p(6),p(7)) - F)./F;
p5 = lsqnonlin(r5, p0_5, lb5, ub5, fitopts);
results = add_result(results, 'Log-linear (Huang)', 7, p5, r5, ...
    @(lam,p) model5fun(lam,p(1),p(2),p(3),p(4),p(5),p(6),p(7)), n);

rms5_fresh = 100*sqrt(sum(r5(p5).^2)/(n-7));
fprintf('\nM5 fresh fit: RMS=%.2f%% (ref 0.57%%)', rms5_fresh);
if p5(4) >= ub5(4)-1e-6 || p5(4) <= lb5(4)+1e-6
    fprintf('  [h3 at bound -- known instability, see header]');
end
fprintf('\n');

Jac = jacobian_fd(r5, p5);
sigma2 = sum(r5(p5).^2) / (n - 7);
try
    cov5 = sigma2 * pinv(Jac'*Jac);
    se5 = sqrt(diag(cov5));
    fprintf('Std. errors (this run):  A %.0f  T %.1f  h2 %.4f  h3 %.2f  h4 %.4f  h5 %.4f  h6 %.4f\n', se5);
catch
    fprintf('Std. error computation skipped (near-singular Jacobian this run).\n');
end

%% ========================================================================
%  7. MODEL 6 — GAUSSIAN PROCESS REGRESSION
%  ========================================================================
% NOTE ON MULTI-START (changelog item 3, CORRECTED after a scale bug):
% fitrgp's default initialization is deterministic given the data, so
% calling it N_GPR_STARTS times with no other variation previously
% produced N_GPR_STARTS *identical* fits.
%
% BUG FOUND AND FIXED HERE: with 'Standardize', true, any user-supplied
% initial 'KernelParameters'/'Sigma' must be given on the STANDARDIZED
% predictor scale (length scale ~ O(1)), not the raw data scale. An
% earlier version of this jitter used std(gp_x) (~600-700, the raw
% nm-scale std of wl) directly as the initial length scale -- 2-3 orders
% of magnitude too large in the space the optimizer actually works in.
% All N_GPR_STARTS then converged to the same badly-scaled, near-
% degenerate fit (LogLikelihood collapsed from ~0.3 to about -98;
% predictions diverged to hundreds-of-thousands of percent RMS). The
% per-fold LOOCV fallback (item 11) silently absorbed the damage in the
% LOOCV column (every fold's "fixed hyperparameter" refit failed and fell
% back to a free refit), so Model 5's selection was never at risk, but
% the full-data GPR curve/statistics in Sections 8 and 10 were wrong.
%
% FIX: fit ONE baseline model with fitrgp's OWN default initialization
% (no custom KernelParameters/Sigma) -- guaranteed correctly scaled,
% since it comes from fitrgp itself -- and jitter subsequent multi-start
% attempts AROUND that baseline instead of around a hand-derived guess.
gp_x = wl;  gp_y = L;
rng(1);

gpr_baseline = fitrgp(gp_x, gp_y, 'KernelFunction','matern52', 'Standardize', true, ...
    'SigmaLowerBound', GPR_SIGMA_MIN, 'FitMethod','exact','PredictMethod','exact');
ell_base  = gpr_baseline.KernelInformation.KernelParameters(1);
sigF_base = gpr_baseline.KernelInformation.KernelParameters(2);
sigN_base = gpr_baseline.Sigma;

best_gp = gpr_baseline; best_ll = gpr_baseline.LogLikelihood;
n_gpr_converged = 1;   % the baseline fit itself counts as the first attempt
for k = 2:N_GPR_STARTS
    ell0  = ell_base  * exp(0.3*randn);
    sigF0 = sigF_base * exp(0.3*randn);
    sigN0 = max(sigN_base * exp(0.3*randn), GPR_SIGMA_MIN);
    try
        gpr_k = fitrgp(gp_x, gp_y, 'KernelFunction','matern52', 'Standardize', true, ...
            'SigmaLowerBound', GPR_SIGMA_MIN, 'FitMethod','exact','PredictMethod','exact', ...
            'KernelParameters', [ell0; sigF0], 'Sigma', sigN0);
        n_gpr_converged = n_gpr_converged + 1;
        ll = gpr_k.LogLikelihood;
        if ll > best_ll, best_ll = ll; best_gp = gpr_k; end
    catch
        continue
    end
end
gprMdl = best_gp;

fprintf('\nM6 (GPR): LogLikelihood=%.2f, %d/%d starts succeeded\n', gprMdl.LogLikelihood, n_gpr_converged, N_GPR_STARTS);
[Lpred, ~] = predict(gprMdl, wl);
pred6 = exp(Lpred) ./ wl.^5;
resid6 = (pred6 - F)./F;
results = add_result_direct(results, 'GPR', 3, resid6, n);

r6T = @(p) (p(1).*planck(wl_um,p(2)) - pred6)./pred6;
p6T = lsqnonlin(r6T, [50,3000], [0,1000],[1e6,6000], fitopts);
fprintf('M6 effective colour temperature (post-hoc): %.0f K (ref 3069)\n', p6T(2));

%% ========================================================================
%  8. METRICS TABLE, AIC/AICc/BIC, LOOCV
%  ========================================================================
fprintf('\n=== METRICS (in-sample) ===\n');
fprintf('%-20s %3s %8s %10s %10s\n','model','k','RMS(%)','AIC','AICc');
for i = 1:numel(results)
    fprintf('%-20s %3d %8.2f %10.2f %10.2f\n', results(i).name, results(i).k, results(i).rms, results(i).AIC, results(i).AICc);
end
% Changelog item 7: GPR's k=3 undercounts its EFFECTIVE degrees of freedom
% once it is near-interpolating (see the "WHY LOOCV" note above) -- its
% AIC/AICc are shown for completeness but are not a fair basis for
% comparison against the parametric models 1-5.
fprintf(['NOTE: GPR''s AIC/AICc above use k=3 (its declared hyperparameter count), but a\n' ...
    'near-interpolating Gaussian process has a much higher EFFECTIVE degree of freedom\n' ...
    'than that; these two columns are not directly comparable to the parametric models''\n' ...
    'for GPR. LOOCV below, not this table, is what actually selects the model.\n']);

fprintf('\n=== LEAVE-ONE-OUT CROSS-VALIDATION ===\n');
loocv = zeros(1,6);
loocv(1) = loocv_simple(@(idx) fit_bb(wl_um(idx),F(idx),fitopts,planck), @(p,lam) p(1).*planck(lam/1000,p(2)), wl,F,n);
loocv(2) = loocv_simple(@(idx) fit_lin(wl_um(idx),F(idx),fitopts,planck), @(p,lam) p(1).*(1-p(3).*(lam/1000)).*planck(lam/1000,p(2)), wl,F,n);
loocv(3) = loocv_simple(@(idx) fit_tab(wl_um(idx),F(idx),fitopts,planck,eps_table,T_MAX_TABLE), @(p,lam) p(1).*eps_table(lam/1000,p(2)).*planck(lam/1000,p(2)), wl,F,n);
loocv(4) = loocv_simple(@(idx) fit_hyb(wl_um(idx),F(idx),fitopts,planck,eps_table,T_MAX_TABLE), @(p,lam) p(1).*eps_table(lam/1000,p(2)).*(1-p(3).*(lam/1000)).*planck(lam/1000,p(2)), wl,F,n);
loocv(5) = loocv_model5(wl,F,p5,h4_1,h6_1,model5fun,fitopts,lb5,ub5);
loocv(6) = loocv_gpr(gp_x,gp_y,wl,F,gprMdl,GPR_SIGMA_MIN);

for i = 1:6
    results(i).loocv = loocv(i);
    fprintf('%-20s in-sample %6.2f%%   LOOCV %6.2f%%   ratio %.2fx\n', results(i).name, results(i).rms, results(i).loocv, results(i).loocv/results(i).rms);
end

%% ========================================================================
%  9. CONSOLIDATED TABLE AND SELECTION
%  ========================================================================
[~, best_idx] = min(loocv);
fprintf('\n=== SELECTED MODEL: %s (lowest LOOCV, %.2f%%) ===\n', results(best_idx).name, loocv(best_idx));
if best_idx ~= 5
    warning('LOOCV selected a model other than the reference (Model 5). Re-examine before trusting the lamp file below.');
end

Tname = {results.name}';  Tk = [results.k]';  Trms = [results.rms]';
TAIC = [results.AIC]';  TAICc = [results.AICc]';  Tloocv = [results.loocv]';
Tout = table(Tname, Tk, Trms, TAIC, TAICc, Tloocv, 'VariableNames', {'model','k','RMS_pct','AIC','AICc','LOOCV_pct'});
writetable(Tout, OUTPUT_TABLE);
fprintf('Comparison table written: %s\n', OUTPUT_TABLE);

%% ========================================================================
%  10. FIGURES — publication-ready, color AND true grayscale
%  ========================================================================
% Palette: Okabe & Ito (2008), chosen for simultaneous colorblind- and
% print-safety. Every series gets a UNIQUE (color, line style, marker)
% triple so the color channel is redundant, never load-bearing.
STYLE = struct( ...
    'color',  {[0.00 0.00 0.00], [0.90 0.60 0.00], [0.35 0.70 0.90], ...
               [0.00 0.60 0.50], [0.80 0.40 0.00], [0.00 0.45 0.70]}, ...
    'line',   {'-','--',':','-.','-','--'}, ...
    'marker', {'none','s','^','d','o','v'}, ...
    'name',   {'1 Black body','2 Linear emissivity','3 Tabulated real', ...
               '4 Hybrid','5 Log-linear (Huang)','6 GPR'});
MARKER_SZ = 6;  LINE_W = 1.3;  LINE_W_WIN = 2.0;   % Model 5 drawn heavier: the adopted model

lam_plot = linspace(250,2400,600)';
% Changelog item 6: reuse each model's own ffun handle (stored in
% `results` by add_result) instead of re-deriving the same formula here a
% second time -- keeps the plotted curve and the fitted residual tied to
% a single definition per model. Model 6 (GPR) has no closed-form ffun
% (its "model" is the fitted GP object), so it keeps its own predict()
% call below.
curve1 = results(1).ffun(lam_plot, results(1).p);
curve2 = results(2).ffun(lam_plot, results(2).p);
curve3 = results(3).ffun(lam_plot, results(3).p);
curve4 = results(4).ffun(lam_plot, results(4).p);
curve5 = results(5).ffun(lam_plot, results(5).p);
[Lp,Lsdp] = predict(gprMdl, lam_plot);
curve6  = exp(Lp)./lam_plot.^5;
band_hi = exp(Lp+1.96*Lsdp)./lam_plot.^5;
band_lo = exp(Lp-1.96*Lsdp)./lam_plot.^5;

fig = figure('Units','normalized','OuterPosition',[0.03 0.05 0.9 0.85], 'Color','w');
set(fig, 'DefaultAxesFontSize',11, 'DefaultTextFontSize',11);

% --- top: all six curves + certificate, full width ---
ax1 = subplot(2,2,[1 2]); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
fill(ax1, [lam_plot;flipud(lam_plot)], [band_hi;flipud(band_lo)], STYLE(6).color, ...
     'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
plot(ax1, wl, F, 'ko', 'MarkerFaceColor','w', 'MarkerSize',7, 'LineWidth',1.3, 'DisplayName','Certificate');
curves_top = {curve1,curve2,curve3,curve4};
for m = 1:4
    plot(ax1, lam_plot, curves_top{m}, 'Color',STYLE(m).color, 'LineStyle',STYLE(m).line, ...
         'LineWidth',LINE_W, 'DisplayName',STYLE(m).name);
end
plot(ax1, lam_plot, curve5, 'Color',STYLE(5).color, 'LineStyle',STYLE(5).line, ...
     'LineWidth',LINE_W_WIN, 'DisplayName',STYLE(5).name);
plot(ax1, lam_plot, curve6, 'Color',STYLE(6).color, 'LineStyle',STYLE(6).line, ...
     'LineWidth',LINE_W, 'DisplayName',[STYLE(6).name ' (95% CI)']);
xlabel(ax1,'Wavelength (nm)'); ylabel(ax1,'Spectral irradiance (W cm^{-3})');
title(ax1,'Six candidate models vs. certificate (Eppley EN-66)');
legend(ax1,'Location','northeast','Box','off','FontSize',9);
xlim(ax1,[200 2450]);

% --- bottom-left: models 1-4 (the physical progression) ---
ax2 = subplot(2,2,3); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
yline(ax2,0,'k-','HandleVisibility','off');
resids14 = {100*r1(p1), 100*r2(p2), 100*r3(p3), 100*r4(p4)};
for m = 1:4
    plot(ax2, wl, resids14{m}, 'Color',STYLE(m).color, 'LineStyle',STYLE(m).line, ...
         'Marker',STYLE(m).marker, 'MarkerFaceColor',STYLE(m).color, 'MarkerSize',MARKER_SZ, ...
         'LineWidth',LINE_W, 'DisplayName',STYLE(m).name);
end
xlabel(ax2,'Wavelength (nm)'); ylabel(ax2,'Relative residual (%)');
title(ax2,'Models 1-4: the physical progression');
legend(ax2,'Location','best','Box','off','FontSize',8);

% --- bottom-right: models 5-6 (the final two candidates) ---
ax3 = subplot(2,2,4); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
yline(ax3,0,'k-','HandleVisibility','off');
plot(ax3, wl, 100*r5(p5), 'Color',STYLE(5).color, 'LineStyle',STYLE(5).line, ...
     'Marker',STYLE(5).marker, 'MarkerFaceColor',STYLE(5).color, 'MarkerSize',MARKER_SZ, ...
     'LineWidth',LINE_W_WIN, 'DisplayName',STYLE(5).name);
plot(ax3, wl, 100*resid6, 'Color',STYLE(6).color, 'LineStyle',STYLE(6).line, ...
     'Marker',STYLE(6).marker, 'MarkerFaceColor','none', 'MarkerSize',MARKER_SZ, ...
     'LineWidth',LINE_W, 'DisplayName',STYLE(6).name);
xlabel(ax3,'Wavelength (nm)'); ylabel(ax3,'Relative residual (%)');
title(ax3,'Models 5 vs. 6: the final two candidates');
legend(ax3,'Location','best','Box','off','FontSize',9);

sgtitle(fig, 'Model 5 (log-linear, Huang et al. 1998) is adopted — see LOOCV, Section 8', ...
        'FontSize',10, 'FontAngle','italic');

% ---- export: color (raster + vector) ----
exportgraphics(fig, [OUTPUT_FIG_BASE '_color.png'], 'Resolution',FIG_DPI);
exportgraphics(fig, [OUTPUT_FIG_BASE '_color.pdf'], 'ContentType','vector');
fprintf('\nColor figure written: %s_color.png / .pdf\n', OUTPUT_FIG_BASE);

% ---- export: true grayscale ----
% `print -RGBImage` renders the figure to an actual pixel array (R2020a+);
% converting that with the standard ITU-R BT.601 luminance weights is a
% real grayscale conversion of the artwork, not a color palette chosen to
% "probably print okay" -- it is what a print-only venue actually needs,
% and it is only legible here because every series above already carries
% a distinct line style and marker, not color alone.
img = print(fig, '-RGBImage', ['-r' num2str(FIG_DPI)]);
gray = uint8(0.299*double(img(:,:,1)) + 0.587*double(img(:,:,2)) + 0.114*double(img(:,:,3)));
imwrite(gray, [OUTPUT_FIG_BASE '_bw.png']);
fprintf('Grayscale figure written: %s_bw.png\n', OUTPUT_FIG_BASE);

%% ========================================================================
%  11. LAMP FILE (from the SELECTED / reference Model 5 parameters)
%  ========================================================================
lam_grid = (GRID_MIN_NM:GRID_STEP_NM:GRID_MAX_NM)';
F_wcm3   = model5fun(lam_grid, REF5.A, REF5.T, REF5.h2, REF5.h3, REF5.h4, REF5.h5, REF5.h6);
% Convert W/cm^3 (this script's internal unit, matching the certificate)
% to the µW·cm⁻²·nm⁻¹ unit SpectraSuite's Absolute Irradiance lamp-file
% format expects: 1 cm = 1e7 nm, so 1 W/cm^3 = 1 W/(cm^2 · cm)
% = (1e6 µW) / (cm^2 · 1e7 nm) = 0.1 µW·cm⁻²·nm⁻¹.
F_uwcm2nm = F_wcm3 * 0.1;

fid = fopen(OUTPUT_LAMPFILE, 'w');
if fid == -1
    error('Could not open %s for writing (check permissions, or whether it is open in another program).', OUTPUT_LAMPFILE);
end
for i = 1:numel(lam_grid)
    fprintf(fid, '%.1f\t%.6f\n', lam_grid(i), F_uwcm2nm(i));
end
fclose(fid);
peak900 = model5fun(900,REF5.A,REF5.T,REF5.h2,REF5.h3,REF5.h4,REF5.h5,REF5.h6);
fprintf('\nLamp file written: %s (%d rows). F(900nm)=%.3f W/cm^3 (cert. 208.900, %.2f%% off)\n', ...
        OUTPUT_LAMPFILE, numel(lam_grid), peak900, 100*(peak900-208.900)/208.900);

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function results = add_result(results, name, k, p, rfun, ffun, n)
    r = rfun(p);  sse = sum(r.^2);
    rms = 100*sqrt(sse/(n-k));
    AIC  = n*log(sse/n) + 2*k;
    AICc = AIC + 2*k*(k+1)/(n-k-1);
    % ffun is now stored (changelog item 6), not just accepted and dropped,
    % so Section 10 can reuse it instead of re-deriving each curve formula.
    entry = struct('name',name,'k',k,'p',p,'rms',rms,'AIC',AIC,'AICc',AICc,'loocv',NaN,'ffun',ffun);
    if isempty(results), results = entry; else, results(end+1) = entry; end
end

function results = add_result_direct(results, name, k, resid, n)
    sse = sum(resid.^2);
    rms = 100*sqrt(sse/n);
    AIC  = n*log(sse/n) + 2*k;
    AICc = AIC + 2*k*(k+1)/(n-k-1);
    % 'ffun' kept empty here (GPR has no closed-form curve function, see
    % Section 10) so the struct array's fields stay consistent with the
    % entries produced by add_result above.
    entry = struct('name',name,'k',k,'p',[],'rms',rms,'AIC',AIC,'AICc',AICc,'loocv',NaN,'ffun',[]);
    results(end+1) = entry;
end

function f_test(mA, mB, n, label)
    sseA = (mA.rms/100)^2 * (n-mA.k);
    sseB = (mB.rms/100)^2 * (n-mB.k);
    dfA = mA.k; dfB = mB.k;
    Fstat = ((sseA-sseB)/(dfB-dfA)) / (sseB/(n-dfB));
    p = 1 - fcdf(Fstat, dfB-dfA, n-dfB);
    fprintf('%s: F=%.2f, p=%.2e\n', label, Fstat, p);
end

function s = local_sse_stage1(p, design_matrix, L)
    h4=p(1); h6=p(2);
    X = design_matrix(h4,h6);
    coef = lsqminnorm(X, L);
    r = X*coef - L;
    s = sum(r.^2);
end

function J = jacobian_fd(rfun, p)
% Centered finite-difference Jacobian (changelog item 4): O(step^2)
% truncation error instead of the previous forward-difference O(step),
% for a tighter standard-error estimate on Model 5's parameters. Costs
% twice the residual evaluations per parameter versus forward differences,
% which is negligible here (rfun is a 31-point residual, not an expensive
% simulation).
    r0 = rfun(p);   % only used to size J; the difference formula below does not need it
    np = numel(p);  J = zeros(numel(r0), np);
    for i = 1:np
        step = max(1e-6, 1e-6*abs(p(i)));
        dp_plus  = p; dp_plus(i)  = dp_plus(i)  + step;
        dp_minus = p; dp_minus(i) = dp_minus(i) - step;
        J(:,i) = (rfun(dp_plus) - rfun(dp_minus)) / (2*step);
    end
end

function rmsq = loocv_simple(fit_fun, pred_fun, wl, F, n)
    err = zeros(n,1);
    for i = 1:n
        idx = true(n,1); idx(i) = false;
        p = fit_fun(idx);
        pred = pred_fun(p, wl(i));
        err(i) = (pred - F(i))/F(i);
    end
    rmsq = 100*sqrt(mean(err.^2));
end

% fit_bb/fit_lin/fit_tab/fit_hyb: the `w` (1/F^2 weight) argument these
% functions used to accept has been removed (changelog item 5) -- it was
% never referenced in their bodies, since the weighting is already
% implicit in the relative residual (model-F)/F each of them builds.
% They now also receive the shared `planck` handle instead of re-deriving
% the Planck function with a hardcoded C2 (changelog item 2).

function p = fit_bb(wl_um,F,opts,planck)
    r=@(p)(p(1).*planck(wl_um,p(2))-F)./F;
    p = lsqnonlin(r,[50,3000],[0,1000],[1e6,6000],opts);
end
function p = fit_lin(wl_um,F,opts,planck)
    r=@(p)(p(1).*(1-p(3).*wl_um).*planck(wl_um,p(2))-F)./F;
    p = lsqnonlin(r,[50,3000,0.1],[0,1000,-2],[1e6,6000,2],opts);
end
function p = fit_tab(wl_um,F,opts,planck,eps_table,T_max)
    r=@(p)(p(1).*eps_table(wl_um,p(2)).*planck(wl_um,p(2))-F)./F;
    p = lsqnonlin(r,[50000,3000],[0,1000],[1e7,T_max],opts);
end
function p = fit_hyb(wl_um,F,opts,planck,eps_table,T_max)
    r=@(p)(p(1).*eps_table(wl_um,p(2)).*(1-p(3).*wl_um).*planck(wl_um,p(2))-F)./F;
    p = lsqnonlin(r,[50000,3000,-0.1],[0,1000,-2],[1e7,T_max,2],opts);
end

function rmsq = loocv_model5(wl,F,p5,h4,h6,model5fun,opts,lb,ub)
    n = numel(wl);  err = zeros(n,1);
    % idxset (changelog item 8): which of Model 5's 7 parameters
    % [A,T,h2,h3,h4,h5,h6] get refit in each LOOCV fold. [1 2 3 4 6] selects
    % [A, T, h2, h3, h5] -- h4 and h6 (positions 5 and 7, the curvature
    % exponents) are held FIXED at their full-data values, per the header's
    % "WHY LOOCV" note: refitting only the better-posed remaining
    % parameters each fold is a documented, conservative simplification
    % that can only add variability to this model's LOOCV, never hide an
    % overfitting problem.
    idxset = [1 2 3 4 6];
    for i = 1:n
        keep = true(n,1); keep(i) = false;
        r = @(q) (model5fun(wl(keep),q(1),q(2),q(3),q(4),h4,q(5),h6) - F(keep))./F(keep);
        q0 = p5(idxset);
        q = lsqnonlin(r, q0, lb(idxset), ub(idxset), opts);
        pred = model5fun(wl(i), q(1),q(2),q(3),q(4),h4,q(5),h6);
        err(i) = (pred - F(i))/F(i);
    end
    rmsq = 100*sqrt(mean(err.^2));
end

function rmsq = loocv_gpr(gp_x, gp_y, wl, F, gprMdl, sigmaMin)
    n = numel(gp_x);  err = zeros(n,1);
    kp = gprMdl.KernelInformation.KernelParameters;
    sigma0 = gprMdl.Sigma;
    n_fallback = 0;
    for i = 1:n
        keep = true(n,1); keep(i) = false;
        try
            % FitMethod = 'none' (not Optimizer = 'none', which is not a
            % valid optimizer name and was throwing on every single call --
            % see changelog): 'none' tells fitrgp to use the supplied
            % KernelParameters/Sigma AS THE FINAL MODEL, with no
            % optimization step, which is what "fixed hyperparameters per
            % fold" actually requires.
            m = fitrgp(gp_x(keep), gp_y(keep), 'KernelFunction','matern52', 'Standardize', true, ...
                'SigmaLowerBound', sigmaMin, 'KernelParameters', kp, 'Sigma', sigma0, ...
                'FitMethod','none', 'PredictMethod','exact');
        catch
            % Fallback: this fold's fixed-hyperparameter evaluation failed
            % (e.g. a near-singular covariance for this particular n-1
            % subset), so re-fit with free hyperparameters for this fold
            % only -- a methodology change for one fold, flagged explicitly.
            warning('loocv_gpr: fold %d (holding out %.1f nm) failed with fixed hyperparameters; refitting with free hyperparameters for this fold only.', i, wl(i));
            n_fallback = n_fallback + 1;
            m = fitrgp(gp_x(keep), gp_y(keep), 'KernelFunction','matern52', 'Standardize', true, ...
                'SigmaLowerBound', sigmaMin);
        end
        Lp = predict(m, wl(i));
        pred = exp(Lp)/wl(i)^5;
        err(i) = (pred - F(i))/F(i);
    end
    if n_fallback > 0
        fprintf('loocv_gpr: %d of %d folds used the fallback (free-hyperparameter) refit.\n', n_fallback, n);
    end
    rmsq = 100*sqrt(mean(err.^2));
end