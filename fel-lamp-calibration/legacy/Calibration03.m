%% Comparacao de Modelos: Corpo Cinza vs Emissividade Linear vs Emissividade Tabelada
% Ajuste da irradiancia espectral da lampada padrao FEL (Eppley, EN-66)
% Tres modelos nao-lineares sao comparados por AIC/BIC e teste F (quando aninhados)

clear; clc; close all;

%% 1. Dados (Certificado Eppley, S.O. 52435, lampada EN-66)
wavelength_nm = [250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, ...
    400, 450, 500, 555, 600, 654.60, 700, 800, 900, 1050, ...
    1150, 1200, 1300, 1540, 1600, 1700, 2000, 2100, 2300, 2400]';

irradiance = [0.136, 0.247, 0.411, 0.648, 0.985, 1.450, 2.057, 2.839, ...
    3.816, 5.052, 6.537, 18.150, 37.220, 62.670, 93.900, ...
    119.900, 148.900, 169.100, 198.900, 208.900, 203.300, ...
    189.200, 180.200, 162.200, 121.400, 112.100, 97.800, ...
    66.000, 58.600, 45.000, 40.300]';

x_um = wavelength_nm / 1000;   % micrometros
y    = irradiance;
n    = numel(y);
c2   = 14387.77;               % hc/k em um*K
w    = 1 ./ y.^2;              % mesmos pesos relativos usados no ajuste original

%% 2. Curva de emissividade espectral do tungstenio (valores aproximados)
% ATENCAO: estes sao valores ILUSTRATIVOS que reproduzem a FORMA geral
% (decaimento monotonico, ~0.45-0.48 no UV/visivel ate ~0.16-0.18 no IR
% proximo) reportada na literatura classica de tungstenio incandescente
% (De Vos 1954; Larrabee 1959), NAO uma copia da tabela original.
% Para um resultado definitivo/publicavel, substitua eps_tab pelos
% valores exatos da Tabela II de De Vos (1954, Physica 20, 690) ou de
% uma fonte primaria (NIST, Touloukian - Thermophysical Properties of
% Matter), idealmente na temperatura mais proxima de T~3000K.
lambda_tab_nm = [250 300 350 400 450 500 600 700 800 900 1000 ...
    1200 1400 1600 1800 2000 2200 2400 2600]';
eps_tab       = [0.480 0.470 0.462 0.455 0.448 0.440 0.423 0.405 0.385 0.365 ...
    0.345 0.305 0.270 0.240 0.215 0.195 0.180 0.170 0.163]';

eps_interp = @(lam_nm) interp1(lambda_tab_nm, eps_tab, lam_nm, 'pchip', 'extrap');

%% 3. Definicao dos tres modelos (todos como funcao de x em micrometros)
planck_core = @(x,T) 1 ./ (x.^5 .* (exp(c2./(x.*T)) - 1));

% Modelo 1: corpo cinza puro (2 parametros: A, T)
model1 = @(p,x) p(1) .* planck_core(x,p(2));

% Modelo 2: emissividade linear (3 parametros: A, T, alpha) -- o seu atual
model2 = @(p,x) p(1) .* (1 - p(3).*x) .* planck_core(x,p(2));

% Modelo 3: emissividade tabelada do tungstenio (2 parametros: A, T)
model3 = @(p,x) p(1) .* eps_interp(x*1000) .* planck_core(x,p(2));

%% 4. Ajuste ponderado (minimos quadrados nao-lineares) para os tres modelos
opts = optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10,'MaxFunctionEvaluations',5000);

resid_w = @(p,model) sqrt(w) .* (model(p,x_um) - y);

p1_0 = [1e4, 3000];
[p1,~,res1] = lsqnonlin(@(p) resid_w(p,model1), p1_0, [0 1000], [Inf 6000], opts);

p2_0 = [2.8e4, 3020, 0.13];   % proximo do resultado que voce ja obteve com fit()
[p2,~,res2] = lsqnonlin(@(p) resid_w(p,model2), p2_0, [0 1000 -1], [Inf 6000 1], opts);

p3_0 = [7e4, 3000];
[p3,~,res3] = lsqnonlin(@(p) resid_w(p,model3), p3_0, [0 1000], [Inf 6000], opts);

%% 5. Metricas de comparacao (SSE ponderado, R2, AIC, BIC)
compute_metrics = @(res,k) deal( ...
    sum(res.^2), ...
    1 - sum(res.^2)/sum(w.*(y-mean(y)).^2), ...
    n*log(sum(res.^2)/n) + 2*k, ...
    n*log(sum(res.^2)/n) + k*log(n));

[sse1,r2_1,aic1,bic1] = compute_metrics(res1,2);
[sse2,r2_2,aic2,bic2] = compute_metrics(res2,3);
[sse3,r2_3,aic3,bic3] = compute_metrics(res3,2);

fprintf('\n=== COMPARACAO DE MODELOS ===\n');
fprintf('%-32s %4s %10s %8s %8s %8s %8s\n','Modelo','k','SSE_w','R2','AIC','BIC','T(K)');
fprintf('%-32s %4d %10.4f %8.4f %8.2f %8.2f %8.0f\n','1) Corpo cinza (A,T)',        2,sse1,r2_1,aic1,bic1,p1(2));
fprintf('%-32s %4d %10.4f %8.4f %8.2f %8.2f %8.0f\n','2) Emiss. linear (A,T,alpha)',3,sse2,r2_2,aic2,bic2,p2(2));
fprintf('%-32s %4d %10.4f %8.4f %8.2f %8.2f %8.0f\n','3) Emiss. tabelada (A,T)',    2,sse3,r2_3,aic3,bic3,p3(2));
fprintf('\nalpha (modelo 2) = %.4f\n', p2(3));
fprintf('(Menor AIC/BIC = melhor compromisso ajuste-vs-complexidade)\n');

%% 6. Teste F entre modelos ANINHADOS (modelo 1 = modelo 2 com alpha=0)
dfe1 = n - 2; dfe2 = n - 3;
F_stat  = ((sse1 - sse2)/(dfe1-dfe2)) / (sse2/dfe2);
p_value = 1 - fcdf(F_stat, dfe1-dfe2, dfe2);
fprintf('\nTeste F (modelo 1 vs 2, aninhados): F = %.2f, p = %.4g\n', F_stat, p_value);
fprintf('(p < 0.05 => o termo alpha melhora o ajuste de forma estatisticamente significativa)\n');
fprintf('Nota: modelo 3 NAO e aninhado em relacao a 1 ou 2 -> comparar so via AIC/BIC.\n');

%% 7. Graficos: dados + 3 ajustes, e residuos percentuais dos 3 modelos
x_plot = linspace(min(x_um),max(x_um),500)';

figure('Position',[80 80 1300 500]);

subplot(1,2,1);
plot(x_um*1000, y, 'ob','MarkerFaceColor','b','MarkerSize',5); hold on;
plot(x_plot*1000, model1(p1,x_plot), '-g','LineWidth',1.3);
plot(x_plot*1000, model2(p2,x_plot), '-r','LineWidth',1.3);
plot(x_plot*1000, model3(p3,x_plot), '-m','LineWidth',1.5);
xlabel('Comprimento de Onda (nm)'); ylabel('Irradiancia Espectral (W/cm^3)');
legend('Dados', sprintf('Cinza (T=%.0fK)',p1(2)), sprintf('Linear (T=%.0fK)',p2(2)), ...
    sprintf('Tabelada (T=%.0fK)',p3(2)), 'Location','best');
title('Comparacao dos Tres Modelos'); grid on;

subplot(1,2,2);
r1pct = 100*(y - model1(p1,x_um))./y;
r2pct = 100*(y - model2(p2,x_um))./y;
r3pct = 100*(y - model3(p3,x_um))./y;
plot(x_um*1000, r1pct, '-og','MarkerSize',4); hold on;
plot(x_um*1000, r2pct, '-sr','MarkerSize',4);
plot(x_um*1000, r3pct, '-dm','MarkerSize',4,'LineWidth',1.3);
yline(0,'--k');
xlabel('Comprimento de Onda (nm)'); ylabel('Erro Residual (%)');
legend('Cinza','Linear (atual)','Tabelada','Location','best');
title('Residuos: Cinza vs Linear vs Tabelada'); grid on;