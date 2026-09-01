%% Comparação de Modelos: Corpo Cinza vs Linear vs Tabelada (real) vs Hibrida
% Ajuste da irradiancia espectral da lampada padrao FEL (Eppley, EN-66)
% Emissividade do tungstenio: dados REAIS extraidos da Tabela III do
% NASA TN D-1088 (Branstetter, 1961) - relatorio do governo dos EUA,
% dominio publico - que compila e corrige os dados de De Vos (1954,
% Physica 20, 690) com a correção de luz espalhada de Larrabee (1959,
% JOSA 49, 619).

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
w    = 1 ./ y.^2;

%% 2. Emissividade hemisferica do tungstenio: dados reais (NASA TN D-1088, Tabela III)
% Fonte: J. Robert Branstetter, "Radiant Heat Transfer Between Nongray
% Parallel Plates of Tungsten," NASA TN D-1088, 1961 (dominio publico).
% Baseado em De Vos (1954) corrigido por Larrabee (1959).
lam_tab_nm = [200 300 400 500 600 700 800 900 1000 1100 1250 1400 1600 1800 2000 2200 2500]';
eps_T0     = [0.4996 0.4996 0.4965 0.4856 0.4731 0.4639 0.4652 0.4485 0.4205 0.3934 0.3404 0.2995 0.2477 0.2052 0.1710 0.1439 0.1169]';
eps_T2000  = [0.4632 0.4632 0.4628 0.4518 0.4398 0.4296 0.4153 0.4000 0.3824 0.3656 0.3404 0.3182 0.2910 0.2684 0.2492 0.2321 0.2141]';
eps_T4000  = [0.4218 0.4218 0.4344 0.4221 0.4100 0.3916 0.3745 0.3608 0.3545 0.3470 0.3404 0.3361 0.3322 0.3266 0.3223 0.3174 0.3087]';

T_tab   = [0 2000 4000];
eps_tab = [eps_T0, eps_T2000, eps_T4000];   % matriz [lambda x T]

% Interpolação 2D (wavelength, temperature); 'linear' e o metodo
% justificado pelo proprio relatorio (interpolação linear em T entre
% 2000 e 4000K e "razoavelmente bem justificada" - Branstetter 1961).
eps_interp = @(lam_nm,T) interp2(T_tab, lam_tab_nm, eps_tab, ...
    max(min(T,4000),0), lam_nm, 'linear');

%% 3. Definição dos quatro modelos
planck_core = @(x,T) 1 ./ (x.^5 .* (exp(c2./(x.*T)) - 1));

% Modelo 1: corpo cinza puro (2 parametros: A, T)
model1 = @(p,x) p(1) .* planck_core(x,p(2));

% Modelo 2: emissividade linear empirica (3 parametros: A, T, alpha)
model2 = @(p,x) p(1) .* (1 - p(3).*x) .* planck_core(x,p(2));

% Modelo 3: emissividade tabelada REAL do tungstenio (2 parametros: A, T)
model3 = @(p,x) p(1) .* eps_interp(x*1000, p(2)) .* planck_core(x,p(2));

% Modelo 4: hibrido - tabela real x correção residual linear pequena
% (interpretavel como transmitancia do envelope de quartzo + efeitos
% geometricos residuais nao capturados pela emissividade do tungstenio
% nu) (3 parametros: A, T, beta)
model4 = @(p,x) p(1) .* eps_interp(x*1000, p(2)) .* (1 - p(3).*x) .* planck_core(x,p(2));

%% 4. Ajuste ponderado (minimos quadrados nao-lineares)
opts = optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10,'MaxFunctionEvaluations',5000);

resid_w = @(p,model) sqrt(w) .* (model(p,x_um) - y);

p1_0 = [1e4, 3000];
[p1,~,res1] = lsqnonlin(@(p) resid_w(p,model1), p1_0, [0 1000], [Inf 6000], opts);

p2_0 = [2.8e4, 3020, 0.13];
[p2,~,res2] = lsqnonlin(@(p) resid_w(p,model2), p2_0, [0 1000 -1], [Inf 6000 1], opts);

p3_0 = [7e4, 3000];
[p3,~,res3] = lsqnonlin(@(p) resid_w(p,model3), p3_0, [0 1000], [Inf 6000], opts);

p4_0 = [6e4, 3000, 0];
[p4,~,res4] = lsqnonlin(@(p) resid_w(p,model4), p4_0, [0 1000 -1], [Inf 6000 1], opts);

% Modelo 5: inspirado em Huang, Cebula & Hilsenrath (1998, Metrologia
% 35, 381-386) - "New procedure for interpolating NIST FEL lamp
% irradiances". A ideia central do artigo: a irradiancia da lampada FEL
% e o produto (corpo negro) x (emissividade do filamento) x
% (transmitancia do envelope de quartzo). O produto das duas ultimas
% tem curvatura negativa abaixo de ~300nm, e aproximadamente linear
% entre 300-450nm, e curvatura positiva acima de 450nm - descrita aqui
% por uma "gaussiana modificada com inclinação", continua em lambda0.
% Mantemos o Planck COMPLETO (com o -1), diferente do artigo original
% que usa a aproximação de Wien (valida so ate ~800-1600nm; no nosso
% caso o range vai a 2400nm, onde Wien perderia precisao).
lambda0 = 450;  % nm - ponto de inflexao de curvatura reportado no artigo
emiss_huang = @(lam_nm,c2,c3,c4,c5,c6) ...
    (lam_nm<lambda0) .* exp(c2.*(lam_nm-lambda0)/500 - c3.*abs((lam_nm-lambda0)/500).^c4) + ...
    (lam_nm>=lambda0) .* exp(c2.*(lam_nm-lambda0)/500 + c5.*abs((lam_nm-lambda0)/500).^c6);

model5 = @(p,x) p(1) .* emiss_huang(x*1000,p(3),p(4),p(5),p(6),p(7)) .* planck_core(x,p(2));

p5_0 = [7e4, 3000, -0.2, 1, 1.5, 1, 1.5];
lb5  = [0, 1000, -5,  0, 0.5,  0, 0.5];
ub5  = [Inf, 6000,  5, 50,  10, 50,  10];
[p5,~,res5] = lsqnonlin(@(p) resid_w(p,model5), p5_0, lb5, ub5, opts);
% AVISO: com 31 pontos e 7 parametros, c3 tende a colar no limite
% superior (curvatura da regiao <450nm mal restringida pelos dados).
% Verifique se p5(4) esta perto de ub5(4); se sim, considere fixar c4
% ou reduzir parametros livres.

%% 5. Metricas de comparação
compute_metrics = @(res,k) deal( ...
    sum(res.^2), ...
    1 - sum(res.^2)/sum(w.*(y-mean(y)).^2), ...
    n*log(sum(res.^2)/n) + 2*k, ...
    n*log(sum(res.^2)/n) + k*log(n));

[sse1,r2_1,aic1,bic1] = compute_metrics(res1,2);
[sse2,r2_2,aic2,bic2] = compute_metrics(res2,3);
[sse3,r2_3,aic3,bic3] = compute_metrics(res3,2);
[sse4,r2_4,aic4,bic4] = compute_metrics(res4,3);
[sse5,r2_5,aic5,bic5] = compute_metrics(res5,7);

fprintf('\n=== COMPARAÇÃO DE MODELOS (emissividade real: NASA TN D-1088) ===\n');
fprintf('%-38s %4s %10s %8s %8s %8s %8s\n','Modelo','k','SSE_w','R2','AIC','BIC','T(K)');
fprintf('%-38s %4d %10.4f %8.5f %8.2f %8.2f %8.0f\n','1) Corpo cinza (A,T)',                2,sse1,r2_1,aic1,bic1,p1(2));
fprintf('%-38s %4d %10.4f %8.5f %8.2f %8.2f %8.0f\n','2) Emiss. linear empirica (A,T,alpha)',3,sse2,r2_2,aic2,bic2,p2(2));
fprintf('%-38s %4d %10.4f %8.5f %8.2f %8.2f %8.0f\n','3) Emiss. tabelada real (A,T)',        2,sse3,r2_3,aic3,bic3,p3(2));
fprintf('%-38s %4d %10.4f %8.5f %8.2f %8.2f %8.0f\n','4) Hibrido: tabela + residual (A,T,b)',3,sse4,r2_4,aic4,bic4,p4(2));
fprintf('%-38s %4d %10.5f %8.5f %8.2f %8.2f %8.0f\n','5) Huang-like (A,T,c2..c6)',           7,sse5,r2_5,aic5,bic5,p5(2));
fprintf('\nalpha (modelo 2) = %.4f | beta (modelo 4) = %.4f\n', p2(3), p4(3));
fprintf('c3 (modelo 5) = %.3f  (limite superior = %.3f -> ', p5(4), ub5(4));
if abs(p5(4)-ub5(4)) < 1e-3
    fprintf('BATEU NO LIMITE: mal restringido, considere simplificar)\n');
else
    fprintf('OK, nao bateu no limite)\n');
end
fprintf('(Menor AIC/BIC = melhor compromisso ajuste-vs-complexidade)\n');

%% 6. Teste F entre modelos ANINHADOS
% 1 vs 2 (1 = 2 com alpha=0)
dfe1 = n-2; dfe2 = n-3;
F12 = ((sse1-sse2)/(dfe1-dfe2)) / (sse2/dfe2);
p12 = 1 - fcdf(F12, dfe1-dfe2, dfe2);
fprintf('\nTeste F (1 vs 2): F=%.2f, p=%.4g\n', F12, p12);

% 3 vs 4 (3 = 4 com beta=0) -- este e o teste relevante para justificar
% o termo residual de quartzo/geometria em cima da fisica real
dfe3 = n-2; dfe4 = n-3;
F34 = ((sse3-sse4)/(dfe3-dfe4)) / (sse4/dfe4);
p34 = 1 - fcdf(F34, dfe3-dfe4, dfe4);
fprintf('Teste F (3 vs 4): F=%.2f, p=%.4g\n', F34, p34);
fprintf('(p<0.05 => o termo residual beta melhora significativamente o modelo tabelado)\n');
fprintf('Nota: 2 e 4 nao sao aninhados entre si -> comparar so via AIC/BIC.\n');

%% 7. Graficos
x_plot = linspace(min(x_um),max(x_um),500)';

figure('Position',[60 60 1300 850]);

subplot(2,2,[1 2]);
plot(x_um*1000, y, 'ob','MarkerFaceColor','b','MarkerSize',5); hold on;
plot(x_plot*1000, model1(p1,x_plot), '-g','LineWidth',1.0);
plot(x_plot*1000, model2(p2,x_plot), '-r','LineWidth',1.0);
plot(x_plot*1000, model3(p3,x_plot), '-m','LineWidth',1.0);
plot(x_plot*1000, model4(p4,x_plot), '-k','LineWidth',1.2);
plot(x_plot*1000, model5(p5,x_plot), '-c','LineWidth',1.8);
xlabel('Comprimento de Onda (nm)'); ylabel('Irradiancia Espectral (W/cm^3)');
legend('Dados', sprintf('Cinza (T=%.0fK)',p1(2)), sprintf('Linear (T=%.0fK)',p2(2)), ...
    sprintf('Tabelada real (T=%.0fK)',p3(2)), sprintf('Hibrida (T=%.0fK)',p4(2)), ...
    sprintf('Huang-like (T=%.0fK)',p5(2)), 'Location','best');
title('Comparação dos Cinco Modelos'); grid on;

subplot(2,2,3);
r1pct = 100*(y - model1(p1,x_um))./y;
r2pct = 100*(y - model2(p2,x_um))./y;
plot(x_um*1000, r1pct, '-og','MarkerSize',4); hold on;
plot(x_um*1000, r2pct, '-sr','MarkerSize',4);
yline(0,'--k');
xlabel('Comprimento de Onda (nm)'); ylabel('Erro Residual (%)');
legend('Cinza','Linear (empirico)','Location','best');
title('Residuos: Cinza vs Linear'); grid on;

subplot(2,2,4);
r3pct = 100*(y - model3(p3,x_um))./y;
r4pct = 100*(y - model4(p4,x_um))./y;
r5pct = 100*(y - model5(p5,x_um))./y;
plot(x_um*1000, r3pct, '-dm','MarkerSize',4); hold on;
plot(x_um*1000, r4pct, '-ok','MarkerSize',4,'LineWidth',1.1);
plot(x_um*1000, r5pct, '-^c','MarkerSize',5,'LineWidth',1.6);
yline(0,'--k');
xlabel('Comprimento de Onda (nm)'); ylabel('Erro Residual (%)');
legend('Tabelada real','Hibrida','Huang-like','Location','best');
title('Residuos: Tabelada vs Hibrida vs Huang-like'); grid on;