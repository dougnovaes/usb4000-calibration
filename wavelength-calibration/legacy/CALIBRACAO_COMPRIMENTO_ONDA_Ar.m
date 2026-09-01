%% CALIBRACAO DE COMPRIMENTO DE ONDA -- Ocean Optics USB4000
%  Ajuste do polinomio cubico lambda(pixel) = a0 + a1*p + a2*p^2 + a3*p^3
%  a partir de linhas de emissao catalogadas de lampadas de Hg e de Ar.
%
%  ------------------------------------------------------------------------
%  CONTEXTO
%  ------------------------------------------------------------------------
%  Etapa A do pipeline de calibracao absoluta do USB4000 (treino preliminar
%  antes do mesmo procedimento no Jobin-Yvon THR1000 e no McPherson 207).
%  Aplicacao final: caracterizacao espectral de plasma frio DBD.
%
%  O detector nao mede comprimento de onda -- mede qual dos 3648 pixels
%  recebeu luz. A relacao pixel->lambda e uma propriedade da montagem optica
%  (passo da grade, distancia focal, angulo) e precisa ser MEDIDA apontando
%  fontes de lambda conhecido. Como a equacao da grade relaciona lambda ao
%  ANGULO de difracao (e nao diretamente a posicao no array plano), projetar
%  um leque de angulos sobre um detector linear introduz curvatura residual
%  -- dai o polinomio cubico em vez de uma reta.
%
%  ------------------------------------------------------------------------
%  FONTES DOS DADOS
%  ------------------------------------------------------------------------
%  [1] Comprimentos de onda de referencia:
%      Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024).
%      NIST Atomic Spectra Database (ver. 5.12). NIST, Gaithersburg, MD.
%      https://physics.nist.gov/asd    DOI: 10.18434/T4W30F
%      Consulta pelo formulario interativo "Lines Form", com:
%        Spectrum          : 'Hg I-II' / 'Ar I-II'
%        Wavelength in     : opcao contendo Air (200-2000 nm)  <-- ver NOTA 1
%        Rel. intensity min: EM BRANCO                          <-- ver NOTA 2
%
%  [2] Pixels medidos: identificacao de picos no Ocean Optics SpectraSuite
%      (Show Peak Info), coluna 'Wavelength' (pico de intensidade maxima).
%      Hg: tempo de integracao 300 ms.  Ar: 250 ms.               <-- NOTA 3
%
%  [3] Coeficientes de fabrica: Spectrometer > Spectrometer Features >
%      aba Wavelength (lidos ANTES de qualquer gravacao; ver NOTA 5).
%
%  ------------------------------------------------------------------------
%  NOTAS METODOLOGICAS (lecoes praticas deste trabalho)
%  ------------------------------------------------------------------------
%  NOTA 1 (vacuo vs. ar): a tabela curada "Strong Lines of Argon" do NIST
%    Handbook lista 580-1066 nm em comprimento de onda de VACUO. O USB4000
%    mede no ar. Usar o formulario ASD com saida em ar elimina a conversao
%    manual (fator aprox. 1-2,8e-4 no visivel) e a chance de erro associada.
%
%  NOTA 2 (nao filtrar por intensidade): a tabela "Strong Lines" e um
%    subconjunto pequeno; usa-la como fonte unica deixou a maioria dos picos
%    de Ar sem correspondencia confiavel. Pedir a lista completa e filtrar
%    depois pelo criterio proprio e mais seguro.
%
%  NOTA 3 (qual coluna usar): usar SEMPRE 'Wavelength' (pico de maximo),
%    nunca 'Center Wavelength' nem 'Centroid'. Espectrometros Czerny-Turner
%    compactos exibem coma, que alarga picos assimetricamente; medidas de
%    centro baseadas em media ponderada ou ponto medio de banda sao puxadas
%    sistematicamente na direcao da cauda -- vies que NAO cancela com mais
%    medicoes. O pico de maximo e robusto a isso.
%    O campo '90% Bandwidth' cumpre papel diferente e complementar: e o
%    filtro de QUALIDADE (o pico e real?), nao a fonte do VALOR.
%
%  NOTA 4 (threshold do detector de picos): o SpectraSuite usa uma unica
%    linha de threshold para o grafico inteiro. Com picos de alturas muito
%    diferentes, um threshold bom para o pico mais alto faz os mais baixos
%    reportarem banda de 90% absurda (centenas de nm) mesmo sendo reais.
%    Solucao: dar zoom em janela local e reajustar o threshold ali.
%
%  NOTA 5 (backup): a gravacao de novos coeficientes SOBRESCREVE os de
%    fabrica. Anote-os antes (ja feito -- ver variavel COEF_FABRICA).
%
%  ------------------------------------------------------------------------
%  ESTADO / RESSALVA IMPORTANTE
%  ------------------------------------------------------------------------
%  As 8 linhas de Hg tem identificacao SOLIDA (RMS proprio ~0,08 nm).
%  As 12 linhas de Ar tem confianca MENOR: foram casadas por proximidade a
%  previsao do ajuste de Hg + intensidade relativa alta, e o resíduo delas e
%  ~6x maior que o das linhas de Hg. Um pico limpo adicional (pixel 2319,
%  651,63 nm bruto) permanece SEM correspondencia e foi deixado de fora.
%  Estes coeficientes AINDA NAO FORAM GRAVADOS no instrumento.
%
%  Autor: Douglas Novaes      Data: agosto de 2026
%  ------------------------------------------------------------------------

clear; clc; close all;

%% ========================================================================
%  1. DADOS
%  ========================================================================
%  Colunas: [pixel_medido, lambda_medido_nm, banda90_nm, lambda_NIST_nm]
%  lambda_medido e a leitura sob a calibracao DE FABRICA (nao e o alvo;
%  serve so para diagnostico do desvio). O alvo do ajuste e lambda_NIST.

% --- Mercurio (Hg), IT = 300 ms ------------------------------------------
% Identificacao solida. Todas em convencao de ar (padrao para estas linhas).
% A quarta linha (pixel 889) exigiu desempate: o pico cai sobre o tripleto
% 365,0158 / 365,4842 / 366,2887 nm, nao resolvido nesta banda (0,83 nm).
% A intensidade relativa do NIST resolve: 9000 / 3000 / 500. A componente
% dominante (365,0158) produz desvio coerente com as demais 7 linhas; a mais
% fraca (366,2887) produziria um desvio fora do padrao. Ver diagnostico da
% secao 3, que torna esse padrao visivel.
dados_Hg = [
%  pixel  medido    banda90   NIST(ar)   % Rel.Int.  transicao
     361  255.30    0.85     253.6521   % 900000
     640  314.27    0.84     313.1555   %   3000
     889  366.32    0.83     365.0158   %   9000   (tripleto, ver acima)
    1080  405.82    0.41     404.6565   %  12000
    1232  436.99    0.41     435.8335   %  12000
    1780  547.13    0.40     546.0750   %   6000
    1937  577.99    3.32     576.9610   %   1000
    1948  580.14    3.13     579.0670   %    900
];

% --- Argonio (Ar), IT = 250 ms -------------------------------------------
% CONFIANCA MENOR -- ver ressalva no cabecalho. Predominam linhas de Ar II
% (argonio ionizado) nesta faixa; as linhas fortes de Ar I catalogadas ficam
% acima de ~690 nm.
dados_Ar = [
%  pixel  medido    banda90   NIST(ar)   % especie
    1979  586.18    1.17     586.0310   % Ar I
    2027  595.58    0.78     595.0905   % Ar II
    2129  615.26    1.16     613.8655   % Ar II
    2168  622.77    0.38     621.2503   % Ar I
    2194  627.76    0.38     627.7425   % Ar II
    2229  634.47    0.38     634.8227   % Ar II
    2265  641.35    1.53     641.6307   % Ar I
    2368  660.92    0.38     660.4853   % Ar I
    2410  668.86    0.38     668.4292   % Ar II
    2544  693.99    0.37     693.7664   % Ar I
    2599  704.23    0.37     705.4993   % Ar II
    2715  725.66    0.37     727.2936   % Ar I
%   2319  651.63    1.14        NaN     % SEM correspondencia -- excluido
];

% --- Coeficientes de fabrica (BACKUP -- ver NOTA 5) -----------------------
% Ordem crescente [a0 a1 a2 a3]; convertidos para a ordem do polyval abaixo.
COEF_FABRICA = [178.085770, 2.1517764E-1, -3.3423826E-6, -5.9958050E-10];

% --- Parametros de controle ----------------------------------------------
GRAU        = 3;      % grau do polinomio (3 = padrao Ocean Optics)
BANDA_MAX   = 5.0;    % nm; descarta picos com banda de 90% acima disto
N_PIX       = 3648;   % pixels do USB4000 (so para os graficos)

%% ========================================================================
%  2. PREPARACAO E FILTRO DE QUALIDADE
%  ========================================================================
pixel  = [dados_Hg(:,1); dados_Ar(:,1)];
medido = [dados_Hg(:,2); dados_Ar(:,2)];
banda  = [dados_Hg(:,3); dados_Ar(:,3)];
alvo   = [dados_Hg(:,4); dados_Ar(:,4)];
eh_Hg  = [true(size(dados_Hg,1),1); false(size(dados_Ar,1),1)];

% Filtro de qualidade pela banda de 90% (ver NOTA 3/4)
ok = banda <= BANDA_MAX & isfinite(alvo);
if ~all(ok)
    fprintf('Descartados %d ponto(s) por banda > %.1f nm ou alvo ausente.\n', ...
            sum(~ok), BANDA_MAX);
end
pixel = pixel(ok); medido = medido(ok); alvo = alvo(ok);
banda = banda(ok);  eh_Hg = eh_Hg(ok);
n = numel(pixel);

assert(n >= GRAU+2, ['Pontos insuficientes: com %d pontos e grau %d nao ' ...
    'sobra grau de liberdade para avaliar a qualidade do ajuste.'], n, GRAU);

fprintf('\n=== DADOS ===\n');
fprintf('Pontos usados: %d (%d Hg + %d Ar)\n', n, sum(eh_Hg), sum(~eh_Hg));
fprintf('Faixa: pixel %d a %d  |  lambda %.1f a %.1f nm\n', ...
        min(pixel), max(pixel), min(alvo), max(alvo));
fprintf('Graus de liberdade: %d\n', n-(GRAU+1));

%% ========================================================================
%  3. DIAGNOSTICO PREVIO: desvio da calibracao de fabrica
%  ========================================================================
%  Antes de ajustar, vale olhar o desvio (fabrica - NIST) linha a linha.
%  Se as identificacoes estiverem corretas, esse desvio deve variar de forma
%  SUAVE com o pixel. Um ponto fora do padrao denuncia identificacao errada
%  -- foi assim que se detectou o desempate do tripleto de Hg.
lam_fabrica = polyval(fliplr(COEF_FABRICA), pixel);
desvio_fab  = lam_fabrica - alvo;

fprintf('\n=== DESVIO DA CALIBRACAO DE FABRICA (diagnostico) ===\n');
fonte = repmat("Ar", n, 1); fonte(eh_Hg) = "Hg";   % rotulo por ponto
fprintf('%6s %10s %12s %10s %6s\n','pixel','NIST(nm)','fabrica(nm)','desvio','fonte');
for i = 1:n
    fprintf('%6d %10.4f %12.4f %+10.4f %6s\n', pixel(i), alvo(i), ...
            lam_fabrica(i), desvio_fab(i), fonte(i));
end
fprintf('RMS fabrica = %.4f nm   |  max|desvio| = %.4f nm\n', ...
        rms(desvio_fab), max(abs(desvio_fab)));

%% ========================================================================
%  4. AJUSTE DO POLINOMIO
%  ========================================================================
%  BOA PRATICA (condicionamento numerico): ajustar diretamente em pixel bruto
%  gera uma matriz de Vandermonde mal condicionada -- com pixels ate ~2700 e
%  grau 3, cond(V) ~ 1e11, e o MATLAB emite warning de "poorly conditioned".
%  A forma [p,S,mu] = polyfit(...) centra e escala o preditor internamente
%  (s = (x-mu(1))/mu(2)), levando cond(V) para ~9. Os dois ajustes sao
%  matematicamente equivalentes; o escalado e numericamente robusto.
%
%  Como o INSTRUMENTO precisa dos coeficientes em pixel BRUTO, o script
%  ajusta na variavel escalada e converte de volta ao final -- e verifica
%  explicitamente que as duas rotas concordam.

[coef_esc, S, mu] = polyfit(pixel, alvo, GRAU);   % ajuste centrado/escalado

% --- Conversao dos coeficientes escalados para pixel bruto ---------------
% lambda = sum_k coef_esc(k) * s^(GRAU-k+1),  com s = (p - mu(1))/mu(2).
% Expandir em potencias de p via convolucao dos binomios.
coef_raw = zeros(1, GRAU+1);
for k = 1:GRAU+1
    g = GRAU - k + 1;                       % expoente de s neste termo
    termo = 1;
    for j = 1:g                             % (p - mu1)/mu2 elevado a g
        termo = conv(termo, [1/mu(2), -mu(1)/mu(2)]);
    end
    coef_raw(end-g:end) = coef_raw(end-g:end) + coef_esc(k) * termo;
end

% --- Verificacao independente: ajuste direto em pixel bruto --------------
% Rodado apenas para CONFERIR a conversao; nao e usado adiante. O warning de
% condicionamento aqui e esperado e faz parte da demonstracao.
aviso = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
coef_direto = polyfit(pixel, alvo, GRAU);
warning(aviso);
discrep = max(abs(polyval(coef_raw,pixel) - polyval(coef_direto,pixel)));
fprintf('\nConcordancia escalado->bruto vs. ajuste direto: %.2e nm\n', discrep);
assert(discrep < 1e-6, 'Conversao de coeficientes inconsistente.');

% --- Coeficientes finais, na ordem da tela do SpectraSuite ---------------
a = fliplr(coef_raw);           % [a0 a1 a2 a3], ordem crescente
fprintf('\n=== COEFICIENTES AJUSTADOS (ordem do SpectraSuite) ===\n');
fprintf('Intercept      (a0) = %.6f\n',   a(1));
fprintf('1st Coefficient(a1) = %.10E\n',  a(2));
fprintf('2nd Coefficient(a2) = %.10E\n',  a(3));
fprintf('3rd Coefficient(a3) = %.10E\n',  a(4));

%% ========================================================================
%  5. RESIDUOS E METRICAS
%  ========================================================================
previsto = polyval(coef_raw, pixel);
residuo  = previsto - alvo;

rms_tot = rms(residuo);
rms_Hg  = rms(residuo(eh_Hg));
rms_Ar  = rms(residuo(~eh_Hg));

fprintf('\n=== RESIDUOS ===\n');
fprintf('%6s %10s %11s %10s %6s\n','pixel','NIST(nm)','previsto','residuo','fonte');
for i = 1:n
    fprintf('%6d %10.4f %11.4f %+10.4f %6s\n', pixel(i), alvo(i), ...
            previsto(i), residuo(i), fonte(i));
end
fprintf('\nRMS total = %.4f nm  (max|res| = %.4f nm)\n', rms_tot, max(abs(residuo)));
fprintf('RMS so Hg = %.4f nm  |  RMS so Ar = %.4f nm\n', rms_Hg, rms_Ar);
fprintf('Melhoria sobre a fabrica: fator %.2f\n', rms(desvio_fab)/rms_tot);

%% ========================================================================
%  6. VALIDACAO CRUZADA LEAVE-ONE-OUT
%  ========================================================================
%  O RMS da secao 5 e medido nos MESMOS pontos usados para ajustar -- nao
%  distingue um ajuste que capturou a fisica de um que memorizou os pontos.
%  A LOOCV remove um ponto, reajusta com os demais, preve o removido, e
%  repete. Se LOOCV >> RMS interno, ha sobreajuste.
%  (Mesma logica aplicada ao estudo da lampada FEL, onde ela reprovou o GPR.)
err_loo = zeros(n,1);
for i = 1:n
    m = true(n,1); m(i) = false;
    % forma [c,~,mu]: centra/escala internamente (mesma pratica da secao 4)
    [c_i, ~, mu_i] = polyfit(pixel(m), alvo(m), GRAU);
    s_i = (pixel(i) - mu_i(1)) / mu_i(2);
    err_loo(i) = polyval(c_i, s_i) - alvo(i);
end
rms_loo = rms(err_loo);
fprintf('\n=== VALIDACAO CRUZADA (LOOCV) ===\n');
fprintf('RMS LOOCV = %.4f nm  (dentro da amostra: %.4f nm)\n', rms_loo, rms_tot);
fprintf('Razao LOOCV/interno = %.2f  ', rms_loo/rms_tot);
if rms_loo/rms_tot < 1.5
    fprintf('--> sem sinal de sobreajuste.\n');
else
    fprintf('--> ATENCAO: possivel sobreajuste.\n');
end

%% ========================================================================
%  7. GRAFICOS
%  ========================================================================
figure('Units','normalized','OuterPosition',[0.05 0.1 0.9 0.8]);

% (a) curva de calibracao
subplot(2,2,1); hold on; grid on; box on;
pp = linspace(0, N_PIX, 500);
plot(pp, polyval(coef_raw,pp), 'k-', 'LineWidth',1.4);
plot(pp, polyval(fliplr(COEF_FABRICA),pp), '--', 'Color',[.6 .6 .6],'LineWidth',1.1);
plot(pixel(eh_Hg),  alvo(eh_Hg),  'o','MarkerFaceColor',[0 .45 .74],'MarkerEdgeColor','k','MarkerSize',7);
plot(pixel(~eh_Hg), alvo(~eh_Hg), 's','MarkerFaceColor',[.85 .33 .1],'MarkerEdgeColor','k','MarkerSize',7);
xlabel('Pixel'); ylabel('Comprimento de onda (nm)');
title('Curva de calibracao \lambda(pixel)');
legend({'ajuste','fabrica','Hg','Ar'},'Location','southeast');

% (b) residuos do ajuste
subplot(2,2,2); hold on; grid on; box on;
yline(0,'k-');
stem(pixel(eh_Hg),  residuo(eh_Hg), 'o','Color',[0 .45 .74],'MarkerFaceColor',[0 .45 .74]);
stem(pixel(~eh_Hg), residuo(~eh_Hg),'s','Color',[.85 .33 .1],'MarkerFaceColor',[.85 .33 .1]);
yline( rms_tot,':k'); yline(-rms_tot,':k');
xlabel('Pixel'); ylabel('Residuo (nm)');
title(sprintf('Residuos  (RMS = %.3f nm; Hg %.3f / Ar %.3f)', rms_tot, rms_Hg, rms_Ar));
legend({'','Hg','Ar','\pm RMS'},'Location','best');

% (c) comparacao com a fabrica
subplot(2,2,3); hold on; grid on; box on;
yline(0,'k-');
plot(pixel, desvio_fab, 'v--','Color',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6]);
plot(pixel, residuo,    'o-','Color',[0 .45 .74],'MarkerFaceColor',[0 .45 .74]);
xlabel('Pixel'); ylabel('Erro vs. NIST (nm)');
title(sprintf('Fabrica (RMS %.3f nm) vs. ajuste (RMS %.3f nm)', rms(desvio_fab), rms_tot));
legend({'','fabrica','ajuste novo'},'Location','best');

% (d) dispersao ao longo do array
subplot(2,2,4); grid on; box on;
disp_nm = polyval(polyder(coef_raw), pp);
plot(pp, disp_nm, 'k-','LineWidth',1.4);
xlabel('Pixel'); ylabel('d\lambda/dpixel (nm/pixel)');
title('Dispersao (derivada do polinomio)');

%% ========================================================================
%  8. EXPORTACAO
%  ========================================================================
T = table(pixel, alvo, previsto, residuo, banda, ...
          categorical(eh_Hg, [true false], {'Hg','Ar'}), ...
    'VariableNames', {'pixel','lambda_NIST_nm','previsto_nm','residuo_nm', ...
                      'banda90_nm','fonte'});
writetable(T, 'calibracao_residuos.csv');

fid = fopen('calibracao_coeficientes.txt','w');
fprintf(fid, '# Calibracao de comprimento de onda USB4000 -- %s\n', datestr(now,'yyyy-mm-dd'));
fprintf(fid, '# Ajuste: polinomio grau %d, %d pontos (%d Hg + %d Ar)\n', ...
        GRAU, n, sum(eh_Hg), sum(~eh_Hg));
fprintf(fid, '# RMS interno %.4f nm | LOOCV %.4f nm | fabrica %.4f nm\n', ...
        rms_tot, rms_loo, rms(desvio_fab));
fprintf(fid, '# Referencias: NIST ASD ver. 5.12 (DOI 10.18434/T4W30F)\n');
fprintf(fid, '# ATENCAO: identificacao das linhas de Ar e provisoria.\n');
fprintf(fid, 'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n', a(1),a(2),a(3),a(4));
fprintf(fid, '# COEFICIENTES DE FABRICA (backup):\n');
fprintf(fid, '# Intercept\t%.6f\n# 1st\t%.7E\n# 2nd\t%.7E\n# 3rd\t%.7E\n', COEF_FABRICA);
fclose(fid);

fprintf('\nArquivos gravados: calibracao_residuos.csv, calibracao_coeficientes.txt\n');
fprintf(['\nLEMBRETE: estes coeficientes ainda NAO foram gravados no ' ...
         'instrumento.\nAntes de gravar, resolver a linha de Ar sem ' ...
         'correspondencia (pixel 2319)\ne confirmar o backup de fabrica.\n']);