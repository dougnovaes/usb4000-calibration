%% CALIBRAÇÃO DE COMPRIMENTO DE ONDA — Ocean Optics USB4000
%  Ajuste do polinômio cúbico  λ(pixel) = a0 + a1·p + a2·p² + a3·p³
%  a partir de linhas de emissão catalogadas de lâmpadas de Hg e de Ne,
%  com determinação sub-pixel da posição de cada pico.
%
%  ------------------------------------------------------------------------
%  CONTEXTO
%  ------------------------------------------------------------------------
%  Etapa A do procedimento de calibração absoluta do USB4000 (treino
%  preliminar antes do mesmo trabalho no Jobin-Yvon THR1000 e no McPherson
%  207). Aplicação final: caracterização espectral de plasma frio DBD.
%
%  O detector não mede comprimento de onda — mede apenas qual dos 3648
%  pixels recebeu luz, e quanto. A relação pixel → λ é uma propriedade da
%  montagem óptica daquele instrumento específico (passo da grade, distância
%  focal, ângulo) e precisa ser MEDIDA, apontando fontes de λ conhecido.
%  Como a equação da grade relaciona λ ao ÂNGULO de difração, e não
%  diretamente à posição no detector plano, projetar um leque de ângulos
%  sobre um array linear introduz curvatura residual — daí o polinômio
%  cúbico em vez de uma reta.
%
%  ------------------------------------------------------------------------
%  ATENÇÃO — A SEGUNDA LÂMPADA É DE NEÔNIO, NÃO DE ARGÔNIO
%  ------------------------------------------------------------------------
%  Uma versão anterior deste trabalho tratou a segunda lâmpada como sendo de
%  argônio, o que produziu identificações de linha inconsistentes (o desvio
%  em relação à calibração de fábrica saltava de forma fisicamente
%  impossível entre pixels vizinhos). O espectro é inequivocamente de Ne I:
%    · 18 dos 19 picos casam com linhas fortes de Ne I dentro de 0,25 nm
%      apenas extrapolando a cúbica ajustada às linhas de Hg, sem nenhum
%      parâmetro ajustado para forçar o encaixe;
%    · as linhas mais fortes de Ar I (696,5 / 706,7 / 750,4 / 763,5 /
%      811,5 nm) estão AUSENTES do espectro — teste negativo decisivo.
%  Lâmpadas de calibração Hg-Ne são tão comuns quanto as Hg-Ar.
%
%  ------------------------------------------------------------------------
%  FONTES DOS DADOS
%  ------------------------------------------------------------------------
%  [1] Comprimentos de onda de referência:
%      Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024).
%      NIST Atomic Spectra Database (ver. 5.12). NIST, Gaithersburg, MD.
%      https://physics.nist.gov/asd     DOI: 10.18434/T4W30F
%      Consulta pelo formulário interativo "Lines Form", com:
%        Spectrum            : 'Hg I-II'  /  'Ne'
%        Wavelength in       : opção contendo Air (200–2000 nm)   → NOTA 1
%        Rel. intensity min. : EM BRANCO                          → NOTA 2
%
%  [2] Espectros brutos: exportados do Ocean Optics SpectraSuite
%      (arquivo de duas colunas: λ de fábrica [nm], contagens).
%      Hg: tempo de integração 300 ms.  Ne: 250 ms.
%      Correção eletrônica de escuro (Electric Dark Correction) ativa.
%
%  [3] Coeficientes de fábrica: Spectrometer → Spectrometer Features →
%      aba Wavelength, lidos ANTES de qualquer gravação (ver NOTA 5).
%
%  ------------------------------------------------------------------------
%  NOTAS METODOLÓGICAS (lições práticas deste trabalho)
%  ------------------------------------------------------------------------
%  NOTA 1 (vácuo vs. ar): tabelas de referência rápida frequentemente listam
%    comprimentos de onda de VÁCUO acima de 200 nm sem sinalizar isso. O
%    USB4000 mede no ar. Pedir a saída da ASD já em ar elimina a conversão
%    manual (fator ≈ 1 − 2,8×10⁻⁴ no visível) e o erro associado.
%
%  NOTA 2 (não filtrar por intensidade): a tabela curada "Strong Lines" do
%    NIST Handbook é um subconjunto pequeno; usá-la como fonte única deixa a
%    maioria dos picos sem correspondência confiável. Peça a lista completa
%    e filtre depois pelo seu próprio critério.
%
%  NOTA 3 (POSIÇÃO SUB-PIXEL — o ponto central deste script):
%    Usar o pixel de máximo como posição do pico limita a precisão a cerca
%    de ±0,1–0,2 nm (metade de um pixel), o que é da mesma ordem do
%    espaçamento típico entre linhas catalogadas — tornando a identificação
%    ambígua. A parábola ajustada aos 3 pontos centrais do pico devolve a
%    posição com ~1/10 de pixel (≈0,02 nm).
%
%    Por que PARÁBOLA NO TOPO e não gaussiana no perfil inteiro: os picos
%    do USB4000 têm cauda assimétrica para o vermelho (aberração de coma,
%    típica de Czerny-Turner compacto). Uma gaussiana ajustada ao perfil
%    todo é puxada pela cauda; a parábola nos 3 pontos centrais enxerga
%    apenas o núcleo do pico, onde a assimetria ainda é desprezível.
%    Medida neste conjunto de dados: o centroide fica sistematicamente
%    ~0,38 pixel (≈0,08 nm) mais para o vermelho que a parábola — viés
%    sistemático, que NÃO cancela com mais medições.
%
%    Efeito prático: a dispersão do desvio (fábrica − NIST) nas linhas de Hg
%    cai de ~0,20 nm (pixel inteiro) para 0,057 nm (parábola sub-pixel).
%
%  NOTA 4 (limiar do detector de picos do SpectraSuite): o programa usa uma
%    única linha de limiar para o gráfico inteiro. Com picos de alturas
%    muito diferentes, um limiar bom para o pico mais alto faz os mais
%    baixos reportarem banda de 90 % absurda (centenas de nm) mesmo sendo
%    reais. Este script não depende disso — detecta os picos diretamente do
%    espectro bruto — mas a observação vale ao usar a interface.
%
%  NOTA 5 (backup): gravar novos coeficientes SOBRESCREVE os de fábrica.
%    Anote-os antes (ver a variável COEF_FABRICA).
%
%  NOTA 7 (PROEMINÊNCIA — filtro indispensável): um limiar de altura sozinho
%    NÃO basta para detectar picos. Numa primeira versão deste script, uma
%    ondulação de ruído na subida da linha saturada de 253,65 nm (altura
%    1700 contagens, acima do limiar de 1284) foi tomada como pico e casada
%    com a própria linha de 253,65 nm, corrompendo todo o ajuste — a taxa de
%    identificação do Ne despencou de 18/19 para 4/25.
%    A defesa é a PROEMINÊNCIA: quanto o pico se ergue acima do maior dos
%    dois "colos" que o cercam. Naquela ondulação, a proeminência era de
%    apenas 31 contagens, ou 1,8 % da própria altura; num pico real, a razão
%    proeminência/altura fica perto de 100 %. O critério
%    PROEM_MIN = 0,30 separa os dois casos com folga, e ainda preserva
%    dubletos parcialmente resolvidos (como Ne 638,3/640,2 nm, cuja razão
%    fica em 0,97).
%
%  NOTA 6 (saturação): a linha de Hg em 253,65 nm satura o detector nesta
%    aquisição (vários pixels no mesmo valor de topo). O topo achatado torna
%    a parábola inválida, e por isso essa linha é EXCLUÍDA do ajuste. O
%    script detecta e avisa automaticamente.
%
%  ------------------------------------------------------------------------
%  RESULTADO ESPERADO (para conferência)
%  ------------------------------------------------------------------------
%    26 pontos (7 Hg + 19 Ne)
%    Intercept        =  177,259857
%    1º Coeficiente   =  2,1437262236E-01
%    2º Coeficiente   = -2,7305369510E-06
%    3º Coeficiente   = -7,2870128288E-10
%    RMS do resíduo   =  0,0212 nm   (máx. 0,0575 nm)
%    LOOCV            =  0,0395 nm
%    Fábrica          =  1,0424 nm   → melhoria de ~49×
%
%  ------------------------------------------------------------------------
%  ESTADO
%  ------------------------------------------------------------------------
%  Estes coeficientes AINDA NÃO FORAM GRAVADOS no instrumento. Antes de
%  gravar: confirmar na etiqueta/documentação que a lâmpada é mesmo Hg-Ne, e
%  preservar o backup de fábrica.
%
%  Autor: (preencher)          Data: agosto de 2026
%  ------------------------------------------------------------------------

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURAÇÃO
%  ========================================================================
ARQ_HG   = 'espectro_hg.txt';   % duas colunas: λ de fábrica [nm], contagens
ARQ_NE   = 'espectro_ne.txt';
GRAU     = 3;        % grau do polinômio (3 = padrão Ocean Optics)
N_PARAB  = 3;        % pontos usados na parábola do topo (3 é o mais robusto)
FRAC_LIM = 0.02;     % limiar de altura, como fração do pico máximo
PROEM_MIN= 0.30;     % proeminência mínima, como fração da altura (NOTA 7)
SEP_MIN  = 5;        % separação mínima entre picos, em pixels
DESL_MAX = 1.0;      % deslocamento máx. [pixel] do vértice em relação ao máximo
TOL_HG   = 1.2;      % tolerância [nm] para casar as linhas de Hg
TOL_ID   = 0.30;     % tolerância [nm] para casar as linhas de Ne
TAXA_MIN = 0.60;     % taxa mínima de acerto do padrão de Ne para prosseguir
N_PIX    = 3648;     % pixels do USB4000

% Coeficientes de fábrica — BACKUP (ver NOTA 5), ordem crescente [a0 a1 a2 a3]
COEF_FABRICA = [178.085770, 2.1517764E-1, -3.3423826E-6, -5.9958050E-10];

%% ========================================================================
%  2. TABELAS DE REFERÊNCIA (NIST ASD, comprimento de onda no ar)
%  ========================================================================
% Linhas de Hg I esperadas na faixa do USB4000. A de 253,6521 nm é listada
% para referência mas satura nesta aquisição (ver NOTA 6).
% A linha em 365,0158 exige atenção: o pico observado é o tripleto
% 365,0158 / 365,4842 / 366,2887 não resolvido. A intensidade relativa do
% NIST resolve a ambiguidade — 9000 / 3000 / 500 — de modo que a componente
% dominante é a de 365,0158 nm. Usar a mais fraca produziria um desvio fora
% do padrão suave observado nas demais linhas.
LINHAS_HG = [253.6521; 313.1555; 365.0158; 404.6565; 435.8335; ...
             546.0750; 576.9610; 579.0670];

% Linhas fortes de Ne I (ar). Faixa 585–745 nm.
LINHAS_NE = [585.2488; 588.1895; 594.4834; 597.5534; 603.0000; 607.4338; ...
             609.6163; 614.3062; 616.3594; 621.7281; 626.6495; 630.4789; ...
             633.4428; 638.2991; 640.2248; 650.6528; 653.2882; 659.8953; ...
             667.8276; 671.7043; 692.9467; 703.2413; 717.3938; 724.5167; ...
             743.8899];

%% ========================================================================
%  3. LEITURA DOS ESPECTROS
%  ========================================================================
[pix_hg, cnt_hg] = ler_espectro(ARQ_HG, COEF_FABRICA);
[pix_ne, cnt_ne] = ler_espectro(ARQ_NE, COEF_FABRICA);
fprintf('Espectros carregados: Hg (%d pontos), Ne (%d pontos)\n', ...
        numel(cnt_hg), numel(cnt_ne));

%% ========================================================================
%  4. DETECÇÃO DE PICOS E POSIÇÃO SUB-PIXEL
%  ========================================================================
[p_hg, sat_hg] = picos_subpixel(pix_hg, cnt_hg, FRAC_LIM, PROEM_MIN, SEP_MIN, N_PARAB, DESL_MAX);
[p_ne, sat_ne] = picos_subpixel(pix_ne, cnt_ne, FRAC_LIM, PROEM_MIN, SEP_MIN, N_PARAB, DESL_MAX);

if any(sat_hg)
    fprintf(['AVISO: %d pico(s) rejeitado(s) no espectro de Hg ' ...
             '(saturação ou vértice inválido).\n'], sum(sat_hg));
end
if any(sat_ne)
    fprintf('AVISO: %d pico(s) rejeitado(s) no espectro de Ne.\n', sum(sat_ne));
end
p_hg = p_hg(~sat_hg);
p_ne = p_ne(~sat_ne);
fprintf('Picos utilizáveis: Hg = %d, Ne = %d\n', numel(p_hg), numel(p_ne));

%% ========================================================================
%  5. IDENTIFICAÇÃO DAS LINHAS
%  ========================================================================
%  Estratégia em duas etapas — deliberadamente NÃO é um casamento pico a
%  pico pelo mais próximo, que é ambíguo quando as linhas catalogadas estão
%  mais juntas que a incerteza de posição.
%
%  Etapa 1: as linhas de Hg são identificadas primeiro. São poucas, fortes,
%           bem separadas, e a correspondência é inequívoca.
%  Etapa 2: ajusta-se uma cúbica APENAS ao Hg e extrapola-se para a região
%           do Ne. O padrão INTEIRO de picos de Ne deve então cair sobre
%           linhas catalogadas de uma vez. É essa verificação global —
%           muitos picos simultaneamente, não um por vez — que dá confiança
%           na identificação. Se o padrão não encaixar, a hipótese sobre a
%           espécie da lâmpada está errada (foi exatamente assim que a
%           suposição de argônio foi refutada).

% --- Etapa 1: Hg ---------------------------------------------------------
lam_fab_hg = polyval(fliplr(COEF_FABRICA), p_hg);
[id_hg, ok_hg] = casar(lam_fab_hg, LINHAS_HG, TOL_HG);
p_hg = p_hg(ok_hg);  lam_hg = id_hg(ok_hg);
fprintf('\nLinhas de Hg identificadas: %d\n', numel(p_hg));

% Verificação de consistência: o desvio (fábrica − NIST) tem de variar de
% forma SUAVE com o pixel, porque tanto a calibração de fábrica quanto a
% curva verdadeira são polinômios suaves. Um ponto fora do padrão denuncia
% identificação errada — foi este teste que expôs, num estágio anterior
% deste trabalho, que a segunda lâmpada não era de argônio.
desv_hg = polyval(fliplr(COEF_FABRICA), p_hg) - lam_hg;
fprintf('Desvio fábrica−NIST nas linhas de Hg: %.3f a %.3f nm (dispersão %.3f nm)\n', ...
        min(desv_hg), max(desv_hg), std(desv_hg));
if std(desv_hg) > 0.15
    warning(['A dispersão do desvio nas linhas de Hg (%.3f nm) é grande ' ...
             'demais para ser só quantização. Suspeite de identificação ' ...
             'errada ou de pico falso — confira a tabela acima.'], std(desv_hg));
end

assert(numel(p_hg) >= GRAU+2, ...
    'Linhas de Hg insuficientes (%d) para ancorar o ajuste.', numel(p_hg));

% --- Etapa 2: Ne, por extrapolação da cúbica do Hg -----------------------
[c_hg, ~, mu_hg] = polyfit(p_hg, lam_hg, GRAU);
lam_extrap = polyval(c_hg, (p_ne - mu_hg(1))/mu_hg(2));
[id_ne, ok_ne] = casar(lam_extrap, LINHAS_NE, TOL_ID);
p_ne = p_ne(ok_ne);  lam_ne = id_ne(ok_ne);

fprintf(['Linhas de Ne identificadas: %d de %d picos ' ...
         '(tolerância %.2f nm, SEM ajuste livre)\n'], ...
         numel(p_ne), numel(ok_ne), TOL_ID);
taxa = numel(p_ne)/numel(ok_ne);
fprintf('  → taxa de acerto do padrão: %.0f%%\n', 100*taxa);
if taxa < TAXA_MIN
    warning(['Só %.0f%% dos picos de Ne caíram sobre linhas catalogadas. ' ...
             'Quando a identificação está certa, essa taxa passa de 90%%. ' ...
             'Causas prováveis: (a) pico falso corrompendo a cúbica do Hg ' ...
             '(ver NOTA 7); (b) a lâmpada não é de Ne. NÃO grave estes ' ...
             'coeficientes sem investigar.'], 100*taxa);
end

%% ========================================================================
%  6. AJUSTE FINAL
%  ========================================================================
%  BOA PRÁTICA (condicionamento numérico): ajustar diretamente em pixel
%  bruto gera matriz de Vandermonde mal condicionada — com pixels até ~2700
%  e grau 3, cond(V) ≈ 10¹¹, e o MATLAB emite aviso de "poorly conditioned".
%  A forma [p,S,mu] = polyfit(...) centra e escala o preditor internamente,
%  levando cond(V) para ~9. Como o INSTRUMENTO precisa dos coeficientes em
%  pixel bruto, o script converte de volta ao final e VERIFICA que as duas
%  rotas concordam.
pixel = [p_hg(:);  p_ne(:)];
alvo  = [lam_hg(:); lam_ne(:)];
eh_Hg = [true(numel(p_hg),1); false(numel(p_ne),1)];
n     = numel(pixel);

[coef_esc, ~, mu] = polyfit(pixel, alvo, GRAU);
coef_raw = escalado_para_bruto(coef_esc, mu, GRAU);

aviso = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
coef_direto = polyfit(pixel, alvo, GRAU);
warning(aviso);
discrep = max(abs(polyval(coef_raw,pixel) - polyval(coef_direto,pixel)));
fprintf('\nConcordância escalado→bruto vs. ajuste direto: %.2e nm\n', discrep);
assert(discrep < 1e-6, 'Conversão de coeficientes inconsistente.');

a = fliplr(coef_raw);          % [a0 a1 a2 a3], ordem crescente
fprintf('\n=== COEFICIENTES AJUSTADOS (ordem da tela do SpectraSuite) ===\n');
fprintf('Intercept          (a0) = %.6f\n',  a(1));
fprintf('1st Coefficient    (a1) = %.10E\n', a(2));
fprintf('2nd Coefficient    (a2) = %.10E\n', a(3));
fprintf('3rd Coefficient    (a3) = %.10E\n', a(4));

%% ========================================================================
%  7. RESÍDUOS, VALIDAÇÃO E COMPARAÇÃO COM A FÁBRICA
%  ========================================================================
previsto = polyval(coef_raw, pixel);
residuo  = previsto - alvo;
lam_fab  = polyval(fliplr(COEF_FABRICA), pixel);
res_fab  = lam_fab - alvo;

fprintf('\n=== RESÍDUOS ===\n');
fprintf('%10s %12s %12s %10s  %s\n', 'pixel','NIST (nm)','previsto','resíduo','fonte');
for i = 1:n
    if eh_Hg(i), fonte = 'Hg'; else, fonte = 'Ne'; end
    fprintf('%10.3f %12.4f %12.4f %+10.4f  %s\n', ...
            pixel(i), alvo(i), previsto(i), residuo(i), fonte);
end
fprintf('\nRMS total = %.4f nm   (máx |resíduo| = %.4f nm)\n', ...
        rms(residuo), max(abs(residuo)));
fprintf('RMS só Hg = %.4f nm   |   RMS só Ne = %.4f nm\n', ...
        rms(residuo(eh_Hg)), rms(residuo(~eh_Hg)));
fprintf('RMS da calibração de fábrica = %.4f nm\n', rms(res_fab));
fprintf('Melhoria sobre a fábrica: fator %.1f\n', rms(res_fab)/rms(residuo));

% --- Validação cruzada leave-one-out -------------------------------------
%  O RMS acima é medido nos MESMOS pontos usados no ajuste — não distingue
%  um ajuste que capturou a física de um que memorizou os pontos. A LOOCV
%  remove um ponto, reajusta com os demais, prevê o removido e repete.
err_loo = zeros(n,1);
for i = 1:n
    m = true(n,1); m(i) = false;
    [c_i, ~, mu_i] = polyfit(pixel(m), alvo(m), GRAU);
    err_loo(i) = polyval(c_i, (pixel(i)-mu_i(1))/mu_i(2)) - alvo(i);
end
fprintf('\n=== VALIDAÇÃO CRUZADA (LOOCV) ===\n');
fprintf('RMS LOOCV = %.4f nm   (dentro da amostra: %.4f nm)\n', ...
        rms(err_loo), rms(residuo));
razao = rms(err_loo)/rms(residuo);
fprintf('Razão LOOCV/interno = %.2f  ', razao);
if razao < 1.5
    fprintf('→ sem sinal de sobreajuste.\n');
else
    fprintf('→ ATENÇÃO: possível sobreajuste.\n');
end

%% ========================================================================
%  8. GRÁFICOS
%  ========================================================================
figure('Units','normalized','OuterPosition',[0.03 0.08 0.94 0.84]);
cHg = [0 0.45 0.74];  cNe = [0.85 0.33 0.10];  cinza = [0.6 0.6 0.6];

% (a) curva de calibração
subplot(2,3,1); hold on; grid on; box on;
pp = linspace(0, N_PIX, 500);
plot(pp, polyval(coef_raw,pp), 'k-', 'LineWidth',1.4);
plot(pp, polyval(fliplr(COEF_FABRICA),pp), '--', 'Color',cinza,'LineWidth',1.1);
plot(pixel(eh_Hg),  alvo(eh_Hg),  'o','MarkerFaceColor',cHg,'MarkerEdgeColor','k','MarkerSize',6);
plot(pixel(~eh_Hg), alvo(~eh_Hg), 's','MarkerFaceColor',cNe,'MarkerEdgeColor','k','MarkerSize',6);
xlabel('Pixel'); ylabel('Comprimento de onda (nm)');
title('Curva de calibração \lambda(pixel)');
legend({'ajuste','fábrica','Hg','Ne'},'Location','southeast');

% (b) resíduos
subplot(2,3,2); hold on; grid on; box on;
yline(0,'k-');
stem(pixel(eh_Hg),  residuo(eh_Hg), 'o','Color',cHg,'MarkerFaceColor',cHg);
stem(pixel(~eh_Hg), residuo(~eh_Hg),'s','Color',cNe,'MarkerFaceColor',cNe);
yline( rms(residuo),':k'); yline(-rms(residuo),':k');
xlabel('Pixel'); ylabel('Resíduo (nm)');
title(sprintf('Resíduos (RMS = %.4f nm)', rms(residuo)));

% (c) comparação com a fábrica
subplot(2,3,3); hold on; grid on; box on;
yline(0,'k-');
plot(pixel, res_fab, 'v--','Color',cinza,'MarkerFaceColor',cinza);
plot(pixel, residuo, 'o-','Color',cHg,'MarkerFaceColor',cHg,'MarkerSize',4);
xlabel('Pixel'); ylabel('Erro vs. NIST (nm)');
title(sprintf('Fábrica (%.3f nm) vs. ajuste (%.4f nm)', rms(res_fab), rms(residuo)));
legend({'','fábrica','ajuste novo'},'Location','best');

% (d) espectro de Hg com as linhas identificadas
subplot(2,3,4); hold on; grid on; box on;
plot(pix_hg, max(cnt_hg,1), '-', 'Color',cHg, 'LineWidth',0.6);  % max(.,1): log ignora ≤0
for i = find(eh_Hg)'
    xline(pixel(i), ':', sprintf('%.1f', alvo(i)), 'Color','k', ...
          'FontSize',7,'LabelOrientation','horizontal');
end
xlabel('Pixel'); ylabel('Contagens'); title('Espectro de Hg e linhas usadas');
set(gca,'YScale','log');

% (e) espectro de Ne com as linhas identificadas
subplot(2,3,5); hold on; grid on; box on;
plot(pix_ne, max(cnt_ne,1), '-', 'Color',cNe, 'LineWidth',0.6);
for i = find(~eh_Hg)'
    xline(pixel(i), ':', 'Color','k');
end
xlabel('Pixel'); ylabel('Contagens'); title('Espectro de Ne e linhas usadas');
set(gca,'YScale','log');

% (f) dispersão ao longo do array
subplot(2,3,6); grid on; box on;
plot(pp, polyval(polyder(coef_raw), pp), 'k-','LineWidth',1.4);
xlabel('Pixel'); ylabel('d\lambda/dpixel (nm/pixel)');
title('Dispersão (derivada do polinômio)');

%% ========================================================================
%  9. EXPORTAÇÃO
%  ========================================================================
fonte_col = repmat("Ne", n, 1);  fonte_col(eh_Hg) = "Hg";
T = table(pixel, alvo, previsto, residuo, fonte_col, ...
    'VariableNames', {'pixel','lambda_NIST_nm','previsto_nm','residuo_nm','fonte'});
writetable(T, 'calibracao_residuos.csv');

fid = fopen('calibracao_coeficientes.txt','w','n','UTF-8');
fprintf(fid, '# Calibração de comprimento de onda — USB4000 — %s\n', ...
        datestr(now,'yyyy-mm-dd'));
fprintf(fid, '# Ajuste: polinômio grau %d, %d pontos (%d Hg + %d Ne)\n', ...
        GRAU, n, sum(eh_Hg), sum(~eh_Hg));
fprintf(fid, '# RMS interno %.4f nm | LOOCV %.4f nm | fábrica %.4f nm\n', ...
        rms(residuo), rms(err_loo), rms(res_fab));
fprintf(fid, '# Posição dos picos: parábola nos %d pontos centrais (sub-pixel)\n', N_PARAB);
fprintf(fid, '# Referências: NIST ASD ver. 5.12 (DOI 10.18434/T4W30F)\n');
fprintf(fid, 'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n', a);
fprintf(fid, '# COEFICIENTES DE FÁBRICA (backup):\n');
fprintf(fid, '# Intercept\t%.6f\n# 1st\t%.7E\n# 2nd\t%.7E\n# 3rd\t%.7E\n', COEF_FABRICA);
fclose(fid);

fprintf('\nArquivos gravados: calibracao_residuos.csv, calibracao_coeficientes.txt\n');
fprintf(['\nLEMBRETE: estes coeficientes ainda NÃO foram gravados no ' ...
         'instrumento.\nAntes de gravar, confirmar que a lâmpada é ' ...
         'mesmo Hg-Ne e preservar o backup de fábrica.\n']);

%% ========================================================================
%  FUNÇÕES LOCAIS
%  ========================================================================

function [pixel, contagem] = ler_espectro(arquivo, coef_fabrica)
% Lê o arquivo de duas colunas exportado pelo SpectraSuite e devolve o
% índice de pixel (0 a N-1) junto com as contagens.
% O eixo λ do arquivo é a calibração de FÁBRICA avaliada em pixel inteiro —
% por isso o índice da linha é o próprio pixel. A função confere isso.
    M = readmatrix(arquivo);
    contagem = M(:,2);
    pixel = (0:numel(contagem)-1)';
    lam_esperado = polyval(fliplr(coef_fabrica), pixel);
    desvio = max(abs(M(:,1) - lam_esperado));
    if desvio > 0.01
        warning(['O eixo λ do arquivo %s não corresponde à calibração de ' ...
                 'fábrica informada (desvio máx. %.4f nm). Verifique se os ' ...
                 'coeficientes COEF_FABRICA são os do instrumento que ' ...
                 'gerou este espectro.'], arquivo, desvio);
    end
end

function [pos, rejeitado] = picos_subpixel(pixel, contagem, frac_lim, proem_min, sep_min, n_parab, desl_max)
% Detecta picos e devolve a posição sub-pixel de cada um pelo vértice da
% parábola ajustada aos n_parab pontos centrais (ver NOTA 3).
% Rejeita: picos de baixa proeminência (NOTA 7), picos saturados (NOTA 6) e
% picos cujo vértice caia longe demais do pixel de máximo.
    limiar = frac_lim * max(contagem);
    N = numel(contagem);
    h = floor(n_parab/2);

    % --- máximos locais acima do limiar de altura ---
    cand = find(contagem(2:end-1) > limiar & ...
                contagem(2:end-1) >= contagem(1:end-2) & ...
                contagem(2:end-1) >= contagem(3:end)) + 1;

    % --- filtro de PROEMINÊNCIA (NOTA 7) ---
    manter = false(size(cand));
    for k = 1:numel(cand)
        i = cand(k);  hp = contagem(i);
        j = i-1;  minE = hp;
        while j >= 1 && contagem(j) <= hp
            minE = min(minE, contagem(j));  j = j-1;
        end
        j = i+1;  minD = hp;
        while j <= N && contagem(j) <= hp
            minD = min(minD, contagem(j));  j = j+1;
        end
        proem = hp - max(minE, minD);
        manter(k) = (proem >= proem_min * hp);
    end
    cand = cand(manter);

    % --- suprime candidatos vizinhos, mantendo o mais alto de cada grupo ---
    cand = sort(cand);
    manter = true(size(cand));
    for i = 1:numel(cand)-1
        if ~manter(i), continue; end
        j = i+1;
        while j <= numel(cand) && cand(j)-cand(i) < sep_min
            if contagem(cand(j)) > contagem(cand(i))
                manter(i) = false;
            else
                manter(j) = false;
            end
            j = j+1;
        end
    end
    cand = cand(manter);

    % --- posição sub-pixel + rejeições ---
    pos = zeros(numel(cand),1);
    rejeitado = false(numel(cand),1);
    for k = 1:numel(cand)
        i = cand(k);
        % saturação/ceifamento: dois ou mais pixels consecutivos idênticos no topo
        viz = contagem(max(1,i-3):min(N,i+3));
        if sum(abs(viz - contagem(i)) < 1e-6) > 1
            rejeitado(k) = true;
        end
        if i-h < 1 || i+h > N
            pos(k) = pixel(i);  rejeitado(k) = true;  continue;
        end
        x = pixel(i-h:i+h);  y = contagem(i-h:i+h);
        c = polyfit(x, y, 2);
        if c(1) >= 0                       % concavidade errada → inválido
            pos(k) = pixel(i);  rejeitado(k) = true;
        else
            pos(k) = -c(2)/(2*c(1));       % vértice da parábola
            if abs(pos(k) - pixel(i)) > desl_max   % vértice longe do máximo
                rejeitado(k) = true;
            end
        end
    end
end

function [identificado, ok] = casar(lambda_medido, tabela, tolerancia)
% Casa cada λ medido com a linha catalogada mais próxima, dentro da
% tolerância. Devolve o λ catalogado e uma máscara dos que casaram.
    identificado = nan(size(lambda_medido));
    ok = false(size(lambda_medido));
    for i = 1:numel(lambda_medido)
        [d, j] = min(abs(tabela - lambda_medido(i)));
        if d <= tolerancia
            identificado(i) = tabela(j);
            ok(i) = true;
        end
    end
    % impede que duas medidas reivindiquem a mesma linha catalogada:
    % mantém apenas a de menor distância
    [u, ~, g] = unique(identificado(ok));
    idx_ok = find(ok);
    for k = 1:numel(u)
        grupo = idx_ok(g == k);
        if numel(grupo) > 1
            [~, melhor] = min(abs(lambda_medido(grupo) - u(k)));
            descartar = grupo;  descartar(melhor) = [];
            ok(descartar) = false;
            identificado(descartar) = NaN;
        end
    end
end

function coef_raw = escalado_para_bruto(coef_esc, mu, grau)
% Converte os coeficientes do ajuste centrado/escalado (variável
% s = (p - mu(1))/mu(2)) para coeficientes em pixel bruto, expandindo cada
% termo s^g em potências de p por convolução dos binômios.
    coef_raw = zeros(1, grau+1);
    for k = 1:grau+1
        g = grau - k + 1;
        termo = 1;
        for j = 1:g
            termo = conv(termo, [1/mu(2), -mu(1)/mu(2)]);
        end
        coef_raw(end-g:end) = coef_raw(end-g:end) + coef_esc(k) * termo;
    end
end