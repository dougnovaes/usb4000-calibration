%% CALIBRACAO_COMPRIMENTO_DE_ONDA
% Calibração do eixo de comprimento de onda de um espectrômetro Ocean Optics
% USB4000 a partir de lâmpadas de Hg e de Ne, com determinação sub-pixel da
% posição dos picos e ajuste ponderado pela largura de banda medida.
%
% =========================================================================
% SINOPSE
% =========================================================================
%   Ajusta  λ(pixel) = a0 + a1·p + a2·p² + a3·p³  (o polinômio que o
%   firmware do USB4000 aceita) a partir de linhas de emissão catalogadas.
%
%   O que distingue este script de um ajuste ingênuo:
%     1. posição dos picos por parábola sub-pixel, e não pelo pixel de
%        máximo — ganho de ~10× em precisão de posição (NOTA 3);
%     2. filtro de proeminência, que rejeita ondulações de ruído no flanco
%        de linhas fortes (NOTA 7);
%     3. correção do comprimento de onda efetivo de blends não resolvidos
%        (NOTA 8);
%     4. ajuste PONDERADO pela largura de banda de 90 % medida em cada
%        linha, que é simultaneamente uma barra de erro e um detector de
%        blends (NOTA 9).
%
% =========================================================================
% CONTEXTO
% =========================================================================
%   Etapa A do procedimento de calibração absoluta do USB4000, usado como
%   treino antes do mesmo trabalho nos monocromadores Jobin-Yvon THR1000 e
%   McPherson 207. Aplicação final: caracterização espectral de um plasma
%   frio DBD (descarga de barreira dielétrica).
%
%   O detector não mede comprimento de onda — mede apenas qual dos 3648
%   pixels recebeu luz, e quanto. A relação pixel → λ é uma propriedade da
%   montagem óptica daquele instrumento específico (passo da grade,
%   distância focal, ângulo) e precisa ser MEDIDA, apontando fontes de λ
%   conhecido. Como a equação da grade relaciona λ ao ÂNGULO de difração, e
%   não diretamente à posição no detector plano, projetar um leque de
%   ângulos sobre um array linear introduz curvatura residual — daí o
%   polinômio cúbico em vez de uma reta.
%
% =========================================================================
% ENTRADAS
% =========================================================================
%   espectro_hg.txt, espectro_ne.txt
%       Duas colunas separadas por tabulação, exportadas do SpectraSuite:
%       coluna 1 = comprimento de onda segundo a calibração DE FÁBRICA [nm]
%       coluna 2 = contagens (com Electric Dark Correction ativa)
%       Uma linha por pixel, 3648 linhas.
%
%   O script confere que o eixo λ do arquivo corresponde à calibração de
%   fábrica declarada em COEF_FABRICA; se não corresponder, avisa.
%
% =========================================================================
% SAÍDAS
% =========================================================================
%   calibracao_coeficientes.txt  — coeficientes nos três esquemas de ajuste
%   calibracao_residuos.csv      — tabela por linha (pixel, λ, resíduo, banda)
%   figura com 6 painéis de diagnóstico
%
% =========================================================================
% RESULTADO DE REFERÊNCIA (para conferência)
% =========================================================================
%   Ajuste PONDERADO pela banda, 31 linhas (253–744 nm):
%       a0 =  1.7635567179E+02
%       a1 =  2.1606649029E-01
%       a2 = -3.6985346142E-06
%       a3 = -5.5381922989E-10
%       RMS ponderado ....... 0,0339 nm
%       RMS nas linhas limpas 0,0185 nm
%       calibração de fábrica 1,0419 nm   (melhoria de ~31×)
%
%   Ajuste só com as 11 linhas LIMPAS (404–724 nm):
%       a0 =  1.7701734777E+02
%       a1 =  2.1485460168E-01
%       a2 = -3.0222362698E-06
%       a3 = -6.7253730018E-10
%       RMS ................. 0,0133 nm
%
% =========================================================================
% NOTAS METODOLÓGICAS
% =========================================================================
% Cada NOTA abaixo registra um erro cometido durante o desenvolvimento e a
% defesa adotada. Estão aqui porque são o tipo de armadilha que não aparece
% em manual nenhum e custou várias iterações para diagnosticar.
%
% -------------------------------------------------------------------------
% NOTA 1 — Vácuo vs. ar
% -------------------------------------------------------------------------
%   Tabelas de referência frequentemente listam comprimentos de onda de
%   VÁCUO acima de 200 nm sem sinalizar. O USB4000 mede no ar. Ao consultar
%   a NIST ASD, selecione a saída que já vem em ar (opção "Air (200–2000
%   nm)"), eliminando a conversão manual (fator ≈ 1 − 2,8×10⁻⁴ no visível).
%
% -------------------------------------------------------------------------
% NOTA 2 — Não filtrar por intensidade na consulta ao catálogo
% -------------------------------------------------------------------------
%   A tabela curada "Strong Lines" do NIST Handbook é um subconjunto
%   pequeno; usá-la como fonte única deixa a maioria dos picos sem
%   correspondência confiável. Peça a lista completa (campo "Relative
%   intensity minimum" em branco) e filtre depois pelo seu próprio critério.
%
% -------------------------------------------------------------------------
% NOTA 3 — Posição sub-pixel: parábola no topo, não gaussiana no perfil
% -------------------------------------------------------------------------
%   Usar o pixel de máximo limita a precisão a ~±0,1 nm (meio pixel), que é
%   da ordem do espaçamento entre linhas catalogadas — a identificação fica
%   ambígua por construção. A parábola ajustada aos 3 pontos centrais dá a
%   posição com ~1/10 de pixel (≈0,02 nm).
%
%   Por que parábola no topo e NÃO gaussiana no perfil inteiro: os picos do
%   USB4000 têm cauda assimétrica para o vermelho (aberração de coma, típica
%   de Czerny-Turner compacto). Uma gaussiana ajustada ao perfil todo é
%   puxada pela cauda; a parábola nos 3 pontos centrais enxerga apenas o
%   núcleo, onde a assimetria ainda é desprezível.
%
%   Medido neste conjunto de dados: o centroide fica sistematicamente ~0,38
%   pixel (≈0,08 nm) mais para o vermelho que a parábola. É viés
%   sistemático — não cancela com mais medições.
%
%   Efeito prático: a dispersão do desvio (fábrica − NIST) nas linhas de Hg
%   caiu de ~0,20 nm (pixel inteiro) para 0,057 nm (parábola sub-pixel).
%
% -------------------------------------------------------------------------
% NOTA 4 — O limiar do detector de picos do SpectraSuite é global
% -------------------------------------------------------------------------
%   O programa usa uma única linha de limiar para o gráfico inteiro. Com
%   picos de alturas muito diferentes, um limiar bom para o pico mais alto
%   faz os mais baixos reportarem banda de 90 % absurda (centenas de nm)
%   mesmo sendo reais e estreitos. Ao usar a interface, dê zoom numa janela
%   local e reajuste o limiar ali. Este script não depende disso — detecta
%   os picos diretamente do espectro bruto.
%
% -------------------------------------------------------------------------
% NOTA 5 — Backup dos coeficientes de fábrica
% -------------------------------------------------------------------------
%   Gravar novos coeficientes SOBRESCREVE os de fábrica. Anote-os antes
%   (ver COEF_FABRICA).
%
% -------------------------------------------------------------------------
% NOTA 6 — Saturação
% -------------------------------------------------------------------------
%   Um pico ceifado tem topo achatado; a parábola de 3 pontos devolve um
%   vértice sem sentido. O script detecta pixels consecutivos idênticos no
%   topo e rejeita o pico. Se a linha mais forte saturar, reduza o tempo de
%   integração e readquira: perder a linha mais intensa custa caro,
%   sobretudo no UV.
%
% -------------------------------------------------------------------------
% NOTA 7 — Proeminência (filtro indispensável)
% -------------------------------------------------------------------------
%   Um limiar de altura sozinho NÃO basta. Numa versão anterior, uma
%   ondulação de ruído no flanco de subida de uma linha forte (altura 1700
%   contagens, acima do limiar de 1284) foi tomada como pico e casada com a
%   própria linha, corrompendo o ajuste inteiro — a taxa de identificação
%   caiu de 18/19 para 4/25.
%
%   A defesa é a PROEMINÊNCIA: quanto o pico se ergue acima do maior dos
%   dois "colos" que o cercam. Naquela ondulação a proeminência era 31
%   contagens, ou 1,8 % da própria altura; num pico real a razão
%   proeminência/altura fica perto de 100 %. O critério PROEM_MIN = 0,30
%   separa os dois casos com folga e ainda preserva dubletos parcialmente
%   resolvidos (Ne 638,3/640,2 nm tem razão 0,97).
%
% -------------------------------------------------------------------------
% NOTA 8 — Blends não resolvidos: use o λ EFETIVO, não o de catálogo
% -------------------------------------------------------------------------
%   Quando duas ou mais linhas catalogadas caem dentro da largura
%   instrumental (~0,95 nm de FWHM aqui), o detector vê um pico só, e a
%   posição desse pico NÃO é a da componente dominante nem o centroide
%   ponderado pelas intensidades — é o máximo da soma dos perfis.
%
%   Duas linhas de Hg deste conjunto são blends:
%     313,1555 (int. 3000) + 313,1844 (4000)              → efetivo 313,1716
%     365,0158 (9000) + 365,4842 (3000) + 366,2887 (500)  → efetivo 365,1205
%
%   A correção na linha de 365 nm é de +0,105 nm — maior que o RMS inteiro
%   do ajuste. Os valores efetivos foram obtidos simulando a soma dos
%   perfis com a largura instrumental medida (ver blend_efetivo, abaixo).
%
% -------------------------------------------------------------------------
% NOTA 9 — A banda de 90 % como peso E como detector de blends
% -------------------------------------------------------------------------
%   A largura de banda de 90 % (largura do pico a 90 % da altura, acima da
%   linha de base local) cumpre dois papéis de uma vez:
%
%   (a) DETECTOR DE BLENDS, sem consultar catálogo. A menor banda observada
%       no conjunto é 0,37 nm — essa é a resposta instrumental pura, de uma
%       linha isolada. Toda linha com banda acima disso está alargada por
%       alguma coisa. Classificação adotada:
%           razão = banda / banda_mínima
%           ≤1,3×  limpa        1,3–2,5×  alargada
%           2,5–6× blend        >6×       blend severo
%       Esse critério reencontrou, de forma independente, os dois blends de
%       Hg identificados via catálogo na NOTA 8 (313 e 365 aparecem em
%       2,2–2,3×) e sinalizou também a linha de 253,65 nm (2,3×).
%
%   (b) PESO no ajuste. Usa-se w_i = (banda_mínima / banda_i)², ou seja,
%       incerteza de posição proporcional à largura. Isso substitui a
%       decisão binária de incluir/excluir linhas suspeitas por uma
%       gradação contínua, que é estatisticamente mais defensável e
%       preserva a cobertura espectral.
%
%   Caso concreto: a linha de 253,65 nm é a única âncora abaixo de 313 nm,
%   mas seu desvio destoa (+1,56 nm contra 0,93–1,15 nm das demais),
%   provavelmente por autoabsorção — é a linha de ressonância do Hg, e em
%   lâmpadas de baixa pressão o próprio vapor reabsorve o centro do perfil,
%   distorcendo-o. Com peso uniforme ela degrada o RMS de 0,027 para 0,044
%   nm; com peso pela banda ela recebe 19 % do peso máximo, continua
%   ancorando o UV e o RMS ponderado fica em 0,034 nm.
%
% -------------------------------------------------------------------------
% NOTA 10 — Identificação em duas etapas, nunca pico a pico
% -------------------------------------------------------------------------
%   As linhas de Hg são identificadas primeiro (poucas, fortes, bem
%   separadas, correspondência inequívoca). Ajusta-se uma cúbica APENAS ao
%   Hg e extrapola-se para a região do Ne; o padrão INTEIRO de picos de Ne
%   deve então cair sobre linhas catalogadas de uma vez.
%
%   É essa verificação global — muitos picos simultaneamente, não um por
%   vez — que dá confiança na identificação. Casamento pico a pico pelo mais
%   próximo SEMPRE encontra alguma linha (na faixa 580–730 nm há uma linha
%   de Ar catalogada a cada 0,78 nm em média) e por isso mascara erros.
%
%   Foi exatamente assim que se descobriu que a segunda lâmpada era de
%   NEÔNIO e não de argônio: sob a hipótese de argônio, nenhum modelo suave
%   explicava mais que 7 dos 13 picos; sob neônio, 18 de 19 casaram dentro
%   de 0,25 nm sem nenhum parâmetro ajustado. O teste negativo foi igualmente
%   decisivo — as linhas mais fortes de Ar I (696,5 / 706,7 / 750,4 / 763,5
%   / 811,5 nm) estão AUSENTES do espectro.
%
% =========================================================================
% REFERÊNCIAS
% =========================================================================
%   Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2024).
%   NIST Atomic Spectra Database (ver. 5.12). National Institute of
%   Standards and Technology, Gaithersburg, MD.
%   https://physics.nist.gov/asd     DOI: 10.18434/T4W30F
%
%   Ocean Optics. SpectraSuite Installation and Operation Manual.
%
% =========================================================================
% ESTADO
% =========================================================================
%   Os coeficientes calculados aqui AINDA NÃO FORAM GRAVADOS no instrumento.
%   Antes de gravar: confirmar na etiqueta/documentação que a lâmpada é
%   mesmo Hg-Ne, e preservar o backup de fábrica.
%
% Autor: (preencher)                              Licença: MIT
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURAÇÃO
%  ========================================================================
ARQ_HG    = 'espectro_hg.txt';
ARQ_NE    = 'espectro_ne.txt';

GRAU      = 3;      % grau do polinômio (3 = o que o firmware aceita)
N_PARAB   = 3;      % pontos na parábola do topo (NOTA 3)
FRAC_LIM  = 0.02;   % limiar de altura, como fração do pico máximo
PROEM_MIN = 0.30;   % proeminência mínima / altura (NOTA 7)
SEP_MIN   = 5;      % separação mínima entre picos, em pixels
DESL_MAX  = 1.0;    % deslocamento máx. [px] do vértice em relação ao máximo
TOL_HG    = 1.5;    % tolerância [nm] para casar as linhas de Hg
TOL_NE    = 0.30;   % tolerância [nm] para casar as linhas de Ne (NOTA 10)
TAXA_MIN  = 0.60;   % taxa mínima de acerto do padrão de Ne para prosseguir
RAZAO_LIMPA = 1.3;  % banda/banda_mín abaixo da qual a linha é "limpa"
N_PIX     = 3648;

% Coeficientes de fábrica — BACKUP (NOTA 5), ordem crescente [a0 a1 a2 a3]
COEF_FABRICA = [178.085770, 2.1517764E-1, -3.3423826E-6, -5.9958050E-10];

%% ========================================================================
%  2. TABELAS DE REFERÊNCIA (NIST ASD, comprimento de onda no ar)
%  ========================================================================
% Coluna 1: λ a usar no ajuste (já corrigido para blends, NOTA 8)
% Coluna 2: λ de catálogo da componente dominante, só para documentação
LINHAS_HG = [ ...
    253.6521, 253.6521; ...   % ressonância; sujeita a autoabsorção (NOTA 9)
    313.1716, 313.1555; ...   % BLEND 313,1555 + 313,1844
    365.1205, 365.0158; ...   % BLEND 365,0158 + 365,4842 + 366,2887
    404.6565, 404.6565; ...
    435.8335, 435.8335; ...
    546.0750, 546.0750; ...
    576.9610, 576.9610; ...
    579.0670, 579.0670];

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
%  4. DETECÇÃO DE PICOS, POSIÇÃO SUB-PIXEL E LARGURA DE BANDA
%  ========================================================================
[p_hg, bw_hg, rej_hg] = picos(pix_hg, cnt_hg, COEF_FABRICA, ...
                    FRAC_LIM, PROEM_MIN, SEP_MIN, N_PARAB, DESL_MAX);
[p_ne, bw_ne, rej_ne] = picos(pix_ne, cnt_ne, COEF_FABRICA, ...
                    FRAC_LIM, PROEM_MIN, SEP_MIN, N_PARAB, DESL_MAX);

if any(rej_hg), fprintf('AVISO: %d pico(s) rejeitado(s) no Hg.\n', sum(rej_hg)); end
if any(rej_ne), fprintf('AVISO: %d pico(s) rejeitado(s) no Ne.\n', sum(rej_ne)); end
p_hg = p_hg(~rej_hg);  bw_hg = bw_hg(~rej_hg);
p_ne = p_ne(~rej_ne);  bw_ne = bw_ne(~rej_ne);
fprintf('Picos utilizáveis: Hg = %d, Ne = %d\n', numel(p_hg), numel(p_ne));

%% ========================================================================
%  5. IDENTIFICAÇÃO EM DUAS ETAPAS (NOTA 10)
%  ========================================================================
% --- Etapa 1: Hg -------------------------------------------------------
lam_fab = polyval(fliplr(COEF_FABRICA), p_hg);
[id_hg, ok_hg] = casar(lam_fab, LINHAS_HG(:,1), TOL_HG);
p_hg = p_hg(ok_hg);  bw_hg = bw_hg(ok_hg);  lam_hg = id_hg(ok_hg);
fprintf('\nLinhas de Hg identificadas: %d\n', numel(p_hg));

assert(numel(p_hg) >= GRAU+2, ...
    'Linhas de Hg insuficientes (%d) para ancorar o ajuste.', numel(p_hg));

% Verificação de consistência: o desvio (fábrica − NIST) tem de variar de
% forma SUAVE com o pixel, porque tanto a calibração de fábrica quanto a
% curva verdadeira são polinômios suaves. Um ponto fora do padrão denuncia
% identificação errada.
desv = polyval(fliplr(COEF_FABRICA), p_hg) - lam_hg;
fprintf('Desvio fábrica−NIST no Hg: %.3f a %.3f nm (dispersão %.3f nm)\n', ...
        min(desv), max(desv), std(desv));

% --- Etapa 2: Ne, por extrapolação da cúbica do Hg ---------------------
[c_hg, ~, mu_hg] = polyfit(p_hg, lam_hg, GRAU);
lam_extrap = polyval(c_hg, (p_ne - mu_hg(1))/mu_hg(2));
[id_ne, ok_ne] = casar(lam_extrap, LINHAS_NE, TOL_NE);
p_ne = p_ne(ok_ne);  bw_ne = bw_ne(ok_ne);  lam_ne = id_ne(ok_ne);

taxa = numel(p_ne)/numel(ok_ne);
fprintf('Linhas de Ne identificadas: %d de %d picos (%.0f%%, SEM ajuste livre)\n', ...
        numel(p_ne), numel(ok_ne), 100*taxa);
if taxa < TAXA_MIN
    warning(['Só %.0f%% dos picos de Ne caíram sobre linhas catalogadas. ' ...
        'Quando a identificação está certa essa taxa passa de 90%%. Causas ' ...
        'prováveis: (a) pico falso corrompendo a cúbica do Hg (NOTA 7); ' ...
        '(b) a lâmpada não é de Ne (NOTA 10). NÃO grave estes coeficientes ' ...
        'sem investigar.'], 100*taxa);
end

%% ========================================================================
%  6. CLASSIFICAÇÃO PELA BANDA E PESOS (NOTA 9)
%  ========================================================================
pixel = [p_hg(:);  p_ne(:)];
alvo  = [lam_hg(:); lam_ne(:)];
banda = [bw_hg(:); bw_ne(:)];
eh_Hg = [true(numel(p_hg),1); false(numel(p_ne),1)];
[pixel, ord] = sort(pixel);  alvo = alvo(ord);  banda = banda(ord);  eh_Hg = eh_Hg(ord);
n = numel(pixel);

banda_min = min(banda);
razao = banda / banda_min;
peso  = (banda_min ./ banda).^2;         % w = (σ_mín/σ)²
limpa = razao <= RAZAO_LIMPA;

fprintf('\n=== CLASSIFICAÇÃO PELA BANDA DE 90%% ===\n');
fprintf('Banda mínima observada = %.3f nm  →  FWHM instrumental ≈ %.2f nm\n', ...
        banda_min, banda_min/0.39);
fprintf('%10s %8s %7s %8s  %s\n','λ (nm)','banda','razão','peso','classe');
for i = 1:n
    if     razao(i) <= RAZAO_LIMPA, cls = 'limpa';
    elseif razao(i) <= 2.5,         cls = 'alargada';
    elseif razao(i) <= 6.0,         cls = 'blend';
    else,                           cls = 'BLEND SEVERO';
    end
    fprintf('%10.4f %8.3f %6.1fx %8.3f  %s\n', alvo(i), banda(i), razao(i), peso(i), cls);
end
fprintf('Linhas limpas: %d de %d\n', sum(limpa), n);

%% ========================================================================
%  7. AJUSTE — TRÊS ESQUEMAS, PARA COMPARAÇÃO
%  ========================================================================
% BOA PRÁTICA (condicionamento): ajustar em pixel bruto gera matriz de
% Vandermonde mal condicionada (cond ≈ 10¹¹ com grau 3 e pixels até 2800).
% Centra-se e escala-se o preditor internamente (cond ≈ 9) e converte-se de
% volta ao final, porque o instrumento precisa dos coeficientes em pixel
% bruto. A conversão é verificada por asserção.
[a_pond,  r_pond,  rms_pond]  = ajuste(pixel, alvo, peso,          GRAU);
[a_unif,  r_unif,  rms_unif]  = ajuste(pixel, alvo, ones(n,1),     GRAU);
[a_limpa, r_limpa, rms_limpa] = ajuste(pixel(limpa), alvo(limpa), ...
                                       ones(sum(limpa),1), GRAU);

lam_fab_todos = polyval(fliplr(COEF_FABRICA), pixel);
rms_fab = rms(lam_fab_todos - alvo);

fprintf('\n=== COMPARAÇÃO DOS ESQUEMAS DE AJUSTE ===\n');
fprintf('%-34s %4s %12s %12s\n','esquema','n','RMS (nm)','faixa (nm)');
fprintf('%-34s %4d %12.4f %6.0f–%.0f\n','ponderado pela banda (NOTA 9)', ...
        n, rms_pond, min(alvo), max(alvo));
fprintf('%-34s %4d %12.4f %6.0f–%.0f\n','uniforme', ...
        n, rms_unif, min(alvo), max(alvo));
fprintf('%-34s %4d %12.4f %6.0f–%.0f\n','só as linhas limpas', ...
        sum(limpa), rms_limpa, min(alvo(limpa)), max(alvo(limpa)));
fprintf('%-34s %4d %12.4f\n','calibração de fábrica', n, rms_fab);
fprintf('\nRMS do ajuste ponderado avaliado só nas limpas: %.4f nm\n', ...
        rms(r_pond(limpa)));

% Escolha recomendada: ponderado (cobertura completa, cada ponto pesando o
% que vale). Trocar aqui se preferir outro esquema.
a = a_pond;  residuo = r_pond;
fprintf('\n=== COEFICIENTES ADOTADOS (ponderado) ===\n');
fprintf('Intercept          (a0) = %.6f\n',  a(1));
fprintf('1st Coefficient    (a1) = %.10E\n', a(2));
fprintf('2nd Coefficient    (a2) = %.10E\n', a(3));
fprintf('3rd Coefficient    (a3) = %.10E\n', a(4));

%% ========================================================================
%  8. VALIDAÇÃO CRUZADA E ALAVANCAGEM
%  ========================================================================
% O RMS acima é medido nos MESMOS pontos usados no ajuste. A LOOCV remove um
% ponto, reajusta e prevê o removido. ATENÇÃO ao interpretar: uma razão
% LOOCV/interno alta nem sempre é sobreajuste — pode ser ALAVANCAGEM, quando
% poucos pontos isolados ancoram uma extremidade da faixa. Compare a razão
% com a alavancagem máxima antes de concluir.
err_loo = zeros(n,1);
for i = 1:n
    m = true(n,1); m(i) = false;
    ai = ajuste(pixel(m), alvo(m), peso(m), GRAU);
    err_loo(i) = polyval(fliplr(ai), pixel(i)) - alvo(i);
end
mu = mean(pixel); sd = std(pixel);
X  = vander_norm((pixel-mu)/sd, GRAU);
H  = diag(X*pinv(X'*X)*X');

fprintf('\n=== VALIDAÇÃO CRUZADA ===\n');
fprintf('RMS LOOCV = %.4f nm (dentro da amostra: %.4f nm), razão = %.2f\n', ...
        rms(err_loo), rms(residuo), rms(err_loo)/rms(residuo));
fprintf('Alavancagem: máxima = %.3f (λ = %.1f nm), média esperada = %.3f\n', ...
        max(H), alvo(H==max(H)), (GRAU+1)/n);
if rms(err_loo)/rms(residuo) > 1.5 && max(H) > 3*(GRAU+1)/n
    fprintf(['  → razão alta EXPLICADA por alavancagem: a amostragem é\n' ...
             '    desequilibrada (poucos pontos ancorando uma ponta da\n' ...
             '    faixa). Não é sobreajuste. Para reduzir, acrescente\n' ...
             '    linhas de referência na região pouco povoada.\n']);
end

%% ========================================================================
%  9. GRÁFICOS
%  ========================================================================
figure('Units','normalized','OuterPosition',[0.03 0.08 0.94 0.84]);
cHg=[0 0.45 0.74]; cNe=[0.85 0.33 0.10]; cinza=[0.6 0.6 0.6];
pp = linspace(0, N_PIX, 500);

subplot(2,3,1); hold on; grid on; box on;
plot(pp, polyval(fliplr(a),pp), 'k-','LineWidth',1.4);
plot(pp, polyval(fliplr(COEF_FABRICA),pp), '--','Color',cinza,'LineWidth',1.1);
plot(pixel(eh_Hg), alvo(eh_Hg), 'o','MarkerFaceColor',cHg,'MarkerEdgeColor','k','MarkerSize',6);
plot(pixel(~eh_Hg),alvo(~eh_Hg),'s','MarkerFaceColor',cNe,'MarkerEdgeColor','k','MarkerSize',6);
xlabel('Pixel'); ylabel('Comprimento de onda (nm)');
title('Curva de calibração \lambda(pixel)');
legend({'ajuste','fábrica','Hg','Ne'},'Location','southeast');

subplot(2,3,2); hold on; grid on; box on;
yline(0,'k-');
sz = 20 + 120*peso;                       % tamanho ∝ peso
scatter(pixel(eh_Hg),  residuo(eh_Hg),  sz(eh_Hg),  cHg,'filled','MarkerEdgeColor','k');
scatter(pixel(~eh_Hg), residuo(~eh_Hg), sz(~eh_Hg), cNe,'filled','MarkerEdgeColor','k');
yline( rms_pond,':k'); yline(-rms_pond,':k');
xlabel('Pixel'); ylabel('Resíduo (nm)');
title(sprintf('Resíduos (RMS pond. = %.4f nm); área \\propto peso', rms_pond));

subplot(2,3,3); hold on; grid on; box on;
plot(pixel, razao, 'ko-','MarkerFaceColor','k','MarkerSize',4);
yline(RAZAO_LIMPA,'--','limpa','Color',[0 .6 0]);
yline(2.5,'--','blend','Color',[.8 .5 0]);
yline(6.0,'--','severo','Color',[.8 0 0]);
set(gca,'YScale','log'); xlabel('Pixel'); ylabel('banda / banda_{mín}');
title('Diagnóstico de blend pela banda de 90%');

subplot(2,3,4); hold on; grid on; box on;
plot(pix_hg, max(cnt_hg,1), '-','Color',cHg,'LineWidth',0.6);
for i = find(eh_Hg)', xline(pixel(i),':','Color','k'); end
set(gca,'YScale','log'); xlabel('Pixel'); ylabel('Contagens');
title('Espectro de Hg e linhas usadas');

subplot(2,3,5); hold on; grid on; box on;
plot(pix_ne, max(cnt_ne,1), '-','Color',cNe,'LineWidth',0.6);
for i = find(~eh_Hg)', xline(pixel(i),':','Color','k'); end
set(gca,'YScale','log'); xlabel('Pixel'); ylabel('Contagens');
title('Espectro de Ne e linhas usadas');

subplot(2,3,6); grid on; box on;
plot(pp, polyval(polyder(fliplr(a)), pp), 'k-','LineWidth',1.4);
xlabel('Pixel'); ylabel('d\lambda/dpixel (nm/pixel)');
title('Dispersão (derivada do polinômio)');

%% ========================================================================
%  10. EXPORTAÇÃO
%  ========================================================================
fonte = repmat("Ne", n, 1); fonte(eh_Hg) = "Hg";
classe = repmat("BLEND SEVERO", n, 1);
classe(razao<=6.0)         = "blend";
classe(razao<=2.5)         = "alargada";
classe(razao<=RAZAO_LIMPA) = "limpa";
T = table(pixel, alvo, polyval(fliplr(a),pixel), residuo, banda, razao, peso, fonte, classe, ...
    'VariableNames', {'pixel','lambda_NIST_nm','previsto_nm','residuo_nm', ...
                      'banda90_nm','razao_banda','peso','fonte','classe'});
writetable(T,'calibracao_residuos.csv');

fid = fopen('calibracao_coeficientes.txt','w','n','UTF-8');
fprintf(fid,'# Calibração de comprimento de onda — Ocean Optics USB4000\n');
fprintf(fid,'# Gerado em %s\n#\n', datestr(now,'yyyy-mm-dd HH:MM'));
fprintf(fid,'# Fontes: lâmpadas de Hg e de Ne; referências NIST ASD ver. 5.12\n');
fprintf(fid,'# (DOI 10.18434/T4W30F), comprimento de onda no ar.\n');
fprintf(fid,'# Posição dos picos: parábola sub-pixel em %d pontos.\n', N_PARAB);
fprintf(fid,'# Pesos: w = (banda_mín/banda_90%%)^2\n#\n');
fprintf(fid,'[PONDERADO]  n=%d  RMS=%.4f nm  faixa=%.1f-%.1f nm\n', ...
        n, rms_pond, min(alvo), max(alvo));
fprintf(fid,'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n\n', a_pond);
fprintf(fid,'[SO_LIMPAS]  n=%d  RMS=%.4f nm  faixa=%.1f-%.1f nm\n', ...
        sum(limpa), rms_limpa, min(alvo(limpa)), max(alvo(limpa)));
fprintf(fid,'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n\n', a_limpa);
fprintf(fid,'[UNIFORME]   n=%d  RMS=%.4f nm\n', n, rms_unif);
fprintf(fid,'Intercept\t%.10f\n1st\t%.10E\n2nd\t%.10E\n3rd\t%.10E\n\n', a_unif);
fprintf(fid,'[FABRICA_BACKUP]  RMS=%.4f nm\n', rms_fab);
fprintf(fid,'Intercept\t%.6f\n1st\t%.7E\n2nd\t%.7E\n3rd\t%.7E\n', COEF_FABRICA);
fclose(fid);

fprintf('\nArquivos gravados: calibracao_residuos.csv, calibracao_coeficientes.txt\n');
fprintf(['\nLEMBRETE: estes coeficientes ainda NÃO foram gravados no ' ...
    'instrumento.\nAntes de gravar, confirmar que a lâmpada é Hg-Ne e ' ...
    'preservar o backup de fábrica.\n']);

%% ========================================================================
%  FUNÇÕES LOCAIS
%  ========================================================================

function [pixel, contagem] = ler_espectro(arquivo, coef_fabrica)
% Lê o arquivo de duas colunas do SpectraSuite. O eixo λ do arquivo é a
% calibração de FÁBRICA avaliada em pixel inteiro, de modo que o índice da
% linha é o próprio pixel. A função confere isso e avisa se não bater.
    M = readmatrix(arquivo);
    contagem = M(:,2);
    pixel = (0:numel(contagem)-1)';
    desvio = max(abs(M(:,1) - polyval(fliplr(coef_fabrica), pixel)));
    if desvio > 0.01
        warning(['O eixo λ de %s não corresponde à calibração de fábrica ' ...
            'informada (desvio máx. %.4f nm). Verifique se COEF_FABRICA é ' ...
            'do instrumento que gerou este espectro.'], arquivo, desvio);
    end
end

function [pos, banda, rejeitado] = picos(pixel, contagem, coef_fab, ...
                              frac_lim, proem_min, sep_min, n_parab, desl_max)
% Detecta picos e devolve, para cada um: a posição sub-pixel (vértice da
% parábola nos n_parab pontos centrais, NOTA 3) e a largura de banda de 90 %
% em nm (largura a 90 % da altura acima da base local, NOTA 9).
% Rejeita picos de baixa proeminência (NOTA 7), saturados (NOTA 6) e cujo
% vértice caia longe do máximo.
    limiar = frac_lim * max(contagem);
    N = numel(contagem);
    h = floor(n_parab/2);

    cand = find(contagem(2:end-1) > limiar & ...
                contagem(2:end-1) >= contagem(1:end-2) & ...
                contagem(2:end-1) >= contagem(3:end)) + 1;

    % --- filtro de PROEMINÊNCIA (NOTA 7) ---
    manter = false(size(cand));
    for k = 1:numel(cand)
        i = cand(k);  hp = contagem(i);
        j = i-1; minE = hp;
        while j >= 1 && contagem(j) <= hp, minE = min(minE,contagem(j)); j = j-1; end
        j = i+1; minD = hp;
        while j <= N && contagem(j) <= hp, minD = min(minD,contagem(j)); j = j+1; end
        manter(k) = (hp - max(minE,minD)) >= proem_min*hp;
    end
    cand = cand(manter);

    % --- suprime vizinhos, mantendo o mais alto de cada grupo ---
    cand = sort(cand);  manter = true(size(cand));
    for i = 1:numel(cand)-1
        if ~manter(i), continue; end
        j = i+1;
        while j <= numel(cand) && cand(j)-cand(i) < sep_min
            if contagem(cand(j)) > contagem(cand(i)), manter(i) = false;
            else,                                     manter(j) = false; end
            j = j+1;
        end
    end
    cand = cand(manter);

    pos = zeros(numel(cand),1);  banda = nan(numel(cand),1);
    rejeitado = false(numel(cand),1);
    for k = 1:numel(cand)
        i = cand(k);
        % saturação: dois ou mais pixels idênticos no topo (NOTA 6)
        viz = contagem(max(1,i-3):min(N,i+3));
        if sum(abs(viz - contagem(i)) < 1e-6) > 1, rejeitado(k) = true; end
        if i-h < 1 || i+h > N
            pos(k) = pixel(i); rejeitado(k) = true; continue;
        end
        c = polyfit(pixel(i-h:i+h), contagem(i-h:i+h), 2);
        if c(1) >= 0
            pos(k) = pixel(i); rejeitado(k) = true;
        else
            pos(k) = -c(2)/(2*c(1));
            if abs(pos(k)-pixel(i)) > desl_max, rejeitado(k) = true; end
        end
        banda(k) = banda90(pixel, contagem, i, coef_fab);
    end
end

function b = banda90(pixel, contagem, i, coef_fab)
% Largura de banda de 90 %: largura do pico, em nm, no nível
% base + 0,90·(topo − base). A base é o mínimo local de cada lado, o que
% torna a medida robusta a fundo inclinado. Interpola linearmente entre
% pixels para não quantizar a largura.
    N = numel(contagem);  topo = contagem(i);
    j = i; while j > 1  && contagem(j-1) <= contagem(j), j = j-1; end,  baseE = contagem(j);
    j = i; while j < N  && contagem(j+1) <= contagem(j), j = j+1; end,  baseD = contagem(j);
    base = max(baseE, baseD);
    nivel = base + 0.90*(topo - base);
    if topo <= base, b = NaN; return; end

    j = i;  while j > 1 && contagem(j) > nivel, j = j-1; end
    if j < 1 || contagem(j+1) == contagem(j), xE = pixel(max(j,1));
    else, xE = pixel(j) + (nivel-contagem(j))/(contagem(j+1)-contagem(j)); end

    j = i;  while j < N && contagem(j) > nivel, j = j+1; end
    if j > N || contagem(j) == contagem(j-1), xD = pixel(min(j,N));
    else, xD = pixel(j-1) + (contagem(j-1)-nivel)/(contagem(j-1)-contagem(j)); end

    b = abs(polyval(fliplr(coef_fab), xD) - polyval(fliplr(coef_fab), xE));
end

function [identificado, ok] = casar(lambda_medido, tabela, tolerancia)
% Casa cada λ medido com a linha catalogada mais próxima dentro da
% tolerância, impedindo que duas medidas reivindiquem a mesma linha.
    identificado = nan(size(lambda_medido));
    ok = false(size(lambda_medido));
    for i = 1:numel(lambda_medido)
        [d,j] = min(abs(tabela - lambda_medido(i)));
        if d <= tolerancia, identificado(i) = tabela(j); ok(i) = true; end
    end
    idx = find(ok);
    [u,~,g] = unique(identificado(ok));
    for k = 1:numel(u)
        grupo = idx(g==k);
        if numel(grupo) > 1
            [~,melhor] = min(abs(lambda_medido(grupo)-u(k)));
            fora = grupo; fora(melhor) = [];
            ok(fora) = false;  identificado(fora) = NaN;
        end
    end
end

function [a, r, rms_p] = ajuste(pixel, alvo, peso, grau)
% Mínimos quadrados PONDERADOS com centragem/escalamento do preditor, e
% conversão dos coeficientes de volta para pixel bruto (verificada).
    mu = mean(pixel);  sd = std(pixel);
    s  = (pixel-mu)/sd;
    X  = vander_norm(s, grau);
    W  = diag(peso(:));
    c  = (X'*W*X) \ (X'*W*alvo(:));      % c em ordem decrescente
    r  = X*c - alvo(:);
    rms_p = sqrt(sum(peso(:).*r.^2)/sum(peso(:)));

    % converte para pixel bruto expandindo s^g em potências de p
    coef_raw = zeros(1, grau+1);
    for k = 1:grau+1
        g = grau-k+1;  termo = 1;
        for j = 1:g, termo = conv(termo, [1/sd, -mu/sd]); end
        coef_raw(end-g:end) = coef_raw(end-g:end) + c(k)*termo;
    end
    a = fliplr(coef_raw);                % [a0 a1 a2 a3]
    assert(max(abs(polyval(coef_raw,pixel) - X*c)) < 1e-6, ...
           'Conversão escalado→bruto inconsistente.');
end

function X = vander_norm(s, grau)
% Matriz de Vandermonde em ordem decrescente de potências.
    X = zeros(numel(s), grau+1);
    for k = 0:grau, X(:,k+1) = s(:).^(grau-k); end
end