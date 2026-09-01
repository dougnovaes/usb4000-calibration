% 1. Dados Iniciais (Vetorização em Coluna)
wavelength_nm = [250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, ...
    400, 450, 500, 555, 600, 654.60, 700, 800, 900, 1050, ...
    1150, 1200, 1300, 1540, 1600, 1700, 2000, 2100, 2300, 2400]';

irradiance = [0.136, 0.247, 0.411, 0.648, 0.985, 1.450, 2.057, 2.839, ...
    3.816, 5.052, 6.537, 18.150, 37.220, 62.670, 93.900, ...
    119.900, 148.900, 169.100, 198.900, 208.900, 203.300, ...
    189.200, 180.200, 162.200, 121.400, 112.100, 97.800, ...
    66.000, 58.600, 45.000, 40.300]';

% 2. Mudança de Escala para Estabilidade Numérica
% x: comprimento de onda em micrômetros (um)
% y: irradiância original
x_um = wavelength_nm / 1000;
y_data = irradiance;

% Constante C2 = (h*c)/k expressa em (um * K)
c2 = 14387.77; 

% 3. Definição do Modelo com Correção Empírica de Emissividade
% A emissividade do tungstênio decai com o comprimento de onda.
% O fator (1 - alpha * x) aproxima esse decaimento fisicamente não-cinza.
% T é a temperatura (K) e A é a constante geométrica de escala.
eq_modelo = 'A * (1 - alpha*x) / (x^5 * (exp(14387.77 / (x * T)) - 1))';

% 4. Configuração dos Pesos Estatísticos
% Peso inverso ao quadrado da leitura, minimizando desvios relativos em todo o espectro.
pesos = 1 ./ (y_data.^2);

% 5. Parâmetros de Ajuste e Condições Iniciais
% [A, T, alpha]
chute_inicial = [1e4, 3000, 0.1];
limite_inferior = [0, 1000, -1];
limite_superior = [Inf, 6000, 1];

opcoes_ajuste = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Weights', pesos, ...
    'StartPoint', chute_inicial, ...
    'Lower', limite_inferior, ...
    'Upper', limite_superior);

modelo_planck = fittype(eq_modelo, 'independent', 'x', 'dependent', 'y', 'options', opcoes_ajuste);

% 6. Execução do Ajuste
[ajuste_final, estatisticas_gof] = fit(x_um, y_data, modelo_planck);

% Exibição dos Parâmetros Extraídos
disp('--- Resultados do Ajuste ---');
disp(ajuste_final);
disp('Estatísticas (Goodness of Fit):');
disp(estatisticas_gof);

% 7. Geração do Gráfico e Análise de Resíduos
figure('Position', [100, 100, 900, 450]);

% Subplot Principal: O modelo e os dados
subplot(1, 2, 1);
x_plot = linspace(min(x_um), max(x_um), 500)';
y_plot = feval(ajuste_final, x_plot);

plot(x_um * 1000, y_data, 'ob', 'MarkerFaceColor', 'b', 'MarkerSize', 5); hold on;
plot(x_plot * 1000, y_plot, '-r', 'LineWidth', 1.5);
xlabel('Comprimento de Onda (nm)');
ylabel('Irradiância Espectral (W/cm^3)');
legend('Dados Experimentais', sprintf('Ajuste Planck\nT = %.0f K', ajuste_final.T));
title('Curva de Emissão Ajustada');
grid on;

% Subplot Secundário: Resíduos (Diferença Percentual)
subplot(1, 2, 2);
residuos_pct = 100 * (y_data - feval(ajuste_final, x_um)) ./ y_data;
plot(x_um * 1000, residuos_pct, '-sk', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
xlabel('Comprimento de Onda (nm)');
ylabel('Erro Residual (%)');
yline(0, '--r');
title('Análise de Resíduos do Ajuste');
grid on;

% % 1. Dados Iniciais
% wavelength_nm = [250.00, 260.00, 270.00, 280.00, 290.00, 300.00, 310.00, ...
%     320.00, 330.00, 340.00, 350.00, 400.00, 450.00, 500.00, ...
%     555.00, 600.00, 654.60, 700.00, 800.00, 900.00, 1050.00, ...
%     1150.00, 1200.00, 1300.00, 1540.00, 1600.00, 1700.00, ...
%     2000.00, 2100.00, 2300.00, 2400.00]';
% 
% irradiance_w_cm3 = [0.136, 0.247, 0.411, 0.648, 0.985, 1.450, 2.057, 2.839, ...
%     3.816, 5.052, 6.537, 18.150, 37.220, 62.670, 93.900, ...
%     119.900, 148.900, 169.100, 198.900, 208.900, 203.300, ...
%     189.200, 180.200, 162.200, 121.400, 112.100, 97.800, ...
%     66.000, 58.600, 45.000, 40.300]';
% 
% % 2. Conversões para o Sistema Internacional (SI)
% lambda_m = wavelength_nm * 1e-9;        % de nm para metros
% E_m3 = irradiance_w_cm3 * 1e6;          % de W/cm^3 para W/m^3
% 
% % 3. Constantes Físicas
% h = 6.626e-34; % Constante de Planck (J.s)
% c = 3.00e8;    % Velocidade da luz (m/s)
% k = 1.38e-23;  % Constante de Boltzmann (J/K)
% B = (h * c) / k;
% 
% % 4. Definição do Modelo de Planck para Ajuste
% % params(1) = Fator de amplitude (A)
% % params(2) = Temperatura (T) em Kelvin
% planck_model = @(params, lambda) (params(1) ./ (lambda.^5)) ./ (exp(B ./ (lambda .* params(2))) - 1);
% 
% % 5. Função Custo com Pesos (Resíduos Relativos)
% % Usamos lsqnonlin para minimizar a diferença relativa em vez da diferença absoluta
% cost_func = @(params) (planck_model(params, lambda_m) - E_m3) ./ E_m3;
% 
% % 6. Chutes Iniciais (Initial Guesses) e Limites
% % O pico em ~900nm sugere T entre 3000K e 3200K (típico para 7.90 A D.C.)
% T_guess = 3100;
% A_guess = 1e-15; 
% initial_guess = [A_guess, T_guess];
% 
% lower_bounds = [0, 1000];
% upper_bounds = [Inf, 6000];
% 
% % 7. Execução do Algoritmo de Otimização
% options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-8);
% fitted_params = lsqnonlin(cost_func, initial_guess, lower_bounds, upper_bounds, options);
% 
% A_fit = fitted_params(1);
% T_fit = fitted_params(2);
% 
% % 8. Plotagem dos Resultados
% lambda_plot = linspace(min(lambda_m), max(lambda_m), 500)';
% E_fit_plot = planck_model(fitted_params, lambda_plot);
% 
% figure;
% plot(lambda_m * 1e9, E_m3 / 1e6, 'ob', 'MarkerFaceColor', 'b'); hold on;
% plot(lambda_plot * 1e9, E_fit_plot / 1e6, '-r', 'LineWidth', 2);
% xlabel('Comprimento de Onda (nm)');
% ylabel('Irradiância Espectral (W/cm^3)');
% legend('Dados Experimentais', sprintf('Ajuste Planck (T = %.0f K)', T_fit));
% title('Ajuste da Irradiância Espectral à Lei de Planck');
% grid on;