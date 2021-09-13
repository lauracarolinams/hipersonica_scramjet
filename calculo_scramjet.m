% Laura Monteiro
% Trabalho de conclus�o de curso
clc
clear all

%% Dados da atmosfera - altitude = 25 km

MM = 0.0289644;                                 % Massa molecular do ar [kg/kmol]
gamma = 1.4;
const_univ_gas = 8.31432;                       % Constante universal dos gases
const_gas_ar = const_univ_gas/MM;               % Constante dos gases para o ar [J/kg�C]
cp_ar = (gamma*const_gas_ar)/(gamma-1);         % Cp do ar

Temp_station_0 = 221.5520647;                   % Temperatura do escoamento livre [K]
rho_station_0 = 0.040083787;                    % Densidade do escoamento livre [kg/m�]
Pres_station_0 = 2549.216627;                   % Press�o do escoamento livre [Pa]
velocidade_som_station_0 = (gamma*const_gas_ar*Temp_station_0)^(1/2);            % Velocidade do som [m/s]
Mach_station_0 = 6.3;                           % Mach do escoamento livre
velocidade_station_0 = Mach_station_0*velocidade_som_station_0;                         % Velocidade do escoamento livre [m/s]

%% Entrada da c�mara de combust�o

f_st = 3/103;                                  % f para combust�o estequiom�trica
const_gas_H2 = 4124;                           % Constante dos gases para o H2 [J/kg�C]
gamma_H2 = 1.405;                    
cp_H2 = (gamma_H2*const_gas_H2)/(gamma_H2-1);  % Cp do H2
Temperatura_injecao_H2 = 300;                  % Temperatura de inje��o do H2 [K]
Temperatura_ignicao_H2 = 845.15;               % Temperatura de igni��o do H2 [K]

% A temperatura ideal � a temperatura necess�ria para que haja a combust�o
% do H2
Temperatura_station_3_ideal = f_st*(cp_H2/cp_ar)*(Temperatura_ignicao_H2-Temperatura_injecao_H2)+Temperatura_ignicao_H2; % Temperatura na entrada da c�mara de combust�o
Mach_station_3_ideal = ((2/(gamma-1))*(((Temp_station_0/Temperatura_station_3_ideal)*(1+((gamma-1)/2)*Mach_station_0^2))-1))^(1/2); % Mach na entrada da c�mara de combust�o
velocidade_som_station_3_ideal = (gamma*const_gas_ar*Temperatura_station_3_ideal)^(1/2);                  % Velocidade do som [m/s]
velocidade_station_3_ideal = Mach_station_3_ideal*velocidade_som_station_3_ideal;                         % Velocidade do escoamento na esta��o 3 [m/s]


%% Rampas
num_rampas = 3;                              % N�mero de rampas na se��o de compress�o

erro_max_Mach_station_3 = 10^-5;
erro_max_Temperatura_station_3 = 10^-5;

% Intervalo inferior
a1 = 10;
beta_rampas(1) = a1;
rampas(num_rampas, beta_rampas, gamma, const_gas_ar, Mach_station_0, Pres_station_0, rho_station_0, Temp_station_0);
[~, ~, ~, ~, ~, beta_rampas, ~, ~, ~, ~, ~, ~, ~, Mach_station_3, ~, ~, Temp_station_3, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = rampas(num_rampas, beta_rampas, gamma, const_gas_ar, Mach_station_0, Pres_station_0, rho_station_0, Temp_station_0);
erro_1_Mach_station_3 = Mach_station_3 - Mach_station_3_ideal;
erro_1_Temperatura_station_3 = Temp_station_3 - Temperatura_station_3_ideal;

% Intervalo superior
a2 = 15;
beta_rampas(1) = a2;
rampas(num_rampas, beta_rampas, gamma, const_gas_ar, Mach_station_0, Pres_station_0, rho_station_0, Temp_station_0);
[~, ~, ~, ~, ~, beta_rampas, ~, ~, ~, ~, ~, ~, ~, Mach_station_3, ~, ~, Temp_station_3, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = rampas(num_rampas, beta_rampas, gamma, const_gas_ar, Mach_station_0, Pres_station_0, rho_station_0, Temp_station_0);
erro_2_Mach_station_3 = Mach_station_3 - Mach_station_3_ideal;
erro_2_Temperatura_station_3 = Temp_station_3 - Temperatura_station_3_ideal;

erro_Mach_station_3 = max(erro_1_Mach_station_3, erro_2_Mach_station_3);
erro_Temperatura_station_3 = max(erro_1_Temperatura_station_3, erro_2_Temperatura_station_3);


if erro_1_Mach_station_3*erro_2_Mach_station_3 > 0
    disp('Necess�rio alterar valores m�ximo ou m�nimo para beta da primeira rampa')
    
else
    while abs(erro_Mach_station_3) > erro_max_Mach_station_3 && abs(erro_Temperatura_station_3) > erro_max_Temperatura_station_3
        x =(a1 + a2)/2;
        beta_rampas(1) = x;
        rampas(num_rampas, beta_rampas, gamma, const_gas_ar, Mach_station_0, Pres_station_0, rho_station_0, Temp_station_0);
        [Mach_normal_rampas, Temp_total_station_3, Temp_total_comp, Mach_rampas, teta_rampas, beta_rampas, Pres_comp, rho_comp, Temp_comp, velocidade_escoamento, teta_refletido, tg_beta_refletido, Mach_normal_station_3, Mach_station_3, Pres_station_3, rho_station_3, Temp_station_3, velocidade_escoamento_station_3, razao_pressao, razao_rho, razao_temperatura, velocidade_som, razao_pressao_station_3, razao_rho_station_3, razao_temperatura_station_3, velocidade_som_station_3, beta_refletido] = rampas(num_rampas, beta_rampas, gamma, const_gas_ar, Mach_station_0, Pres_station_0, rho_station_0, Temp_station_0);
        erro_Mach_station_3 = Mach_station_3 - Mach_station_3_ideal;
        erro_Temperatura_station_3 = Temp_station_3 - Temperatura_station_3_ideal;

        if erro_Mach_station_3*erro_1_Mach_station_3 > 0
            a1 = x;
        else
            a2 = x;
        end
    end
 beta_rampas(1) = x;
end


%% Dimens�es 
altura_station_0 = 0.38;         % Altura total
altura_station_3 = (rho_station_0*velocidade_station_0*altura_station_0)/(rho_station_3*velocidade_escoamento_station_3);
largura_scramjet = 0.1028;       % Largura do scramjet [m]


%% Se��o de expans�o power-off

Mach_station_4 = Mach_station_3;
Pres_station_4 = Pres_station_3;
Temp_station_4 = Temp_station_3;
rho_station_4 = rho_station_3;
velocidade_station_4 = velocidade_escoamento_station_3;

erro_max_altura_station_10 = 10^-5;

% Intervalo inferior
a1 = 0;
teta_expansao = a1;
expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
[~, ~, ~, ~, altura_station_10, altura_station_10_ideal, ~, ~, ~, ~] = expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
erro_1_altura_station_10 = altura_station_10 - altura_station_10_ideal;

% Intervalo superior
a2 = 20;
teta_expansao = a2;
expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
[~, ~, ~, ~, altura_station_10, altura_station_10_ideal, ~, ~, ~, ~] = expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
erro_2_altura_station_10 = altura_station_10 - altura_station_10_ideal;

erro_altura_station_10 = max(erro_1_altura_station_10, erro_2_altura_station_10);


if erro_1_altura_station_10*erro_2_altura_station_10 > 0
    disp('Necess�rio alterar valores m�ximo ou m�nimo para teta da expans�o power-off')
    
else
    while abs(erro_altura_station_10) > erro_max_altura_station_10
        y =(a1 + a2)/2;
        teta_expansao = y;
        expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
        [Mach_prandtl_meyer, Temp_station_10, Pres_station_10, Mach_station_10, altura_station_10, altura_station_10_ideal, rho_station_10, velocidade_som_station_10, velocidade_station_10, altura_prandtl_meyer] = expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
        erro_altura_station_10 = altura_station_10 - altura_station_10_ideal;
        
        if erro_altura_station_10*erro_1_altura_station_10 > 0
            a1 = y;
        else
            a2 = y;
        end
    end
teta_expansao = y;
end

Mach_prandtl_meyer_power_off = Mach_prandtl_meyer;
Temp_station_10_power_off = Temp_station_10;
Mach_station_10_power_off = Mach_station_10;
altura_station_10_power_off = altura_station_10;
rho_station_10_power_off = rho_station_10;
velocidade_som_station_10_power_off = velocidade_som_station_10;
velocidade_escoamento_station_10_power_off = velocidade_station_10;
altura_prandtl_meyer_power_off = altura_prandtl_meyer;
teta_expansao_power_off = teta_expansao;


%% C�mara de combust�o

Mach_station_4 = 1.25;
vazao_massica_ar = rho_station_3*velocidade_escoamento_station_3*(altura_station_3/2)*largura_scramjet;
vazao_massica_H2 = f_st*vazao_massica_ar;        
calor_reacao_H2 = 119954*(10^3);                       % Calor da rea��o [J/kg]
adicao_calor_estequiometrica = calor_reacao_H2*vazao_massica_H2;       % [J/s]

razao_pressao_station_4 = (1 + gamma*Mach_station_3^2)/(1 + gamma*Mach_station_4^2);
Pres_station_4 = razao_pressao_station_4*Pres_station_3;

razao_rho_station_4 = ((1 + gamma*Mach_station_4^2)/(1 + gamma*Mach_station_3^2))*((Mach_station_3/Mach_station_4)^2);
rho_station_4 = razao_rho_station_4*rho_station_3;

razao_temperatura_station_4 = (((1 + gamma*Mach_station_3^2)/(1 + gamma*Mach_station_4^2))^2)*((Mach_station_4/Mach_station_3)^2);
Temp_station_4 = razao_temperatura_station_4*Temp_station_3;

velocidade_som_station_4 = sqrt(gamma*const_gas_ar*Temp_station_4);
velocidade_station_4 = velocidade_som_station_4*Mach_station_4;

Temp_total_station_4 = (1 + ((gamma-1)/2)*Mach_station_4^2)*Temp_station_4;
adicao_calor = cp_ar*(Temp_total_station_4 - Temp_total_station_3);


%% Se��o de expans�o power-on

erro_max_altura_station_10 = 10^-5;

% Intervalo inferior
a1 = 0;
teta_expansao = a1;
expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
[~, ~, ~, ~, altura_station_10, altura_station_10_ideal, ~, ~, ~, ~] = expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
erro_1_altura_station_10 = altura_station_10 - altura_station_10_ideal;

% Intervalo superior
a2 = 40;
teta_expansao = a2;
expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
[~, ~, ~, ~, altura_station_10, altura_station_10_ideal, ~, ~, ~, ~] = expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
erro_2_altura_station_10 = altura_station_10 - altura_station_10_ideal;

erro_altura_station_10 = max(erro_1_altura_station_10, erro_2_altura_station_10);


if erro_1_altura_station_10*erro_2_altura_station_10 > 0
    disp('Necess�rio alterar valores m�ximo ou m�nimo para teta da expans�o')
    
else
    while abs(erro_altura_station_10) > erro_max_altura_station_10
        y =(a1 + a2)/2;
        teta_expansao = y;
        expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
        [Mach_prandtl_meyer, Temp_station_10, Pres_station_10, Mach_station_10, altura_station_10, altura_station_10_ideal, rho_station_10, velocidade_som_station_10, velocidade_station_10, altura_prandtl_meyer] = expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4);
        erro_altura_station_10 = altura_station_10 - altura_station_10_ideal;
        
        if erro_altura_station_10*erro_1_altura_station_10 > 0
            a1 = y;
        else
            a2 = y;
        end
    end
teta_expansao = y;
end


%% Empuxo n�o instalado

empuxo_power_on = vazao_massica_ar*(velocidade_station_10 - velocidade_station_0) + (Pres_station_10 - Pres_station_0)*altura_station_10*largura_scramjet;