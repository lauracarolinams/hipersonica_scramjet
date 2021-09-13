function[Mach_normal_rampas, Temp_total_station_3, Temp_total_comp, Mach_rampas, teta_rampas, beta_rampas, Pres_comp, rho_comp, Temp_comp, velocidade_escoamento, teta_refletido, tg_beta_refletido, Mach_normal_station_3, Mach_station_3, Pres_station_3, rho_station_3, Temp_station_3, velocidade_escoamento_station_3, razao_pressao, razao_rho, razao_temperatura, velocidade_som, razao_pressao_station_3, razao_rho_station_3, razao_temperatura_station_3, velocidade_som_station_3, beta_refletido] = rampas(num_rampas, beta_rampas, gamma, const_gas_ar, Mach_station_0, Pres_station_0, rho_station_0, Temp_station_0)

%% Para a primeira rampa

Mach_normal_rampas(1) = Mach_station_0*sind(beta_rampas(1));               % Intensidade do primeiro choque
tg_teta(1) = 2*cotd(beta_rampas(1))*((Mach_normal_rampas(1)^2-1)/(Mach_station_0^2*(gamma+cosd(2*beta_rampas(1)))+2));
teta_rampas(1) = atand(tg_teta(1));           % Ângulo da primeira rampa
Mach_rampas(1) = (sqrt((Mach_normal_rampas(1)^2+(2/(gamma-1)))/(((2*gamma)/(gamma-1))*Mach_normal_rampas(1)^2-1)))/sind(beta_rampas(1)-teta_rampas(1));
sen_beta(1) = sind(beta_rampas(1));

razao_pressao(1) = 1+((2*gamma)/(gamma+1))*(Mach_normal_rampas(1)^2-1);
Pres_comp(1) = razao_pressao(1)*Pres_station_0;

razao_rho(1) = ((gamma+1)*Mach_normal_rampas(1)^2)/((gamma-1)*Mach_normal_rampas(1)^2+2);
rho_comp(1) = razao_rho(1)*rho_station_0;

razao_temperatura(1) = razao_pressao(1)/razao_rho(1);
Temp_comp(1) = razao_temperatura(1)*Temp_station_0;

velocidade_som(1) = sqrt(gamma*const_gas_ar*Temp_comp(1));
velocidade_escoamento(1) = velocidade_som(1)*Mach_rampas(1);

Temp_total_comp(1) = (1 + ((gamma-1)/2)*Mach_rampas(1)^2)*Temp_comp(1);

num = 1;                                      % Índice das rampas

for num = 2: num_rampas
    
    sen_beta(num) = Mach_normal_rampas(num-1)/Mach_rampas(num-1);
    beta_rampas(num) = asind(sen_beta(num));
    Mach_normal_rampas(num) = Mach_rampas(num-1)*sind(beta_rampas(num));
    tg_teta(num)= 2*cotd(beta_rampas(num))*((Mach_normal_rampas(num)^2-1)/(Mach_rampas(num-1)^2*(gamma+cosd(2*beta_rampas(num)))+2));
    teta_rampas(num) = atand(tg_teta(num));
    Mach_rampas(num) = (sqrt((Mach_normal_rampas(num)^2+(2/(gamma-1)))/(((2*gamma)/(gamma-1))*Mach_normal_rampas(num)^2-1)))/sind(beta_rampas(num)-teta_rampas(num));
    
    razao_pressao(num) = 1+((2*gamma)/(gamma+1))*(Mach_normal_rampas(num)^2-1);
    Pres_comp(num) = razao_pressao(num)*Pres_comp(num-1);
    
    razao_rho(num) = ((gamma+1)*Mach_normal_rampas(num)^2)/((gamma-1)*Mach_normal_rampas(num)^2+2);
    rho_comp(num) = razao_rho(num)*rho_comp(num-1);
    
    razao_temperatura(num) = razao_pressao(num)/razao_rho(num);
    Temp_comp(num) = razao_temperatura(num)*Temp_comp(num-1);
    
    velocidade_som(num) = sqrt(gamma*const_gas_ar*Temp_comp(num));
    velocidade_escoamento(num) = velocidade_som(num)*Mach_rampas(num);
    
    Temp_total_comp(num) = (1 + ((gamma-1)/2)*Mach_rampas(num)^2)*Temp_comp(num);
    
    num = num + 1;
    
end

%% Choque refletido

teta_refletido = sum(teta_rampas);
lambda = sqrt((Mach_rampas(num_rampas)^2-1)^2-3*(1+((gamma-1)/2)*Mach_rampas(num_rampas)^2)*(1+((gamma+1)/2)*Mach_rampas(num_rampas)^2)*(tand(teta_refletido))^2);
X = ((Mach_rampas(num_rampas)^2-1)^3-9*(1+((gamma-1)/2)*Mach_rampas(num_rampas)^2)*(1+((gamma-1)/2)*Mach_rampas(num_rampas)^2+((gamma+1)/4)*Mach_rampas(num_rampas)^4)*(tand(teta_refletido))^2)/lambda^3;
tg_beta_refletido = (Mach_rampas(num_rampas)^2-1+2*lambda*cos((4*pi+(acos(X)))/3))/(3*(1+((gamma-1)/2)*Mach_rampas(num_rampas)^2)*tand(teta_refletido));
beta_refletido = atand(tg_beta_refletido);
Mach_normal_station_3 = Mach_rampas(num_rampas)*sind(beta_refletido);
Mach_station_3 = (sqrt((Mach_normal_station_3^2+(2/(gamma-1)))/(((2*gamma)/(gamma-1))*Mach_normal_station_3^2-1)))/sind(beta_refletido-teta_refletido);

razao_pressao_station_3 = 1+((2*gamma)/(gamma+1))*(Mach_normal_station_3^2-1);
Pres_station_3 = razao_pressao_station_3*Pres_comp(num_rampas);
    
razao_rho_station_3 = ((gamma+1)*(Mach_normal_station_3)^2)/((gamma-1)*Mach_normal_station_3^2+2);
rho_station_3 = razao_rho_station_3*rho_comp(num_rampas);
    
razao_temperatura_station_3 = razao_pressao_station_3/razao_rho_station_3;
Temp_station_3 = razao_temperatura_station_3*Temp_comp(num_rampas);

velocidade_som_station_3 = sqrt(gamma*const_gas_ar*Temp_station_3);
velocidade_escoamento_station_3 = velocidade_som_station_3*Mach_station_3;

Temp_total_station_3 = (1 + ((gamma-1)/2)*Mach_station_3^2)*Temp_station_3;
