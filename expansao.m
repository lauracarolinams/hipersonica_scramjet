function[Mach_prandtl_meyer, Temp_station_10, Pres_station_10, Mach_station_10, altura_station_10, altura_station_10_ideal, rho_station_10, velocidade_som_station_10, velocidade_station_10, altura_prandtl_meyer] = expansao(Mach_station_4, altura_station_3, teta_expansao, gamma, Pres_station_0, Pres_station_4, Temp_station_4, const_gas_ar, altura_station_0, rho_station_0, velocidade_station_0, rho_station_4, velocidade_station_4)

mu_frente_onda = asind(1/Mach_station_4);
nu_Mach_station_4 = rad2deg(sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(Mach_station_4^2-1)))-atan(sqrt(Mach_station_4^2-1)));
nu_Mach_prandtl_meyer_ideal = teta_expansao + nu_Mach_station_4;

erro_max_Mach_prandtl_meyer = 10^-5;

% Intervalo inferior
a1 = 0;
nu_Mach_prandtl_meyer = rad2deg(sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(a1^2-1)))-atan(sqrt(a1^2-1)));
erro_1_Mach_prandtl_meyer = nu_Mach_prandtl_meyer - nu_Mach_prandtl_meyer_ideal;

% Intervalo superior
a2 = 10;
nu_Mach_prandtl_meyer = rad2deg(sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(a2^2-1)))-atan(sqrt(a2^2-1)));
erro_2_Mach_prandtl_meyer = nu_Mach_prandtl_meyer - nu_Mach_prandtl_meyer_ideal;

erro_Mach_prandtl_meyer = max(erro_1_Mach_prandtl_meyer, erro_2_Mach_prandtl_meyer);


if erro_1_Mach_prandtl_meyer*erro_2_Mach_prandtl_meyer > 0
    disp('Necessário alterar valores máximo ou mínimo para Mach de Prandtl-Meyer')
    
else
    while abs(erro_Mach_prandtl_meyer) > erro_max_Mach_prandtl_meyer
        z =(a1 + a2)/2;
        nu_Mach_prandtl_meyer = rad2deg(sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(z^2-1)))-atan(sqrt(z^2-1)));
        erro_Mach_prandtl_meyer = nu_Mach_prandtl_meyer - nu_Mach_prandtl_meyer_ideal;

        if erro_Mach_prandtl_meyer*erro_1_Mach_prandtl_meyer > 0
            a1 = z;
        else
            a2 = z;
        end
    end
 Mach_prandtl_meyer = z;
end

mu_cauda_onda = asind(1/Mach_prandtl_meyer);

AD = (altura_station_3/2)/(2*sind(mu_frente_onda));
BD = (sind(teta_expansao + mu_frente_onda)*AD);
AB = cosd(teta_expansao + mu_frente_onda)*AD; 
BC = BD/tand(mu_frente_onda - teta_expansao);
h_estrela = (AB + BC)*sind(teta_expansao);
altura_prandtl_meyer = 2*h_estrela + altura_station_3/2;

%% Propriedades Prandtl-Meyer

Temp_prandtl_meyer = ((1+((gamma-1)/2)*Mach_station_4^2)/(1+((gamma-1)/2)*Mach_prandtl_meyer^2))*Temp_station_4;
Pres_prandtl_meyer = ((Temp_prandtl_meyer/Temp_station_4)^(gamma/(gamma-1)))*Pres_station_4;
rho_prandtl_meyer = Pres_prandtl_meyer/(const_gas_ar*Temp_prandtl_meyer);
velocidade_som_prandtl_meyer = sqrt(gamma*const_gas_ar*Temp_prandtl_meyer);
velocidade_escoamento_prandtl_meyer = velocidade_som_prandtl_meyer*Mach_prandtl_meyer;


%% Propriedades da seção 10

Pres_station_10 = Pres_station_0;
Temp_station_10 = ((Pres_station_10/Pres_prandtl_meyer)^((gamma-1)/gamma))*Temp_prandtl_meyer;
Mach_station_10 = sqrt((((1+((gamma-1)/2)*Mach_prandtl_meyer^2)*(Temp_prandtl_meyer/Temp_station_10))-1)/((gamma-1)/2));
rho_station_10 = Pres_station_10/(const_gas_ar*Temp_station_10);
velocidade_som_station_10 = sqrt(gamma*const_gas_ar*Temp_station_10);
velocidade_station_10 = velocidade_som_station_10*Mach_station_10;
altura_station_10_ideal = 2*(altura_station_0/2 - altura_station_3/2) + altura_station_3/2;
% altura_station_10_ideal = (rho_station_4*velocidade_station_4*(altura_station_3/2))/(rho_station_10*velocidade_escoamento_station_10*(sind(teta_expansao)/teta_expansao));


altura_station_10 = (rho_prandtl_meyer*velocidade_escoamento_prandtl_meyer*altura_prandtl_meyer)/(rho_station_10*velocidade_station_10);

