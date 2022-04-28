%% Parametros iniciais 
g = 9.8;
mb = 0.064;
mv = 0.65;
R = 0.0254;
L = 0.425;
d = 0.12;
delta = 0.2;
Km = 0.00767;
Ki = 0.00767;
Kg = 14 
Rm = 2.6 
N_motor = 0.69 
N_gearbox = 0.85
N_total = N_motor + N_gearbox;
Jv =  (1/2) * mv * L^2;

torque_v = (Kg*Ki*N_total/Rm)/(L/d) 