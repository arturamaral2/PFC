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
%%  representação espaço estados 
a22 = -(  ((Kg^2)*Ki*Km*N_total)/(Rm *(Jv + mb * delta^2 ) ) *(L^2/d^2))
a23 = - ( ( mb * Jv *g + (mb^2 )*g * delta^2 )  / ((Jv + mb*delta^2)^2))
a41 = -5*g/7;
b21 = ((Kg*Ki*N_total)/(Rm* (Jv + (mb*(delta^2) )))) * (L/d)

A = [ 0 1 0  0 ; 
      0  a22  a23 0; 
      0 0 0 1; 
      a41 0 0 0 ]

B = [0 ;
     b21;
     0 ;
     0]

C = [0 0 1 0 ]

D = [0]

%% sistema 
sys = ss(A,B,C,D)
[a,b] = ss2tf(A,B,C,D)
G  = tf(a,b)