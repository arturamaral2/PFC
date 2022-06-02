clear
clc
syms  k1 k2 k0 k3 s;


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
Kg = 14 ;
Rm = 2.6 ;
N_motor = 0.69 ;
N_gearbox = 0.85;
N_total = N_motor + N_gearbox;
Jv =  (1/2) * mv * L^2;

%%  representação espaço estados 

a22 = -(((Kg^2)*Ki*Km*N_total)/(Rm *(Jv + mb * delta^2 ) ) *(L^2/d^2));
a23 = - (( mb * Jv *g + (mb^2 )*g * delta^2 )  / ((Jv + mb*delta^2)^2));
a41 = -5*g/7;
b21 = ((Kg*Ki*N_total)/(Rm* (Jv + (mb*(delta^2) )))) * (L/d);

A = [ 0 1 0  0 ; 
      0  a22  a23 0; 
      0 0 0 1; 
      a41 0 0 0 ];

B = [0 ;
     b21;
     0 ;
     0];

C = [0 0 1 0];

D = [0];

E = eig(A)
%% sistema em função de transferencia

sys = ss(A,B,C,D)
[a,b] = ss2tf(A,B,C,D)
F  = tf(a,b)
F0 = a(5)/b(5) % funcao de transferencia com S avaliado em zero
ganho_entrada = (Kg * Ki *N_total * L)/(Rm * d) 

%% Projeto de controlador LQR

Q  = eye(4);
R = 0.7;
% equação de riccati
P = icare(A,B,Q,R)
K = (R^-1)*transpose(B)*P
K2 = lqr(A,B,Q,R)

M = (F0^-1)*transpose(K)*(K*transpose(K))^-1


%% Projeto controlador alocação de polos

p= [-0.5 -1.5 -3 -5]
K_alocacao = place(A,B,p)
l = place(A',C',p)
M_alocacao = (F0^-1)*transpose(K_alocacao)*(K_alocacao*transpose(K_alocacao))^-1

%K_alocacao = [k0 k1 k2 k3]
%I = eye(4)
%matriz_alocacao = (s*I - (A - B*K_alocacao))
%polinomio_alocacao = det(matriz_alocacao)
%eqn = polinomio_alocacao == 0
%polos_alocao = solve(eqn,s)
