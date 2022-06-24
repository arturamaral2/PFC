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


amplitude_entrada = 1;
bool_degrau = 1;
bool_sine = 0;
bool_linear = 0; 


% tensões no ponto de operacao
V = ((L*mv*g/2) + mb*g*delta ) * ((Rm * d)/ (L * Kg * Ki * N_total) )

V2 = ((L*mv*g/2) + mb*g*delta ) * ((Rm * L)/ (d * Kg * Ki * N_total) )

%%  representação espaço estados 

a22 = -(((Kg^2)*Ki*Km*N_total)/(Rm *(Jv + mb * delta^2 ) ) *(L^2/d^2));
a23 = -(( mb * Jv *g + (mb^2 )*g * delta^2 )/((Jv + mb*delta^2)^2));
a41 = -5*g/7;
b21 = ((Kg*Ki*N_total)/(Rm* (Jv + (mb*(delta^2) )))) * (L/d);

% modelagem pela tensão 
ganho_entrada_antigo = ((Kg*Ki*N_total)/(Rm)) *(L/d)
ganho_angulo_antigo = (((Kg^2)*Km*Ki*N_total)/(Rm)) * ((L^2)/(d^2))

ganho_entrada = Ki/Rm 
ganho_angulo = (Ki*Km/Rm)

a22 = - ganho_angulo_antigo/(Jv + mb * delta^2 )
b21 = ganho_entrada_antigo /(Jv + mb * delta^2 )


a22_novo = - ganho_angulo /(Jv + mb * delta^2 ) 
b21_novo =  ganho_entrada/(Jv + mb * delta^2 )


A = [ 0 1 0  0 ; 
      0  0  a23 0; 
      0 0 0 1; 
      a41 0 0 0 ];

B = [0 ;
     1;
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
 

%% Projeto de controlador LQR

Q  = eye(4);
R = 0.5;
% equação de riccati
P = icare(A,B,Q,R)
K = (R^-1)*transpose(B)*P
K2 = lqr(A,B,Q,R)

M = (F0^-1)*transpose(K)*(K*transpose(K))^-1

%% Projeto controlador alocação de polos

p= [-0.3 -2.5 -5 -10]
K_alocacao = place(A,B,p)
l = place(A',C',p)
M_alocacao = (F0^-1)*transpose(K_alocacao)*(K_alocacao*transpose(K_alocacao))^-1

%K_alocacao = [k0 k1 k2 k3]
%I = eye(4)
%matriz_alocacao = (s*I - (A - B*K_alocacao))
%polinomio_alocacao = det(matriz_alocacao)
%eqn = polinomio_alocacao == 0
%polos_alocao = solve(eqn,s)


[a2,b2] = ss2tf(A-B*K, B*K*M ,C-D*K,D*K*M)
F_lqr_malha_fechada  = tf(a2,b2)


%% seguimento de referencia erro nulo  LQR

A_nulo = [A; 
          C]
coluna_nula = zeros(5,1)
A_nulo = [A_nulo coluna_nula]

B_nulo = [B; 
     0]
C_nulo = [C 1]
r_nulo = [0;0;0;0;-1]
z_nulo = C_nulo - r_nulo

Q_nulo = eye(5);
R = 0.5;
% equação de riccati
P_nulo = icare(A_nulo,B_nulo,Q_nulo,R)
K_nulo = (R^-1)*transpose(B_nulo)*P_nulo
M_nulo = (F0^-1)*transpose(K_nulo)*(K_nulo*transpose(K_nulo))^-1


k_seguimento_lqr = K_nulo(1:4)
ki_seguimento_lqr = K_nulo(5)

%% seguimento nulo alocacao 
nulvec = [0; 0; 0; 0]

pb = [-0.3 -2.5 -5 -10 -11]
kb_seguimento_alocacao = acker(A_nulo, B_nulo, pb)
k_seguimento_alocacao = kb_seguimento_alocacao(1:4)
ki_seguimento_alocacao = kb_seguimento_alocacao(5)



Az = [A-B*k_seguimento_lqr B*ki_seguimento_lqr; -C 0];
bz = [nulvec ;1];
cz = [C 0 ];
dz = 0 ;

%% Observador de estados 
polo_mais_rapido_do_sistema = min(real(eig(A)))

polo_de_f = 1.5 * polo_mais_rapido_do_sistema
sistema_f = (s - polo_de_f)^length(A)
dem = sym2poly(sistema_f)

tf_sistema = tf(1,dem)
[num, dem] = tfdata(tf_sistema, 'v')
[F,lixo1,lixo2,lixo3] = tf2ss(num,dem)


Lchp = [0; 0; 0; 1]

estados_incontrolaveis = length(A) - rank(ctrb(F,Lchp))


T = lyap(F, -A,Lchp*C)
L_observador = inv(T)*Lchp
 

