%% 1.1
V = 60; % mph
W = 3000; % lb
d = 50; % in
h1 = 25; % in
h2 = 20; % in
Fd = 60; % lb
Ft = Fd;
syms Wf Wr
[Wf Wr] = vpasolve(Wr+Wf==W,Fd*(h1-h2)+Wf*d+Ft*h2==Wr*d)
P = Ft*V*5280/60/33000
r = 12.5; %in
T = Ft/2*r/12
w = V*63360/60/r
%% 1.4
P = 1;
n1 = 1800;
T12 = P*5252/n1;

n3 = 6000;
T23 = P*5252/n3;

syms F1 F2;
G = 2;
[F1 F2] = vpasolve(F1+F2==G, F1*5/12+T23==F2*5/12+T12)
%% 1.5
W = 800; % N
r1 = 330;
r2 = 40;
r3 = 100;
d1 = 400;
d2 = 600;
d3 = 160;

T = W*d3/r3;
Ft = W*(d3)/r3*r2/r1;
syms Fr Ff;
[Fr Ff] = solve(Fr+Ff==200,Ff*330+Fr*1330-W*770-Ft*330==0)
%% 2.1
w = 10/100; %m
h = 30/100; %m
P = 750; %N
Mz = 60; %Nm

yA = h/2;
yB = 0;
yC = h*0.25;

I = w*h^3/12;

Normal_A = (P/(w*h) + Mz*yA/I)/1000
Normal_B = P/(w*h)/1000
Normal_C = (P/(w*h) - Mz*yC/I)/1000
%% 2.2
L = 400/1000; %m
d = 40/1000; %m
P = 1600; %N
L1 = 15/1000; %m

Mx = P*L;
Tz = P*L1;

I = pi*d^4/64;
J = pi*d^4/32;
A = pi*d^2/4;

Normal_z = Mx*d/2/I*10^-6
Shear_xz = Tz*d/2/J*10^-6
Shear_max = sqrt(Shear_xz^2 + Normal_z^2/4)
%% 3.1
A = 20*35/100/100; % m^2
L = 2.9; % m

Le = 2*L;
I = 0.35*0.2^3/12;
rho = sqrt(I/A);
Le/rho

A = 35*35/10000;
L = 2;

Le = 2*L;
I = 0.35*0.35^3/12;
rho = sqrt(I/A);
Le/rho
%% 3.2
L = 4; %m
E = 200*10^9; %Pa
A = 20*60/1000/1000; %m^2
Sfc = 2;

I = 1/12*0.06*0.02^3;
Le = 2*L;
Pcr = pi^2*E*I/Le^2;
P = Pcr/Sfc;
W = P/9.81
%% 3.3
A = 1.64; %in^2
I11 = 2.5; %in^4
I22 = 0.46; %in^4
E = 30*10^6; %psi
Sy = 42; %ksi
Le = 10*12; %in
Sfc = 2;

P = pi^2*E*I22/Le^2/Sfc
%% 4.1
d = 30;    % mm
b = 200;    % mm
h = 25;    % mm
P = 65000;    % N

Kt = 2.5;

stress_nom = P/((b-d-d)*h)

stress_max = stress_nom*Kt

Kt = 2.45;

stress_nom = (P/2)/(h*(b/2-d))

stress_max = stress_nom*Kt
%% 4.2
b = 60;    % mm
d = 10;    % mm
h = 20;    % mm
P = 400*1000;    % Pa

stress_nom = P/((b-d)*h);

Kt = 2.5;

stress_max = stress_nom*Kt
%% 5.1
P = 50000; %lb
K_IC = 55; %ksi*in^0.5
d = 3; %inch
w = 5; %inch
c = 1; %inch
t = 1; %inch
Y = 1.2;

stress_normal = P/(w*t);
K0 = stress_normal*sqrt(pi*c)/1000;
KI_axial = K0*Y

M = P*0.5;
Y = 1.05;
K0 = 6*M*sqrt(pi*c)/w^2/t/1000;
KI_bend = K0*Y

KI_total = KI_axial + KI_bend

K_IC/KI_total
%% 5.2
shear_xy = 100; %MPa
normal_x = 50; %MPa
Sy = 500; %MPa

normal_VM = sqrt(normal_x^2+3*shear_xy^2)
Sy/normal_VM


S_sy = 100;
shear_max = sqrt(50^2+200^2)/2;

S_sy/shear_max
%% 6.1
W = 10000; % lb
A = 2.5; % in^2
E = 12*10^6; % psi
v = 400*12/60; % in/s
L = 70*12; % in
g = 32*12; % in/s^2

k = E*A/L;
d_st = W/k;
IF = 1 + sqrt(1+v^2/g/d_st);
d = IF*d_st % maximum elongation, inch
F_max = d*k;
Ts_max = F_max/A/1000
%% 6.2
k = 5000*1000; % N/m
L = 5; % m
v = 4000/60/60; % m/s
m = 1400; % kg

U = m*v^2/2
d = sqrt(2*U/k); % m
F = sqrt(2*U*k); % N
d*1000
F/1000
%% 6.3
W = 1400*9.81;
v = 4/3.6;
L = 5;
k = 5*10^6;
g = 9.81;

d_st = W/k
d = d_st*(1+sqrt(1+v^2/g/d_st))
Fe = W*(1+sqrt(1+v^2/g/d_st))
%% 7.2
clc;clear;
Su = 708*10^6; %Pa
Gc = 0.8;
Gs = 0.75;
GT = 1;
GR = 1;
%% 7.3 (Unknown)
clc;clear;
Su = 800; % MPa
Sy = 500; %MPa
Cg = 0.7;
Kt = 1.9;
q = 0.86;
CL = 0.58;
CG = 1;
CS = 0.9;
CT = 1;
CR = 0.868;
S= 0.5*Su*CL*CG*CS*CT*CR
SF = Sy/S
%% 7.4
n1 = 1;
N1 = 1.6*10^4;
n2 = 2;
N2 = 3.8*10^4;
n3 = 5;
N3 = 10^5;
syms SF ncycles;
[SF, ncycles] = vpasolve(ncycles==1/SF/(n1/N1+n2/N2+n3/N3),20*ncycles == 24*60*60)
%% 8.1
clc;clear;
dm = 22*10^-3; %m
pitch = 0.5*10^-3; %m
F = 100; %N
r = 0.01; %m
f = 0.125;
T = F*r;
an = 14.5/360*2*pi;
L = 1*10^-3;
dc = 40*10^-3;
syms W;
[W] = vpasolve(T==W*dm/2*(f*pi*dm+L*cos(an))/(pi*dm*cos(an)-f*L)+W*f*dc/2)
syms f;
[f] = vpasolve(1.25*f==L*cos(an)/(pi*dm))
%% 8.2
clc;clear;
Fbolt = 6000; % lb
Sf = 3;
Sp = 85*10^3;
syms At;
[At] = vpasolve(Fbolt*Sf==0.9*At*Sp)
syms t;
d = 5/8;
At = 0.256;
Sy_bolt = 36;
Sy_nut = 92;
[t] = vpasolve(pi*d*(0.75*t)*Sy_bolt*0.58==At*Sy_nut)
N = 18*t
%% 8.3
clc;clear;
SF = 2;
F = SF*9*1000/2;
syms At;
Sp = 380*10^6;
[At] = vpasolve(F==At*0.9*Sp);
At*10^6
At = 28.9*10^-6; %m^2
F = 0.9*At*380*10^6; %N
d =  7*10^-3; %m
T = 0.2*F*d
S = 380/2*0.9/3.8/2*10^6;
F = S*At
%% 8.4
clc;clear;
F = 15;
D = 10*10^-3; %m
d = 1*10^-3; %m
C = D/d;
Ks = 1+0.5/C;
tmax = 8*F/(pi*d^2)*C*Ks;
Su = 1300*10^6; %Pa
clash = 1.1;
syms SF;
[SF] = vpasolve(tmax==0.45*Su/(clash*SF))
syms N;
delta = 150*10^-3;
G = 79*10^9; %Pa
[N] = vpasolve(delta/N==8*F*D^3/(d^4*G))
Lf = (clash*SF*N*d+0.15)*1000
%% 8.5
clc;clear;
Su = 1300*10^6; %Pa
d = 5*10^-3; %m
C = 4;
clash = 1.1;
tmax = 0.43*Su/clash^2;
Kw = (4*C-1)/(4*C-4)+0.615/C
F = tmax*pi*d^2/C/Kw
%% 9.1
T = 600; % Nm
r_out = 150*10^-3/2; % m
p_max = 0.8*10^6; %Pa
f = 0.15; 

r_in = sqrt(1/3)*r_out;
N = T/pi/p_max/r_in/f/(r_out^2-r_in^2)
N = 10;
F = 2*T/f/(r_out+r_in)/N*10^-3
%% 9.2
T = 275; % Nm
N = 2;
f = 0.25;
p_max = 350*1000; %Pa

syms ri ro;
[ri ro] = vpasolve(ri == sqrt(1/3)*ro, T==pi*p_max*ri*f*(ro^2-ri^2)*N,ri,ro)
ri = ri(3,1);
ro = ro(3,1);
d = ro*2*1000
F = 2*T/f/(ro+ri)/N*10^-3

%Even number rounded up - 向上取偶数
%% 9.3
clc;clear;
w = 35*10^-3; %m
f = 0.25;
p_max = 700*1000; %Pa
c = (320+250)*10^-3; %m
b = 250*10^-3; %m
a = 40*10^-3; %m
r = 200*10^-3; %m

A = w*0.15;
N = A*p_max

F = N*(b-f*a)/c
T = N*f*r
%% 9.4
clc;clear;
w = 30*10^-3; %m
f = 0.35;
T = 800; %Nm
F = 300; %N
r = 0.5; %m

syms P1 phi;
[P1, phi] = vpasolve(T==(P1-F)*r,P1/exp(f*phi)==F,P1,phi)
phi = vpa(phi/2/pi*360)
