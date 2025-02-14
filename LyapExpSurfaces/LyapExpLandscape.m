%% Lyapunov exponent landscape for a three nonidentical laser network
% This code analyses the stability of a three nonidentical lasers by
% calculating the Lyapunov exponent for each heterogeneous combinations of
% a given parameter (one of the 3):
%  - alpha (linewidth enhancement factor)
%  - kappa (coupling strength)
%  - omega (frequency detuning)

clear all; close all; clc;

%%  Defining the parameters  
N_0 = 1.5;     param.N_0 = N_0;         % Carrier number at transparency   (1.5e08, rescaled N_0 -> N_0*1e-8)
g = 1.5e03;    param.g = g;             % Differential gain coefficient    (1.5e-05, rescaled g -> g*1e8)
s = 1e-03;     param.s = s;             % Gain saturation coefficient      (1e-07, rescaled s -> s*1e4)
gamma = 500;   param.gamma = gamma;     % Cavity loss
alpha0 = 5;    param.alpha0 = alpha0;   % Linewidth enhancement factor
gamma_n = 0.5; param.gamma_n = gamma_n; % Carrier loss rate
sigma0 = 3;    param.sigma0 = sigma0;   % Variance of detuning
omega0 = 0;    param.omega0 = omega0;   % Frequency detuning   
a = 2.5455;    param.a = a;             % Pump factor   
tau = 0.15;    param.tau = tau;         % Time delay 
J_0 = a*gamma_n*(N_0+(gamma/g)); param.J_0 = J_0; % Pump current
M = 10;        param.M = M;             % Number of lasers
dx = 0.95;     param.dx = dx;           % Decaying factor
neig = M;      param.neig = neig;       % Number of neighbors 

% Defining topology of coupling (all-2-all)
adj = ones(M,M);   % Adjacency matrix
coup_str = 4.75;   % Coupling strength
deg = 12.5;        % Degree 
neig = 10;         % Number of neighbors

% Coupling parameters     
param.adj = adj;
param.deg = deg; 
param.neig = neig;

%% Finding the sync solution (transcendental equations)
x00 = [5.3413, 0];
x= fminsearch(@(x)sync(x,param),x00);

r = x(1);
Omega = x(2);
N =((((g*r^2)/(1+s*r^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r^2)/(1+s*r^2))*N_0*1e-4);
r0 = repmat(r,M,1); N00 = repmat(N,M,1); Omega0 = Omega; delta = zeros(M,1);


%%%%%--- Calculating Matrices ---%%%%%
M1 = zeros(3,3);
M2 = zeros(3,3);
M3 = zeros(3,3);

% Defining Df matrix
M1(1,1) = (1/2)*((g*(N-N_0)/((1+s*r^2)^2))*(1-s*r^2) - gamma);
M1(1,3) = (g/2)*(r/(1+s*r^2));
M1(2,1) = -(alpha0*g*s*r)*((N-N_0)/((1+s*r^2)^2));
M1(2,3) = (alpha0*g/2)*(1/(1+s*r^2));
M1(3,1) = -2*g*r*((N-N_0)/((1+s*r^2)^2))*1e-4;
M1(3,3) = -(gamma_n + 1e-4*g*((r^2)/(1+s*r^2)));

% Defining D(t)h
M2(1,2) = -r*sin(Omega*tau);
M2(2,1) = (1/r)*sin(Omega*tau);
M2(2,2) = -cos(Omega*tau);

% Defining D(t-tau)h
M3(1,1) = cos(Omega*tau);
M3(1,2) = r*sin(Omega*tau);
M3(2,1) = -(1/r)*sin(Omega*tau);
M3(2,2) = cos(Omega*tau);

J1 = M1 + (coup_str/neig)*deg*M2;

m = 200;
re_fin = 25;
im_fin = 10;
Re_alpha = linspace(-re_fin,re_fin,m);
Im_beta = linspace(-im_fin,im_fin,m);
mtle = zeros(m,m);

for ii = 1:m
    ii
    for jj = 1:m
        J2 = (Re_alpha(ii) + 1i*Im_beta(jj))*(coup_str/neig)*M3;
        mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau);
    end
end



