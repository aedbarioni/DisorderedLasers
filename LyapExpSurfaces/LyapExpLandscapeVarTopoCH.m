%% Lyapunov exponent landscape for a three nonidentical laser network for multiple topologies
% This code analyses the stability of a three nonidentical lasers by
% calculating the Lyapunov exponent for the constrained heterogeneity (CH)
% (where the synchronous solution remain invariant)

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

%% Defining topology of coupling 
topo = input('Topology: ','s');         % Choice of topology
if topo == "A"
    adj = [1 1 1;1 1 1;1 1 1];  
    kappa = 4.82;
elseif topo == "B"
    adj = [1 0 1;1 1 0;0 1 1];  
    kappa = 4.82;
elseif topo == "C"
    adj = [1 0 0;1 1 1;0 0 1];
    kappa = 2.89;
elseif topo == "D"
    adj = [1 0 0;1 1 0;0 1 1]; 
    kappa = 4.7;
elseif topo == "E"
    adj = [1 1 0;1 1 1;0 1 1];
    kappa = 4.0;
elseif topo == "F"
    adj = [1 0 1;1 1 0;1 1 1];  
    kappa = 4.45;
elseif topo == "G"
    adj = [1 0 1;1 1 1;1 0 1];
    kappa = 4.82;
elseif topo == "H"
    adj = [1 1 1;1 1 1;1 0 1];
    kappa = 5.82;
elseif topo == "I"
    adj = [1 1 0;0 1 0;0 1 1];
    kappa = 7.33;
elseif topo == "J"
    adj = [1 0 1;1 1 1;0 0 1];
    kappa = 5.15;
elseif topo == "K"
    adj = [1 1 1;0 1 0;1 1 1];
    kappa = 5.86;
elseif topo == "L"
    adj = [1 1 0;1 1 1;0 0 1];
    kappa = 4.98;
elseif topo == "M"
    adj = [1 1 0;1 1 0;0 1 1];
    kappa = 4.82;
elseif topo == "N"
    adj = [1 0 0;1 1 0;1 1 1];
    kappa = 2.89;
end


sum_deg = 0;
for jjj = 1:M
    sum_deg = sum_deg + adj(1,jjj);
    if adj(1,jjj)==0
        neig = neig-1;
    end
end
deg = sum_deg;

% Coupling parameters
param.adj = adj;     % Adjacency matrix
param.deg = deg;     % Degree (row-sum)
param.kappa = kappa; % Coupling strength

%% Finding the sync solution (transcendental equations)
x0 = [5.3401, -5.0386];
x= fminsearch(@(x)sync(x,param),x0);

r0 = x(1);
omega = x(2);
N_h =((((g*r0^2)/(1+s*r0^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r0^2)/(1+s*r0^2))*N_0*1e-4) ;

param_h = alpha0*ones(M,1);
het_choice = "alpha";
r_h = repmat(r0,M,1); delta = zeros(M,1); N = repmat(N_h,M,1);
[J1,J2] = Jac_freq_het_gen(param,param_h,het_choice,r_h,omega,delta,N,tau);
mtle_ini = dde_rightmost_eig(J1,J2,tau,M);
disp(mtle_ini)


het_choice = "alpha";
param_h0 = alpha0;
m = 501;

delta_ini = -param_h0/2;
delta_end = param_h0/2;

if het_choice == "omega"
    delta_ini = -1.5;
    delta_end = 1.5;
end

delta1 = linspace(delta_ini,delta_end,m);
delta2 = linspace(delta_ini,delta_end,m);
mtle = ones(m);
order = ones(m);

for ii = 1:m
	ii
	parfor jj = 1:m
        param_h = param_h0*ones(M,1) + [delta1(ii);delta2(jj);-delta1(ii)-delta2(jj)];
        r = r0*ones(M,1);
        Omega = omega;
        N = N_h*ones(M,1);
        delta0 = zeros(M,1);
        found =1;
		if found ==1 
            [J1,J2] = Jac_het(param,param_h,het_choice,r,Omega,delta0,N,tau);
            mtle(ii,jj) = dde_rightmost_eig(J1,J2,tau);
		else
			mtle(ii,jj) = Inf;
		end
	
	end
end

dlmwrite(strcat('filename.dat'),mtle);

%%
figure();
imagesc(mtle(end:-1:1,1:1:end)) 
xlabel('\delta_1')
ylabel('\delta_2')
title('Lyapunov Exponent')
colormap(bluewhitered),colorbar

