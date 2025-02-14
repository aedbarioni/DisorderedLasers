%% MSF dependence of coupling strength and degree
% This code calculates the MSF and determine the deepness and size of the
% stable region in the complex plane for varying values of coupling
% strength and degree

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
neig = 10; param.neig = neig;           % Number of neighbors


dif = 0; 
m = 200;
re_fin = 10;
im_fin = 5;
deg_v = linspace(0,15,m);
kappa_v = linspace(0,15,m);
mtle = zeros(m,m);
sizeXi_m = zeros(m,m);
deepXi_m = zeros(m,m);


for ii = 1:m
    ii
    parfor jj = 1:m
        deg = deg_v(ii);
        coup_str = kappa_v(jj);
        x = deg-1;
        adj = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0; x, 1, 0, 0, 0, 0, 0, 0, 0, 0; 0, x, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, x, 1, 0, 0, 0, 0, 0, 0; 0, 0, 0, x, 1, 0, 0, 0, 0, 0; 0, 0, 0, 0, x, 1, 0, 0, 0, 0; x, 0, 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 0, x, 1, 0, 0; x, 0, 0, 0, 0, 0, 0, 0, 1, 0; 0, 0, 0, 0, 0, 0, 0, 0, x, 1];        
        r0_ini1 = [5.3269;5.3269;5.3269;5.3269;5.3269;5.3269;5.3269;5.3269;5.3269;5.3269;-4.5090;0;0;0;0;0;0;0;0;0];
        r0_ini = ObjAdjSol(coup_str,dif,adj,r0_ini1,param);
        r = r0_ini(1:2);
        Omega = r0_ini(11);
        delta = r0_ini(12:13);
        N = zeros(2,1);
        N(1) =((((g*r(1)^2)/(1+s*r(1)^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r(1)^2)/(1+s*r(1)^2))*N_0*1e-4);
        N(2) =((((g*r(2)^2)/(1+s*r(2)^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r(2)^2)/(1+s*r(2)^2))*N_0*1e-4);
        [J1,J2] = Jac_red(deg,coup_str,neig,tau,r,Omega,delta,N,param);
        mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau);
    end
end

filename = "filename.mat";
save(filename)

%% Calculating the Jacobians
function [J1,J2] = Jac_red(deg,coup_str,neig,tau,r,Omega,delta,N,param)
% Parameters 
N_0 = param.N_0;       
g = param.g;       
s = param.s;        
gamma = param.gamma;      
alpha0 = param.alpha0;       
gamma_n = param.gamma_n;    
kappa = coup_str;    
epsilon = param.epsilon;   
tau_tilde = param.tau_tilde;
tau_prime = tau - epsilon*tau_tilde; 
phase_lag = delta;
M=2;

J1 = zeros(3*M);
J2 = zeros(3*M);
adj = [1 0; (deg-1) 1];
for kk = 1:M
    J1(kk,kk) = (1/2)*(g*((N(kk)-N_0)/(1+s*(r(kk))^2)^2)*(1-s*(r(kk))^2) - gamma);          % (\partial F1)/(\partial r)
    J1(kk,kk+M) = (kappa/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau_prime,kk,1);                % (\partial F1)/(\partial phi)
    J1(kk,kk+2*M) = (1/2)*((g*r(kk))/(1+s*(r(kk))^2));                                      % (\partial F1)/(\partial N)
    
    J1(kk+M,kk) = -alpha0*g*s*((N(kk)-N_0)/((1+s*(r(kk))^2)^2))*r(kk) - (kappa/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau_prime,kk,2);        % (\partial F2)/(\partial r)     
    J1(kk+M,kk+M) = -(kappa/neig)*sum_cos(M,adj,r,phase_lag,Omega,tau_prime,kk);                                                                  % (\partial F2)/(\partial phi)
    J1(kk+M,kk+2*M) = (alpha0/2)*(g/(1+s*(r(kk))^2));                                                                                     % (\partial F2)/(\partial N)
    
    J1(kk+2*M,kk) = -2*g*r(kk)*((N(kk)-N_0)/((1+s*r(kk)^2)^2))*(1e-4);                 % (\partial F3)/(\partial r)
    J1(kk+2*M,kk+2*M) = -gamma_n - g*1e-4*((r(kk))^2/(1+s*r(kk)^2));                   % (\partial F3)/(\partial N)

    for kkk=1:M
        J2(kk,kkk) = (kappa/neig)*adj(kk,kkk)*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau_prime);                       % (\partial F1)/(\partial r)
        J2(kk,kkk+M) = -(kappa/neig)*adj(kk,kkk)*r(kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau_prime);             % (\partial F1)/(\partial phi)

        J2(kk+M,kkk) = (kappa/(neig*r(kk)))*adj(kk,kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau_prime);                       % (\partial F2)/(\partial r)
        J2(kk+M,kkk+M) = (kappa/neig)*adj(kk,kkk)*(r(kkk)/r(kk))*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau_prime);              % (\partial F2)/(\partial phi)
    end            
end

end

%% Calculating sync solution
function r0_new = ObjAdjSol(kappa,dif,adj0,r0,param)
% Inputs:
%        dif ---- vector with nondiagonal perturbation
%        adj0   ---- initial (unperturbed) matrix
%        r0     ---- sync sol for the initial matrix
% Outputs:
%        mtle   ---- Lyapunov exponent of solution of perturbed matrix
%        r0     ---- sync sol for the perturbed matrix

M = param.M;
g = param.g;
s = param.s;
gamma_n = param.gamma_n;
J_0 = param.J_0;
N_0 = param.N_0;
tau = param.tau;

adj_h = assemble_dif(dif,adj0,M);

r = r0(1:M); Omega = r0(M+1); Delta = r0(M+2:2*M);

found = 0;
options = optimoptions('fsolve','Display','off','MaxIterations',1000,'MaxFunctionEvaluations',100000);
for count = 1:10
    		[sol,~,~,~] = fsolve(@(x)sync_het_net(x,param,adj_h,kappa),[r;Omega;Delta],options);
            found = 1;
            r = sol(1:M); Omega = sol(M+1); delta = [0;sol(M+2:2*M)]; 
end

N = zeros(1,M);
for ind = 1:M
    N(ind) =((((g*(r(ind))^2)/(1+s*(r(ind))^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*(r(ind))^2)/(1+s*(r(ind))^2))*N_0*1e-4);
end

r0_new = zeros(2*M,1);
r0_new(1:M) = r; r0_new(M+1) = Omega; r0_new(M+2:2*M) = Delta;

end

%% Transcendental equations
function y = sync_het_net(x,param,adj,coup_str)
% Parameters 
N_0 = param.N_0;       
g = param.g;       
s = param.s;        
gamma = param.gamma;      
alpha0 = param.alpha0;       
gamma_n = param.gamma_n;    
kappa = coup_str;      
sigma0 = param.sigma0;       
omega0 = param.omega0;      
a = param.a;            
tau = param.tau;  
M = param.M; 

neig_ = param.neig;
kappa = (M/neig_)*kappa;

neig = M*ones(M,1);

% Finding the synchronous solution
y = zeros(2*M,1);
for j = 1:M
	if j == 1
		y(j) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma)*x(j) + (kappa/neig(j))*adj(j,:)*(x(1:M).*cos([0;x(M+2:2*M)] - x(M+1)*tau)); %
		y(M+j) = - x(M+1) + (alpha0/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega0 + (kappa/neig(j))*adj(j,:)*((x(1:M)/x(j)).*sin([0;x(M+2:2*M)] - x(M+1)*tau));
	else
		y(j) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma)*x(j) + (kappa/neig(j))*adj(j,:)*(x(1:M).*cos([0;x(M+2:2*M)] - x(M+j) - x(M+1)*tau));
		y(M+j) = - x(M+1) + (alpha0/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega0 + (kappa/neig(j))*adj(j,:)*((x(1:M)/x(j)).*sin([0;x(M+2:2*M)] - x(M+j) - x(M+1)*tau));
	end
end
end

%% Calculating adjacency matrix
function [adj_h] = assemble_dif(dif,adj0,M)
adjdif = zeros(M,M);
idx = 1;

pairs = [2 1;3 2; 4 3; 5 4; 6 5; 7 1; 8 7; 9 1; 10 9];

for ii = 1:9
    i = pairs(ii,1);
    j = pairs(ii,2);
    adjdif(i,j) = dif;
    idx = idx+1;

end

adj_h = adj0 + adjdif;
end