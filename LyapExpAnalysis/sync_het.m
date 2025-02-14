function y = sync_het(x,param,param_h,het_choice)
% This code evaluates how close a heterogeneous state 'x' is from being a
% solution for the transcendental equations. The more negative 'y' the
% closest the state 'x' is to solving the transcendental equations
%  
% x           - Homogeneous state [r;Omega;deltas]
% param       - System parameters
% param_h     - Heterogeneous parameter
% het_choice  - Choice of heterogeneous parameter

% Parameters 
N_0 = param.N_0;       
g = param.g;       
s = param.s;        
gamma = param.gamma;      
alpha0 = param.alpha0;       
gamma_n = param.gamma_n;    
kappa = param.kappa;      
sigma0 = param.sigma0;       
omega0 = param.omega0;      
a = param.a;            
tau = param.tau; 
M = param.M; 
adj = param.adj;
neig = param.neig;

alpha_v = alpha0*ones(M,1);
omega_v = omega0*ones(M,1);
kappa_v = kappa*ones(M,1);

if het_choice == "alpha"
    alpha_v = param_h;
elseif het_choice == "omega"
    omega_v = param_h;
elseif het_choice == "kappa"
    kappa_v = param_h;    
end

% Finding the synchronous solution 
y = zeros(2*M,1);
for j = 1:M
	if j == 1
		y(j) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma)*x(j) + (kappa_v(j)/neig)*adj(j,:)*(x(1:M).*cos([0;x(M+2:2*M)] - x(M+1)*tau)); %
		y(M+j) = - x(M+1) + (alpha_v(j)/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega_v(j) + (kappa_v(j)/neig)*adj(j,:)*((x(1:M)/x(j)).*sin([0;x(M+2:2*M)] - x(M+1)*tau));
	else
		y(j) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma)*x(j) + (kappa_v(j)/neig)*adj(j,:)*(x(1:M).*cos([0;x(M+2:2*M)] - x(M+j) - x(M+1)*tau));
		y(M+j) = - x(M+1) + (alpha_v(j)/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g(j)*N_0+ a*gamma)-gamma) + sigma0*omega_v(j) + (kappa_v(j)/neig)*adj(j,:)*((x(1:M)/x(j)).*sin([0;x(M+2:2*M)] - x(M+j) - x(M+1)*tau));
	end
end
end
