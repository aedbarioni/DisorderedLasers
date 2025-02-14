function y = sync_het_multi(x,param,param_h1,het_choice1,param_h2,het_choice2)
% This code evaluates how close a heterogeneous state 'x' is from being a
% solution for the transcendental equations. The more negative 'y' the
% closest the state 'x' is to solving the transcendental equations
%  
% x           - Homogeneous state [r;Omega;deltas]
% param       - System parameters
% param_h1     - First heterogeneous parameter
% het_choice1  - Choice of first heterogeneous parameter
% param_h2     - Second heterogeneous parameter
% het_choice2  - Choice of second heterogeneous parameter

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


if het_choice1 == "alpha"
    alpha_v = param_h1;
elseif het_choice1 == "omega"
    omega_v = param_h1;         
elseif het_choice1 == "kappa"
    kappa_v = param_h1;            
end

if het_choice2 == "alpha"
    alpha_v = param_h2;
elseif het_choice2 == "omega"
    omega_v = param_h2;
elseif het_choice2 == "kappa"
    kappa_v = param_h2;            
end

% Finding the synchronous solution 
y = zeros(2*M,1);
for j = 1:M
	if j == 1
		y(j) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma)*x(j) + (kappa_v(j)/neig)*adj(j,:)*(x(1:M).*cos([0;x(M+2:2*M)] - x(M+1)*tau)); %
		y(M+j) = - x(M+1) + (alpha_v(j)/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega_v(j) + (kappa_v(j)/neig)*adj(j,:)*((x(1:M)/x(j)).*sin([0;x(M+2:2*M)] - x(M+1)*tau));
	else
		y(j) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma)*x(j) + (kappa_v(j)/neig)*adj(j,:)*(x(1:M).*cos([0;x(M+2:2*M)] - x(M+j) - x(M+1)*tau));
		y(M+j) = - x(M+1) + (alpha_v(j)/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(j))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega_v(j) + (kappa_v(j)/neig)*adj(j,:)*((x(1:M)/x(j)).*sin([0;x(M+2:2*M)] - x(M+j) - x(M+1)*tau));
	end
end   
end