function y = sync_het_net(x,param,adj)
% This code evaluates how close a heterogeneous state 'x' is from being a
% solution for the transcendental equations when the adjacency matrix is 
% perturbed. The more negative 'y' the closest the state 'x' is to solving 
% the transcendental equations
%  
% x           - Homogeneous state [r;Omega;deltas]
% param       - System parameters
% adj         - Adjacency matrix

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