function y = sync(x,param)
% This code evaluates how close a homogeneous state 'x' is from being a
% solution for the transcendental equations. The more negative 'y' the
% closest the state 'x' is to solving the transcendental equations
%  
% x           - homogeneous state [r;Omega]
% param       - System parameters

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
deg = param.deg;
neig = param.neig;

F = zeros(2,1);
    
    F(2) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(1))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + (kappa/neig)*deg*cos(x(2)*tau);
    F(1) = -x(2) + (alpha0/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(1))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega0 - (kappa/neig)*deg*sin(x(2)*tau);   

    y = log10(norm(F));
end