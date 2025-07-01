function y = sync_hom_var(x,param,param_var,var_choice)
% This code evaluates how close a homogeneous state 'x' is from being a
% solution for the transcendental equations. The more negative 'y' the
% closest the state 'x' is to solving the transcendental equations
%  
% x           - homogeneous state [r;Omega]
% param       - System parameters
% param_var   - Homogeneous parameter that is varied 
% var_choice  - Choice of homogeneous parameter that is varied

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
deg = param.deg;
neig = param.neig;

if var_choice == "alpha"
    alpha0 = param_var;
elseif var_choice == "N_0"
    N_0 = param_var;
elseif var_choice == "g"
    g = param_var;
elseif var_choice == "s"
    s = param_var;
elseif var_choice == "gamma"
    gamma = param_var;
elseif var_choice == "gamma_n"
    gamma_n = param_var;
elseif var_choice == "tau"
    tau = param_var; 
elseif var_choice == "kappa"
    kappa = param_var;     
elseif var_choice == "deg"
    deg = param_var;     
end

F = zeros(2,1);

    F(2) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(1))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + (kappa/neig)*deg*cos(x(2)*tau);
    F(1) = -x(2) + (alpha0/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(1))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega0 - (kappa/neig)*deg*sin(x(2)*tau);   

    y = log10(norm(F));
end