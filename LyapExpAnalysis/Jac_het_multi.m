function [J1,J2] = Jac_het_multi(param,param_h1,het_choice1,param_h2,het_choice2,r,Omega,delta,N,tau)
% This code calculates the Jacobians J1 and J2 for the heterogeneous sync
% state
%  
% param        - System parameters
% param_h1     - First heterogeneous parameter
% het_choice1  - Choice of first heterogeneous parameter
% param_h2     - Second heterogeneous parameter
% het_choice2  - Choice of second heterogeneous parameter
% r            - vector of amplitudes
% Omega        - Frequency shift
% delta        - Vector of phase shifts
% N            - Carrier number
% tau          - Time delay

%%%%--- Parameters ---%%%%%
N_0 = param.N_0;       
g = param.g;       
s = param.s;        
gamma = param.gamma;      
alpha0 = param.alpha0;       
gamma_n = param.gamma_n;    
kappa = param.kappa;      
M = param.M; 
adj = param.adj;
phase_lag = delta;
neig = param.neig;



alpha_v = alpha0*ones(M,1);
kappa_v = kappa*ones(M,1);

if het_choice1 == "alpha"
    alpha_v = param_h1;           
elseif het_choice1 == "kappa"
    kappa_v = param_h1;            
end

if het_choice2 == "alpha"
    alpha_v = param_h2;
elseif het_choice2 == "kappa"
    kappa_v = param_h2;            
end





J1 = zeros(3*M);
J2 = zeros(3*M);

for kk = 1:M
    J1(kk,kk) = (1/2)*(g*((N(kk)-N_0)/(1+s*(r(kk))^2)^2)*(1-s*(r(kk))^2) - gamma);          % (\partial F1)/(\partial r)
    J1(kk,kk+M) = (kappa_v(kk)/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,1);                % (\partial F1)/(\partial phi)
    J1(kk,kk+2*M) = (1/2)*((g*r(kk))/(1+s*(r(kk))^2));                                      % (\partial F1)/(\partial N)
    
    J1(kk+M,kk) = -alpha_v(kk)*g*s*((N(kk)-N_0)/((1+s*(r(kk))^2)^2))*r(kk) - (kappa_v(kk)/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,2);        % (\partial F2)/(\partial r)     
    J1(kk+M,kk+M) = -(kappa_v(kk)/neig)*sum_cos(M,adj,r,phase_lag,Omega,tau,kk);                                                                  % (\partial F2)/(\partial phi)
    J1(kk+M,kk+2*M) = (alpha_v(kk)/2)*(g/(1+s*(r(kk))^2));                                                                                     % (\partial F2)/(\partial N)
    
    J1(kk+2*M,kk) = -2*g*r(kk)*((N(kk)-N_0)/((1+s*r(kk)^2)^2))*(1e-4);                 % (\partial F3)/(\partial r)
    J1(kk+2*M,kk+2*M) = -gamma_n - g*1e-4*((r(kk))^2/(1+s*r(kk)^2));                   % (\partial F3)/(\partial N)

    for kkk=1:M
        J2(kk,kkk) = (kappa_v(kk)/neig)*adj(kk,kkk)*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                       % (\partial F1)/(\partial r)
        J2(kk,kkk+M) = -(kappa_v(kk)/neig)*adj(kk,kkk)*r(kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);             % (\partial F1)/(\partial phi)

        J2(kk+M,kkk) = (kappa_v(kk)/(neig*r(kk)))*adj(kk,kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                       % (\partial F2)/(\partial r)
        J2(kk+M,kkk+M) = (kappa_v(kk)/neig)*adj(kk,kkk)*(r(kkk)/r(kk))*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);              % (\partial F2)/(\partial phi)
    end            
end

end
        