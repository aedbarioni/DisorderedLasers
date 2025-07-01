function [J1,J2] = Jac_het_net(param,adj,r,Omega,delta,N,tau)

%%%%--- Parameters ---%%%%%
N_0 = param.N_0;       
g = param.g;       
s = param.s;        
gamma = param.gamma;      
alpha0 = param.alpha0;       
gamma_n = param.gamma_n;    
kappa = param.kappa;      
% sigma0 = param.sigma0;       
% omega0 = param.omega0;      
epsilon = param.epsilon;   
% a = param.a;            
% tau = param.tau;  
tau_tilde = param.tau_tilde;
% J_0 = param.J_0;
tau_prime = tau - epsilon*tau_tilde;
M = param.M; 
% dx = param.dx;
% adj = param.adj;
% deg = param.deg; 
phase_lag = delta;
neig = param.neig;

% neig = M*ones(M,1);
% for iii = 1:M
%     for jjj = 1:M
%         if adj(iii,jjj)==0
%             neig(iii) = neig(iii)-1;
%         end
%     end
% end






J1 = zeros(3*M);
J2 = zeros(3*M);

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

        