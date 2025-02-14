function [J1,J2] = Jac_het(param,param_h,het_choice,r,Omega,delta,N,tau)
% This code calculates the Jacobians J1 and J2 for the heterogeneous sync
% state
%  
% param       - System parameters
% param_h     - Heterogeneous parameter
% het_choice  - Choice of heterogeneous parameter
% r           - vector of amplitudes
% Omega       - Frequency shift
% delta       - Vector of phase shifts
% N           - Carrier number
% tau         - Time delay

% Parameters
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


J1 = zeros(3*M);
J2 = zeros(3*M);
if het_choice == "omega" 
    for kk = 1:M
        J1(kk,kk) = (1/2)*(g*((N(kk)-N_0)/(1+s*(r(kk))^2)^2)*(1-s*(r(kk))^2) - gamma);          
        J1(kk,kk+M) = (kappa/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,1);              
        J1(kk,kk+2*M) = (1/2)*((g*r(kk))/(1+s*(r(kk))^2));                                     
        
        J1(kk+M,kk) = -alpha0*g*s*((N(kk)-N_0)/((1+s*(r(kk))^2)^2))*r(kk) - (kappa/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,2);       
        J1(kk+M,kk+M) = -(kappa/neig)*sum_cos(M,adj,r,phase_lag,Omega,tau,kk);                                                                  
        J1(kk+M,kk+2*M) = (alpha0/2)*(g/(1+s*(r(kk))^2));                                                                                     
        
        J1(kk+2*M,kk) = -2*g*r(kk)*((N(kk)-N_0)/((1+s*r(kk)^2)^2))*(1e-4);                
        J1(kk+2*M,kk+2*M) = -gamma_n - g*1e-4*((r(kk))^2/(1+s*r(kk)^2));                  
    
        for kkk=1:M
            J2(kk,kkk) = (kappa/neig)*adj(kk,kkk)*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                      
            J2(kk,kkk+M) = -(kappa/neig)*adj(kk,kkk)*r(kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);            
    
            J2(kk+M,kkk) = (kappa/(neig*r(kk)))*adj(kk,kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                     
            J2(kk+M,kkk+M) = (kappa/neig)*adj(kk,kkk)*(r(kkk)/r(kk))*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);            
        end            
    end
elseif het_choice == "alpha"
    for kk = 1:M
        J1(kk,kk) = (1/2)*(g*((N(kk)-N_0)/(1+s*(r(kk))^2)^2)*(1-s*(r(kk))^2) - gamma);          
        J1(kk,kk+M) = (kappa/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,1);             
        J1(kk,kk+2*M) = (1/2)*((g*r(kk))/(1+s*(r(kk))^2));                                      
        
        J1(kk+M,kk) = -param_h(kk)*g*s*((N(kk)-N_0)/((1+s*(r(kk))^2)^2))*r(kk) - (kappa/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,2);       
        J1(kk+M,kk+M) = -(kappa/neig)*sum_cos(M,adj,r,phase_lag,Omega,tau,kk);                                                                 
        J1(kk+M,kk+2*M) = (param_h(kk)/2)*(g/(1+s*(r(kk))^2));                                                                                     
        
        J1(kk+2*M,kk) = -2*g*r(kk)*((N(kk)-N_0)/((1+s*r(kk)^2)^2))*(1e-4);                 
        J1(kk+2*M,kk+2*M) = -gamma_n - g*1e-4*((r(kk))^2/(1+s*r(kk)^2));                   
    
        for kkk=1:M
            J2(kk,kkk) = (kappa/neig)*adj(kk,kkk)*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                       
            J2(kk,kkk+M) = -(kappa/neig)*adj(kk,kkk)*r(kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);             
    
            J2(kk+M,kkk) = (kappa/(neig*r(kk)))*adj(kk,kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                    
            J2(kk+M,kkk+M) = (kappa/neig)*adj(kk,kkk)*(r(kkk)/r(kk))*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);           
        end            
    end
elseif het_choice == "kappa" 
    for kk = 1:M
        J1(kk,kk) = (1/2)*(g*((N(kk)-N_0)/(1+s*(r(kk))^2)^2)*(1-s*(r(kk))^2) - gamma);       
        J1(kk,kk+M) = (param_h(kk)/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,1);    
        J1(kk,kk+2*M) = (1/2)*((g*r(kk))/(1+s*(r(kk))^2));                                   
        
        J1(kk+M,kk) = -alpha0*g*s*((N(kk)-N_0)/((1+s*(r(kk))^2)^2))*r(kk) - (param_h(kk)/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,2);     
        J1(kk+M,kk+M) = -(param_h(kk)/neig)*sum_cos(M,adj,r,phase_lag,Omega,tau,kk);                                                          
        J1(kk+M,kk+2*M) = (alpha0/2)*(g/(1+s*(r(kk))^2));                                                                                
        
        J1(kk+2*M,kk) = -2*g*r(kk)*((N(kk)-N_0)/((1+s*r(kk)^2)^2))*(1e-4);                
        J1(kk+2*M,kk+2*M) = -gamma_n - g*1e-4*((r(kk))^2/(1+s*r(kk)^2));                  
    
        for kkk=1:M
            J2(kk,kkk) = (param_h(kk)/neig)*adj(kk,kkk)*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                       
            J2(kk,kkk+M) = -(param_h(kk)/neig)*adj(kk,kkk)*r(kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);             
    
            J2(kk+M,kkk) = (param_h(kk)/(neig*r(kk)))*adj(kk,kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                      
            J2(kk+M,kkk+M) = (param_h(kk)/neig)*adj(kk,kkk)*(r(kkk)/r(kk))*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);             
        end            
    end
end        