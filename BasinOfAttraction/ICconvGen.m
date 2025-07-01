function conv = ICconvGen(adj,r_pert,r_sol,Omega_sol,Omega_pert,delta_sol,alpha_h,omega_h,kappa_h,param)

% Parameters 
N_0 = param.N_0;       
g = param.g;       
s = param.s;        
gamma_n = param.gamma_n;         
J_0 = param.J_0;
M = param.M; 
tau = param.tau;

N_pert = zeros(1,M);
for ind = 1:M
    N_pert(ind) =((((g*(r_pert(ind))^2)/(1+s*(r_pert(ind))^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*(r_pert(ind))^2)/(1+s*(r_pert(ind))^2))*N_0*1e-4);
end


% Simulating the system for T time units
T = 500;
time = [0 T];
opts = ddeset('RelTol',1e-3,'AbsTol',1e-6, 'MaxStep', 1e-2);

rhs = @(t,x,Z)rhs_aux(t,x,Z,param,adj,alpha_h,omega_h,kappa_h);
rhshist = @(t)rhshist_aux(t,r_pert,Omega_pert,N_pert,delta_sol);

% Solving the DDE
sol = dde23(rhs,tau,rhshist,time,opts);

% Find the first index where the entry is larger than the threshold
iniTime = find(sol.x > T -100, 1);

freq = sol.yp(M+1:2*M,iniTime:end);
rs = sol.y(1:M,iniTime:end);
phi = sol.y(M+1:2*M,iniTime:end);
phases = wrapToPi(phi);
szd = size(phi,2);
deltas = zeros(M-1,szd);


conv = 1;
for i = 1:M
    max_val = max(rs(i,:));
    min_val = min(rs(i,:));
    if abs(max_val - min_val) < 1e-4 && abs(max_val - r_sol(i)) < 1e-4 && abs(min_val -r_sol(i)) < 1e-4
        continue
    else
        conv = 0;
    end
end


% Calculate the max and min of freq
for i = 1:M
    max_val = max(freq(i,:));
    min_val = min(freq(i,:));
    
    % Check if the absolute difference is less than 10^-4
    if abs(max_val - min_val) < 1e-4 && abs(max_val - Omega_sol) < 1e-4 && abs(min_val - Omega_sol) < 1e-4
    else
        conv = 0;
    end
end

end

%% Defining the DDEs 
function dxdt = rhs_aux(t,x,Z,param,adj,alpha_h,omega_h,kappa_h)  
    %%%%--- Parameters ---%%%%%
    N_0 = param.N_0;       
    g = param.g;       
    s = param.s;        
    gamma = param.gamma;        
    gamma_n = param.gamma_n;      
    sigma0 = param.sigma0;       
    J_0 = param.J_0;
    M = param.M; 
    neig = param.neig;

	dxdt = zeros(3*(M),1);
	xlag = Z(:,1);

    % Noise term (Gaussian noise)
    noise_std = 0; % Standard deviation of the noise
    noise = noise_std * randn(3*M, 1);

	for iii = 1:M
        sum1=0;
        sum2=0;
        for jjj = 1:M
            sum1 = sum1 + adj(iii,jjj)*(xlag(jjj)*cos(xlag(M+jjj)-x(M+iii)));
            sum2 = sum2 +adj(iii,jjj)*((xlag(jjj)/x(iii))*sin(xlag(M+jjj)-x(M+iii)));
        end
            
		dxdt(iii) = (1/2)*(g*((x(2*M+iii)-N_0)/(1+s*x(iii)^2))-gamma)*x(iii) + (kappa_h(iii)/neig)*sum1; 
		dxdt(M+iii) = (alpha_h(iii)/2)*(g*((x(2*M+iii)-N_0)/(1+s*x(iii)^2))-gamma) + sigma0*omega_h(iii) + (kappa_h(iii)/neig)*sum2; 
        dxdt(2*M+iii) = J_0 - gamma_n*x(2*M+iii) - 1e-4*g*((x(2*M+iii)-N_0)/(1+s*x(iii)^2))*x(iii)^2; 
    end
        

    dxdt = dxdt + noise;

end

%% Defining the initial conditions 
function S = rhshist_aux(t,r,Omega,N,delta)	
	histr = r;
	histphi = Omega*t + delta;
    histN = N';
	S = [histr;histphi;histN];
end