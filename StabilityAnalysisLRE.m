%% Stability analysis of the laser rate equations 
% This code analyses the stability of coupled lasers as parameter disorder 
% is continually introduced in the frequency detuning (omega)

%% Defining the parameters
gamma = 1.7;           param.gamma = gamma;  
tau_c = 13.3;          param.tau_c = tau_c;   
omega = 0;             param.omega = omega;  
tau_f = tau_c*10^3;    param.tau_f = tau_f; 
J_0 = 19.8;            param.J_0 = J_0;
s = 0.25;              param.s = s;
kappa = 0.13;          param.kappa = kappa;
M = 10;                param.M = M;
tau = input('tau: ');  param.tau = tau;
dx = 0.95;             param.dx = dx;

%% Defining topology of coupling 
topo = input('Topology: ','s') % Choice of topology
%       'topo = A' -> All-to-all topology
%       'topo = D' -> Decaying topology
%       'topo = R' -> Ring topology

if topo == "D"
    % Defining decaying coupling network
    adj = zeros(1,M);
    if mod(M,2)==0
        for ii = 1:(M/2)+1
            adj(ii) = dx^(ii-1);
        end
        for ii = (M/2)+2:M
            adj(ii) = dx^(M-ii+1);
        end
    end
    
    ring = adj;
    for ii = 1:M-1
      adj = [adj;circshift(ring,[0,ii])];
    end
    if M ==2
        adj = [1 dx; dx 1];
    end
    
    sum_deg = 0;
    for jjj = 1:M
        sum_deg = sum_deg + adj(1,jjj);
        if adj(1,jjj)==0
            neig = neig-1;
        end
    end
    deg = sum_deg;
elseif topo == "R"
    % Defining ring network 
    adj = zeros(1,M);
    adj(1,2) = 1; adj(1,end) = 0;
    ring = adj;
    for ii = 1:M-1
      adj = [adj;circshift(ring,[0,ii])];
    end
    adj = dx*(adj + transpose(adj)) + eye(M);
    if M ==2
        adj = [1 dx; dx 1];
    end
    
    sum_deg = 0;
    for jjj = 1:M
        sum_deg = sum_deg + adj(1,jjj);
        if adj(1,jjj)==0
            neig = neig-1;
        end
    end
    deg = sum_deg;
elseif topo == "A"
    % Defining all-2-all network
    adj = ones(M,M);
    sum_deg = 0;
    for jjj = 1:M
        sum_deg = sum_deg + adj(1,jjj);
    end
    deg = sum_deg;
end

dist = 10000;
for i = 1:100
    r_0 = 5 + + 10*(rand -0.5);
    Omega_0 = 200*(rand -0.5); 
    G_0 = 200*(rand -0.5);
    x0_trial = [r_0,Omega_0,G_0];
    x_cand = fminsearch(@(x)syncLRE(x,param),x0_trial);
    z_cand = syncLRE(x_cand,param);

    if z_cand <-10
            x_now = x_cand;
            if abs(x_now(2))<dist
                x_min = x_now;
            end
    end
end

%% Finding the sync solution (transcendental equations)
x= fminsearch(@(x)syncLRE(x,param),x0);

r0 = x_min(1);
Omega0 = x_min(2);
G0 =x_min(3);

r = repmat(r0,M,1); Omega = Omega0; delta_ = zeros(M,1); G = repmat(G0,M,1);
[J1,J2] = JacLRE_het(param,r,delta_,G,Omega);
mtle_0 = dde_rightmost_eig(J1,J2,tau,M);
disp(mtle_0)

%% Defining the perturbation variables 
trial = 100;                       % Number of realizations of random heterogeneity
sigma_end = input('sigma_end: ');  % Maximum value of sigma_p
m = 1000*sigma_end;                % Step-size for analytical continuation
sigma = linspace(0,sigma_end,m);   % level of heterogeneity measured by standard deviation;

%% Initialize variables 
mtle = Inf(trial,m);              % Initialize maximal transverse Lyapunov exponent
negativity = zeros(trial,m);      % Initialize negativity matrix (to determine fraction of stabilized systems) 
Trial_neg = zeros(trial,1);       % Initialize vector for trials that successfully achieve stability
Sol_m = cell(trial,m);            % Initialize cell that saves solution at each step 
Param_het = cell(trial,m);        % Initialize cell that saves heterogeneous parameter at each step 
seeds = randi([0 100000],1,trial);

%% Main loop
parfor jj = 1:trial
    disp(['jj= ',num2str(jj)])
    % Random seed
    rng(seeds(jj))

    % Random realization of heterogeneity
    sample = normrnd(0,1,M,1);  % unifrnd(0,0,M); %
    sample = (sample - mean(sample))/std(sample);

    % Initialize sync state
	r = repmat(r0,M,1); Omega = Omega0; Delta = zeros(M-1,1); G = repmat(G0,M,1);
    found = 0;
    negPar = 0;
    
    for ii = 1:m
        if found == 2
            disp("Algorithm didn't find limit cycle solution")
            break
        end
        disp(['ii= ',num2str(sigma(ii))])
        hetero = sigma(ii); % Level of heterogeneity
        
        % Adding heterogeneity to the specified parameter(s)
        dif = sample*hetero; % heterogeneity profile scaled by heterogeneity level
        param_h = omega*ones(M,1) + dif;
        
        % Finding the synchronous solution 
        tol = 1e-15;
		option = optimset('Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolFun',tol,'TolX',tol);
		options = optimoptions('fsolve','Display','off','MaxIterations',1000,'MaxFunctionEvaluations',100000);
        for count = 1:10
    		[sol,fval,exitflag,output] = fsolve(@(x)syncLRE_het(x,param,param_h),[r;Omega;Delta],options);

            if abs(sol(M+1)- Omega)<1e-4
                found = 1;
                r = sol(1:M); Omega = sol(M+1); delta = [0;sol(M+2:2*M)]; Delta = sol(M+2:2*M); 
                break
            else
                found = 2;              
            end
        end

        G = zeros(M,1);
        for ind = 1:M
            G(ind) = J_0/(1+(s*r(ind)^2));
        end

        if found == 1 % Sync state identified successfully
            [J1,J2] = JacLRE_het(param,r,delta,G,Omega);
            mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau,M);
            if  mtle(jj,ii)<0
                negativity(jj,ii) = 1; 
                Trial_neg(jj) = 1;
            end

        else % No sync state found
            mtle(jj,ii) = Inf;                       
            disp(param_h)
        end
    
    end
end

disp(mean(Trial_neg));
filename = "filename.mat";
save(filename)

%%
function y = syncLRE(x,param)
% This code evaluates how close a homogeneous state 'x' is from being a
% solution for the transcendental equations. The more negative 'y' the
% closest the state 'x' is to solving the transcendental equations
%  
% x           - homogeneous state [r;Omega]
% param       - System parameters

% Parameters 
gamma = param.gamma;       
tau_c = param.tau_c;       
omega = param.omega;        
tau_f = param.tau_f;      
J_0 = param.J_0;       
s = param.s;
kappa = param.kappa;     
deg = param.deg;

F = zeros(3,1);
    
    F(1) = ((x(3)-gamma)/tau_c) +(1/tau_c)*kappa*deg ;
    F(2) = -x(2)+omega;   
    F(3) = (1/tau_f)*(J_0- x(3)*((s*x(1)^2)+1));

    y = log10(norm(F));
end

%%
function y = syncLRE_het(x,param,param_h)
% This code evaluates how close a heterogeneous state 'x' is from being a
% solution for the transcendental equations. The more negative 'y' the
% closest the state 'x' is to solving the transcendental equations
%  
% x         - Homogeneous state [r;Omega;deltas]
% param     - System parameters
% param_h   - Heterogeneous parameter

% Parameters
gamma = param.gamma;       
tau_c = param.tau_c;        
J_0 = param.J_0;       
s = param.s;
kappa = param.kappa;         
M = param.num_nodes;
A = param.A;
tau = param.tau;
omega_v = param_h;

% Finding the synchronous solution 
y = zeros(2*M,1);
for j = 1:M
    Gj = J_0/(1+(s*x(j)^2));
	if j == 1
		y(j) = ((Gj-gamma)/tau_c)*x(j) + (kappa/tau_c)*A(j,:)*(x(1:M).*cos([0;x(M+2:2*M)] - x(M+1)*tau)); 
		y(M+j) = - x(M+1) + omega_v(j) +  (kappa/tau_c)*A(j,:)*((x(1:M)/x(j)).*sin([0;x(M+2:2*M)] - x(M+1)*tau));
	else
		y(j) = ((Gj-gamma)/tau_c)*x(j) + (kappa/tau_c)*A(j,:)*(x(1:M).*cos([0;x(M+2:2*M)] - x(M+j) - x(M+1)*tau));
		y(M+j) = - x(M+1) + omega_v(j) +  (kappa/tau_c)*A(j,:)*((x(1:M)/x(j)).*sin([0;x(M+2:2*M)] - x(M+j) - x(M+1)*tau));
    end
end
    
end


%%
function [J1,J2] = JacLRE_het(param,r,delta,G,Omega)
% This code calculates the Jacobians J1 and J2 for the heterogeneous sync
% state
%  
% param    - System parameters
% r        - vector of amplitudes
% delta    - Vector of phase shifts
% G        - Vector of gains
% Omega    - Frequency shift

% Parameters 
gamma = param.gamma;       
tau_c = param.tau_c;         
tau_f = param.tau_f;        
s = param.s;
kappa = param.kappa;       
A = param.A;
M = param.num_nodes;
phase_lag = delta;
tau = param.tau;


J1 = zeros(3*M);
J2 = zeros(3*M);

for kk = 1:M
    J1(kk,kk) = (G(kk)-gamma)/tau_c;          
    J1(kk,kk+M) = -(1/tau_c)*sum_sinLRE(M,A,r,phase_lag,Omega,tau,kk,1);
    J1(kk,kk+2*M) = r(kk)/tau_c;                                      
    
    J1(kk+M,kk) = -(1/tau_c)*sum_sinLRE(M,A,r,phase_lag,Omega,tau,kk,2);
    J1(kk+M,kk) = (1/tau_c)*sum_cosLRE(M,A,r,phase_lag,Omega,tau,kk);

    J1(kk+2*M,kk) = -(2/(tau_f))*s*r(kk)*G(kk);                 
    J1(kk+2*M,kk+2*M) = -(1/tau_f)*(s*(r(kk))^2 + 1);                   

    for kkk=1:M
        J2(kk,kkk) = (kappa/tau_c)*A(kk,kkk)*cos(Omega*tau + phase_lag(kk)-phase_lag(kkk));                      
        J2(kk,kkk+M) =(kappa/tau_c)*A(kk,kkk)*r(kkk)*sin(-Omega*tau+ phase_lag(kkk)-phase_lag(kk));             

        J2(kk+M,kkk) = (kappa/(tau_c*r(kk)))*A(kk,kkk)*sin(-Omega*tau + phase_lag(kkk)-phase_lag(kk));                       
        J2(kk+M,kkk+M) = (kappa/(tau_c*r(kk)))*A(kk,kkk)*r(kkk)*cos(Omega*tau + phase_lag(kk)-phase_lag(kkk));              
    end            
end
end

%%
function sum = sum_sinLRE(M,adj,r,phase_lag,Omega,tau,i,type)
    sum = 0;
    if type == 1
        for j=1:M
            sum = sum + adj(i,j)*r(j)*sin(-Omega*tau - phase_lag(i) +phase_lag(j));
        end
    elseif type == 2
        for j=1:M
            sum = sum + adj(i,j)*(r(j)/(r(i))^2)*sin(-Omega*tau - phase_lag(i) +phase_lag(j));
        end
    end
end

%%
function sum = sum_cosLRE(M,adj,r,phase_lag,Omega,tau,i)
    sum = 0;
    
    for j=1:M
    sum = sum + adj(i,j)*(r(j)/r(i))*cos(-Omega*tau - phase_lag(i) +phase_lag(j));
    end    
end

