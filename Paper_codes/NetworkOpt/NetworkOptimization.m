%% Optimization of Adjacency matrix to minimize the Lyapunov exponent
% This code optimizes the adjacency matrix with the constraint that the the
% sum of all the entries of the matrix remains constant

clear all; close all; clc;

%%  Defining the parameters  
N_0 = 1.5;     param.N_0 = N_0;         % Carrier number at transparency   (1.5e08, rescaled N_0 -> N_0*1e-8)
g = 1.5e03;    param.g = g;             % Differential gain coefficient    (1.5e-05, rescaled g -> g*1e8)
s = 1e-03;     param.s = s;             % Gain saturation coefficient      (1e-07, rescaled s -> s*1e4)
gamma = 500;   param.gamma = gamma;     % Cavity loss
alpha0 = 5;    param.alpha0 = alpha0;   % Linewidth enhancement factor
gamma_n = 0.5; param.gamma_n = gamma_n; % Carrier loss rate
sigma0 = 3;    param.sigma0 = sigma0;   % Variance of detuning
omega0 = 0;    param.omega0 = omega0;   % Frequency detuning   
a = 2.5455;    param.a = a;             % Pump factor   
tau = 0.15;    param.tau = tau;         % Time delay 
J_0 = a*gamma_n*(N_0+(gamma/g)); param.J_0 = J_0; % Pump current
M = 10;        param.M = M;             % Number of lasers
dx = 0.95;     param.dx = dx;           % Decaying factor
neig = 3;      param.neig = neig;       % Number of neighbors 

kappa = 3*0.475;

%% Defining topology of coupling 
topo = input('Topology: ','s') % Choice of topology
%       'topo = A' -> All-to-all topology
%       'topo = R' -> Ring topology
if topo == "R"
    %%%%--- Defining the ring network ---%%%%%
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
    x00 = [5.3393, 0];
elseif topo == "A"
    %%%%--- Defining the all-2-all network ---%%%%%
    adj = ones(M,M);
    sum_deg = 0;
    for jjj = 1:M
        sum_deg = sum_deg + adj(1,jjj);
    end
    deg = sum_deg;
    x00 = [5.3413, 0]; 
end

% Coupling parameters
param.adj = adj;     % Adjacency matrix
param.deg = deg;     % Degree (row-sum)
param.kappa = kappa; % Coupling strength

%% Finding the sync solution (transcendental equations)
x= fminsearch(@(x)sync(x,param),x00);

r00 = x(1);
omega00 = x(2);
N00 =((((g*r00^2)/(1+s*r00^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r00^2)/(1+s*r00^2))*N_0*1e-4);

r0_ini = zeros(2*M,1);
r0_ini(1:M) = r00*ones(M,1); r0_ini(M+1) = omega00; r0_ini(M+2:2*M) = zeros(M-1,1); 
r = r0_ini(1:M); Omega = omega00; delta = zeros(M,1); N = N00*ones(M,1);

% Stability of the initial configuration
[J1,J2] = Jac_freq_het_net(param,adj,r,Omega,delta,N,tau);
mtle_0 = dde_rightmost_eig(J1,J2,tau);

%% Defining the perturbation variables 
trial = 100; % Number of realizations of random heterogeneity
m = 500; % Maximum number of steps
sstep = 1; % Size of step 

%% Initialize variables 
mtle = Inf(trial,m);
numSteps = zeros(trial,1);

adj_best_v = cell(trial, m);
mtle_best_v = Inf(trial,m);

mtle_best_v(:,1) = mtle_0*ones(trial,1);
for jj = 1:trial
    adj_best_v{jj,1} = adj;
end



%% Main loop 
parfor jj = 1:trial
    disp(['jj= ',num2str(jj)])
    adj0 = adj;
    mtle_best = mtle_0;
    r0 = r0_ini;
    for ii = 2:m
        % Define the maximum number of retries and initial sstep value
        max_retries = 10; 
        retry_count = 0;
        success = false; % Flag to check if optimization is successful

        % a random realization of heterogeneity
        szhet = M*(M-1);
        sample = normrnd(0,1,szhet,1); 
        sample = (sample - mean(sample))/std(sample);
        sstep_current = sstep; % Start with the initial sstep value


        % Set optimization options
        optionss = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter','MaxIterations', 5);
    
        % Define the objective function
        obj_fun = @(x) ObjAdj(x,adj0,r0,param);

        % Loop to attempt optimization with adjusted sstep
        while ~success && retry_count < max_retries
            x0 = (sstep_current / (5 * norm(sample))) * sample; % Initial guess

            % Define the nonlinear constraints function
            nonlcon = @(x) constraints(x, sstep_current);

            % Run the optimization
            [x_opt, fval, exitflag] = fmincon(obj_fun, x0, [], [], [], [], -ones(M*(M-1), 1), [], nonlcon, optionss);
            success = true;
        end

        if ~success
            disp('Optimization failed after maximum retries. Moving to next iteration.');
            continue; % Skip to the next iteration of the loop if unsuccessful
        end 
          
        % Run the optimization
        [mtle,r0_new,adj_h] = ObjAdj_full(x_opt,adj0,r0,param);
        
        if abs(mtle - mtle_best)<1e-4
            disp("stopping because the change is too small")
            adj0 = adj_h;
            r0 = r0_new;
            mtle_best = mtle;
    
            mtle_best_v(jj,ii) = mtle_best;
            adj_best_v{jj,ii} = adj0;
            numSteps(jj) = ii;
            break
        end

        adj0 = adj_h;
        r0 = r0_new;
        mtle_best = mtle;

        mtle_best_v(jj,ii) = mtle_best;
        adj_best_v{jj,ii} = adj0;
        numSteps(jj) = ii;        
    end    
end

[mtle_best_tot,indmin] = min(mtle_best_v(:));
[jm, im] = ind2sub(size(mtle_best_v), indmin);
adj_best_tot = adj_best_v{jm,im};

% 
filename = "filename.mat";
save(filename)
runtime = toc

%% Defining optimization constraints
function [c, ceq] = constraints(x, sstep)
    % Inequality constraints
    % The norm of the vector x should be less than or equal to eps
    c = norm(x) - sstep;

    % Equality constraints (ceq)
    % The sum of the entries of the vector x should be zero
    ceq = sum(x);
end

%% Calculating Lyapunov exponent and new sync state for perturbed matrix
function [mtle,r0_new,adj_h] = ObjAdj_full(dif,adj0,r0,param)
% Inputs:
%        dif    ---- vector with nondiagonal perturbation
%        adj0   ---- initial (unperturbed) matrix
%        r0     ---- sync sol for the initial matrix
% Outputs:
%        mtle   ---- Lyapunov exponent of solution of perturbed matrix
%        r0_new ---- sync sol for the perturbed matrix
%        adj_h  ---- perturbed matrix

M = param.M;
g = param.g;
s = param.s;
gamma_n = param.gamma_n;
J_0 = param.J_0;
N_0 = param.N_0;
tau = param.tau;

adj_h = assemble_dif(dif,adj0,M);

r = r0(1:M); Omega = r0(M+1); Delta = r0(M+2:2*M);

found = 0;
options = optimoptions('fsolve','Display','off','MaxIterations',1000,'MaxFunctionEvaluations',100000);
for count = 1:10
    		[sol,~,~,~] = fsolve(@(x)sync_het_net(x,param,adj_h),[r;Omega;Delta],options);

            if norm(sol(1:M)-r)<sqrt(M)*1e-2 && abs(sol(M+1)- Omega)<1e-2
                found = 1;
                r = sol(1:M); Omega = sol(M+1); delta = [0;sol(M+2:2*M)]; Delta = sol(M+2:2*M);     
                break
            else
                found = 2;
            end
end

N = zeros(1,M);
for ind = 1:M
    N(ind) =((((g*(r(ind))^2)/(1+s*(r(ind))^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*(r(ind))^2)/(1+s*(r(ind))^2))*N_0*1e-4);
end

if found == 1 % Sync state identified successfully
    [J1,J2] = Jac_het_net(param,adj_h,r,Omega,delta,N,tau);
    mtle = dde_rightmost_eig(J1,J2,tau);
else % No sync state found
    mtle = Inf;
end

r0_new = zeros(2*M,1);
r0_new(1:M) = r; r0_new(M+1) = Omega; r0_new(M+2:2*M) = Delta;

end

%% Calculating perturbed matrix
function [adj_h] = assemble_dif(dif,adj0,M)

adjdif = zeros(M,M);
idx = 1;
for i = 1:M
    for j = 1:M
        if i == j
            continue
        end
        adjdif(i,j) = dif(idx);
        idx = idx+1;
    end
end

adj_h = adj0 + adjdif;
end

