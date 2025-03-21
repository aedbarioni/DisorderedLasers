%% Stability analysis of the Lang-Kobayashi laser model for disorder in multiple parameters
% This code analyses the stability of coupled lasers as parameter disorder 
% is continually introduced in a combination of the parameters:
%  - alpha (linewidth enhancement factor)
%  - kappa (coupling strength)
%  - omega (frequency detuning)

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
neig = M;      param.neig = neig;       % Number of 
disp("updated_paper3")

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
    x0 = [5.3412, -7.1611]; 
    kappa = 5.38; 
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
    x0 = [5.3393, -6.9461];     
    kappa = 4.33;
elseif topo == "A"
    % Defining all-2-all network 
    adj = ones(M,M);
    sum_deg = 0;
    for jjj = 1:M
        sum_deg = sum_deg + adj(1,jjj);
    end
    deg = sum_deg;
    x0 = [5.3413, -7.1731]; 
    kappa = 4.75; 
end

% Coupling parameters
param.adj = adj;     % Adjacency matrix
param.deg = deg;     % Degree (row-sum)
param.kappa = kappa; % Coupling strength

%% Finding the sync solution (transcendental equations)
x= fminsearch(@(x)sync(x,param),x0);

r0 = x(1);
omega = x(2);
N =((((g*r0^2)/(1+s*r0^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r0^2)/(1+s*r0^2))*N_0*1e-4) ;


%% Parameter to be disordered
het_choice = input('Heterogeneous parameter: ','s');
het_choice2 = input('Heterogeneous parameter2: ','s');
    % 'alpha' - heterogeneity will be introduced in the linewidth factor 
    % 'kappa' - heterogeneity will be introduced in the coupling strength
    % 'omega' - heterogeneity will be introduced in the frequency detuning

%% Defining the perturbation variables 
trial = 10000;                    % Number of realizations of random heterogeneity
sigma_end = input('sigma_end: '); % Maximum value of sigma_p
m = 500*sigma_end;                % Step-size for analytical continuation
sigma = linspace(0,sigma_end,m);  % Level of heterogeneity measured by standard deviation;

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
    sample1 = normrnd(0,1,M,1);
	sample1 = (sample1 - mean(sample1))/std(sample1);
	sample2 = normrnd(0,1,M,1);
	sample2 = (sample2 - mean(sample2))/std(sample2);

    % Initialize sync state
	r = repmat(r0,M,1); Omega = omega; Delta = zeros(M-1,1);
    found = 0;
    negPar = 0;
    
    for ii = 1:m
        if found == 2
            disp("Algorithm didn't find sync solution")            
            break
        end

        disp(['ii= ',num2str(sigma(ii))])
        hetero = sigma(ii); % Level of heterogeneity
        dif = sample*hetero;
		dif1 = sample1*hetero;
		dif2 = sample2*hetero;
        
        % Adding heterogeneity to the specified parameter 1
        if het_choice == "alpha"
            param_h = alpha0*(ones(M,1)+dif);
        elseif het_choice == "omega"            
            param_h = omega0*ones(M,1) + dif;        
        elseif het_choice == "kappa"
            param_h = kappa*(ones(M,1)+dif);            
        end

        if het_choice ~= "omega"
            for j = 1:M
                if param_h(j)<0
                    negPar = 1;
                end
            end  
        end
        
        % Adding heterogeneity to the specified parameter 2
        if het_choice2 == "alpha"
            param_h2 = alpha0*(ones(M,1)+dif1);
        elseif het_choice2 == "omega"            
            param_h2 = omega0*ones(M,1) + dif1;
        elseif het_choice2 == "kappa"
            param_h2 = kappa*(ones(M,1)+dif1);            
        end

        if het_choice2 ~= "omega"
            for j = 1:M
                if param_h2(j)<0
                    negPar = 1;
                end
            end  
        end

        if negPar ==1
            disp("One of the parameters is negative")
            break
        end

        % Finding the synchronous solution 
        tol = 1e-15;
		option = optimset('Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolFun',tol,'TolX',tol);
		options = optimoptions('fsolve','Display','off','MaxIterations',1000,'MaxFunctionEvaluations',100000);
        for count = 1:10
    		[sol,fval,exitflag,output] = fsolve(@(x)sync_het_multi(x,param,param_h,het_choice,param_h2,het_choice2),[r;Omega;Delta],options);

            if norm(sol(1:M)-r)<sqrt(M)*1e-2 && abs(sol(M+1)- Omega)<1e-2
                found = 1;
                r = sol(1:M); Omega = sol(M+1); delta = [0;sol(M+2:2*M)]; Delta = sol(M+2:2*M);              
                break
            else
                found = 2;
            end
        end

        alpha_v = alpha0*ones(M,1);
        omega_v = omega0*ones(M,1);        
        kappa_v = kappa*ones(M,1);

        if het_choice == "alpha"
            alpha_v = param_h;
        elseif het_choice == "omega"
            omega_v = param_h;                   
        elseif het_choice == "kappa"
            kappa_v = param_h;            
        end

        if het_choice2 == "alpha"
            alpha_v = param_h2;
        elseif het_choice2 == "omega"
            omega_v = param_h2;        
        elseif het_choice2 == "kappa"
            kappa_v = param_h2;            
        end

        N = zeros(1,M);
        for ind = 1:M
            N(ind) =((((g*(r(ind))^2)/(1+s*(r(ind))^2))*1e-4 + gamma_n)^(-1))*(a*gamma_n*(N_0+(gamma/g)) + ((g*(r(ind))^2)/(1+s(ind)*(r(ind))^2))*N_0*1e-4);       
        end
        
        if found == 1 % Sync state identified successfully
            [J1,J2] = Jac_het_multi(param,param_h,het_choice,param_h2,het_choice2,r,Omega,delta,N,tau);
            mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau);
            if  mtle(jj,ii)<0
                negativity(jj,ii) = 1; 
                Trial_neg(jj) = 1;
            end
        else % No sync state found
            mtle(jj,ii) = Inf;
        end
    
    end
end

disp(mean(Trial_neg));
filename = "filename.mat";
save(filename)

%% Plotting all curves 
figure();
plot(sigma,mtle(1,:))
hold on
for iik= 2:trial
    plot(sigma,mtle(iik,:))
end
xlabel('\sigma_h'); ylabel('\Lambda')
hold off

%% Plotting median 
figure();
Median = median(mtle);
plot(sigma,Median,'LineWidth',2.0)
xlabel('\sigma_h'); ylabel('\Lambda')

%% Plotting Fraction of stabilized systems 
figure();
Negativ = mean(negativity);
plot(sigma,Negativ,'LineWidth',2.0)
xlabel('\sigma_h'); ylabel('fraction of trials with negative \Lambda')
