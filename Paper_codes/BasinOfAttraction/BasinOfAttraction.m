%% Stability analysis of the Lang-Kobayashi laser model
% This code analyses the stability of coupled lasers as parameter disorder 
% is continually introduced in one of the 3:
%  - alpha (linewidth enhancement factor)
%  - kappa (coupling strength)
%  - omega (frequency detuning)
% This code requires loading a .mat file with the solution for the sync
% state and the associated heterogeneous parameter at each step and for
% each trial

M = param.M; 
%% Initialize variables 
Bfcoeff = Inf(trial,m);      % Initialize matrix for fraction of converging initial conditions
conv_class = cell(trial,m);  % Initialize cell for convergence of each trial of initial conditions
r_class = cell(trial,m);     % Initialize cell for amplitudes of each trial of initial conditions
Omega_class = cell(trial,m); % Initialize cell for frequency shifts of each trial of initial conditions

%% Main loop 
for jj = 1:trial
    disp(['jj= ',num2str(jj)])
    for ii = 1:m
        sol_now = Sol_m{jj,ii};
        if isempty(sol_now)
            continue
        end
        param_h_now = Param_het{jj,ii};
        r = sol_now(1:M); Omega = sol_now(M+1); delta = [0;sol_now(M+2:2*M)]; 
        if het_choice == "alpha"
            omega_h = param.omega0*ones(M,1);
            alpha_h = param_h_now;
            kappa_h = param.kappa*ones(M,1);
        elseif het_choice == "omega"
            omega_h = param_h_now;
            alpha_h = param.alpha0*ones(M,1);
            kappa_h = param.kappa*ones(M,1);
        elseif het_choice == "kappa"
            omega_h = param.omega0*ones(M,1);
            alpha_h = param.alpha0*ones(M,1);
            kappa_h = param_h_now;
        end
        
        trials = 1000; 
        seeds = randi([0 100000],1,trials);
        conv_m = zeros(trials,1);
        r_m = zeros(trials,M);
        Omega_m = zeros(trials,M);
        tic    
        parfor i = 1:trials
            % random seed
            rng(seeds(i))
            r_now = 5.2 + (5.4 - 5.2) * rand(M,1);  % Element-wise multiplication to apply the perturbation  
            Omega_now = -10 + (30 - (-10)) * rand(M,1); 
            delta = -pi + 2*pi*rand(M,1); 
            conv_m(i) = ICconvGen(adj,r_now,r,Omega,Omega_now,delta,alpha_h,omega_h,kappa_h,param);
            r_m(i,:) = r_now;
            Omega_m(i,:) = Omega_now;
        end 
        
        
        Bf = sum(conv_m)/trials;
        Bfcoeff(jj,ii) = Bf;
        conv_class{jj,ii} = conv_m;
        r_class{jj,ii} = r_m;
        Omega_class{jj,ii} = Omega_m;

    end       
end

filename = "filename.mat";
save(filename)

