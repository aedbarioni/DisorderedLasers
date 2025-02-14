function [] = Multistability(tau,kappa, topo)
% This code calculates all the solutions of the transcendental equations  
% and determines it's stability (Lyapunov exponent)
%  
% tau     - Time delay
% kappa   - Coupling strength
% topo    - Topology
%       'topo = A' -> All-to-all topology
%       'topo = D' -> Decaying topology
%       'topo = R' -> Ring topology

% Defining the parameters 
N_0 = 1.5;     param.N_0 = N_0;         % Carrier number at transparency   (1.5e08, rescaled N_0 -> N_0*1e-8)
g = 1.5e03;    param.g = g;             % Differential gain coefficient    (1.5e-05, rescaled g -> g*1e8)
s = 1e-03;     param.s = s;             % Gain saturation coefficient      (1e-07, rescaled s -> s*1e4)
gamma = 500;   param.gamma = gamma;     % Cavity loss
alpha0 = 5;    param.alpha0 = alpha0;   % Linewidth enhancement factor
gamma_n = 0.5; param.gamma_n = gamma_n; % Carrier loss rate
sigma0 = 3;    param.sigma0 = sigma0;   % Variance of detuning
omega0 = 0;    param.omega0 = omega0;   % Frequency detuning   
a = 2.5455;    param.a = a;             % Pump factor   
J_0 = a*gamma_n*(N_0+(gamma/g)); param.J_0 = J_0; % Pump current
M = 10;        param.M = M;             % Number of lasers
dx = 0.95;     param.dx = dx;           % Decaying factor
neig = M;      param.neig = neig;       % Number of neighbors 



% Defining topology of coupling
if topo == "D"
    % Defining Decaying coupling network 
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

% Coupling parameters
param.adj = adj;     % Adjacency matrix
param.deg = deg;     % Degree (row-sum)

% Defining the perturbation variasyncbles 
m = 510; 
Omega_c = cell(m,1); % Cell containing the values of frequency shift
r_c = cell(m,1);     % Cell containing the values of amplitude

% Main loop 
r_v= zeros(1000,1);
Omega_v= zeros(1000,1);
for i = 1:1000
    r_0 = 5.4255 + 10*(rand -0.5);  % Random initial amplitude
    Omega_0 = 200*(rand -0.5); % Random initial frequency shift
    x0_trial = [r_0,Omega_0]; % Initial condition
    x_cand = fminsearch(@(x)sync_hom_var(x,param,kappa,"kappa"),x0_trial);
    z_cand = sync_hom_var(x_cand,param,kappa,var_choice);

    if z_cand <-12 % Threshold to accept solution
        x_now = x_cand;
        r_v(i) = x_now(1); 
        Omega_v(i) = x_now(2);             
    end
end

% Combine the two vectors into a matrix where rows represent the pairs
combined = [r_v(:), Omega_v(:)];

% Find unique rows while preserving order ('stable')
[unique_combined, idx] = unique(combined, 'rows', 'stable');

% Extract the unique r_v and Omega_v values
r_c{ii,1} = unique_combined(:, 1);
Omega_c{ii,1} = unique_combined(:, 2);

% Extract current frame's data
r_plot = r_c{ii};
Omega_plot = Omega_c{ii};

% Exclude (0,0) points for the current frame
non_zero_idx = ~(r_plot == 0 & Omega_plot == 0);
r_non_zero = r_plot(non_zero_idx);
Omega_non_zero = Omega_plot(non_zero_idx);

% Calculate stability of each point
for i = 1:length(r_non_zero)
    r_now = r_non_zero(i);     
    Omega_now = Omega_non_zero(i); 
    N_now = ((((g*r_now^2)/(1+s*r_now^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r_now^2)/(1+s*r_now^2))*N_0*1e-4);
    r = repmat(r_now,M,1); Omega = Omega_now; N = repmat(N_now,M,1); delta = zeros(M,1);
    [J1,J2] = Jac_hom(param,kappa,var_choice,r,Omega,delta,N,tau);
    mtle_fin(i) = dde_rightmost_eig(J1,J2,tau);
end

% Plot points with Z > 0 in red
idx_positive = mtle_fin > 0;

% Plot points with Z <= 0 in blue
idx_negative = mtle_fin <= 0;

filename = "filename.mat";
save(filename)

figure();
plot(r_non_zero(idx_positive), Omega_non_zero(idx_positive),'o', 'MarkerFaceColor',[0.6350 0.0780 0.1840], 'MarkerEdgeColor','none','MarkerSize', 5);
hold on;
plot(r_non_zero(idx_negative), Omega_non_zero(idx_negative), 'o', 'MarkerFaceColor',[0 0.4470 0.7410], 'MarkerEdgeColor','none', 'MarkerSize', 5);
yline(0, 'k', 'LineWidth', 0.5); 
hold off;