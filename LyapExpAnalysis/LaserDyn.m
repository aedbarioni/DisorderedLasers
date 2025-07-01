function [stdFFT,Ptot,OrderCoh] = LaserDyn(r,Omega,delta,N,param_h,het_choice,param)
% This code solves the laser dynamics and calculates the following
% quantities
%   - The standard deviation of the power spectrum 'stdFFT' 
%   - The power of the combined beam 'Ptot'
%   - The order coherence of the combined beam 'OrderCoh'
%  
% r           - Vector of amplitudes of the synchronized state
% Omega       - Frequency shift of the synchronized state
% delta       - Phase mismatch of the synchronized state
% N           - Carrier number of the synchronized state
% param_h     - heterogeneous parameter
% het_choice  - Parameter chosen to be heterogeneous
% param       - System parameters

M = param.M;
tau = param.tau;

% Introduce heterogeneity into chosen parameter 
if het_choice == "alpha"
    omega_h = param.omega0*ones(M,1);
    alpha_h = param_h;
    kappa_h = param.kappa*ones(M,1);
elseif het_choice == "omega"
    omega_h = param_h;
    alpha_h = param.alpha0*ones(M,1);
    kappa_h = param.kappa*ones(M,1);
elseif het_choice == "kappa"
    omega_h = param.omega0*ones(M,1);
    alpha_h = param.alpha0*ones(M,1);
    kappa_h = param_h;
end

% Simulating the system for T time units
T = 300;
time = [0 T];
opts = ddeset('RelTol',1e-3,'AbsTol',1e-6, 'MaxStep', 1e-2);

% Perturbation added to the sync state at t=0
perturb = 1e-5*rand(3*M,1); param.perturb = perturb;

% Initial conditions to solve dde
rhs = @(t,x,Z)rhs_aux(t,x,Z,alpha_h,omega_h,kappa_h,param);
rhshist = @(t)rhshist_aux(t,r,Omega,N,delta,param);

% solving the DDE
sol = dde23(rhs,tau,rhshist,time,opts);

rs = sol.y(1:M,:);
phi = sol.y(M+1:2*M,:);
Ns = sol.y(2*M+1:3*M,:);
freq = sol.yp(M+1:2*M,:);
phases = wrapToPi(phi);

% End
endind = size(phases,2);
iniind = endind-702;
Freq = zeros(10,702);
Efield_end = zeros(10,702);
ComplexAngle = zeros(10,702);

for ii=1:M
    for jj = 1:length(Freq)
        real = rs(ii,jj+iniind)*1e2*cos(phases(ii,jj+iniind));
        imag = rs(ii,jj+iniind)*1e2*sin(phases(ii,jj+iniind));
        Efield_end(ii,jj) = complex(real,imag);
        ComplexAngle(ii,jj) = complex(cos(phases(ii,jj+iniind)),sin(phases(ii,jj+iniind)));
    end
end

x_end = sol.x(iniind+1:endind);
EfieldSum_end = sum(Efield_end, 1);
fft_end = fft(EfieldSum_end');

% std FFT 
[p, f] = pspectrum(EfieldSum_end', 'TwoSided', true); % To plot power spectrum, just call the function (remove "[p, f] = ")
p_normalized = p / sum(p);
mean_freq = sum(f .* p_normalized);
stdFFT = sqrt(sum((f - mean_freq).^2 .* p_normalized));


% Total Power
Ptot = (abs(EfieldSum_end(end,1)))^2;

% Order Coherence 
OrderCoh = (abs(sum(ComplexAngle(end,1))))^2;


end

%% Defining the DDEs 
function dxdt = rhs_aux(t,x,Z,alpha_h,omega_h,kappa_h,param)  
    % Parameters
    N_0 = param.N_0;       
    g = param.g;       
    s = param.s;        
    gamma = param.gamma;        
    gamma_n = param.gamma_n;     
    sigma0 = param.sigma0;       
    J_0 = param.J_0;
    M = param.M; 
    adj = param.adj; 
    neig = param.neig;

	dxdt = zeros(3*(M),1);
	xlag = Z(:,1);

    % Noise term (Gaussian noise)
    noise_std = 0; % standard deviation of the noise
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
function S = rhshist_aux(t,r,Omega,N,delta,param)	
    % Parameters 
    M = param.M; 
    perturb = param.perturb;

	histr = r;
	histphi = Omega*t*ones(M,1) + delta;
    histN = N';
	S = [histr;histphi;histN] + perturb*heaviside(t);
end
