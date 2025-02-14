%% MSF dependence of coupling strength and degree
% This code calculates the MSF and determine the deepness and size of the
% stable region in the complex plane for varying values of coupling
% strength and degree

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

% Defining topology of coupling (all-2-all)
adj = ones(M,M);  % Adjacency matrix
coup_str = 4.75;  % Coupling strength
deg = 10;         % Degree 
neig = 10;        % Number of neighbors

% Coupling parameters     
param.adj = adj;
param.deg = deg; 
param.neig = neig;

m = 101;
re_fin = 50;
im_fin = 25;
deg_v = linspace(0,30,m);
kappa_v = linspace(0,30,m);
mtle = zeros(m,m);
sizeXi_m = zeros(m,m);
deepXi_m = zeros(m,m);
x00 = [5.3413, 0];

for ii = 1:m
    ii
    x0 = x00;
    for jj = 1:m
        deg = deg_v(ii);
        coup_str = kappa_v(jj);

        x= fminsearch(@(x)sync(x,deg,coup_str,param),x00);
        r = x(1);
        Omega = x(2);
        N =((((g*r^2)/(1+s*r^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r^2)/(1+s*r^2))*N_0*1e-4);
        
        M1 = zeros(3,3);
        M2 = zeros(3,3);
        M3 = zeros(3,3);

        % Defining Df matrix
        M1(1,1) = (1/2)*((g*(N-N_0)/((1+s*r^2)^2))*(1-s*r^2) - gamma);
        M1(1,3) = (g/2)*(r/(1+s*r^2));
        M1(2,1) = -(alpha0*g*s*r)*((N-N_0)/((1+s*r^2)^2));
        M1(2,3) = (alpha0*g/2)*(1/(1+s*r^2));
        M1(3,1) = -2*g*r*((N-N_0)/((1+s*r^2)^2))*1e-4;
        M1(3,3) = -(gamma_n + 1e-4*g*((r^2)/(1+s*r^2)));
        
        % Defining D(t)h
        M2(1,2) = -r*sin(Omega*tau);
        M2(2,1) = (1/r)*sin(Omega*tau);
        M2(2,2) = -cos(Omega*tau);
        
        % Defining D(t-tau)h
        M3(1,1) = cos(Omega*tau);
        M3(1,2) = r*sin(Omega*tau);
        M3(2,1) = -(1/r)*sin(Omega*tau);
        M3(2,2) = cos(Omega*tau);

        J1 = M1 + (coup_str/neig)*deg*M2;
        J2 = deg*(coup_str/neig)*M3;
        mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau);
        [sizeXi,deepXi] = XiChar(deg,coup_str,neig,tau,M,M1,M2,M3);
        sizeXi_m(jj,ii) = sizeXi;
        deepXi_m(jj,ii) = deepXi;
        x0 = [r,Omega];   
    end
    filename = "filename.mat";
    save(filename)
end


%% Plotting stability for varying kappa and indegree
for i = 1:m
    for j = 1:m
        if mtle(i,j)<0
            mtle(i,j)=2;
        elseif mtle(i,j)>0
            mtle(i,j)=1;
        elseif abs(mtle(i,j))<1e-5
            mtle(i,j)=2;
        end
    end
end

mtle(1,:) = ones(1,m);
mtle(:,14) = ones(m,1);

c2 = [184,0,0]/255;
c1 = [0,85,212]/255;
custom_colormap = [c2; c1];

figure();
imagesc(mtle(end:-1:1,14:1:end)) 
xlabel('indegree')
ylabel('\kappa')
title('Lyapunov Exponent')
colormap(custom_colormap);
xticks(1:(m/(15)):m); xticklabels(1:1:15);
yticks(1:(m/(15)):m); yticklabels(15:-1:0);

%% Plotting figure for deepness and size of stability region 
figure();
plot(kappa_v,10*sizeXi_m(:,1),'LineWidth',2.0,'color',[0.4660 0.6740 0.1880]);
hold on
plot(kappa_v,deepXi_smooth(:,1),'LineWidth',2.0,'color',[0.8500 0.3250 0.0980]);
yline(0, 'k', 'LineWidth', 0.5); 
xline(0.5, 'k', 'LineWidth', 0.5); 
xline(4.75, 'k', 'LineWidth', 0.5); 
xline(7, 'k', 'LineWidth', 0.5); 
xlabel('\kappa'); ylabel('\Lambda')
xlim([0,10]) ; ylim([-10,10])
hold off 

%% Calculating the deepness and size of the stable region of the MSF in the complex plane
function [sizeXi,deepXi] = XiChar(deg,coup_str,neig,tau,M,M1,M2,M3)
    J1 = M1 + (coup_str/neig)*deg*M2;

    m = 51;
    re_fin = 25;
    im_fin = 10;
    Re_alpha = linspace(-re_fin,re_fin,m);
    Im_beta = linspace(-im_fin,im_fin,m);
    mtle = zeros(m,m);
    
    for ii = 1:m    
        for jj = 1:m
            J2 = (Re_alpha(ii) + 1i*Im_beta(jj))*(coup_str/neig)*M3;
            if (ii == 36 && jj ==27)
                disp("stop")
            end
            mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau,M);
        end
    end
    deepXi = min(min(mtle));
    sizeXi = sum(mtle(:) < 0)/(m^2);

end

%% Finding the synchronous solution 
function y = sync(x,kappa,deg,param)
N_0 = param.N_0;        
g = param.g;       
s = param.s;        
gamma = param.gamma;      
alpha0 = param.alpha0;       
gamma_n = param.gamma_n;     
sigma0 = param.sigma0;       
omega0 = param.omega0;       
a = param.a;            
tau = param.tau;          
neig = param.neig;

F = zeros(2,1);
    
    F(2) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(1))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + (kappa/neig)*deg*cos(x(2)*tau);
    F(1) = -x(2) + (alpha0/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(1))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega0 - (kappa/neig)*deg*sin(x(2)*tau);   

    y = log10(norm(F));
end

