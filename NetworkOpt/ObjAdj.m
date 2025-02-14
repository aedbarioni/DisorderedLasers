function [mtle] = ObjAdj(dif,adj0,r0,param)
% This code calculates the Lyapunov exponent for a matrix perturbed by
% 'dif'
% dif    - vector with nondiagonal perturbation
% adj0   - initial (unperturbed) matrix
% r0     - sync sol for the initial matrix


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
    		[sol,~,~,~] = fsolve(@(x)sync_freq_hetnet(x,param,adj_h),[r;Omega;Delta],options);

            if norm(sol(1:M)-r)<sqrt(M)*1e-2 && abs(sol(M+1)- Omega)<0.2
                found = 1;
                r = sol(1:M); Omega = sol(M+1); delta = [0;sol(M+2:2*M)];     
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
    [J1,J2] = Jac_freq_het_net(param,adj_h,r,Omega,delta,N,tau);
    mtle = dde_rightmost_eig(J1,J2,tau);
else % No sync state found
    mtle = Inf;
end

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