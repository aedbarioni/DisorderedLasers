function mtle = dde_rightmost_eig(J1,J2,tau)
% This code calculate the rightmost eigenvalue of the generalized
% characteristic equation. Notice that this calls a code that is a modified
% version of one of the codes in the DDE-BIFTOOL (stst_stabil_mod)
%  
% J1,J2    - Jacobian matrices
% tau      - Time delay

if tau~=0
    flag_newhheur = 1;
    method = df_mthod('stst',flag_newhheur);
    method.stability.minimal_real_part = -10;
    
    stability = stst_stabil_mod(J1,J2,tau,method);
    
    mtle = unique(max(real(stability.l1)));
    if any(size(mtle) == 0)
        mtle = unique(max(real(stability.l1)));
        if any(size(mtle) == 0)
            mtle = 0;
        end
    end
end

if tau == 0 
    e = eig(J1+J2);
    e_max = max(real(e));
    % This part of the code disconsider a 0 eigenvalue (only keep for
    % systems where this is needed)
    % e_max_index = find(real(e)==e_max);
    % if abs(e_max)<1e-4
    %     e(e_max_index) = [];
    %     e_max = max(real(e));
    %     % e_max_index = find(abs(real(e)-e_max)<1e-4);
    % end
    % e_lambda = e(e_max_index);
    % mtle = unique(real(e_lambda));
    mtle = e_max;
end