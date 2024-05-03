function [beta,rho] = update_PMM_parameters(res_p,res_d,new_res_p,new_res_d,beta,rho,reg_limit)
% ================================================================================================ 
% This function takes as input the current and the last step's primal and dual 
% infeasibilities as well as the current penalty parameters of PMM (alongside the penalty limit)
% and updates the penalty parameters.
%
% If the overall primal and dual residual error is decreased, 
% we increase the penalty parameters aggressively.
% If not, we continue increasing the penalty parameters slower,  
% limiting the increase to the value of the regularization threshold.
% ------------------------------------------------------------------------------------------------ 
    cond_p = 0.95*res_p > new_res_p;
    cond_d = 0.95*res_d > new_res_d;
    if (cond_p || cond_d)
        % Standard PMM (heuristic)  
        beta = min(reg_limit,beta*1.2);              
    else
        beta = min(reg_limit,beta*1.1);
    end
    if (cond_d || cond_p)
        rho = min(1e2*reg_limit,rho*1.4);
    else
        rho = min(1e2*reg_limit,rho*1.1);
    end
% ________________________________________________________________________________________________ 
end

