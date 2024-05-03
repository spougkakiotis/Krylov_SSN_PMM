function [res_p,res_d,compl] = compute_residual(Q,C,A,b,d,c,lb,ub,x,y1,y2,z)
% ============================================================================================== 
% This function takes the problem data as well as the current iterate as input, and outputs the 
% l2-norm of the scaled primal, dual residuals, as well as complementarity.
% ---------------------------------------------------------------------------------------------- 
    % Dual residual
    res_d = norm(c + Q*x - A'*y1 + C'*y2 + z)/(norm(c)+1);

    % Primal residual (Projection onto M due to the max{,0} terms in the objective)
    temp_res_p = y2 + C*x + d;
    temp_idx_u = (temp_res_p > 1);
    temp_idx_l = (temp_res_p < 0);
    temp_res_p(temp_idx_u) = 1;
    temp_res_p(temp_idx_l) = 0;
    res_p = norm([A*x-b; y2 - temp_res_p])/(norm([b;d])+1);
    
    % Complementarity
    temp_compl = x + z;
    temp_lb = (temp_compl < lb);
    temp_ub = (temp_compl > ub);
    temp_compl(temp_lb) = lb(temp_lb);
    temp_compl(temp_ub) = ub(temp_ub);
    compl = norm(x -  temp_compl);                                     
% ______________________________________________________________________________________________ 
end
% ********************************************************************************************** 
% END OF FILE.
% ********************************************************************************************** 
