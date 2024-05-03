function [dx,iter,Krylov_stat] = Newton_Iterative_Solver(NS,rhs,PS,maxit,tol,pl)
% ==========================================================================================================
% Newton Iterative Solver 
% ----------------------------------------------------------------------------------------------------------
% [dx,iter,instability,drop_direction] = Newton_Iterative_Solver(NS,rhs,PS,maxit,tol)
%
%                              INPUT: takes as an input all relevant 
%                              information needed  to build the semismooth Newton system at 
%                              the j-th iteration of the SSN solver as well as its preconditioner. 
%                              Then, it employs Preconditioned MINRES
%                              to solve the associated linear system.
%
%                              OUTPUT: it returns the SSN direction dx, the number of inner 
%                              iterations required to find it, as well as some flag variables 
%                              indicating ill-conditioning or whether the direction is inaccurate
%                              and should be dropped.
% 
% Author: Spyridon Pougkakiotis.
% __________________________________________________________________________________________________________
    accuracy_bound = 1e-3;
    % Krylov info: Krylov_stat(1) -> instability, Krylov_stat(2) -> drop direction, Krylov_stat(3) -> flag
    Krylov_stat = [false;false;0;0];
    tol = min(tol,1e-3);
    dx = zeros(NS.n,1);
    warn_stat = warning;
    warning('off','all');
    [lhs, flag, res, iter] = minres(@(x) AS_multiplier(x,NS), rhs, tol, maxit, @(x) Precond_Operator(x,PS));
    warning(warn_stat);
    if (pl >= 3)
        fprintf(NS.fid,['-------------------------------***Krylov method: MINRES***',...
                        '-----------------------------------------\n']);
        fprintf(NS.fid,['Krylov Iterations                     Krylov Flag          ',...
                       '   Residual                    \n']);
        fprintf(NS.fid,['%4d                                %4d                     ',...
                        '%9.2e              \n'],iter,flag,res);
        fprintf(NS.fid,['-----------------------------------------------------------',...
                        '----------------------------------------\n']);
    end
    Krylov_stat(3,1) = flag;
    Krylov_stat(4,1) = res;

    % If something went wrong, assume that the preconditioner is not good enough -> increase quality.
    if (flag > 0) 
        if (flag ~= 3)
            iter = maxit;
        end
        if (res > accuracy_bound)
            Krylov_stat(2,1) = true;
            return;
        elseif (flag == 2 || flag == 4 || flag == 5)
            Krylov_stat(1) = true;
            fprintf('Instability detected during the iterative method. flag = %d.\n',flag);
            return;
        end
    end
    % Check for ill-conditioning.
    if (nnz(isnan(lhs)) > 0 || nnz(isinf(lhs)) > 0 || (max(lhs) == 0 && min(lhs) == 0)) 
        Krylov_stat(1) = true;
        iter = maxit;
        fprintf('Instability detected during the iterative method.\n');
        return;
    end
    dx = lhs(1:NS.n,1);
end

