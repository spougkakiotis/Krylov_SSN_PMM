function [x,iter,tol_achieved,total_Krylov_iters,num_of_factor] = SSN(NS,tol,maxit,pl)
% ======================================================================================================= 
% Semismooth_Newton: 
% ------------------------------------------------------------------------------------------------------- 
% x = Semismooth_Newton(NS,tol,maxit,method)
%     
% Takes as an input a MATLAB struct NS that containing relevant information needed 
% to build the semismooth Newton system corresponding to minimizing the augmented 
% Lagrangian with respect to x, the minimum accuracy tolerance 
% and the maximum number of SSN iterations, the version of SSN-PMM used, as 
% well as a structure for potential low-rank updates. It employs a semismooth 
% Newton iteration (given the current iterate (x_k,y_k,z_k)) and returns the  
% accepted "optimal" solution x alongside some performance metrics.
%
% 
% Author: Spyridon Pougkakiotis.
% _______________________________________________________________________________________________________ 
    n = NS.n;
    m = NS.m;
    l = NS.l;
    % Starting point for SSN -> x_k0 = x_k.
    x = NS.x;               
    % =================================================================================================== 
    % Set the semismooth Newton parameters
    % --------------------------------------------------------------------------------------------------- 
    % Maximum tolerance allowed when solving the corresponding linear systems.
    eta_1 = (1e-1)*tol;
    % Determines the rate of convergence of SSN (that is, the rate is (1+eta_2)).
    % The trade-off: for larger eta_2, SSN is faster but CG is slower. (eta_2 in (0,1].)
    eta_2 = 0.1;        
    % Fraction of the decrease in Lagrangian predicted by linear extrapolation that we accept.
    mu = (0.4995/2);  
    % Maximum possible step-length used within the backtracking line-search.
    delta = 0.995; 
    % ___________________________________________________________________________________________________ 
    
    % =================================================================================================== 
    % Initialize metric and preconditioning struct.
    % --------------------------------------------------------------------------------------------------- 
    % count overall Krylov iterates
    total_Krylov_iters = 0;    
    % keep track of SSN iterations
    iter = 0;       
    % keep track of the number of factorizations performed
    num_of_factor = 0; 
    % maximum number of outliers allowed before recomputing the preconditioner
    max_outliers = 0;  
    % maximum number of Krylov iterations
    maxit_Krylov = 150;    
    % A struct to store all information concerning the preconditioner
    PS = struct();
    PS.n = n;    PS.m = m;    PS.l = l;  
    prec_enabled = NS.prec_enabled;
    % ___________________________________________________________________________________________________ 
 
    while (iter < maxit)
    % --------------------------------------------------------------------------------------------------- 
    % SSN Main Loop structure:
    % Let L(x) be the proximal augmented Lagrangian associated with the subproblem of interest
    % Until (|| \nabla L(x_{k_j}) || <= tol) do
    %   Compute a Clarke subgradient M of \nabla L(x_{k_j}) and solve
    %               M dx = - \nabla L(x_{k_j})
    %   using a preconditioned Krylov subspace method.
    %   Perform Line-search with Backtracking
    %   j = j + 1;
    % End
    % --------------------------------------------------------------------------------------------------- 
        % =============================================================================================== 
        % Compute and store an element in the Clarke subdifferential of Proj_{K}(w)
        % evaluated at (beta^{-1} z_k + x_j). 
        % ----------------------------------------------------------------------------------------------- 
        w = (1/NS.beta).*(NS.z) + x;  
        w_lb = (w <= NS.lb);
        w_ub = (w >= NS.ub);
        B_delta = ones(n,1);
        B_delta(w_lb) = zeros(nnz(w_lb),1);
        w(w_lb) = NS.lb(w_lb);
        B_delta(w_ub) = zeros(nnz(w_ub),1);
        w(w_ub) = NS.ub(w_ub);
        NS.B_delta = B_delta;
        % _______________________________________________________________________________________________ 

        % =============================================================================================== 
        % Compute and store an element in the Clarke subdifferential of Proj_{C}(u)
        % evaluated at (y2(1:l) + beta(Cx+d)) (Nonseparable max{,0} terms)
        % ----------------------------------------------------------------------------------------------- 
        u = NS.y2(1:l,1) + NS.beta.*(NS.C*x + NS.d);
        u_lb = (u <= 0);
        u_ub = (u >= 1);
        B_h_1 = ones(l,1);
        B_h_1(u_lb) = zeros(nnz(u_lb),1);
        u(u_lb) = 0;
        B_h_1(u_ub) = zeros(nnz(u_ub),1);
        u(u_ub) = 1;
        NS.B_h_1 = B_h_1;
        active_C = (B_h_1 == 1);
        inactive_C = (B_h_1 == 0);
        % _______________________________________________________________________________________________ 

        % =============================================================================================== 
        % Compute and store an element in the Clarke subdifferential of Proj_{C}(v)
        % evaluated at (y2(l+1:l+n) + beta(Dx)) (Separable max{,0} terms)
        % ----------------------------------------------------------------------------------------------- 
        v = NS.y2(l+1:end,1) + NS.beta.*(NS.D.*x);
        v_lb = (v < 0);
        v_ub = (v > 1);
        B_h_2 = ones(n,1);
        B_h_2(v_lb) = zeros(nnz(v_lb),1);
        v(v_lb) = 0;
        B_h_2(v_ub) = zeros(nnz(v_ub),1);
        v(v_ub) = 1;
        NS.B_h_2 = B_h_2;
        % _______________________________________________________________________________________________ 

        % Coefficient matrix in the (2,1) block of the augmented system
        NS.G = [NS.A; -(NS.C_tr(:,active_C))'];
        % =============================================================================================== 
        % Compute the right-hand-side and check the termination criterion
        % -----------------------------------------------------------------------------------------------
        rhs = [NS.c + NS.Q*x + NS.D.*(v) + (NS.C_tr(:,inactive_C)*u(inactive_C))...
                + NS.z + NS.beta.*x - NS.beta.*w + (1/NS.rho).*(x-NS.x);
               (1/NS.beta).*NS.y1 - NS.A*x + NS.b;
               (1/NS.beta).*u(active_C)];          
        NS.Schur_size = size(rhs,1) - n;
        [~,residual] = Evaluate_Prox_Lagrangian(NS,x,true);
        res_error = norm(residual);
        % Stop if the desired accuracy is reached.
        if (res_error < tol)     
            break;
        end

        iter = iter + 1;
        if (pl > 1)
            if (iter == 1)
                fprintf(NS.fid,['__________________________________________________________',...
                                '_________________________________________\n']);
                fprintf(NS.fid,['___________________________________Semismooth Newton method',...
                               '________________________________________\n']);
            end
            fprintf(NS.fid,['SSN Iteration                                   ',...
                            '  Residual Infeasibility                    \n']);
            fprintf(NS.fid,['%4d                                              ',...
                            '%9.2e                      \n'],iter,res_error);
        end
        % _______________________________________________________________________________________________ 

        % =============================================================================================== 
        % Check if we need to re-compute the preconditioner. If so compute and store its factors in PS
        % ----------------------------------------------------------------------------------------------- 
        if (iter == 1)
            % true if the preconditioner must be recomputed
            update_preconditioner = true; 
        elseif ((nnz(PS.B_delta - B_delta) + nnz(PS.B_h_2 - B_h_2) > max_outliers) ...
                || nnz(PS.active_C - active_C))
            update_preconditioner = true;
        else
            update_preconditioner = false;
        end
        % update the (1,1) diagonal block anyway!
        PS.H_tilde = NS.Q_diag + (NS.beta + (1/NS.rho)).*ones(n,1)  ...
                     - NS.beta.*B_delta + NS.beta.*(((NS.D).^2).*B_h_2); 
        if (update_preconditioner)
            PS.B_delta = B_delta;   PS.B_h_2 = B_h_2;   PS.active_C = active_C;
            idx_set = ((B_delta == 1) & ((B_h_2 == 0 | NS.D == 0)));
            card_idx_set = nnz(idx_set);
            if (card_idx_set && prec_enabled)
               num_of_factor = num_of_factor + 1;
               Schur_tilde = NS.G(:,idx_set)*...
                             (spdiags((1./PS.H_tilde(idx_set)), 0, card_idx_set, card_idx_set)*...
                             (NS.G(:,idx_set)')) + (1/NS.beta).*speye(NS.Schur_size);
            else
                Schur_tilde = (1/NS.beta).*speye(NS.Schur_size);
            end
            % Cholesky factorization of the approximated Schur complement
            [PS.L_S,chol_flag,PS.Perm] = chol(Schur_tilde,'lower','vector');
            if (chol_flag)
                fprintf("Numerical instability in the Cholesky decomposition of the preconditioner.\n");
                NS.beta = NS.beta*0.5;
                NS.rho = NS.rho*0.5;
                continue;
            end
            PS.PermInv = 1:NS.Schur_size;
            PS.PermInv(PS.Perm) = 1:NS.Schur_size;
        end
        % _______________________________________________________________________________________________ 
        
        % =============================================================================================== 
        % Call MINRES using the previous preconditioner to approximately solve the SSN sub-problem.
        % ----------------------------------------------------------------------------------------------- 
        Krylov_tol = max(min(min(min(eta_1,res_error^(1+eta_2)),1),1e-10),1e-12);
        % Call preconditioned MINRES to solve the system
        [dx,Krylov_iter,Krylov_stat] = Newton_Iterative_Solver(NS,rhs,PS,maxit_Krylov,Krylov_tol,pl);
        instability = Krylov_stat(1,1);
        drop_direction = Krylov_stat(2,1);
        Krylov_flag = Krylov_stat(3,1);
        Krylov_res = Krylov_stat(4,1);
        if (Krylov_iter > 100 || (Krylov_flag && Krylov_res > 1e-8))
            prec_enabled = true;
        end
        total_Krylov_iters = Krylov_iter + total_Krylov_iters;
        % If something went wrong, decraese the penalty parameters, drop the direction and re-solve.
        if (drop_direction || instability) 
            NS.beta = NS.beta*0.1;
            NS.rho = NS.rho*0.1;
            continue;
        end
        % _______________________________________________________________________________________________ 
        if (iter == 1)
          alpha = 0.995;
        else
            alpha = Backtracking_Line_Search(NS,x,dx,mu,delta);
        end
        x = x + alpha.*dx;
    end
    tol_achieved = norm(residual);
    if (pl > 1)
        fprintf(NS.fid,['______________________________________________________',...
                        '_____________________________________________\n']);
    end
end 
% ******************************************************************************************************* 
% END OF FILE.
% ******************************************************************************************************* 
