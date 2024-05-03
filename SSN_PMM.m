function [solution_struct] = SSN_PMM(pb_struct)
% ============================================================================================================ 
% This function is a primal-dual Semismooth Newton Proximal Method of Multipliers, suitable for solving convex 
% optimization problems of the following form
%
%                            min   (1/2) x^t Q x + c^t x + \sum_{i=1}^l max((Cx+d)_i,0) + \|Dx\|_1,
%    (P)                     s.t.  A x = b,
%                                  lb <= x <= ub,
%
% where D is a diagonal weight matrix. The method returns an approximate primal-dual solution,
% or a message indicating that the optimal solution was not found. 
% 
% The dual of the problem reads
%
%
%    (D)                    max    b^t y1 + d^t y2 - (1/2) x^t Q x - delta*_K(z) - delta_M(y2),
%                           s.t.   c + Qx - A^t y1 + C^t y2 + z = 0.
%    
% where M = {x in R^n | x_i \in [0,1], i \in {1,...,n}}. 
% ============================================================================================================
% INPUT: pb_struct
% ------------------------------------------------------------------------------------------------------------
% A MATLAB struct as an input which contains the following fields (see also (P) above):
%       .Q     -> the coefficient matrix of the quadratic in the objective 
%                 *** if empty, provide as sparse(n,n) ***
%       .A     -> the linear equalities coefficient matrix 
%                 *** if empty, provide as sparse(0,n) ***
%       .C     -> the coefficient matrix appearing within the max{,0} term in the objective 
%                 *** if empty, provide as sparse(0,n) ***
%       .c     -> the linear coefficients of the objective function (Default: all zeros)
%       .b     -> the right hand side of the linear equalities (Default: all zeros)
%       .d     -> the constant displacement in max terms in the objective (Default: all zeros)
%       .lb    -> the lower bound vector on the primal variables x (Default: all -Inf)
%       .ub    -> the upper bound vector on the dual variables x (Default: all Inf)
%       .D     -> weight vector for possible ell-1 norm in the objective (Default: all zeros)
%       .tol   -> the specified tolerance for termination (Default: 1e-4)
%       .maxit -> the maximum allowed number of PMM iterations (Default: 200)
%       .pl    -> 0 for not printing intermediate iterates
%                 1 for printing only PMM iterates (Default)
%                 2 for additionally printing SNN iterates
%                 3 for additionally printing Krylov iterate info
%       .p_fid -> file ID to print output (Default: 1 (prints on workspace)).
% ------------------------------------------------------------------------------------------------------------
%
% ============================================================================================================
% OUTPUT: solution_struct
% ------------------------------------------------------------------------------------------------------------
% A MATLAB struct containing the following fields:
%       .opt         -> a integer variable indicating the termination status:
%                         status = 0 <=> "optimal solution found"
%                         status = 1 <=> "maximum number of iterations reached"
%                         status = 2 <=> "termination due to numerical errors"
%       .x           -> Optimal primal solution
%       .y1          -> Lagrange multiplier vector corresponding to equality constraints
%       .y2          -> Lagrange multiplier vector corresponding to the max terms (see (D) above)
%       .z           -> Lagrange multiplier vector corresponding to box constraints on x
%       .PMM_iter    -> number of PMM iterations to termination
%       .SSN_iter    -> number of SSN iterations to termination
%       .Krylov_iter -> number of Krylov iterations to termination
%       .num_fact    -> the total number of factorizations performed (preconditioner factorizations).
% ------------------------------------------------------------------------------------------------------------
%
% Author: Spyridon Pougkakiotis, March 2024, Dundee, Scotland, UK.
% ____________________________________________________________________________________________________________ 

    % ======================================================================================================== 
    % Initialize problem data and test input correctness
    % -------------------------------------------------------------------------------------------------------- 
    % Set the default values if not provided
    if (~isfield(pb_struct,'A') || ~isfield(pb_struct,'Q') || ~isfield(pb_struct,'C'))
        error("Incorrect input: Missing data. For more information, type ""Help SSN_PMM""");
    else
        A = pb_struct.A;    Q = pb_struct.Q;    C = pb_struct.C;
        [m,n] = size(A);    [l,~] = size(C);
    end
    pb_struct = rmfield(pb_struct,["A","Q","C"]);

    % Make sure that A, C, Q are sparse.
    if (~issparse(A))    A = sparse(A);    end
    if (~issparse(C))    C = sparse(C);    end
    if (~issparse(Q))    Q = sparse(Q);    end


    if (~isfield(pb_struct,'c') || isempty(pb_struct.c))
        c = zeros(n,1);
    else
        c = pb_struct.c;
    end
    if (~isfield(pb_struct,'b') || isempty(pb_struct.b))
        b = zeros(m,1);
    else
        b = pb_struct.b;
    end
    if (~isfield(pb_struct,'d') || isempty(pb_struct.d))
        d = zeros(l,1);
    else
        d = pb_struct.d;
    end
    if (~isfield(pb_struct,'lb') || isempty(pb_struct.lb))
        lb = -Inf.*ones(n,1);
    else
        lb = pb_struct.lb;
    end
    if (~isfield(pb_struct,'ub') || isempty(pb_struct.ub))
        ub = Inf.*ones(n,1);
    else
        ub = pb_struct.ub;
    end
    if (~isfield(pb_struct,'D') || isempty(pb_struct.D))
        D = zeros(n,1);
    else
        D = pb_struct.D;
    end
    if (~isfield(pb_struct,'tol'))
        tol = 1e-4;
    else
        tol = pb_struct.tol;
    end
    if (~isfield(pb_struct,'maxit'))
        maxit = 200;
    else
        maxit = pb_struct.maxit;
    end
    if (~isfield(pb_struct,'pl') || isempty(pb_struct.pl))
        pl = 1;
    else
        pl = pb_struct.pl;
    end
    if (~isfield(pb_struct,'p_fid') || isempty(pb_struct.p_fid))
        p_fid = 1;
    else
        p_fid = pb_struct.p_fid;
    end
    clear pb_struct;   
  
    % Make sure that b, d, c are full.
    if (issparse(b))    b = full(b);       end  
    if (issparse(c))    c = full(c);       end
    if (issparse(d))    d = full(d);       end

    % Make sure that b, c and d are column vectors of dimension m, n and l, repsectively.
    if (size(b,2) > 1)   b = (b)';     end
    if (size(c,2) > 1)   c = (c)';     end
    if (size(d,2) > 1)   d = (d)';     end

    if (~isequal(size(c,1),n) || ~isequal(size(b,1),m) || ~isequal(size(d,1),l))
        error('Incorrect problem dimensions.');
    end
    if (size(D) == size(Q))
        D = spdiags(D,0);
    elseif (size(D) ~= size(c))
        error('Vector D representing the ell-1 weights has wrong dimensions.');
    end
    % ________________________________________________________________________________________________________ 


    % ======================================================================================================== 
    % Transform the ell-1 term \|Dx\|_1 into a max term by adjusting c, d, and C.
    % --------------------------------------------------------------------------------------------------------
    c = c -(ones(n,1).*D);
    C = [C; spdiags(2.*D,0,n,n);];
    d = [d; zeros(n,1)];
    l = l+n;
    % ________________________________________________________________________________________________________
    
    % ======================================================================================================== 
    % Initialization - Starting Point
    % --------------------------------------------------------------------------------------------------------- 
    % Initial primal and dual regularization values.
    beta = 5e1;   rho = 1e2;                             
    x = zeros(n,1);  y1 = zeros(m,1);  y2 = zeros(l,1);  z = zeros(n,1);
    % ________________________________________________________________________________________________________ 

    % ======================================================================================================== 
    % Initialize parameters
    % --------------------------------------------------------------------------------------------------------
    % Iteration counters (and max iterations) for PMM, SSN and counter for Krylov solver iterates
    max_SSN_iters = 4000; 
    SSN_maxit = 40;   
    PMM_iter = 0;                                                       
    SSN_iter = 0;                                                       
    Krylov_iter = 0;
    % Number of factorizations (of preconditioner) performed
    num_fact = 0;         
    % Inner tolerance for Semismooth Newton method.
    in_tol_thr = tol; 
    % Variable monitoring the optimality.
    opt = 0;        
    % Set the printing choice.
    PMM_header(p_fid,pl);   
    % Maximum value for the penalty parameters.
    reg_limit = 1e+6;     
    % Struct to contain output information.
    solution_struct = struct();       
    % Compute the resudual
    [res_p,res_d,compl] = compute_residual(Q,C,A,b,d,c,lb,ub,x,y1,y2,z);
    % ________________________________________________________________________________________________________ 

    while (PMM_iter < maxit)
    % -------------------------------------------------------------------------------------------------------- 
    % SSN-PMM Main Loop structure:
    % Until (primal infeasibility < tol && dual infeasibility < tol && complementarity < tol) do
    %   Call Semismooth Newton method to approximately minimize the augmented Lagrangian w.r.t. x;
    %   update y1,y2,z;
    %   update the reuglarization paramters;
    %   k = k + 1;
    % End
    % -------------------------------------------------------------------------------------------------------- 
        % ==================================================================================================== 
        % Check termination criteria
        % ---------------------------------------------------------------------------------------------------- 
        if (res_p < tol && res_d < tol && compl < tol)
            fprintf('optimal solution found\n');
            opt = 1;
            break;
        end
        PMM_iter = PMM_iter+1;
        % ____________________________________________________________________________________________________ 
                
        % Build or update the Newton structure
        if (PMM_iter == 1) 
            NS = build_Newton_structure(Q,C,A,b,d,c,x,y1,y2,z,beta,rho,lb,ub,PMM_iter,p_fid);
        else
            NS.x = x; NS.y1 = y1; NS.y2 = y2; NS.z = z; 
            NS.beta = beta; NS.rho = rho; NS.PMM_iter = PMM_iter; 
        end
        
        % ==================================================================================================== 
        % Call semismooth Newton method to find the x-update.
        % ---------------------------------------------------------------------------------------------------- 
        res_vec = [0.1*res_p; 0.1*res_d; 1];
        in_tol = max(min(res_vec),in_tol_thr);
        SSN_tol_achieved = 2*max(res_vec);
        counter = 0; 
        while (SSN_tol_achieved > max(1e-1*max(res_vec),min(res_vec))) 
            counter = counter + 1;
            % Call Semismooth Newton method 
            [x, SSN_in_iters, SSN_tol_achieved, Krylov_in_iters, num_of_fact] = SSN(NS,in_tol,SSN_maxit,pl);
            num_fact = num_of_fact + num_fact;
            if (num_of_fact)
                NS.prec_enabled = true;
            end
            SSN_iter = SSN_iter + SSN_in_iters;
            Krylov_iter = Krylov_iter + Krylov_in_iters;
            % Update current x  in the Newton structure in case of a restart
            NS.x = x;
            if (SSN_iter >= max_SSN_iters)
                break;
            end
        end     
        % ____________________________________________________________________________________________________ 
     
        % ==================================================================================================== 
        % Perform the dual update (y1).
        % ---------------------------------------------------------------------------------------------------- 
        y1 = y1 - beta.*(A*x-b);
        % ____________________________________________________________________________________________________ 

        % ==================================================================================================== 
        % Perform the dual update (y2).
        % ---------------------------------------------------------------------------------------------------- 
        tmp_y2 = y2 + beta.*(C*x+d);
        tmp_y2_lb = (tmp_y2 < 0);
        tmp_y2_ub = (tmp_y2 > 1);
        tmp_y2(tmp_y2_lb) = 0;
        tmp_y2(tmp_y2_ub) = 1;
        y2 = tmp_y2;
        % ____________________________________________________________________________________________________ 

        % ==================================================================================================== 
        % Perform the dual update (z).
        % ---------------------------------------------------------------------------------------------------- 
        tmp_vec = (1/beta).*(NS.z) + x;
        tmp_lb = (tmp_vec < lb);
        tmp_ub = (tmp_vec > ub);
        tmp_vec(tmp_lb) = lb(tmp_lb);
        tmp_vec(tmp_ub) = ub(tmp_ub);
        z = NS.z + beta.*x - beta.*tmp_vec;
        % ____________________________________________________________________________________________________ 

        % Compute the new residual norms
        [new_res_p,new_res_d,compl] = compute_residual(Q,C,A,b,d,c,lb,ub,x,y1,y2,z);

        % Check if we need to exit due to excessive SSN iterations
        if (SSN_iter >= max_SSN_iters)
            fprintf('Maximum number of inner iterations is reached. Terminating without optimality.\n');
            break;
        end

        % Update primal-dual PMM penalty parameters (heuristic)
        [beta,rho] = update_PMM_parameters(res_p,res_d,new_res_p,new_res_d,beta,rho,reg_limit);
        res_p = new_res_p;
        res_d = new_res_d;

        % Print iteration output.  
        PMM_output(p_fid,pl,PMM_iter,SSN_iter,res_p,res_d,compl,SSN_tol_achieved,beta,rho);
    end % while (iter < maxit)

    % ========================================================================================================   
    % The PMM has terminated. Print results, and prepare output.
    % -------------------------------------------------------------------------------------------------------- 
    fprintf(p_fid,'outer iterations: %5d\n', PMM_iter);
    fprintf(p_fid,'inner iterations: %5d\n', SSN_iter);
    [res_p,res_d,compl] = compute_residual(Q,C,A,b,d,c,lb,ub,x,y1,y2,z);
    fprintf(p_fid,'primal feasibility: %8.2e\n', res_p);
    fprintf(p_fid,'dual feasibility: %8.2e\n', res_d);
    fprintf(p_fid,'complementarity: %8.2e\n', compl);  
    fprintf(p_fid,'total number of factorizations: %5d\n', num_fact);
    solution_struct.x = x;  solution_struct.y1 = y1; solution_struct.y2 = y2;  solution_struct.z = z;
    solution_struct.opt = opt;  
    solution_struct.PMM_iter = PMM_iter;
    solution_struct.SSN_iter = SSN_iter;   
    solution_struct.Krylov_iter = Krylov_iter;
    solution_struct.num_fact = num_fact;
    solution_struct.obj_val = c'*x + (1/2)*(x'*(Q*x)) + sum(max(C*x+d,zeros(l,1)));
    % ________________________________________________________________________________________________________ 
end
% ************************************************************************************************************ 
% END OF FILE
% ************************************************************************************************************ 