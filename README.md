# An active-set method for Convex QP with piecewise linear terms

This is a MATLAB implementation of an active-set method suitable for convex quadratic programming problems with piecewise linear terms of the following form:

$$ \min_{x \in \mathbb{R}^n}\  c^\top x + \frac{1}{2} x^\top Q x + \sum_{i=1}^l \max((Cx+d)_i,0),\qquad \text{s.t., }Ax = b,\ x \in [a_l,a_u],$$

where $Q \in \mathbb{R}^{n\times n}$ is a positive semidefinite matrix, and $C \in \mathbb{R}^{l\times n}$, $d \in \mathbb{R}^l$.



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
