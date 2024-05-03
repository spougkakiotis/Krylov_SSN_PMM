function w = Precond_Operator(x,PS)
% ===================================================================================================================== %
% Preconditioner Operator:
% --------------------------------------------------------------------------------------------------------------------- %
% [w] = Precond_Operator(x,PS) takes as an input the struct containing the preconditioner blocks, as well as 
% a vector of size n or m, and returns the matrix-vector product of the inverse preconditioner by this vector.
% _____________________________________________________________________________________________________________________ %
    warn_stat = warning;
    warning('off','all');
    u_1 = x(1:PS.n,1);
    u_1 = (1./PS.H_tilde).*u_1;
    u_2 = x(PS.n+1:end,1);
    u_2 = u_2(PS.Perm,1);
    u_2 = PS.L_S'\(PS.L_S\u_2);
    u_2 = u_2(PS.PermInv,1);
    w = [u_1; u_2];
    warning(warn_stat);
end

