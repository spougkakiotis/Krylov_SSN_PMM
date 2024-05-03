function [w] = AS_multiplier(x,NS)
% ===================================================================================================================== %
% Augmented System Operator:
% --------------------------------------------------------------------------------------------------------------------- %
% [w] = AS_multiplier(NS,x) takes as an input the struct containing the Newton blocks
% and returns the matrix-vector product of the augmented system matrix by this vector.
% _____________________________________________________________________________________________________________________ %
    n = size(NS.x,1);
    m = NS.Schur_size;
    w = zeros(n+m,1);
    x_1 = x(1:n,1); x_2 = x(n+1:end,1);
    w(1:n,1) = - (NS.Q*x_1 + (NS.beta+ (1/NS.rho)).*x_1 ...
                  - NS.beta.*(NS.B_delta.*x_1) + NS.beta.*((((NS.D).^2).*NS.B_h_2).*x_1))  + NS.G'*x_2;
    w(n+1:end,1) = NS.G*x_1 + (1/NS.beta).*x_2;
end

