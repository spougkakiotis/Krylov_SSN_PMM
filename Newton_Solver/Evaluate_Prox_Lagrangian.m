function [phi,nabla_phi] = Evaluate_Prox_Lagrangian(NS,x,grad)
% ========================================================================================= 
% [phi,nabla_phi] = Evaluate_Prox_Lagrangian(NS,x)
% ----------------------------------------------------------------------------------------- 
% This function takes as an input the Newton structure 
% containing information necessary to build the proximal 
% augmented Lagrangian, which it evaluates at the current point x.
%
% The third argument is a logical variable. If true, returns the gradient as well
% otherwise, it returns Inf for the gradient.
% _________________________________________________________________________________________ 
    if (nargin < 3 || isempty(grad))
        grad = false;
    end
    % ===================================================================================== 
    % Compute the projection onto K.
    % ------------------------------------------------------------------------------------- 
    w = (1/NS.beta).*(NS.z) + x;  
    tmp_lb = (w < NS.lb);
    tmp_ub = (w > NS.ub);
    w(tmp_lb) = NS.lb(tmp_lb);
    w(tmp_ub) = NS.ub(tmp_ub);
    % _____________________________________________________________________________________ 

    % ===================================================================================== 
    % Evaluate Proj_{C}(y2(1:l) + beta(Cx+d)) (Nonseparable max{,0} terms)
    % -------------------------------------------------------------------------------------
    max_res = NS.C*x + NS.d;
    u = NS.y2(1:NS.l,1) + NS.beta.*(max_res);
    u_lb = (u < 0);
    u_ub = (u > 1);
    u(u_lb) = 0;
    u(u_ub) = 1;
    % _____________________________________________________________________________________
   
    % ===================================================================================== 
    % Evaluate Proj_{C}(y2(l+1:l+n) + beta(Dx)) (Separable max{,0} terms)
    % ------------------------------------------------------------------------------------- 
    v = NS.y2(NS.l+1:end,1) + NS.beta.*(NS.D.*x);
    v_lb = (v < 0);
    v_ub = (v > 1);
    v(v_lb) = 0;
    v(v_ub) = 1;
    % _____________________________________________________________________________________

    % ===================================================================================== 
    % Evaluate the proximal augmented Lagrangian at x
    % ------------------------------------------------------------------------------------- 
    pr_res = NS.A*x - NS.b;
    phi = NS.c'*x + (1/2).*(x'*(NS.Q*x)) - NS.y1'*(pr_res) ...
          + (NS.beta/2)*(norm(pr_res)^2) + (max_res)'*u + (NS.D.*x)'*v ...
          - (1/(2*NS.beta))*(norm(NS.y2-[u;v])^2) - (1/(2*NS.beta))*norm(NS.z)^2 ...
          + (1/(2*NS.beta))*(norm(NS.z + NS.beta.*x - NS.beta.*w)^2)...
          + (1/(2*NS.rho))*(norm(x-NS.x))^2;
    % _____________________________________________________________________________________

    % ===================================================================================== 
    % Evaluate the gradient of the proximal augmented Lagrangian at x (if asked to)
    % ------------------------------------------------------------------------------------- 
    if (grad)
        nabla_phi = NS.c + NS.Q*x - NS.A'*NS.y1 + NS.beta.*(NS.A'*(pr_res))...
                    + NS.C'*u + NS.D.*v + NS.z + NS.beta.*x - NS.beta.*w ...
                    + (1/NS.rho).*(x-NS.x);
    else
        nabla_phi = Inf;
    end
    % _____________________________________________________________________________________
end

