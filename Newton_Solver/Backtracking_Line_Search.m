function alpha = Backtracking_Line_Search(NS,x,dx,mu,delta)
% =================================================================================================
% Backtracking Line Search
% -------------------------------------------------------------------------------------------------
% alpha = Backtracking_Line_Search(NS,start_point,direction,mu,delta)
%
% This function performs backtracking linesearch based on the proximal augmented Lagrangain penalty
%                                                           
% Author: Spyridon Pougkakiotis.
% _________________________________________________________________________________________________

    % =============================================================================================
    % Evaluate the proximal AL penalty and its gradient
    % ---------------------------------------------------------------------------------------------
    [phi,nabla_phi] = Evaluate_Prox_Lagrangian(NS,x,true);
    % _____________________________________________________________________________________________
    
    % =============================================================================================
    % Let alpha = delta, and evaluate phi(x + alpha.*dx).
    % ---------------------------------------------------------------------------------------------
    alpha = delta;
    x_new = x + alpha.*dx;
    [phi_new,~] = Evaluate_Prox_Lagrangian(NS,x_new,false);
    % _____________________________________________________________________________________________
    
    % =============================================================================================
    % Iterate until you find a step-length satisfying the Armijo-Goldstein condition.
    % ---------------------------------------------------------------------------------------------
    counter = 1;
    while (phi_new > phi + mu*(alpha)*(nabla_phi'*dx))
        counter = counter + 200;
        alpha = delta^counter;
        x_new = x + alpha.*dx;
        [phi_new,~] = Evaluate_Prox_Lagrangian(NS,x_new,false);
        if (alpha < 1e-3)
            break;
        end
    end
    % _____________________________________________________________________________________________
end
% *************************************************************************************************
% END OF FILE.
% *************************************************************************************************

