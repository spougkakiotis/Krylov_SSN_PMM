function NS = build_Newton_structure(Q,C,A,b,d,c,x,y1,y2,z,beta,rho,lb,ub,PMM_iter,fid)
% ==================================================================================================================== %
% build_Newton_structure: Store all relevant information about the Newton system.
% -------------------------------------------------------------------------------------------------------------------- %
% NS = build_Newton_structure(Q,D,C,A,b,c_1,c_2,x,w,v,z,beta,rho,zeta,lb,ub,PMM_iter) returns a MATLAB 
%      struct that holds the relevant information of the Newton system, required for solving the step equations in
%      the SSN-PMM.
% 
% Author: Spyridon Pougkakiotis.
% ____________________________________________________________________________________________________________________ %
    % ================================================================================================================ %
    % Store all the relevant information required from the semismooth Newton's method.
    % ---------------------------------------------------------------------------------------------------------------- %
    NS = struct();
    NS.x = x;
    NS.y1 = y1;
    NS.y2 = y2;
    NS.z = z;
    NS.m = size(b,1);
    NS.n = size(c,1);
    % We only treat the nonseparable terms
    NS.l = size(d,1) - NS.n;
    NS.b = b;
    NS.d = d(1:NS.l,1);
    NS.c = c;
    NS.beta = beta;
    NS.rho = rho;
    NS.PMM_iter = PMM_iter;
    NS.lb = lb;
    NS.ub = ub;
    NS.A = A;
    NS.A_tr = A';
    % Separate the nonseparable terms from the separable ones.
    NS.C = C(1:NS.l,:);
    NS.D = spdiags(C(NS.l+1:end,:),0);
    NS.C_tr = NS.C';
    NS.Q = Q;
    NS.Q_diag = spdiags(Q,0);
    NS.fid = fid;
    NS.prec_enabled = false;
    % ________________________________________________________________________________________________________________ %
end
% ******************************************************************************************************************** %
% END OF FILE.
% ******************************************************************************************************************** %
