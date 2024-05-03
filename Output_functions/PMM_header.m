function PMM_header(fid,pl)
% =================================================================================================================== 
% header function: (for output on a given file, indicated by fid)
% pl >= 1: primal-dual infeasibility, complementarity as well as number of PMM and SSN iteration is printed.
% ------------------------------------------------------------------------------------------------------------------- 
    if (pl >= 1)
        fprintf(fid,' ');
        fprintf(fid,'%5s    ', 'PMM iter');
        fprintf(fid,'%4s   ', 'SSN iter');
        fprintf(fid,'%1s   ','SSN tol achieved');
        fprintf(fid,'%5s   ',   'pr feas');
        fprintf(fid,'%4s  ',   'dl feas');
        fprintf(fid,'%7s',   'compl.');
        fprintf(fid,'%8s ',   'beta');
        fprintf(fid,'%8s  ',   'rho');
    end
    if (pl >= 1)
        fprintf(fid,'\n ========    ========   ================   ========  ========  ========  ========  ========');
    end
    if (pl >= 1) fprintf(fid,'\n'); end
% ___________________________________________________________________________________________________________________ 
end
