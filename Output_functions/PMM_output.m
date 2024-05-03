function PMM_output(fid,pl,o_it,i_it,xinf,sinf,compl,SSN_tol_achieved,beta,rho)
% ============================================================================= 
% This function outputs the residual infeasibilities
% and other statistics of SNN-PMM, to a given FID.
% ----------------------------------------------------------------------------- 
    if (pl >= 1)
        fprintf(fid,' ');
        fprintf(fid,'%5d    ', o_it);
        fprintf(fid,'%8d    ', i_it);
        fprintf(fid,'%15e   ', SSN_tol_achieved);
        fprintf(fid,'%11.2e  ', xinf);
        fprintf(fid,'%8.2e  ', sinf);
        fprintf(fid,'%8.2e  ', compl);
        fprintf(fid,'%8.2e  ', beta);
        fprintf(fid,'%8.2e  ', rho);            
    end
    if (pl >= 1) fprintf('\n'); end
% _____________________________________________________________________________
end
