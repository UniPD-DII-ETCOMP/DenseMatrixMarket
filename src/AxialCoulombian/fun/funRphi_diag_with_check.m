function [Rcdiag] = funRphi_diag_with_check(ex,N,rho_c,Area_c,Matrix_P0_c,F1_c)
N_GAUSS=4;
if ex.vol_con 
    tic
    [Rcdiag] = funRphi(N.face_con,rho_c,Area_c,N_GAUSS,Matrix_P0_c,F1_c);
    toc
else
    Rcdiag=zeros(N.face_con,1);
end
end

