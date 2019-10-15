N_GAUSS=4;
N_GAUSS1=1; N_GAUSS2=2;
%% Matrix_Rc
disp('... Rc ...')
if ex.vol_con
    tic
    [Rcdiag] = funRphi(N.face_con,rho_c,Area_c,N_GAUSS,Matrix_P0_c,F1_c);
    Matrix_Rc=sparse(diag(Rcdiag));
    toc
else
    Matrix_Rc=sparse(N.face_con,N.face_con);
end
disp('... done!')
%% Matrix_Rm
disp('... Rm ...')
if ex.vol_mag
    tic
    [Matrix_Rm] = funRrz_face2_acc(N.face_mag,(rho_m).',N_GAUSS,Matrix_P0_m,F1_m,N.edge_mag,C1_m,Fac_Ed_loc_m,Aed_m.',Led_m.',Area_m);
    Matrix_Rm=sparse(Matrix_Rm);
    toc  
else
    Matrix_Rm=sparse(N.edge_mag,N.edge_mag);  
end
disp('... done!')
%% Matrix_Rm Lcc
disp('... Lcc...')
if ex.vol_con 
tic
    [Matrix_Lcet] = mu_0*funLphiphi3_acc(N.face_con+N.face_ext,...
                                        N_GAUSS1,...
                                        N_GAUSS2,...
                                       [Matrix_P0_c;Matrix_P0_e],...
                                       [F1_c,N.node_con+F1_e]);
    indc=1:N.face_con;
    inde=N.face_con+1:N.face_con+N.face_ext;
    Matrix_Lcc=    Matrix_Lcet(indc,indc);
    Matrix_Lce=    Matrix_Lcet(indc,inde);
    Matrix_Lee=    Matrix_Lcet(inde,inde);
    clear Matrix_Lcet
toc
else
    Matrix_Lcc=    sparse(N.face_con,N.face_con);
    Matrix_Lce=    sparse(N.face_con,N.face_ext);
    Matrix_Lee=    sparse(N.face_ext,N.face_ext);
end
disp('... done!')
%% Matrix_Pm_ss
disp('... Pm_ss ...')
if ex.vol_mag
    tic
    [Matrix_Pm_ss] = funPphiphi3_ss_acc(N.edge_free_mag,...
                                    N_GAUSS1,...
                                    N_GAUSS2,...                               
                                    Matrix_P0_m,...
                                    G1_m(:,ind_edge_free_mag),...
                                    Aed_m(ind_edge_free_mag).');
    Matrix_Pm_ss=(1/(mu_0))*Matrix_Pm_ss;
    toc    
else
    Matrix_Pm_ss=sparse(N.edge_free_mag,N.edge_free_mag);  
end
disp('... done!')
%% Matrix_Nme 
disp('... Nmc Nme ...')
if ex.vol_mag && (ex.vol_con||ex.vol_ext)
    tic
            [Matrix_Projm] = funRrz_face2_acc(N.face_mag,ones(N.face_mag,1),N_GAUSS,Matrix_P0_m,F1_m,N.edge_mag,C1_m,Fac_Ed_loc_m,Aed_m.',Led_m.',Area_m);
            Matrix_Projm=sparse(Matrix_Projm);
    toc
    tic
            [Matrix_Nme_ce] = funNme_acc(N.face_con+N.face_ext,...
                                         N_GAUSS,...
                                         N.node_mag,...
                                         [Area_c;Area_e].',...
                                         [Matrix_P0_c;Matrix_P0_e],...
                                         [F1_c,N.node_con+F1_e],...
                                         Matrix_P0_m);             
    toc             
    indc=1:N.face_con;
    inde=N.face_con+1:N.face_con+N.face_ext;                  
    Matrix_Nmc=Matrix_Nme_ce(:,indc);
    Matrix_Nme=Matrix_Nme_ce(:,inde);
    clear Matrix_Nme_ce
else
    Matrix_Projm=sparse(N.edge_mag,N.edge_mag);
    Matrix_Nmc=sparse(N.node_mag,N.face_con);
    Matrix_Nme=sparse(N.node_mag,N.face_ext);
end
%%
return

                         