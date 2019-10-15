%% Rhs-term External fields volu
tic
if ext_Field_ext==1
    %%
    if ex.vol_con
        N_GAUSS=3;
        [e_dual_c] = fun_compute_ext_field_phi(Area_c,Matrix_P0_c,N.face_con,F1_c,E_ext,N_GAUSS);
    else
        e_dual_c=[];
    end
    %%
    if ex.vol_mag
         N_GAUSS=4;
         [h_dual_m] = fun_compute_ext_field_rz_face2(N.face_mag,N_GAUSS,Matrix_P0_m,F1_m,N.edge_mag,C1_m,H_ext,Fac_Ed_loc_m,Aed_m,Led_m);
    else
        h_dual_m=[];
    end
    %%
    b(1:N.face_con+N.edge_mag,1)=[e_dual_c;h_dual_m];
    brhs=brhs+b;
end
toc
%%