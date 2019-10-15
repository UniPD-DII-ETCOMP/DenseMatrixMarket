% 
if ex.vol_mag
[sol.M_r,sol.M_i,sol.M_norm_r,sol.M_norm_i,sol.M_bar] = fun_post_J_Quad_face2(NaN,N.face_mag,Matrix_P0_m,F1_m,C1_m,sol.m,Fac_Ed_loc_m,Aed_m,Led_m);
sol.M_r=-sol.M_r;
sol.M_i=-sol.M_i;
end
