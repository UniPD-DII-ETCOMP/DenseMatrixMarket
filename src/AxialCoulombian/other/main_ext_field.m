%% main_ext_field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E0=1;
%% coil J_ext current in the external coil, if any
J_ext=@(r,phi,z) [0, 1, 0];
%% external fields
E_ext=@(r,phi,z) [0, -mu_0*1j*w(1)*r*0.5, 0];
H_ext=@(r,phi,z) [0 ,0, mu_0/mu_0];
B_ext=@(r,phi,z) [0 ,0, mu_0];
%%
