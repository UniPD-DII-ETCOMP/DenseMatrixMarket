%% Main_AXIAL_EH
clear
close all
clc
%% BEGIN USER SETTINGS
% problem definition
test_case_dir.con = 'test1'; % select source directory for conductive media, set test_case_dir.con = '' to exclude conductive media
test_case_dir.mag = 'test1'; % select source directory for magnetic media, set test_case_dir.mag = '' to exclude magnetic media
test_case_dir.ext = 'test1'; % select source directory for exteranl coil, set test_case_dir.ext = '' to exclude exteranl coil
rrho_con=1/2e3; % electric resistivity of conductive media [Ohm/m^2]
mmu_r=5;% permeability of magnetic  media
f=50; % frequency [Hz]
eps_0 = 8.85418781762e-12; % vacuum permittivity
mu_0 = 4*pi*1e-7; % vacuum permeability 
plotflag = 1; %graphics postprocessing flag (1 = yes, 0 = no)
ext_Field_ext=1; %0=external field off ||| 1=extarnal field on
E_ext=@(r,phi,z) [0, -mu_0*1j*2*pi*f*r*0.5, 0]; % external electric fiedl [V/m]
H_ext=@(r,phi,z) [0 ,0, mu_0/mu_0]; % external magnetic field [A/m]
J_ext=@(r,phi,z) [0, 1, 0]; % coil J_ext current in the external coil [A/m^2], if any
% END USER SETTINGS
%% Pre-processing of data
disp('-------------------------------------------------------------------')
disp('pre-processing...')
cd('fun'); addpath(pwd); cd ..
[N,Matrix_P0_c,Matrix_P0_m,Matrix_P0_e,w,rho_c,Area_c,rho_m,Fac_Ed_loc_m,Aed_m,...
      Led_m,Area_m,ind_edge_free_mag,Area_e,ind_edge_free_con,bar_e,Matrix_P0,...
      Matrix_Cb_m,phi_m_n0,ex,F1_c,F1_m,C1_m,F1_e,G1_m,Matrix_G_m,Matrix_C_m]=...
      fun_pre_processing(f,rrho_con,mmu_r,mu_0,test_case_dir);
disp('...pre-processing done!')
disp('-------------------------------------------------------------------')
%% Computing matrices RLP
disp('-------------------------------------------------------------------')
disp('computing matrices RLP...')
[Matrix_Rc,Matrix_Rm,Matrix_Lce,Matrix_Lee,Matrix_Pm_ss,Matrix_Lcc,...
    Matrix_Projm,Matrix_Nmc,Matrix_Nme,brhs] ...
    = fun_compute_matrices(N,ex,rho_c,Area_c,Matrix_P0_c,F1_c,rho_m,...
    Matrix_P0_m,F1_m,C1_m,Fac_Ed_loc_m,Aed_m,Led_m,Area_m,mu_0,...
    Matrix_P0_e,F1_e,G1_m,ind_edge_free_mag,Area_e,ext_Field_ext,E_ext,...
    H_ext,bar_e,J_ext,w,Matrix_G_m);
disp('...RLP done!')
disp('-------------------------------------------------------------------')
%% Assebling System
disp('-------------------------------------------------------------------')
disp('Assembling system...')
SYSTEM=[Matrix_Rc+1j*w(1)*Matrix_Lcc,...
        -(Matrix_Projm*(Matrix_G_m*Matrix_Nmc)).'*Matrix_G_m(:,phi_m_n0);... 
        Matrix_G_m(:,phi_m_n0).'*Matrix_Projm*(Matrix_G_m*Matrix_Nmc),...
        Matrix_G_m(:,phi_m_n0).'*(Matrix_Rm+(1/(1j*w))*(Matrix_Cb_m.'*Matrix_Pm_ss*Matrix_Cb_m))*Matrix_G_m(:,phi_m_n0)];
brhs=blkdiag(speye(N.face_con),Matrix_G_m(:,phi_m_n0)).'*brhs;
disp('... done!')
disp('-------------------------------------------------------------------')
%% SOLVING
disp('-------------------------------------------------------------------')
disp('Solving ...')
tic
yq=SYSTEM\brhs;
x=blkdiag(speye(N.face_con),Matrix_G_m(:,phi_m_n0))*yq;
toc
disp('... done!')
disp('-------------------------------------------------------------------')
%% PLOT Results
if plotflag
disp('-------------------------------------------------------------------')
disp('post-processing...')
cd('fun'); addpath(pwd); cd ..
[sol]=fun_post_processing(x,N,Area_c,w,mu_0,Matrix_C_m,...
    Matrix_Cb_m,ex,Matrix_P0_m,F1_m,C1_m,Fac_Ed_loc_m,Aed_m,...
    Led_m,F1_c,Matrix_P0_c,F1_e,Matrix_P0_e,Matrix_P0);
disp('...post-processing done!')
disp('-------------------------------------------------------------------')
end
%%
return
