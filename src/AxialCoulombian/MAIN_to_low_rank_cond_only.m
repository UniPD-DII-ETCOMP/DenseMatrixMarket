%% Main_AXIAL_EH
clear
close all
clc
%%
if exist('hm-toolbox-master','dir')
    cd('hm-toolbox-master'); addpath(pwd); cd ..
else
    error('Download hm-toolbox from https://github.com/numpi/hm-toolbox and extract it') 
end
%% BEGIN USER SETTINGS
% problem definition
test_case_dir.con = 'test2_small'; % select source directory for conductive media, set test_case_dir.con = '' to exclude conductive media
test_case_dir.mag = ''; % not supported yet
test_case_dir.ext = ''; % not supported yet
rrho_con=1/2e3; % electric resistivity of conductive media [Ohm m]
% mmu_r=5;% permeability of magnetic  media
f=50; % frequency [Hz]
eps_0 = 8.85418781762e-12; % vacuum permittivity
mu_0 = 4*pi*1e-7; % vacuum permeability 
plotflag = 1; %graphics postprocessing flag (1 = yes, 0 = no)
ext_Field_ext=1; %0=external field off ||| 1=extarnal field on
E_ext=@(r,phi,z) [0, -mu_0*1j*2*pi*f*r*0.5, 0]; % external electric fiedl [V/m]
H_ext=@(r,phi,z) [0 ,0, mu_0/mu_0]; % external magnetic field [A/m]
J_ext=@(r,phi,z) [0, 1, 0]; % coil J_ext current in the external coil [A/m^2], if any
N_thread =22; % for mexed fortran function
% END USER SETTINGS
%% Pre-processing of data
disp('-------------------------------------------------------------------')
disp('pre-processing...')
cd('fun'); addpath(pwd); cd ..
if ~strcmp(test_case_dir.mag,''); error('magnetic media not supported yet'); end
if ~strcmp(test_case_dir.ext,''); error('external coil  not supported yet');end
[N,Matrix_P0_c,Matrix_P0_m,Matrix_P0_e,w,rho_c,Area_c,rho_m,Fac_Ed_loc_m,Aed_m,...
      Led_m,Area_m,ind_edge_free_mag,Area_e,ind_edge_free_con,bar_e,Matrix_P0,...
      Matrix_Cb_m,phi_m_n0,ex,F1_c,F1_m,C1_m,F1_e,G1_m,Matrix_G_m,Matrix_C_m]=...
      fun_pre_processing(f,rrho_con,NaN,mu_0,test_case_dir);
disp('...pre-processing done!')
disp('-------------------------------------------------------------------')
%% Rhs-term External fields 
disp('-------------------------------------------------------------------')
disp('rhs...')
[brhs] = fun_compute_ext_field_phi(Area_c,Matrix_P0_c,N.face_con,F1_c,E_ext,3);
disp('...done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp('compute R matrix and handle-L...')
[Rcdiag] = funRphi_diag_with_check(ex,N,rho_c,Area_c,Matrix_P0_c,F1_c);
try 
 Lfun = @(ii,jj)mu_0*fun_L_interface(N.node_con,N.face_con,Area_c,...
                Matrix_P0_c,F1_c,N_thread,ii,jj);
 check=Lfun(1,1); mtlb=0;
catch ME
    warning('FORTAN MEX-FILE is not supported, try to re-mex it in /MEXfortran and replace it in /fun')
    warnin(['MESSAGE ERROR:' ME.message]);
    warning('Slow Matlab function is used instead');
    mtlb=1;
     [PPP1,PPP2,NPP,NP1,NP2,WW] = funLphiphi3_acc_PrePro(N.face_con,1,2,Matrix_P0_c,F1_c);
     Lfun = @(ii,jj) mu_0*funLphiphi3_acc_low_rank(PPP1,PPP2,NPP,NP1,NP2,WW,ii,jj);   
end
if mtlb==0; disp('...FORTAN MEX-FUNCTION supported...'); end
disp('... done!')
disp('-------------------------------------------------------------------')
%% compress HSS/HODLR
hodlroption('block-size',100);
hodlroption('compression','qr');
hodlroption('threshold',1e-4);
hssoption('block-size',100);
hssoption('compression','qr');
hssoption('threshold',1e-4);
disp('R matrix to hodlr/hss...')
tic
HR=hodlr('diagonal',Rcdiag);
toc
disp('L matrix to hodlr/hss...')
tic
HL=hodlr('handle',Lfun,N.face_con,N.face_con);
toc
compr=100*getSize(HL)/(N.face_con*getSize(Lfun(1,1:N.face_con)));
disp(['Compression ratio perc. =',num2str(compr),'%'])
%% SOLVING
disp('-------------------------------------------------------------------')
disp('Solving ...')
tic
x=(HR+1j*HL)\brhs;
toc
disp('... done!')
disp('-------------------------------------------------------------------')
%%
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

