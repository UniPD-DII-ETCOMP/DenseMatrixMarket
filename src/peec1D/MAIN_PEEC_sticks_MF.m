clc
clear
close all
%% BEGIN USER SETTINGS
% problem definition
test_case_dir = 'pier';
plotflag = 1; %graphics postprocessing flag (1 = yes, 0 = no)
freq = 3*3e0; % set the value of the frequency [Hz]
radius = 0.1e-3; % set the radius of the wires  [m]
rho = 1/57e6;% set the value of the resistivity [Ohm m]
% algorithmic settings
ret = 1; %consider retarded potentials (1 = yes, 0 = no)
Wint = 1; %consider internal energy of wires for the self inductance computation (1 = yes, 0 = no)
np_gauss = 2; % number of gauss points for the numerical integration
E_external = 0; % external electric field (1 = yes, 0 = no)
% External field (only active if E_external==1)
E_ext=@(x,y,z) [(y-0.5),-(x-0.5),0]; % set the external electric field distribution
% END USER SETTINGS
%% data input (load data.mat examples from "test_case_dir" directory in "data_generation" directory)
cd test_cases
cd(test_case_dir)
load data.mat
type description.txt;
disp(' ')
cd ../..
%% Time
gtic = tic; %global tic
%% Pre-processing of data
addpath([pwd,'\fun'])
[N,nSticks,PPgh,ll_h,ut_h,gauss_W,ne_cap_max,...
    Cap_Elem,NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,R,omega,Aobj,...
    Geapp,A_eapp,R_app,L_app,Cinv_app,A_napp,...
    beta,radiusvec,b,Is,nNodes,lgt,radius_cap,bar_stick,barDoFs]...
    =fun_PEEC_1D_preprocessing(Voltage_Source,Current_Source,...
    NN,G,radius,freq,np_gauss,rho,E_external,ret,E_ext);
%% Matrix free system
disp('------------------------------------')
disp('Creating funtion handled (matrix free system)...') % sys_MF is a function handled which generates a generic matrix block of the system of equation
sys_MF=@(ii,jj) fun_system_MF(NN,G,radiusvec,np_gauss,Wint,PPgh,ll_h,ut_h,gauss_W,...
                               ne_cap_max,Cap_Elem,...
                               NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,...
                               nSticks,N.vol_sou,N.node_app,R,omega,Aobj,Geapp,A_eapp,R_app,L_app,Cinv_app,A_napp,...
                               ret,beta,...
                               ii,... % ii=[valMin_ii:valMax_ii]
                               jj);   % jj=[valMin_jj:valMax_jj]
disp('...done!') % N.B.: system M is not forced to be symmetric 
disp('------------------------------------')
%%
disp('Build system M from matrix-free function "sys_MF"...')
tic 
   M=sys_MF(1:nSticks+N.vol_sou+N.node_app,1:nSticks+N.vol_sou+N.node_app);                        
toc
disp('...done!')
disp('------------------------------------')
%% solve
disp('------------------------------------')
disp('solving system...')
tic
x = M\b;
toc
disp('...done!')
disp('------------------------------------')
%% post-processing
if plotflag
[I_obj,I_app,Q,U_obj,rho_charge,Jrenorm,Jimnorm,Jre,Jim] = fun_PEEC_1D_postprocessing(nSticks,x,N,Aobj,A_eapp,Is,...
                   omega,ret,NN,G,radiusvec,np_gauss,ne_cap_max,Cap_Elem,...
                   NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,gauss_W,...
                   nNodes,beta,lgt,radius_cap);
%% plot
[~] = fun_PEEC_1D_plot(U_obj,NN,rho_charge,I_obj,bar_stick,Jrenorm,Jimnorm,G,Jre,Jim);
end
%%
gtoc = toc(gtic) %global TOC
%%
 

                           
                           