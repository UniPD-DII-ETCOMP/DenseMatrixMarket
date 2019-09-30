clc
clear
close all
%%
if exist('hm-toolbox-master','dir')
    cd('hm-toolbox-master'); addpath(pwd); cd ..
    texthss=textread('hss.m', '%s','delimiter', '\n'); %#ok<DTXTRD>
    texthss=char(texthss(2,:));
    if strcmp(texthss,'% mod_by_RT')
        disp('modified hss.m file detected')
    else
        error('copy lowrank/hss.m file to hm-toolbox-master/@hss ')
    end
else
    error('Download hm-toolbox from https://github.com/numpi/hm-toolbox and extract it') 
end
%% BEGIN USER SETTINGS
% problem definition
test_case_dir = 'lineHSS';
plotflag = 1; %graphics postprocessing flag (1 = yes, 0 = no)
freq = 3*3e4; % set the value of the frequency [Hz]
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
cd('fun'); addpath(pwd); cd ..
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
disp('------------------------------------')
hssoption('block-size',floor((nSticks+N.vol_sou+N.node_app)/10));
hssoption('compression','svd');
hssoption('threshold',1e-6);
hodlroption('block-size',floor((nSticks+N.vol_sou+N.node_app)/10));
hodlroption('compression','svd'); %hodlr with 'handle' uses ACA
hodlroption('threshold',1e-6);
%
disp('build HSS/HODLR matrix')
tic
H=hss('function',sys_MF,nSticks+N.vol_sou+N.node_app,nSticks+N.vol_sou+N.node_app);
% H=hodlr('handle',sys_MF,nSticks+N.vol_sou+N.node_app,nSticks+N.vol_sou+N.node_app);
toc
%
compr=100*getSize(H)/((nSticks+N.vol_sou+N.node_app)^2*16);
disp(['Compression ratio HSS/HODLR =',num2str(compr),'%'])
%
figure
spy(H)
axis equal
title('rank of off-diagonal blocks')
drawnow
%
disp('Solution HSS/HODLR')
tic
x_fromHSS = H\b; 
toc
disp('------------------------------------')
%%
disp('------------------------------------')
disp('Build full system M from matrix-free function "sys_MF"...')
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
%%
err=100*norm(x-x_fromHSS)/norm(x);
disp(['Solution perc. error =',num2str(err),'%'])
disp(['Relative residual =',num2str(100*norm(H*x-b)/norm(b)),'%']);
%%
gtoc = toc(gtic); %global TOC
%%
 

                           
                           
