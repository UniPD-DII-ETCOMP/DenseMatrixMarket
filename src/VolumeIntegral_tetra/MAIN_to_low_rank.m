clear
close all
clc
%% BEGIN USER SETTINGS 
% N.B. non-simply connected domains are not supported, cavities are supported 
test_case_dir = 'test2';
f = 50;% Hz
rrho_c=1/56e3; % Ohm m
N.thread = 22; 
reord=1; %reordering of unknowns flag (1 = yes, 0 = no)
E_ext=@(x,y,z) [-1j*y, 1j*x, 0];
plotflag = 1; %graphics postprocessing flag (1 = yes, 0 = no)
thr=0.3; %threshold for close to far field
% END USER SETTINGS
%% 
if exist('hm-toolbox-master','dir')
    cd('hm-toolbox-master'); addpath(pwd); cd ..
else
    error('Download hm-toolbox from https://github.com/numpi/hm-toolbox and extract it') 
end
%% load data
cd test_cases; cd(test_case_dir); load data.mat; cd ..; cd ..
%% pre_procesing
cd('fun'); addpath(pwd); cd ..
[rho_c,N,ind,EV,VE,EV_ind,VE_ind,Q,brhs,norm_eps,w,bared,...
    xmin,xmax,ymin,ymax,zmin,zmax] = ...
    fun_pre_processing(F1,D1,C1,...
    G1,VP,Matrix_D,Matrix_C,Matrix_G,Matrix_P0,E_ext,rrho_c,N,f);
%% reordering
if reord==1
    dims = 2; % 1=1D geometry, 2=2D geometry, 3=3D geometry
    base=12; lev=4; % base=how many sons? lev=how many levels?
    [map,invmap] = fun_reo_DMM(bared,N,ind,dims,base,lev,...
        xmin,xmax,ymin,ymax,zmin,zmax);
else
    map=1:N.cotree;  invmap=1:N.cotree;
end
%% R
disp('-------------------------------------------------------------------')
disp('MAKING R ...')
try 
    tic
    [Ree,ne] =fun_const_rot_edge_for90_sp2sp(N.edge,Matrix_P0,N.node,N.volu,...
        N.thread,VE(:,:),VE_ind(:,:),VP(:,:),rho_c(:),N.n_EV_max,N.nnzR);
    toc
    R=sparse(Ree(1,1:ne),Ree(2,1:ne),Ree(3,1:ne),N.edge,N.edge); clear Ree;
    R=R(ind.cotree,ind.cotree);
    mtlb=0;
catch ME
    warning('FORTAN MEX-FILE is not supported, try to re-mex it in /MEXfortran/R and replace it in /fun')
    warning(['MESSAGE ERROR:' ME.message]);
    warning('Slow Matlab function is used instead');    
    mtlb=1;
    tic
    [R] = funRmatlab(N,VP,Matrix_P0,VE,VE_ind,rho_c,ind);
    toc
end
if mtlb==0; disp('...FORTAN MEX-FUNCTION supported...'); end
disp('... done!')
disp('-------------------------------------------------------------------')
%% L
disp('-------------------------------------------------------------------')
disp('MAKING RL...')
global nncoeff; 
nncoeff=0; %#ok<NASGU>
try 
    [RLfun] = @(ii,jj) full(R(map(ii),map(jj)))+1j*w*...
        fun_L_interface1(N,EV,Matrix_P0,VP,EV_ind,...
        norm_eps,ind.cotree(map(ii)),ind.cotree(map(jj)),bared,thr);    
    check=RLfun(1,1); mtlb=0;
catch ME
    warning('FORTAN MEX-FILE is not supported, try to re-mex it in /MEXfortran/L and replace it in /fun')
    warning(['MESSAGE ERROR:' ME.message]);
    warning('Slow Matlab function is used instead');    
    mtlb=1;
    tic 
    [we,X8e,Y8e,Z8e,W8e,X11e,Y11e,Z11e,W11e] = ...
        fun_L_curledge_curledge_PrePro(N.edge,...
        N.n_EV_max,EV,Matrix_P0.',VP,EV_ind,N.n_EV);
    toc
    [RLfun] = @(ii,jj) full(R(map(ii),map(jj)))+1j*w*...
        fun_L_curledge_curledge_pp(EV,norm_eps,N.n_EV,...
        X8e,Y8e,Z8e,...
        W8e,X11e,Y11e,Z11e,W11e,we,length(ii),length(jj),...
        ind.cotree(map(ii)),ind.cotree(map(jj)));
end
disp('... done!')
disp('-------------------------------------------------------------------')
%%  HSS/HODLR SETTINGS 
hodlroption('block-size',min(floor(N.cotree/20),100));
hodlroption('compression','qr');
hodlroption('threshold',1e-6);
hssoption('block-size',min(floor(N.cotree/20),100));
hssoption('compression','qr');
hssoption('threshold',1e-6);
%% SYSTEM in HODLR format
disp('-------------------------------------------------------------------')
disp('SYS=R+1j*w*L in HODLR format...')
nncoeff=0;
tic
SYS=hodlr('handle',RLfun,N.cotree,N.cotree); % (in this case ACA is used)
toc
compr=100*getSize(SYS)/(N.cotree*getSize(RLfun(1,1:N.cotree)));
disp(['Compression ratio =',num2str(compr),'%'])
disp(['Coefficients evaluation ratio =',num2str(100*nncoeff/(N.cotree^2)),'%'])
disp('... done!')
disp('-------------------------------------------------------------------')
%%
figure
spy(SYS)
drawnow
axis equal
title('ranks')
%%
disp('-------------------------------------------------------------------')
disp('solving...')
tic
x=SYS\brhs(map); x=x(invmap);
toc
disp('... done!')
disp('-------------------------------------------------------------------')
%%
if plotflag
disp('-------------------------------------------------------------------')
disp('post-processing...')
x_f=Matrix_C(:,ind.cotree)*x;
[J_r,J_i,J_norm_r,J_norm_i,J_bar] = fun_post_J(N.face,N.volu,...
    N.volu,Matrix_P0,VP,D1,x_f,1:N.face,1:N.volu);
fun_plot_J(F1,ind,Matrix_P0,J_bar,J_r,J_norm_r,J_norm_i,J_i,...
    xmin,xmax,ymin,ymax,zmin,zmax);
disp('... done!')
disp('-------------------------------------------------------------------')
end
%%
return
