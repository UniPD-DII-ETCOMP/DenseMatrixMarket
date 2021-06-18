clear
close all
clc
restoredefaultpath
%% BEGIN USER SETTINGS
test_case_dir = 'test_spherical_shell_from_user';
f = 50;% Hz
rrho_c=1/56e3;% Ohm m
N.thread = 22; 
E_ext=@(x,y,z) [-1j*y, 1j*x, 0];
plotflag = 1; %graphics postprocessing flag (1 = yes, 0 = no)
thr=1.3; %threshold for close to far field, set this value greater than the max distance between two mesh elements to improve accuracy
% END USER SETTINGS
%% load data
cd test_cases; cd(test_case_dir); load data.mat; cd ..; cd ..
%% pre_procesing
cd('fun'); addpath(pwd); cd ..
[rho_c,N,ind,EV,VE,EV_ind,VE_ind,Q,brhs,norm_eps,w,bared,...
    xmin,xmax,ymin,ymax,zmin,zmax] = fun_pre_processing(F1,D1,C1,...
    G1,VP,Matrix_D,Matrix_C,Matrix_G,Matrix_P0,E_ext,rrho_c,N,f);
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
hhout=ind.cotree;
kkout=ind.cotree;
disp('MAKING L...')
try 
tic
[L] = fun_L_curledge_curledge_for90_st1(N.volu,N.node,N.edge,N.n_EV_max,...
    [EV;zeros(1,N.edge)],Matrix_P0.',VP,[EV_ind;zeros(1,N.edge)],0,...
    N.n_EV,N.thread,length(kkout),length(hhout),hhout,kkout,bared,thr); L=0.5*(L+L.');
toc
mtlb=0;
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
    tic
    [L] = fun_L_curledge_curledge_pp(EV,norm_eps,N.n_EV,X8e,Y8e,Z8e,...
        W8e,X11e,Y11e,Z11e,W11e,we,length(hhout),length(kkout),hhout,kkout); L=0.5*(L+L.');
    toc
end
disp('... done!')
disp('-------------------------------------------------------------------')
%% SYS
disp('-------------------------------------------------------------------')
disp('R+1j*w*L...')
tic
SYS=R+1j*w*L;
toc
disp('... done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp('solving...')
tic
x=SYS\brhs;
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
fun_plot_normJ(J_norm_r,J_norm_i,Matrix_D,ind,xmin,xmax,ymin,ymax,zmin,zmax,F1,Matrix_P0)
disp('... done!')
disp('-------------------------------------------------------------------')
end
%%
return
