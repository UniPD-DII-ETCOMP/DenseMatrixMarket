function [L] = fun_L_interface(N,EV,Matrix_P0,VP,EV_ind,norm_eps,ii,jj)

if size(ii,1)~=1
    ii=ii.';
end
if size(jj,1)~=1
    jj=jj.';
end
    L=fun_L_curledge_curledge_for90_st(N.volu,N.node,N.edge,N.n_EV_max,...
        [EV;zeros(1,N.edge)],Matrix_P0.',VP,[EV_ind;zeros(1,N.edge)],norm_eps,...
        N.n_EV,N.thread,length(jj),length(ii),ii,jj); 
    
global nncoeff;
nncoeff=nncoeff+length(jj)*length(ii);
end

