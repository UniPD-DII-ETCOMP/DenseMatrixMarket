function [R] = funRmatlab(N,VP,Matrix_P0,VE,VE_ind,rho,ind)
RR=zeros(3,6*N.n_EV_max*6*N.n_EV_max*N.edge);
qq=1;
for ii=1:N.volu
   pvt=Matrix_P0(VP(:,ii),1:3).';
    [nnorm,area,vol] = tetareavol(pvt);
   for jj=1:6
      ind1=abs(VE(jj,ii));
      w1 = fun_whitney_tetra_rot_edge(nnorm,area,vol,VE_ind(jj,ii));
      for hh = 1:6
        ind2=abs(VE(hh,ii)); 
        w2 = fun_whitney_tetra_rot_edge(nnorm,area,vol,VE_ind(hh,ii));
        RR(1,qq)=ind1;
        RR(2,qq)=ind2;
        RR(3,qq)=rho(ii)*vol*sign(VE(jj,ii))*sign(VE(hh,ii))*fun_my_dot(w1,w2);
        qq=qq+1;
      end
   end
end
R=sparse(RR(1,1:qq-1),RR(2,1:qq-1),RR(3,1:qq-1),N.edge,N.edge);
R=R(ind.cotree,ind.cotree);
end

