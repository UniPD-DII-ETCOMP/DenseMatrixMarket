function [dist] = dist_edge_edge_MF(NN,EP,nout,mout)
dist=zeros(length(nout),length(mout));
for ii = 1:length(nout)
    n=nout(ii);
    centre_n=0.5*(NN(:,EP(1,n))+NN(:,EP(2,n)));
    m=mout;
    centre_m=0.5*(NN(:,EP(1,m))+NN(:,EP(2,m))); 
    dist(ii,:)=fun_my_norm(centre_n-centre_m);
end 
end  