function [dist] = dist_node_node_MF(NN,kkout,jjout)
dist=zeros(length(kkout),length(jjout));
ii = 1;
for kk=kkout(1):kkout(end)
jj=jjout;
dist(ii,:)=fun_my_norm(NN(1:3,kk)-NN(1:3,jj));
ii = ii+1;
end 
end 