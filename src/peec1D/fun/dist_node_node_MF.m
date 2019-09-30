function [dist] = dist_node_node_MF(NN,kkout,jjout)
dist=zeros(length(kkout),length(jjout));
for kkk=1:length(kkout)
kk=kkout(kkk);
jj=jjout;
dist(kkk,:)=fun_my_norm(NN(1:3,kk)-NN(1:3,jj));
end 
end 