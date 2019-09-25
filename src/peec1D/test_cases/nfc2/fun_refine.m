function [G,NN,Matrix_G] = fun_refine(NN,G)
nNode=size(NN,2);
nSticks=size(G,2);
bar_stick = 0.5*(NN(:,G(1,:))+NN(:,G(2,:)));
NN=[NN,bar_stick];
G2=[G;nNode+1:nNode+nSticks];
G=[G2([1,3],:),G2([2,3],:)];
nn=size(NN,2);
ne=size(G,2);
Matrix_G = sparse(1:ne,G(1,1:ne),-ones(ne,1),ne,nn);
Matrix_G = Matrix_G+sparse(1:ne,G(2,1:ne),ones(ne,1),ne,nn);
end

