function [rotwi]=fun_whitney_tetra_rot_edge(nnorm,area,vol,k)
locedg = [1,2,3,4,4,4;2,3,1,1,2,3].';
[gradlambda]= x2gradlambda(nnorm,area,vol);
nd1loc=locedg(k,1);
nd2loc=locedg(k,2);
rotwi = vecprod(gradlambda(:,nd1loc),gradlambda(:,nd2loc));
rotwi=2.0d0*rotwi;
end 