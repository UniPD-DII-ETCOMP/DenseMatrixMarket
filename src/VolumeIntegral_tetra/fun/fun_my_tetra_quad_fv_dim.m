function [X,Y,Z,W] = fun_my_tetra_quad_fv_dim(tetra,V)
F1=(tetra(1:3,1)+tetra(1:3,2)+tetra(1:3,3))/3.0;%sum(tetra(1:3,1:3))/3;
F2=(tetra(1:3,2)+tetra(1:3,3)+tetra(1:3,4))/3.0;%sum(tetra(2:4,1:3))/3;
F3=(tetra(1:3,1)+tetra(1:3,3)+tetra(1:3,4))/3.0;%sum(tetra([1,3,4],:))/3;
F4=(tetra(1:3,1)+tetra(1:3,2)+tetra(1:3,4))/3.0;%sum(tetra([1,2,4],:))/3;
X = [tetra(1,1:4),F1(1),F2(1),F3(1),F4(1)];
Y = [tetra(2,1:4),F1(2),F2(2),F3(2),F4(2)];
Z = [tetra(3,1:4),F1(3),F2(3),F3(3),F4(3)];
W  = [V/40,V/40,V/40,V/40,9*V/40,9*V/40,9*V/40,9*V/40]; %[V/40*ones(4,1);9*V/40*ones(4,1)]
end 