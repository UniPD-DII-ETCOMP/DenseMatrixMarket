%%
clear 
close all
clc
%%
C = textread('data_comsol.mphtxt', '%s','delimiter', '\n'); %#ok<DTXTRD>
[Matrix_P0,VP] = fun_comsol_extract(C);
%%
[G1,C1]=gcd_quad(VP,size(VP,2));
%% number of nodes edges faces 
nn = max(max(VP))
ne = size(G1,2)
nf = size(C1,2)
%% G edges xnodes
Matrix_G = sparse(1:ne,G1(1,1:ne),-ones(ne,1),ne,nn);
Matrix_G = Matrix_G+sparse(1:ne,G1(2,1:ne),ones(ne,1),ne,nn);
%% C faces x edges
Cpos = (C1+abs(C1))/2;
Cneg = (C1-abs(C1))/2;
[r_pos_c,c_pos_c,val_pos_c] = find(Cpos);
[r_neg_c,c_neg_c,val_neg_c] = find(Cneg);
Matrix_C = sparse(c_pos_c,val_pos_c,ones(size(c_pos_c,1),1),nf,ne);
Matrix_C = Matrix_C+sparse(c_neg_c,abs(val_neg_c),-ones(size(c_neg_c,1),1),nf,ne);
%% check rot(grad())
disp('check rot(grad()) = 0...')
CxG=Matrix_C*Matrix_G;
find(CxG)
%% save
clearvars -except keepVariables G1 C1 VP Matrix_P0 Matrix_C Matrix_G 
%% barycenter
N.face=size(C1,2);
barf=zeros(3,N.face);
for ii = 1:N.face
    barf(1:3,ii)=(sum(Matrix_P0(VP(1:4,ii),:))/4).';
end
%%
F1=VP;
%%
N.edge=size(G1,2);
N.face=size(C1,2);
%%
figure
patch('Faces',F1(1:4,:).','Vertices',Matrix_P0,'Facecolor','b','FaceAlpha',0.2)
axis equal
view(0,0)
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
%%
G2=fun_G2(F1,C1,N,G1);
%%
figure
patch('Faces',F1(1:4,:).','Vertices',Matrix_P0,'Facecolor','g','FaceAlpha',0.2)
axis equal
view(0,0)
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
%% check normal vector
vec=zeros(N.face,3);
for ii = 1:N.face
   vec(ii,1:3)=cross(Matrix_P0(F1(2,ii),:)-Matrix_P0(F1(1,ii),:),...
                     Matrix_P0(F1(3,ii),:)-Matrix_P0(F1(2,ii),:));
   if vec(ii,2)>0
      F1(:,ii)=flip(F1(:,ii)); 
      vec(ii,2)=-vec(ii,2);
   end
end
%%
C1_m=C1;
G1_m=G1;
Matrix_G_m=Matrix_G;
Matrix_C_m=Matrix_C;
Matrix_P0_m=Matrix_P0;
F1_m=F1;
G2_m=G2;
%%
figure
hold on
patch('Faces',F1(1:4,:).','Vertices',Matrix_P0,'Facecolor','g','FaceAlpha',0.2)
quiver3(barf(1,:),barf(2,:),barf(3,:),...
        vec(:,1).',vec(:,2).',vec(:,3).')
axis equal
view(0,0)
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
%%
save data_mag.mat Matrix_G_m G2_m C1_m  F1_m G1_m Matrix_C_m  Matrix_P0_m  




