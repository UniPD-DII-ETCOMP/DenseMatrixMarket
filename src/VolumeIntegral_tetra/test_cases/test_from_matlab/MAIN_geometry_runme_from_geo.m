%%
clear 
close all
clc
%%
model = createpde(1);
importGeometry(model,'Plate10x10x1.stl');
msh=generateMesh(model,'GeometricOrder','linear');
pdeplot3D(model)
Matrix_P0=msh.Nodes.';
VP=msh.Elements;

%%
dad=pwd;
cd ..
[G1,C1,D1,F1,Matrix_C,Matrix_G,Matrix_D] = fun_create_data(VP);
cd(dad);
%%
ind_face_free=find(F1(5,:)==0);
figure
patch('Faces',F1(1:3,ind_face_free).','Vertices',Matrix_P0,'Facecolor','b','FaceAlpha',0.2)
axis equal
view(3)
%%
save data.mat Matrix_G  C1 D1 F1 G1 Matrix_C Matrix_D Matrix_P0 VP 