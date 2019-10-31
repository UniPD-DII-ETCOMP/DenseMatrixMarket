%%
clear 
close all
clc
%%
n=8;
outer_radius=1.00;
inner_radius=0.9;
%%
[Matrix_P0,VPq] = cubedsphere(n); VPq=VPq.';
figure
patch('Faces',VPq(1:4,:).','Vertices',Matrix_P0,'Facecolor','b','FaceAlpha',0.2)
axis equal
view(3)
title('quadrilateral')
drawnow
%%
Npq=size(Matrix_P0,1);
Nq=size(VPq,2);
for ii = 1:Nq % reo quad face
    p1=Matrix_P0(VPq(1,ii),:);
    p2=Matrix_P0(VPq(2,ii),:);
    p3=Matrix_P0(VPq(3,ii),:);
    p4=Matrix_P0(VPq(4,ii),:);
    cent=(p1+p2+p3+p4)*0.25;
    v21=p2-p1;
    v32=p3-p2;
    nout=cross(v21,v32);
    val=dot(nout,cent);
    if sign(val)<0
        VPq(:,ii)=VPq([1,4,3,2],ii);
    end
end
%% from quad to hexa
Matrix_P0=[Matrix_P0*outer_radius;Matrix_P0*inner_radius];
Nhexa=Nq;
VPh=[VPq;VPq+Npq];
figure
hexa_mesh2(VPh,Matrix_P0,0.2,'r')
axis equal
view(3)
title('hexahedra')
drawnow
%% from hexa to tetra
Ntetra=Nhexa*6
VP=zeros(4,Ntetra);
hh=0;
for ii = 1:Nhexa/6*2
    VP(1:4,hh+1)=VPh([1,4,6,5],ii);
    VP(1:4,hh+2)=VPh([1,4,2,6],ii);
    VP(1:4,hh+3)=VPh([2,3,6,4],ii);
    VP(1:4,hh+4)=VPh([4,8,6,5],ii);
    VP(1:4,hh+5)=VPh([4,8,7,6],ii);
    VP(1:4,hh+6)=VPh([4,3,7,6],ii);
    hh=hh+6;
end
for ii = Nhexa/6*2+1:Nhexa/6*4
    map=[4,3,2,1];
    h=VPh([map,map+4],ii);
    VP(1:4,hh+1)=h([1,2,8,4]);
    VP(1:4,hh+2)=h([1,2,5,8]);
    VP(1:4,hh+3)=h([5,6,2,8]);
    VP(1:4,hh+4)=h([2,3,4,8]);
    VP(1:4,hh+5)=h([2,7,3,8]);
    VP(1:4,hh+6)=h([6,7,2,8]);
    hh=hh+6;
end
for ii = Nhexa/6*4+1:Nhexa/6*5
    map=[3,4,1,2];
    h=VPh([map,map+4],ii);
    VP(1:4,hh+1)=h([1,2,8,4]);
    VP(1:4,hh+2)=h([1,2,5,8]);
    VP(1:4,hh+3)=h([5,6,2,8]);
    VP(1:4,hh+4)=h([2,3,4,8]);
    VP(1:4,hh+5)=h([2,7,3,8]);
    VP(1:4,hh+6)=h([6,7,2,8]);
    hh=hh+6;
end
for ii = Nhexa/6*5+1:Nhexa/6*6
    map=[4,1,2,3];
    h=VPh([map,map+4],ii);
    VP(1:4,hh+1)=h([1,2,8,4]);
    VP(1:4,hh+2)=h([1,2,5,8]);
    VP(1:4,hh+3)=h([5,6,2,8]);
    VP(1:4,hh+4)=h([2,3,4,8]);
    VP(1:4,hh+5)=h([2,7,3,8]);
    VP(1:4,hh+6)=h([6,7,2,8]);
    hh=hh+6;
end
%% reordering
if 1
a = zeros(size(VP,2),1);
for ii = 1:size(VP,2)
    e12 = Matrix_P0(VP(2,ii),:)-Matrix_P0(VP(1,ii),:);
    e13 = Matrix_P0(VP(3,ii),:)-Matrix_P0(VP(1,ii),:);
    e14 = Matrix_P0(VP(4,ii),:)-Matrix_P0(VP(1,ii),:);
    vec1=cross(e12,e13);
    a(ii,1) = dot(vec1,e14);
    if a(ii,1) < 0
        VP(:,ii) = VP([1,2,4,3],ii);
    end
end
end
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
title('tetrahedra')
figure
patch('Faces',F1(1:3,ind_face_free).','Vertices',Matrix_P0,'Facecolor',[0.7,0.7,0.7],'FaceAlpha',1)
axis equal
view(3)
title('tetrahedra')

%% check
N_face_free=length(ind_face_free)
%%
save data.mat Matrix_G  C1 D1 F1 G1 Matrix_C Matrix_D Matrix_P0 VP 