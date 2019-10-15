function [wf] = fun_whitney_face_hexa_6_quad(Hexa,P,NP)
%% descrizione
% compute shape function hexa 
% face node
% face(1,1:4)=[1 2 3 4];
% face(2,1:4)=[4 3 7 8];
% face(3,1:4)=[8 7 6 5];
% face(4,1:4)=[5 6 2 1];
% face(5,1:4)=[2 6 7 3];
% face(6,1:4)=[4 8 5 1];
% edge_node
% edge(1 ,1:2)=[1 2];
% edge(2 ,1:2)=[2 3];
% edge(3 ,1:2)=[3 4];
% edge(4 ,1:2)=[4 1];
% edge(5 ,1:2)=[5 6];
% edge(6 ,1:2)=[6 7];
% edge(7 ,1:2)=[7 8];
% edge(8 ,1:2)=[8 5];
% edge(9 ,1:2)=[1 5];
% edge(10,1:2)=[2 6];
% edge(11,1:2)=[3 7];
% edge(12,1:2)=[4 8];
% edge_face
% edge_face(1 ,1:2)=[6 5];
% edge_face(2 ,1:2)=[4 2];
% edge_face(3 ,1:2)=[5 6];
% edge_face(4 ,1:2)=[2 4];
% edge_face(5 ,1:2)=[6 5];
% edge_face(6 ,1:2)=[4 2];
% edge_face(7 ,1:2)=[5 6];
% edge_face(8 ,1:2)=[2 4];
% edge_face(9 ,1:2)=[1 3];
% edge_face(10,1:2)=[1 3];
% edge_face(11,1:2)=[1 3];
% edge_face(12,1:2)=[1 3];
%% Input 
% Hexa =  xyz hexa (3,8)
% NP = N target points 
% P = target points in local coordinates
%% Output
% w = shape face (3,NP)
%% local coordinates
xi =P(1,1:NP);
eta=P(2,1:NP);
zet=P(3,1:NP);
%% Hexa master
% Hexa_master=[-1,-1, 1;...  %
%              -1,-1,-1;...
%               1,-1,-1;...
%               1,-1, 1;...
%              -1, 1, 1;...
%              -1, 1,-1;...
%               1, 1,-1;...
%               1, 1, 1];
%% grad_face (xi, eta,zeta)
grad_snF=zeros(6,3);
grad_snF(1,1:3)=[0.0d0,-0.5d0,0.0d0];% fac 1 -> node(1 2 3 4) 
grad_snF(2,1:3)=[0.5d0,0.0d0,0.0d0];% fac 2 -> node(4 3 7 8)
grad_snF(3,1:3)=[0.0d0,0.5d0,0.0d0];% fac 3 -> node(8 7 6 5)
grad_snF(4,1:3)=[-0.5d0,0.0d0,0.0d0];% fac 4 -> node(5 6 2 1)
grad_snF(5,1:3)=[0.0d0,0.0d0,-0.5d0];% fac 5 -> node(2 6 7 3)
grad_snF(6,1:3)=[0.0d0,0.0d0,0.5d0];% fac 6 -> node(4 8 5 1)
%% shape 
wf=zeros(3,NP,4);
%% 
for ii = 1:NP
%  [--+][---][+--][+-+][-++][-+-][++-][+++]
J(1:3,1)=(-(1.0d0-eta(ii))*(1.0d0+zet(ii))*Hexa(1,1:3)-(1.0d0-eta(ii))*(1.0d0-zet(ii))*Hexa(2,1:3) ...
		  +(1.0d0-eta(ii))*(1.0d0-zet(ii))*Hexa(3,1:3)+(1.0d0-eta(ii))*(1.0d0+zet(ii))*Hexa(4,1:3) ...
		  -(1.0d0+eta(ii))*(1.0d0+zet(ii))*Hexa(5,1:3)-(1.0d0+eta(ii))*(1.0d0-zet(ii))*Hexa(6,1:3) ...
		  +(1.0d0+eta(ii))*(1.0d0-zet(ii))*Hexa(7,1:3)+(1.0d0+eta(ii))*(1.0d0+zet(ii))*Hexa(8,1:3))*0.125d0;

J(1:3,2)=(-(1.0d0- xi(ii))*(1.0d0+zet(ii))*Hexa(1,1:3)-(1.0d0- xi(ii))*(1.0d0-zet(ii))*Hexa(2,1:3) ...
		  -(1.0d0+ xi(ii))*(1.0d0-zet(ii))*Hexa(3,1:3)-(1.0d0+ xi(ii))*(1.0d0+zet(ii))*Hexa(4,1:3) ...
		  +(1.0d0- xi(ii))*(1.0d0+zet(ii))*Hexa(5,1:3)+(1.0d0- xi(ii))*(1.0d0-zet(ii))*Hexa(6,1:3) ...
		  +(1.0d0+ xi(ii))*(1.0d0-zet(ii))*Hexa(7,1:3)+(1.0d0+ xi(ii))*(1.0d0+zet(ii))*Hexa(8,1:3))*0.125d0;

J(1:3,3)=(+(1.0d0- xi(ii))*(1.0d0-eta(ii))*Hexa(1,1:3)-(1.0d0- xi(ii))*(1.0d0-eta(ii))*Hexa(2,1:3) ...
		  -(1.0d0+ xi(ii))*(1.0d0-eta(ii))*Hexa(3,1:3)+(1.0d0+ xi(ii))*(1.0d0-eta(ii))*Hexa(4,1:3) ...
		  +(1.0d0- xi(ii))*(1.0d0+eta(ii))*Hexa(5,1:3)-(1.0d0- xi(ii))*(1.0d0+eta(ii))*Hexa(6,1:3) ...
		  -(1.0d0+ xi(ii))*(1.0d0+eta(ii))*Hexa(7,1:3)+(1.0d0+ xi(ii))*(1.0d0+eta(ii))*Hexa(8,1:3))*0.125d0;
%% interpolanti nodali 
sn(1)=1.0d0/8.0d0*(1.0d0-xi(ii))*(1.0d0-eta(ii))*(1.0d0+zet(ii));
sn(2)=1.0d0/8.0d0*(1.0d0-xi(ii))*(1.0d0-eta(ii))*(1.0d0-zet(ii));
sn(3)=1.0d0/8.0d0*(1.0d0+xi(ii))*(1.0d0-eta(ii))*(1.0d0-zet(ii));
sn(4)=1.0d0/8.0d0*(1.0d0+xi(ii))*(1.0d0-eta(ii))*(1.0d0+zet(ii));
sn(5)=1.0d0/8.0d0*(1.0d0-xi(ii))*(1.0d0+eta(ii))*(1.0d0+zet(ii));
sn(6)=1.0d0/8.0d0*(1.0d0-xi(ii))*(1.0d0+eta(ii))*(1.0d0-zet(ii));
sn(7)=1.0d0/8.0d0*(1.0d0+xi(ii))*(1.0d0+eta(ii))*(1.0d0-zet(ii));
sn(8)=1.0d0/8.0d0*(1.0d0+xi(ii))*(1.0d0+eta(ii))*(1.0d0+zet(ii)); 
%% Fun face 
% Fsn(1)=sn(1)+sn(2)+sn(3)+sn(4);
% Fsn(2)=sn(4)+sn(3)+sn(7)+sn(8);
% Fsn(3)=sn(8)+sn(7)+sn(6)+sn(5);
% Fsn(4)=sn(5)+sn(6)+sn(2)+sn(1);
% Fsn(5)=sn(2)+sn(6)+sn(7)+sn(3);
% Fsn(6)=sn(4)+sn(8)+sn(5)+sn(1);
%% shape 
wf_loc(1:3,3)=sn(4)*cross(grad_snF(6,1:3),grad_snF(1,1:3))+sn(3)*cross(grad_snF(1,1:3),grad_snF(5,1:3))+sn(7)*cross(grad_snF(5,1:3),grad_snF(3,1:3))+sn(8)*cross(grad_snF(3,1:3),grad_snF(6,1:3));
wf_loc(1:3,1)=sn(5)*cross(grad_snF(6,1:3),grad_snF(3,1:3))+sn(6)*cross(grad_snF(3,1:3),grad_snF(5,1:3))+sn(2)*cross(grad_snF(5,1:3),grad_snF(1,1:3))+sn(1)*cross(grad_snF(1,1:3),grad_snF(6,1:3));
wf_loc(1:3,2)=sn(2)*cross(grad_snF(1,1:3),grad_snF(4,1:3))+sn(6)*cross(grad_snF(4,1:3),grad_snF(3,1:3))+sn(7)*cross(grad_snF(3,1:3),grad_snF(2,1:3))+sn(3)*cross(grad_snF(2,1:3),grad_snF(1,1:3));
wf_loc(1:3,4)=sn(4)*cross(grad_snF(1,1:3),grad_snF(2,1:3))+sn(8)*cross(grad_snF(2,1:3),grad_snF(3,1:3))+sn(5)*cross(grad_snF(3,1:3),grad_snF(4,1:3))+sn(1)*cross(grad_snF(4,1:3),grad_snF(1,1:3));
%% Piola Transfrom div-conforming 
for jj = 1:4
wf(1:3,ii,jj)=(1/det(J))*(J(1:3,1:3))*wf_loc(1:3,jj);
end
%% Piola Transfrom di curl-conforming 
% for jj = 1:12
% ww(1:3,ii,jj)=INV_T(J(1:3,1:3))*we_loc(1:3,jj);
% end
end
%%
end

