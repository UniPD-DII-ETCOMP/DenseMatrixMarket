function [PP]=funTrilinear2(csi,eta,zeta,Face)
%% orientazione Paolo
% kcsi=[-1,+1,+1,-1,-1,+1,+1,-1];
% keta=[-1,-1,+1,+1,-1,-1,+1,+1];
% kzeta=[-1,-1,-1,-1,+1,+1,+1,+1];
% sn=0.125*(1+kcsi*csi).*(1+keta*eta).*(1+kzeta*zeta);
% PP(1)=sn*QQ(:,1);
% PP(2)=sn*QQ(:,2);
% PP(3)=sn*QQ(:,3);
%%
QQ(1:4,[1,3])=Face(1:4,[1,3]);
QQ(5:8,[1,3])=Face(1:4,[1,3]);
QQ(1:4,2)=-1.0;
QQ(5:8,2)=1.0;
%% orientazione Riccardo
kcsi= [-1,-1,+1,+1,-1,-1,+1,+1];
keta= [-1,-1,-1,-1,+1,+1,+1,+1];
kzeta=[+1,-1,-1,+1,+1,-1,-1,+1];
sn=0.125*(1+kcsi.*csi).*(1+keta.*eta).*(1+kzeta.*zeta);
PP(1)=sn*QQ(:,1);
PP(2)=sn*QQ(:,2);
PP(3)=sn*QQ(:,3);
end
