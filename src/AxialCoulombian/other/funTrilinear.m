function [PP]=funTrilinear(csi,eta,zeta,QQ)
%% orientazione Paolo
% kcsi=[-1,+1,+1,-1,-1,+1,+1,-1];
% keta=[-1,-1,+1,+1,-1,-1,+1,+1];
% kzeta=[-1,-1,-1,-1,+1,+1,+1,+1];
% sn=0.125*(1+kcsi*csi).*(1+keta*eta).*(1+kzeta*zeta);
% PP(1)=sn*QQ(:,1);
% PP(2)=sn*QQ(:,2);
% PP(3)=sn*QQ(:,3);
%% orientazione Riccardo
kcsi= [-1,-1,+1,+1,-1,-1,+1,+1];
keta= [-1,-1,-1,-1,+1,+1,+1,+1];
kzeta=[+1,-1,-1,+1,+1,-1,-1,+1];
sn=0.125*(1+kcsi*csi).*(1+keta*eta).*(1+kzeta*zeta);
PP(1)=sn*QQ(:,1);
PP(2)=sn*QQ(:,2);
PP(3)=sn*QQ(:,3);
end
