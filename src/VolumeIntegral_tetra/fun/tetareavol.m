function [nnorm,area,vol] = tetareavol(v)
nnorm(1:3,1:4)=0;
area(1:4,1)=0;
permut =[ 1,2,3,4; 4,3,4,1; 2,1,1,3; 3,4,2,2 ].';
for i=1:4
r(:)=v(:,permut(i,3))-v(:,permut(i,2));
s(:)=v(:,permut(i,4))-v(:,permut(i,2));
vp(1)=r(2)*s(3)-r(3)*s(2);
vp(2)=-(r(1)*s(3)-r(3)*s(1));
vp(3)=r(1)*s(2)-r(2)*s(1);
area(i)=sqrt(vp(1)^2+vp(2)^2+vp(3)^2);
nnorm(:,i)=vp(:)/area(i);
end
area(:)=0.5d0*area(:);
t(:)=v(:,permut(4,2))-v(:,permut(4,1));
sp=fun_my_dot(vp,t);
vol=abs(sp)/6.0d0;
end  