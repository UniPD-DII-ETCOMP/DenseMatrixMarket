function [vp] = vecprod(v1,v2)
vp(1:3)=0;
vp(1)=   v1(2)*v2(3)-v1(3)*v2(2);
vp(2)= -(v1(1)*v2(3)-v1(3)*v2(1));
vp(3)=   v1(1)*v2(2)-v1(2)*v2(1);
end   
