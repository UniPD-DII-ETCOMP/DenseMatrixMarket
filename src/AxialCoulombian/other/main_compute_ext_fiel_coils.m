%%
jext=zeros(N.face_ext,1);
for ii = 1:N.face_ext
    jextii=J_ext(bar_e(1:3,ii));
    jext(ii)=jextii(2)*Area_e(ii);
end
%%
if ex.vol_con
    b_e=-1j*w(1)*Matrix_Lce*jext;
else
b_e=[];
end
%%
if ex.vol_mag
        b_h=-Matrix_Projm*Matrix_G_m*Matrix_Nme*jext;
else
b_h=[];
end
%%
brhs=[b_e;b_h];