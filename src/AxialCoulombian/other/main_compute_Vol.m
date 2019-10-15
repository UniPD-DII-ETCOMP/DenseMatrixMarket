%%  con
Vol_c=zeros(N.face_con,1);
for ii = 1:N.face_con
   Vol_c(ii,1)=Area_c(ii,1)*2*pi*bar_c(1,ii);
end
%% mag
Vol_m=zeros(N.face_mag,1);
for ii = 1:N.face_mag
   Vol_m(ii,1)=Area_m(ii,1)*2*pi*bar_m(1,ii);
end
