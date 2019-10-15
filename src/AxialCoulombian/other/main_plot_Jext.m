%%
sol.Jext=jext./Area_e;
if ex.vol_ext
figure
hold on
% subplot(2,1,1)
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'CData',NaN*ones(N.face_con,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'CData',NaN*ones(N.face_mag,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'CData',NaN*ones(N.face_ext,1),'Facecolor','flat','EdgeColor',grey)

patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e, 'CData',real(sol.Jext),'Facecolor','flat','FaceAlpha',1.0,'Edgecolor','none') 
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
colormap jet 
colorbar
title('J \Re')
cmin=min(real(sol.Jext));
cmax=max(real(sol.Jext));
if cmin==cmax
else
caxis([cmin cmax])
end
figure
hold on
% subplot(2,1,2)
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'CData',NaN*ones(N.face_con,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'CData',NaN*ones(N.face_mag,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'CData',NaN*ones(N.face_ext,1),'Facecolor','flat','EdgeColor',grey)


patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e, 'CData',imag(sol.Jext),'Facecolor','flat','FaceAlpha',1.0,'Edgecolor','none') 
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
colormap jet 
colorbar
title('J \Im')
cmin=min(imag(sol.Jext));
cmax=max(imag(sol.Jext));
if cmin==cmax
else
caxis([cmin cmax])
end
end
%%
