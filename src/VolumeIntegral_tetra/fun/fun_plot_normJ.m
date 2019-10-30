function [] = fun_plot_normJ(J_norm_r,J_norm_i,Matrix_D,ind,xmin,xmax,ymin,ymax,zmin,zmax,F1,Matrix_P0)
mapJv2f=Matrix_D(:,ind.face_free).'*J_norm_r;
cmax=max(max(mapJv2f));
cmin=min(min(mapJv2f));
figure
% patch('Faces',F1(1:4,ind.face_free).','Vertices',Matrix_P0,'CData',mapJv2f,'Facecolor','flat','EdgeColor',[0.721568644046783 0.721568644046783 0.721568644046783])
patch('Faces',F1(1:4,ind.face_free).','Vertices',Matrix_P0,'CData',mapJv2f,'Facecolor','flat','EdgeColor','none')
colormap jet
colorbar;
axis equal
caxis([cmin cmax]);
axis([xmin,xmax,ymin,ymax,zmin,zmax])
view(3)
title('|\Re(J)|')
mapJv2f=Matrix_D(:,ind.face_free).'*J_norm_i;
cmax=max(max(mapJv2f));
cmin=min(min(mapJv2f));
figure
patch('Faces',F1(1:4,ind.face_free).','Vertices',Matrix_P0,'CData',mapJv2f,'Facecolor','flat','EdgeColor','none')
% patch('Faces',F1(1:4,ind.face_free).','Vertices',Matrix_P0,'CData',mapJv2f,'Facecolor','flat','EdgeColor',[0.721568644046783 0.721568644046783 0.721568644046783])
colormap jet
colorbar;
axis equal
caxis([cmin cmax]);
axis([xmin,xmax,ymin,ymax,zmin,zmax])
view(3)
title('|\Im(J)|')
end

