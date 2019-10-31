function [] = fun_plot_geo(F1,ind,Matrix_P0,xmin,xmax,ymin,ymax,zmin,zmax)
tic
grey= [0.721568644046783 0.721568644046783 0.721568644046783];
figure
hold on
plot3(0,0,0,'b')
alpha=0.3;
patch('Faces',F1(1:3,ind.face_free).','Vertices',Matrix_P0,'Facecolor',grey,'FaceAlpha',alpha) 
axis equal
title('mesh')
xlabel('x')
ylabel('y')
zlabel('z')
axis([xmin,xmax,ymin,ymax,zmin,zmax])
view(3)
toc
end