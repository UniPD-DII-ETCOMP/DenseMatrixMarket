tic
green=[0 0.600000023841858 0];
green_surf=[0    0.8  0];
grey= [0.721568644046783 0.721568644046783 0.721568644046783];
blue =[0.5  0.5  1 ];
alpha=0.8;
figure
hold on
patch('Faces',[1:4].','Vertices',zeros(4,3)*NaN,'Facecolor','b','FaceAlpha',alpha) 
patch('Faces',[1:4].','Vertices',zeros(4,3)*NaN,'Facecolor',green,'FaceAlpha',alpha) 
patch('Faces',[1:4].','Vertices',zeros(4,3)*NaN,'Facecolor','r','FaceAlpha',alpha) 
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'Facecolor','b','FaceAlpha',alpha) 
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'Facecolor',green,'FaceAlpha',alpha) 
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'Facecolor','r','FaceAlpha',alpha) 
legend('con','mag','ext-coil')
axis equal
title('Boundary')
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
Matrix_P0=[Matrix_P0_c;Matrix_P0_m;Matrix_P0_e];
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
toc