function [] = fun_plot_J(F1,ind,Matrix_P0,J_bar,J_r,J_norm_r,J_norm_i,J_i)
tic
grey= [0.721568644046783 0.721568644046783 0.721568644046783];
%% real
figure
hold on
patch('Faces',F1(1:3,ind.face_free).','Vertices',Matrix_P0,'Facecolor','b','FaceAlpha',0.01,'EdgeColor',grey) 
cmax=max(max(J_norm_r));
cmin=min(min(J_norm_r));
scal=1;
quiver3_c_scal(J_bar(:,1),J_bar(:,2),J_bar(:,3),J_r(:,1),J_r(:,2),J_r(:,3),J_norm_r,scal);
axis equal
zlim auto
xlabel('x')
ylabel('y')
zlabel('z')
colorbar;
caxis([cmin cmax]);
view(3)
title('\Re(J)')
%% imag
figure
hold on
patch('Faces',F1(1:3,ind.face_free).','Vertices',Matrix_P0,'Facecolor','b','FaceAlpha',0.01,'EdgeColor',grey) 
cmax=max(max(J_norm_i));
cmin=min(min(J_norm_i));
scal=1;
quiver3_c_scal(J_bar(:,1),J_bar(:,2),J_bar(:,3),J_i(:,1),J_i(:,2),J_i(:,3),J_norm_i,scal);
axis equal
zlim auto
xlabel('x')
ylabel('y')
zlabel('z')
colorbar;
caxis([cmin cmax]);
view(3)
title('\Im(J)')
%%
toc