function [sol] = fun_post_processing(x,N,Area_c,w,mu_0,...
    Matrix_C_m,Matrix_Cb_m,ex,Matrix_P0_m,F1_m,C1_m,Fac_Ed_loc_m,...
    Aed_m,Led_m,F1_c,Matrix_P0_c,F1_e,Matrix_P0_e,Matrix_P0)
%% Extracting solution
disp('-------------------------------------------------------------------')
disp('...Extracting solution ...')
sol.jphi=x(1:N.face_con); % j DoFs
sol.Jphi=sol.jphi./Area_c; % J
sol.jm=x(N.face_con+1:N.face_con+N.edge_mag); % jm DoFs
sol.m=(1/(1j*w(1)*mu_0))*sol.jm; % m DoFs
sol.rho_m=mu_0*Matrix_C_m*sol.m; % rhom DoFs
sol.sig_m=mu_0*Matrix_Cb_m*sol.m; % sigm DoFs
disp('... done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp('...Post-processing M ...')
if ex.vol_mag
[sol.M_r,sol.M_i,sol.M_norm_r,sol.M_norm_i,sol.M_bar] = fun_post_J_Quad_face2(NaN,N.face_mag,Matrix_P0_m,F1_m,C1_m,sol.m,Fac_Ed_loc_m,Aed_m,Led_m);
sol.M_r=-sol.M_r;
sol.M_i=-sol.M_i;
end
disp('... done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp('...Plot JM ...')
green=[0 0.600000023841858 0];
green_surf=[0    0.8  0];
grey= [0.721568644046783 0.721568644046783 0.721568644046783];
blue =[0.5  0.5  1 ];
alpha=0.8;
if ex.vol_con
figure
hold on
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'CData',NaN*ones(N.face_con,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'CData',NaN*ones(N.face_mag,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'CData',NaN*ones(N.face_ext,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c, 'CData',real(sol.Jphi),'Facecolor','flat','FaceAlpha',1.0,'Edgecolor','none') 
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
colormap jet 
colorbar
title('J \Re')
cmin=min(real(sol.Jphi));
cmax=max(real(sol.Jphi));
if cmin==cmax
else
caxis([cmin cmax])
end
figure
hold on
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'CData',NaN*ones(N.face_con,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'CData',NaN*ones(N.face_mag,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'CData',NaN*ones(N.face_ext,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c, 'CData',imag(sol.Jphi),'Facecolor','flat','FaceAlpha',1.0,'Edgecolor','none') 
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
colormap jet 
colorbar
title('J \Im')
cmin=min(imag(sol.Jphi));
cmax=max(imag(sol.Jphi));
if cmin==cmax
else
caxis([cmin cmax])
end
end
%%
if ex.vol_mag	
figure
hold on
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'CData',NaN*ones(N.face_con,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'CData',NaN*ones(N.face_mag,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'CData',NaN*ones(N.face_ext,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m, 'CData',sol.M_norm_r,'Facecolor','flat','FaceAlpha',1.0,'Edgecolor','none') 
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
colormap jet 
colorbar
title('M \Re')
cmin=min(sol.M_norm_r);
cmax=max(sol.M_norm_r);
if cmin==cmax
else
caxis([cmin cmax])
end
figure
hold on
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'CData',NaN*ones(N.face_con,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'CData',NaN*ones(N.face_mag,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'CData',NaN*ones(N.face_ext,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m, 'CData',sol.M_norm_i,'Facecolor','flat','FaceAlpha',1.0,'Edgecolor','none') 
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
colormap jet 
colorbar
title('M \Im')
cmin=min(sol.M_norm_i);
cmax=max(sol.M_norm_i);
if cmin==cmax
else
caxis([cmin cmax])
end
end
%%
if ex.vol_mag	
figure
hold on
hold on
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'CData',NaN*ones(N.face_con,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'CData',NaN*ones(N.face_mag,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'CData',NaN*ones(N.face_ext,1),'Facecolor','flat','EdgeColor',grey)
quiver3_c(sol.M_bar(:,1),sol.M_bar(:,2),sol.M_bar(:,3),sol.M_r(:,1),sol.M_r(:,2),sol.M_r(:,3),sol.M_norm_r,1);
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
colormap jet 
colorbar
title('M \Re')
cmin=min(sol.M_norm_r);
cmax=max(sol.M_norm_r);
if cmin==cmax
else
caxis([cmin cmax])
end
figure
hold on
hold on
patch('Faces',F1_c(1:4,:).','Vertices',Matrix_P0_c,'CData',NaN*ones(N.face_con,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'CData',NaN*ones(N.face_mag,1),'Facecolor','flat','EdgeColor',grey)
patch('Faces',F1_e(1:4,:).','Vertices',Matrix_P0_e,'CData',NaN*ones(N.face_ext,1),'Facecolor','flat','EdgeColor',grey)
quiver3_c(sol.M_bar(:,1),sol.M_bar(:,2),sol.M_bar(:,3),sol.M_i(:,1),sol.M_i(:,2),sol.M_i(:,3),sol.M_norm_i,1);
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
colormap jet 
colorbar
title('M \Im')
cmin=min(sol.M_norm_i);
cmax=max(sol.M_norm_i);
if cmin==cmax
else
caxis([cmin cmax])
end
end
disp('... done!')
disp('-------------------------------------------------------------------')
end

