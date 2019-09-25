function [cmin] = fun_PEEC_1D_plot(U_obj,NN,rho_charge,I_obj,bar_stick,Jrenorm,Jimnorm,G,Jre,Jim)
[xmin,xmax,ymin,ymax,zmin,zmax] = fun_domain_limits(NN);
%% plot electric potential
%% Potential 
figure
cmin=min(min(real(U_obj)));
cmax=max(max(real(U_obj)));
scatter3(NN(1,:),NN(2,:),NN(3,:),3,'CData',real(U_obj),'MarkerFaceColor','flat')
axis equal
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title('Phi_e \Re')
caxis([cmin cmax])
%% Potential 
figure
cmin=min(min(imag(U_obj)));
cmax=max(max(imag(U_obj)));
scatter3(NN(1,:),NN(2,:),NN(3,:),3,'CData',imag(U_obj),'MarkerFaceColor','flat')
axis equal
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title('Phi_e \Im')
caxis([cmin cmax])
%% Charge 
figure
cmin=min(min(real(rho_charge)));
cmax=max(max(real(rho_charge)));
scatter3(NN(1,:),NN(2,:),NN(3,:),3,'CData',real(rho_charge),'MarkerFaceColor','flat')
axis equal
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title('Q_e \Re')
caxis([cmin cmax])
%% Charge 
figure
cmin=min(min(imag(rho_charge)));
cmax=max(max(imag(rho_charge)));
scatter3(NN(1,:),NN(2,:),NN(3,:),3,'CData',imag(rho_charge),'MarkerFaceColor','flat')
axis equal
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title('Q_e \Im')
caxis([cmin cmax])
%% Current
figure 
cmax=max(max(abs(real(I_obj))));
cmin=min(min(abs(real(I_obj))));
scatter3(bar_stick(:,1),bar_stick(:,2),bar_stick(:,3),3,'CData',abs(real(I_obj)),'MarkerFaceColor','flat')
axis equal
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title('I \Re')
caxis([cmin cmax])
%%
figure 
cmax=max(max(abs(imag(I_obj))));
cmin=min(min(abs(imag(I_obj))));
scatter3(bar_stick(:,1),bar_stick(:,2),bar_stick(:,3),3,'CData',abs(imag(I_obj)),'MarkerFaceColor','flat')
axis equal
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title('I \Im')
caxis([cmin cmax])
%% J vec re
cmax=max(max(Jrenorm));
cmin=min(min(Jrenorm));
scal=0.7;
figure
hold on
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color',[0.65 0.65 0.65])
quiver3_c_scal(bar_stick(:,1).',bar_stick(:,2).',bar_stick(:,3).',(Jre(1,:)),(Jre(2,:)),(Jre(3,:)),Jrenorm,scal);
axis equal
zlim auto
xlabel('x')
ylabel('y')
zlabel('z')
title('J vec \Re')
colorbar;
caxis([cmin cmax]);
view(3)
xlim([xmin xmax])
ylim([ymin ymax])
zlim([zmin zmax])
%% J vec im
cmax=max(max(Jimnorm));
cmin=min(min(Jimnorm));
figure
hold on
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color',[0.65 0.65 0.65])
quiver3_c_scal(bar_stick(:,1).',bar_stick(:,2).',bar_stick(:,3).',(Jim(1,:)),(Jim(2,:)),(Jim(3,:)),Jimnorm,scal);
axis equal
zlim auto
xlabel('x')
ylabel('y')
zlabel('z')
title('J vec imag')
colorbar;
caxis([cmin cmax]);
view(3)
xlim([xmin xmax])
ylim([ymin ymax])
zlim([zmin zmax])
end

