if 0
[Matrix_Rm_im] = funRrz_for(N.node_mag,N.face_mag,imag(rho_m).',N_GAUSS,Matrix_P0_m,F1_m,N.edge_mag,C1_m,Fac_Ed_loc_m,Aed_m.',Led_m.',N.thread);
[R] = funRrz_face2(N.face_mag,imag(rho_m).',N_GAUSS,Matrix_P0_m,F1_m,N.edge_mag,C1_m,Fac_Ed_loc_m,Aed_m,Led_m);
end
%%
if 0
Nf=2000;
tic
[Matrix_Lcet] = mu_0*funLphiphi3_for(N.node_con+N.node_ext,...
                                    Nf,...
                                   [Area_c(1:Nf)].',...
                                    N_GAUSS1,...
                                    N_GAUSS2,...
                                   [Matrix_P0_c;Matrix_P0_e],...
                                   [F1_c(:,1:Nf),N.node_con+F1_e],...
                                    N.thread);
toc1=toc
tic
[L] = mu_0*funLphiphi3_acc(Nf,N_GAUSS1,N_GAUSS2,Matrix_P0_c,F1_c(:,1:Nf));                                
toc2=toc
toc2/toc1
end
%%
if 0
    tic
    [Pm_for] = funPphiphi3_ss_for(N.node_mag,... %funPphiphi3_ss (matlab)
                                    N.edge_free_mag,...
                                    N_GAUSS1,...
                                    N_GAUSS2,...
                                    N.thread,...
                                    Matrix_P0_m,...
                                    G1_m(:,ind_edge_free_mag),...
                                    Aed_m(ind_edge_free_mag).');
    toc1=toc    
    tic
    [Pm_mat] = funPphiphi3_ss_acc(N.node_mag,... %funPphiphi3_ss (matlab)
                                N.edge_free_mag,...
                                N_GAUSS1,...
                                N_GAUSS2,...
                                N.thread,...
                                Matrix_P0_m,...
                                G1_m(:,ind_edge_free_mag),...
                                Aed_m(ind_edge_free_mag).');
    toc2=toc       
end   
%%
Nf=2000;
tic
Nme_ce_for  = funNme_for(N.node_con+N.node_ext,... % funNme (matlab)
                         Nf,...
                         [Area_c(1:Nf);Area_e].',...
                         N_GAUSS,...
                         N.node_mag,...
                         [Matrix_P0_c;Matrix_P0_e],...
                         [F1_c(:,1:Nf),N.node_con+F1_e],...
                         N.thread,...
                         Matrix_P0_m);              
toc1=toc 
tic
[Nme_ce_mat]= funNme_acc(Nf,...
                    N_GAUSS,...
                    N.node_mag,...
                     [Area_c(1:Nf);Area_e].',...
                     [Matrix_P0_c;Matrix_P0_e],...
                     [F1_c(:,1:Nf),N.node_con+F1_e],...
                     Matrix_P0_m); 
toc2=toc
%%
tic
if 0 
figure
hold on
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'Facecolor',green,'FaceAlpha',alpha) 
for ii = 1:N.edge_mag
    p1=Matrix_P0_m(G1_m(1,ii),:);
    p2=Matrix_P0_m(G1_m(2,ii),:);
    bar=(p1+p2)*0.5;
    text(bar(1),bar(2),bar(3),num2str(ii))
end
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
Matrix_P0=[Matrix_P0_c;Matrix_P0_m];
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
toc
end
%%
if 0
figure
hold on
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'Facecolor',green,'FaceAlpha',alpha) 
for ii = 1:N.edge_free_mag
    p1=Matrix_P0_m(G1_m(1,ind_edge_free_mag(ii)),:);
    p2=Matrix_P0_m(G1_m(2,ind_edge_free_mag(ii)),:);
    bar=(p1+p2)*0.5;
    text(bar(1),bar(2),bar(3),num2str(ii))
end
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
Matrix_P0=[Matrix_P0_c;Matrix_P0_m];
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
toc
end
%%
if 0
figure
hold on
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'Facecolor',green,'FaceAlpha',alpha) 
for ii = 1:N.face_mag
    text(bar_m(1,ii),bar_m(2,ii),bar_m(3,ii),num2str(ii))
end
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
Matrix_P0=[Matrix_P0_c;Matrix_P0_m];
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
toc
end
%%
if 0
figure
hold on
% patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'Facecolor',green,'FaceAlpha',alpha) 
for ii = 1:N.edge_free_mag
    p1=Matrix_P0_m(G1_m(1,ind_edge_free_mag(ii)),:);
    p2=Matrix_P0_m(G1_m(2,ind_edge_free_mag(ii)),:);
    bar=(p1+p2)*0.5;
    vec=(p2-p1)*0.5;    
    text(bar(1),bar(2),bar(3),num2str(ii))
    quiver3(bar(1),bar(2),bar(3),vec(1),vec(2),vec(3))
end
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
Matrix_P0=[Matrix_P0_c;Matrix_P0_m];
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
end
%%
if 0
figure
hold on
% patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'Facecolor',green,'FaceAlpha',alpha) 
for ii = 1:N.edge_mag
    p1=Matrix_P0_m(G1_m(1,(ii)),:);
    p2=Matrix_P0_m(G1_m(2,(ii)),:);
    bar=(p1+p2)*0.5;
    vec=(p2-p1)*0.8;    
    text(bar(1),bar(2),bar(3),num2str(ii))
    quiver3(p1(1),p1(2),p1(3),vec(1),vec(2),vec(3))
end
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
Matrix_P0=[Matrix_P0_c;Matrix_P0_m];
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
toc
end
%%
if 0
[xabsc1,weig1] = lgwt(N_GAUSS1,-1, 1); % calcolo i punti e pesi di Gauss nella linea -[1 1]
NP1=N_GAUSS1; %numero di punti di Gauss nell'esaedro 
xi1=zeros(NP1,1);
% eta1=zeros(NP1,1);
% zet1=zeros(NP1,1);
WG_loc1=zeros(NP1,1);
hh=1;
for ii = 1:N_GAUSS1 % 
   xi1(hh) =xabsc1(ii);
   WG_loc1(hh)= weig1(ii); % peso di Gauss
   hh=hh+1;
end
sum(WG_loc1)
end
%%
if 0 
figure
hold on
patch('Faces',F1_m(1:4,:).','Vertices',Matrix_P0_m,'Facecolor',green,'FaceAlpha',alpha) 
for ii = 1:N.node_mag
    text(Matrix_P0_m(ii,1),Matrix_P0_m(ii,2),Matrix_P0_m(ii,3),num2str(ii))
end
axis equal
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
Matrix_P0=[Matrix_P0_c;Matrix_P0_m];
axis([min(Matrix_P0(:,1)) max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
toc
end
%%
if 0 
%%
En=0.5*(Area_m.'*Matrix_Lm_vv*Area_m)
I=sum(Area_m)
LL=En*2/(I^2)
%%
En=0.5*(Area_c.'*Matrix_Lcc*Area_c)
I=sum(Area_c)
LL=En*2/(I^2)
%%
En=0.5*(Area_c.'*Matrix_Lcc2*Area_c)
I=sum(Area_c)
LL=En*2/(I^2)
%%
cent=sum(bar_c,2)/N.face_con
rm=cent(1);
dr=2*(cent(1)-min(Matrix_P0_c(:,1)));
dz=2*(cent(3)-min(Matrix_P0_c(:,3)))
[L] = L_self_Axial_rect_cross(rm,dr,dz)
L=mu_0*L
%% Lcc
N_GAUSS1=3;
N_GAUSS2=4;
F11_c=[1,2,3,4].';
Matrix_P0=[4.5 0  0.5;...
           4.5 0 -0.5;...
           5.5 0 -0.5;...
           5.5 0  0.5];
%%
tic
[LLL] = mu_0*funLphiphi(1,1,N_GAUSS1,N_GAUSS2,Matrix_P0,F11_c)
toc
%%
tic
[LLL2] = mu_0*funLphiphi2(1,1,N_GAUSS1,N_GAUSS2,Matrix_P0,F11_c)
toc
%%
tic
[LLL3] = mu_0*funLphiphi3(1,1,N_GAUSS1,N_GAUSS2,Matrix_P0,F11_c)
toc
%%
cent=[5 0 0];
rm=cent(1);
dr=2*(cent(1)-min(Matrix_P0(:,1)));
dz=2*(cent(3)-min(Matrix_P0(:,3)))
[L] = L_self_Axial_rect_cross(rm,dr,dz);
L=mu_0*L
%%
LLL2/LLL
%%
tic
[LGa] = mu_0*funLphiphiGarrett(1,1,N_GAUSS1,N_GAUSS2,Matrix_P0,F11_c);
toc    
end
%%
if 0
    N_GAUSS1=2;
    N_GAUSS2=3;
tic
[Matrix_Lcm] = mu_0*funLphiphi3(N.face_con+N.face_mag,...
                                [Area_c;Area_m],...
                                N_GAUSS1,...
                                N_GAUSS2,...
                                [Matrix_P0_c;Matrix_P0_m],...
                                [F1_c,N.node_con+F1_m]);
toc
tic
[Matrix_Lcm_for] = mu_0*funLphiphi3_for(N.node_con+N.node_mag,...
                                        N.face_con+N.face_mag,...
                                        [Area_c;Area_m].',...
                                        N_GAUSS1,...
                                        N_GAUSS2,...
                                        [Matrix_P0_c;Matrix_P0_m],...
                                        [F1_c,N.node_con+F1_m],...
                                        N.thread);
toc
sum(sum(Matrix_Lcm))
sum(sum(Matrix_Lcm_for))
figure
hold on
plot(reshape(Matrix_Lcm,numel(Matrix_Lcm),1),'or')
plot(reshape(Matrix_Lcm_for,numel(Matrix_Lcm_for),1),'xb')
end
%%
if 0 
N_GAUSS1=4;
N_GAUSS2=5;
tic
[Matrix_L_ss] = mu_0*funLphiphi3_ss(N.edge_free_con,...
                                     N_GAUSS1,...
                                     N_GAUSS2,...
                                     Matrix_P0_c,...
                                     G1_c(:,ind_edge_free_con));
toc
tic
[L_ss_for] = mu_0*funLphiphi3_ss_for(N.node_con,...
                                     N.edge_free_con,...
                                     N_GAUSS1,...
                                     N_GAUSS2,...
                                     N.thread,...
                                     Matrix_P0_c,...
                                     G1_c(:,ind_edge_free_con));
toc
sum(sum(Matrix_Lcm))
sum(sum(Matrix_Lcm_for))
figure
hold on
plot(reshape(Matrix_L_ss,numel(Matrix_L_ss),1),'or')
plot(reshape(L_ss_for,numel(L_ss_for),1),'xb')
end
%%
if 0 
N_GAUSS1=2;
N_GAUSS2=3;
tic
[L_mat]         = mu_0*funLphiphi3_vs(N.edge_free_con,...
                                      N.face_con,...
                                      N_GAUSS1,...
                                      N_GAUSS2,...
                                      Matrix_P0_c,...
                                      G1_c(:,ind_edge_free_con),...
                                      Area_c,...
                                      F1_c);
toc
tic
[L_for]     = mu_0*funLphiphi3_vs_for(N.node_con,...
                                      N.edge_free_con,...
                                      N_GAUSS1,...
                                      N_GAUSS2,...
                                      N.thread,...
                                      N.face_con,...                                      
                                      Matrix_P0_c,...
                                      G1_c(:,ind_edge_free_con),...
                                      F1_c,...
                                      Area_c.');
toc
figure
hold on
plot(reshape(L_mat,numel(L_mat),1),'or')
plot(reshape(L_for,numel(L_for),1),'xb')
end
%%
if 0
[Matrix_Lcm] = mu_0*funLphiphi3_for(N.node_ext,...
                                    N.face_ext,...
                                   [Area_e].',...
                                    N_GAUSS1,...
                                    N_GAUSS2,...
                                   [Matrix_P0_e],...
                                   [F1_e],...
                                    N.thread);
end
%%  
if 0 
N_GAUSS1=4;
N_GAUSS2=5;
tic
[Matrix_Lcm] = mu_0*funLphiphi3_for(N.node_con+N.node_mag+N.node_ext,...
                                    N.face_con+N.face_mag+N.face_ext,...
                                   [Area_c;Area_m;Area_e].',...
                                    N_GAUSS1,...
                                    N_GAUSS2,...
                                   [Matrix_P0_c;Matrix_P0_m;Matrix_P0_e],...
                                   [F1_c,N.node_con+F1_m,N.node_con+N.node_mag+F1_e],...
                                    N.thread);
toc
sum(sum(Matrix_Lcm))
%%
N_GAUSS1=4;
N_GAUSS2=5;
[Matrix_Lm_ss] = mu_0*funLphiphi3_ss_for(N.node_mag,...
                                         N.edge_free_mag,...
                                         N_GAUSS1,...
                                         N_GAUSS2,...
                                         N.thread,...
                                         Matrix_P0_m,...
                                         G1_m(:,ind_edge_free_mag));
sum(sum(Matrix_Lm_ss))
%%
N_GAUSS1=4;
N_GAUSS2=5;
tic
[Matrix_Lmc_vs] = mu_0*funLphiphi3_vs_for(N.node_con+N.node_mag+N.node_ext,...
                                      N.edge_free_mag,...
                                      N_GAUSS1,...
                                      N_GAUSS2,...
                                      N.thread,...
                                      N.face_mag+N.face_con+N.face_ext,...                                      
                                      [Matrix_P0_m;Matrix_P0_c;Matrix_P0_e],...
                                      G1_m(:,ind_edge_free_mag),...
                                      [F1_m,N.node_mag+F1_c,N.node_con+N.node_mag+F1_e],...
                                      [Area_m;Area_c;Area_e].');

toc
sum(sum(Matrix_Lmc_vs))
end