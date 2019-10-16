function [N,Matrix_P0_c,Matrix_P0_m,Matrix_P0_e,w,rho_c,Area_c,rho_m,...
    Fac_Ed_loc_m,Aed_m,Led_m,Area_m,ind_edge_free_mag,Area_e,...
    ind_edge_free_con,bar_e,Matrix_P0,Matrix_Cb_m,phi_m_n0,ex,F1_c,F1_m,C1_m,...
    F1_e,G1_m,Matrix_G_m,Matrix_C_m] = ...
    fun_pre_processing(f,rrho_con,mmu_r,mu_0,test_case_dir)
%%
if ~strcmp(test_case_dir.con,'');cd test_cases; cd('cond'); cd(test_case_dir.con); load data_con.mat; ex.vol_con=1; cd ..; cd ..; cd ..; else ex.vol_con=0; end
if ~strcmp(test_case_dir.mag,'');cd test_cases; cd('mag'); cd(test_case_dir.mag); load data_mag.mat;  ex.vol_mag=1; cd ..; cd ..; cd ..; else ex.vol_mag=0; end
if ~strcmp(test_case_dir.ext,'');cd test_cases; cd('ext'); cd(test_case_dir.ext); load data_ext.mat;  ex.vol_ext=1; cd ..; cd ..; cd ..; else ex.vol_ext=0; end
%%
w = 2*pi*f; % omega
%% tol_x
tol_x=1e-4; % to avoid that some edges touch the z axis 
%% convenctions 
loc_edge(1 ,1:2)=[1 2];
loc_edge(2 ,1:2)=[2 3];
loc_edge(3 ,1:2)=[3 4];
loc_edge(4 ,1:2)=[4 1];
loc_edge(5 ,1:2)=[5 6];
loc_edge(6 ,1:2)=[6 7];
loc_edge(7 ,1:2)=[7 8];
loc_edge(8 ,1:2)=[8 5];
loc_edge(9 ,1:2)=[1 5];
loc_edge(10,1:2)=[2 6];
loc_edge(11,1:2)=[3 7];
loc_edge(12,1:2)=[4 8];
%% COND
if ex.vol_con==false
    F1_c=zeros(4,0);
    G2_c=zeros(8,0);
    Matrix_P0_c=zeros(0,3);
    Matrix_C_c=zeros(0,0);
end
N.face_con=size(F1_c,2);
N.edge_con=size(G2_c,2);
N.node_con=size(Matrix_P0_c,1);
Area_c=zeros(N.face_con,1);
vec_c=zeros(N.face_con,3);
bar_c=zeros(3,N.face_con);
for ii = 1:N.face_con
   vec_c(ii,1:3)=0.5*cross(Matrix_P0_c(F1_c(2,ii),:)-Matrix_P0_c(F1_c(1,ii),:),...
                       Matrix_P0_c(F1_c(3,ii),:)-Matrix_P0_c(F1_c(1,ii),:))+...
                 0.5*cross(Matrix_P0_c(F1_c(4,ii),:)-Matrix_P0_c(F1_c(3,ii),:),...
                       Matrix_P0_c(F1_c(1,ii),:)-Matrix_P0_c(F1_c(4,ii),:));
   Area_c(ii,1)=norm(vec_c(ii,1:3));
   bar_c(1:3,ii)=(sum(Matrix_P0_c(F1_c(1:4,ii),:))/4).';
end
ind_edge_free_con=find(G2_c(8,:)==0);
ind_edge_shar_con=find(G2_c(8,:)~=0);
N.edge_free_con=length(ind_edge_free_con);
N.edge_shar_con=length(ind_edge_shar_con);
indn_c=find(Matrix_P0_c(:,1)<tol_x);
Matrix_P0_c(indn_c,1)=Matrix_P0_c(indn_c,1)+tol_x;
%% MAG
rrho_mag=1/(1j*w(1)*mu_0*(mmu_r-1)); % magnetic resistivity of magnetic volume media
if ex.vol_mag==false
    Matrix_G_m=zeros(0,0);
    F1_m=zeros(4,0);
    G2_m=zeros(8,0); 
    Matrix_P0_m=zeros(0,3);
    Matrix_C_m=zeros(0,0);
    G1_m=zeros(2,0);
    C1_m=zeros(4,0);
end
N.face_mag=size(F1_m,2);
N.edge_mag=size(G2_m,2);
N.node_mag=size(Matrix_P0_m,1);
Area_m=zeros(N.face_mag,1);
vec_m=zeros(N.face_mag,3);
bar_m=zeros(3,N.face_mag);
for ii = 1:N.face_mag
   vec_m(ii,1:3)=0.5*cross(Matrix_P0_m(F1_m(2,ii),:)-Matrix_P0_m(F1_m(1,ii),:),...
                       Matrix_P0_m(F1_m(3,ii),:)-Matrix_P0_m(F1_m(1,ii),:))+...
                 0.5*cross(Matrix_P0_m(F1_m(4,ii),:)-Matrix_P0_m(F1_m(3,ii),:),...
                       Matrix_P0_m(F1_m(1,ii),:)-Matrix_P0_m(F1_m(4,ii),:));
   Area_m(ii,1)=norm(vec_m(ii,1:3));
   bar_m(1:3,ii)=(sum(Matrix_P0_m(F1_m(1:4,ii),:))/4).';
end
ind_edge_free_mag=find(G2_m(8,:)==0);
ind_edge_shar_mag=find(G2_m(8,:)~=0);
N.edge_free_mag=length(ind_edge_free_mag);
N.edge_shar_mag=length(ind_edge_shar_mag);
N.node_mag=size(Matrix_P0_m,1);
Fac_Ed_loc_m=zeros(4,N.face_mag); % local index of the edge
for ii = 1:N.face_mag
   for jj=1:4
       ed=sort(G1_m(1:2,abs(C1_m(jj,ii)))); % 
       for hh = 1:4
           loc_ed=sort(F1_m(loc_edge(hh,1:2),ii));
           if loc_ed==ed
               Fac_Ed_loc_m(jj,ii)=hh;
           end
       end
   end
end
indn_m=find(Matrix_P0_m(:,1)<tol_x);
Matrix_P0_m(indn_m,1)=Matrix_P0_m(indn_m,1)+tol_x;
% AREA
N_GAUSS=3;
[xi,WG_loc] = lgwt(N_GAUSS,-1, 1); % 
Aed_m=zeros(N.edge_mag,1);
PP=zeros(N_GAUSS,3);
Led_m=zeros(N.edge_mag,1);
for kk = 1:N.edge_mag
    Ed1=Matrix_P0_m(G1_m(1:2,kk),1:3);
    Cent1=0.5*(Ed1(1,1:3)+Ed1(2,1:3));
    vec1=(Ed1(2,1:3)-Ed1(1,1:3))*0.5;
    lung=norm(vec1)*2.0;  
    for qq = 1:N_GAUSS
        [PP(qq,1:3)]=Cent1+vec1*xi(qq);
    end    
    for qq = 1:N_GAUSS
        Aed_m(kk)=Aed_m(kk)+2*pi*(PP(qq,1))*WG_loc(qq)*lung*0.5;
    end    
    Led_m(kk)=lung;
end
%%
if ex.vol_ext==false
    F1_e=zeros(4,0);
    G2_e=zeros(8,0);
    Matrix_P0_e=zeros(0,3);
    Matrix_C_e=zeros(0,0);
end
N.face_ext=size(F1_e,2);
N.edge_ext=size(G2_e,2);
N.node_ext=size(Matrix_P0_e,1);
Area_e=zeros(N.face_ext,1);
vec_e=zeros(N.face_ext,3);
bar_e=zeros(3,N.face_ext);
for ii = 1:N.face_ext
   vec_e(ii,1:3)=0.5*cross(Matrix_P0_e(F1_e(2,ii),:)-Matrix_P0_e(F1_e(1,ii),:),...
                       Matrix_P0_e(F1_e(3,ii),:)-Matrix_P0_e(F1_e(1,ii),:))+...
                 0.5*cross(Matrix_P0_e(F1_e(4,ii),:)-Matrix_P0_e(F1_e(3,ii),:),...
                       Matrix_P0_e(F1_e(1,ii),:)-Matrix_P0_e(F1_e(4,ii),:));
   Area_e(ii,1)=norm(vec_e(ii,1:3));
   bar_e(1:3,ii)=(sum(Matrix_P0_e(F1_e(1:4,ii),:))/4).';
end
ind_edge_free_ext=find(G2_e(8,:)==0);
ind_edge_shar_ext=find(G2_e(8,:)~=0);
N.edge_free_ext=length(ind_edge_free_ext);
N.edge_shar_ext=length(ind_edge_shar_ext);
indn_e=find(Matrix_P0_e(:,1)<tol_x);
Matrix_P0_e(indn_e,1)=Matrix_P0_e(indn_e,1)+tol_x;
%% disp
disp('...SIMULATION DATA... ')
disp(['...frequency = ' num2str(f) ' [Hz]'])
disp(['...N_face_con = ',num2str(N.face_con),'...'])
disp(['...N_edge_con = ',num2str(N.edge_con),'...'])
disp(['...N_node_con = ',num2str(N.node_con),'...'])
disp(['...N_face_mag = ',num2str(N.face_mag),'...'])
disp(['...N_edge_mag = ',num2str(N.edge_mag),'...'])
disp(['...N_node_mag = ',num2str(N.node_mag),'...'])
%disp('-------------------------------------------------------------------')
%% set material vector
%disp('-------------------------------------------------------------------')
disp('...Set material...')
%% COND
rho_c=zeros(N.face_con,1);
for ii = 1:N.face_con
    rho_c(ii)=rrho_con;
end
%% MAG 
rho_m=zeros(N.face_mag,1);
mu_r_vec=zeros(N.face_mag,1);
for ii = 1:N.face_mag
    rho_m(ii)=rrho_mag;
    mu_r_vec(ii)=mmu_r;
end
%% Boundary plot
%disp('-------------------------------------------------------------------')
disp('...geo-boundary plot...')
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
title('Mesh')
xlabel('r')
ylabel('\phi')
zlabel('z')
view(0,0)
Matrix_P0=[Matrix_P0_c;Matrix_P0_m;Matrix_P0_e];
axis([min(Matrix_P0(:,1))*0 max(Matrix_P0(:,1)) min(Matrix_P0(:,2))-1 max(Matrix_P0(:,2))+1 min(Matrix_P0(:,3)) max(Matrix_P0(:,3))])
toc
pause(0.5)
disp('... done!')
%disp('-------------------------------------------------------------------')
%% Computing matrices incidances
%disp('-------------------------------------------------------------------')
disp('...computing matrices GCD...')
%% Incidance matrices Conductive vol
%disp('-------------------------------------------------------------------')
disp('...conductor...')
tic
Matrix_Cb_c=zeros(N.edge_free_con,N.edge_con);
[~,ind]=intersect(1:N.edge_con,ind_edge_free_con);
Matrix_Cb_c(1:N.edge_free_con,ind.')=eye(N.edge_free_con);
for ii = 1:N.edge_free_con 
   [~,~,v]=find(Matrix_C_c(:,ind(ii)));
   if v == -1
       Matrix_Cb_c(:,ind(ii))=Matrix_Cb_c(:,ind(ii));
   elseif v == 1
       Matrix_Cb_c(:,ind(ii))=-Matrix_Cb_c(:,ind(ii));
   end
end
Matrix_Cb_c=sparse(Matrix_Cb_c);
Matrix_Ca_c=[Matrix_C_c;Matrix_Cb_c];
toc
%disp('-------------------------------------------------------------------')
%% Incidance matrices Magnetic vol
%disp('-------------------------------------------------------------------')
disp('...magnetic...')
tic
Matrix_Cb_m=zeros(N.edge_free_mag,N.edge_mag);
[~,ind]=intersect(1:N.edge_mag,ind_edge_free_mag);
Matrix_Cb_m(1:N.edge_free_mag,ind.')=eye(N.edge_free_mag);
for ii = 1:N.edge_free_mag 
   [~,~,v]=find(Matrix_C_m(:,ind(ii)));
   if v == -1
       Matrix_Cb_m(:,ind(ii))=Matrix_Cb_m(:,ind(ii));
   elseif v == 1
       Matrix_Cb_m(:,ind(ii))=-Matrix_Cb_m(:,ind(ii));
   end
end
Matrix_Cb_m=sparse(Matrix_Cb_m);
Matrix_Ca_m=[Matrix_C_m;Matrix_Cb_m];
%
toc
%disp('-------------------------------------------------------------------')
%%
%disp('-------------------------------------------------------------------')
disp('...Change of variables (jc,m)-->(jc,phim)...')
CON_MAT=sparse(G1_m(1,:),G1_m(2,:),ones(N.edge_mag,1),N.node_mag,N.node_mag);
CON_MAT=CON_MAT+CON_MAT.';
[bins_m, rts] = graph_connected_components(logical(CON_MAT));
% GRAPH_m = graph(G1_m(1,:),G1_m(2,:));
% bins_m = conncomp(GRAPH_m);
N.dom_m=length(unique(bins_m));
phi_m0=zeros(N.dom_m,1);
for ii = 1:N.dom_m
   ind=find(bins_m==ii);
   phi_m0(ii)=ind(1);
end
phi_m_n0=setdiff(1:N.node_mag,phi_m0); 
%disp('-------------------------------------------------------------------')
end

