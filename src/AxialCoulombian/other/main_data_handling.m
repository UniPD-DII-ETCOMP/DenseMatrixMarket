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
if ex.vol_mag==false
    Matrix_G_m=zeros(0,0);
    F1_m=zeros(4,0);
    G2_m=zeros(8,0); 
    Matrix_P0_m=zeros(0,3);
    Matrix_C_m=zeros(0,0);
    G1_m=zeros(2,0);
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