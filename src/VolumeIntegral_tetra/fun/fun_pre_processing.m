function [rho_c,N,ind,EV,VE,EV_ind,VE_ind,Q,brhs,norm_eps,w,bared,...
    xmin,xmax,ymin,ymax,zmin,zmax] = ...
    fun_pre_processing(F1,D1,C1,G1,VP,Matrix_D,...
    Matrix_C,Matrix_G,Matrix_P0,E_ext,rrho_c,N,f)
%%
xmin=min(Matrix_P0(:,1));
xmax=max(Matrix_P0(:,1));
ymin=min(Matrix_P0(:,2));
ymax=max(Matrix_P0(:,2));
zmin=min(Matrix_P0(:,3));
zmax=max(Matrix_P0(:,3));
%% arranging data & conta degli elementi
disp('-------------------------------------------------------------------')
disp('ARRANGING DATA ...')
tic
N.volu = size(Matrix_D,1); % numero di volumi (tetraedri)
N.face = size(Matrix_D,2); % numero di facce 
N.edge = size(Matrix_G,1); % numero di lati (geometrici)
N.node = size(Matrix_G,2); % numero di nodi (geometrici)
N.face_free = size(F1(5,(F1(5,:)==0)),2);
N.face_shar = N.face - N.face_free;
ind.face_free=find(F1(5,:)==0); 
ind.face_shar=find(F1(5,:)~=0); 
ind.edge=unique(abs(reshape(C1(:,:),1,3*N.face))); % lati dei conduttori
N.edge=size(ind.edge,2);                                       % numero di lati dei conduttori
ind.edge_ext=unique(abs(reshape(C1(:,ind.face_free),1,3*N.face_free))); % lati esterni dei conduttori
N.edge_ext=size(ind.edge_ext,2);                                             % numero di lati esterni dei conduttori
ind.edge_int=setdiff(ind.edge,ind.edge_ext);              % lati interni conduttori
N.edge_int=size(ind.edge_ext,2);                               % numero di lati interni dei conduttori
% node
ind.node=unique(abs(reshape(G1(:,ind.edge),1,2*N.edge))); % nodi conduttori
N.node=size(ind.node,2);                                       % numero di nodi conduttori
ind.node_ext=unique(abs(reshape(G1(:,ind.edge_ext),1,2*N.edge_ext))); % nodi esterni conduttori
N.node_ext=size(ind.node_ext,2);                                           % numero di nodi esterni conduttori
ind.node_int=setdiff(ind.node,ind.node_ext);                          % nodi interni conduttori
N.node_int=size(ind.node_int,2);                                           % numero di nodi interni conduttori
% VP_fv
VP_fv = zeros(4,N.volu);   % matrice che comunica con D1; D1(ii,jj)=faccia ii del tetraedro jj -- VP_fv(ii,jj) verice opposto alla faccia ii del tetraedro jj  
for ii = 1:N.volu 
    a = (D1(:,ii));
   for jj = 1:4
       b = abs(a(jj));
       fp = F1(1:3,b); 
       VP_fv(jj,ii) = setdiff(VP(:,ii),fp);
   end
end
% EV e EV_ind
Matrix_DC=abs(Matrix_D)*abs(Matrix_C); % volumexedge
N.n_EV=full(sum(Matrix_DC,1))/2; % per ogni lato dico quanti tetraedri ho attaccati
N.n_EV_max=max(N.n_EV); % massimo numero di tetraedri attaccati ad un lato
N.nnzR=sum(N.n_EV*6);
EV=zeros(N.n_EV_max,N.edge); % matrice che mi dà i tetraedri (con relativo segno) attaccati al lato 
EV_ind=zeros(N.n_EV_max,N.edge); % matrice che mi dà l'indice locale(1,2,3,4,5,6) del lato del tetraedro che coincide con l'indice della colonna della matrice (lato considerarto)
loc_edg=[1,2;2,3;3,1;4,1;4,2;4,3];
for ii = 1:N.edge
   EV(1:N.n_EV(ii),ii)= find(Matrix_DC(:,ii));
   ed=G1(:,ii)'; % nodi (orientati) che compongono il lato
   for jj=1:N.n_EV(ii)
      p_V=VP(:,EV(jj,ii)); 
      ed_V=p_V(loc_edg);
      check=(ed_V(:,1) == ed(1) & ed_V(:,2) == ed(2)); % check=sum(ed_V(:, 1) == ed(1) & ed_V(:, 2) == ed(2));%sum((ed_V(:, 1) == ed(1)).* (ed_V(:, 2) ==ed(2)));%ismember(ed_V,ed,'rows');
      check_s=sum(check);
      if check_s==0
          EV(jj,ii)=-EV(jj,ii);
          check2=(ed_V(:, 1) == ed(2) & ed_V(:, 2) == ed(1));
          EV_ind(jj,ii)=find(check2==1);
      elseif check_s==1
          EV_ind(jj,ii)=find(check==1);
      end
   end
end
% VE e VE_ind
VE=zeros(6,N.volu);
VE_ind=zeros(6,N.volu);
for ii = 1:N.volu
    VE(1:6,ii)=find(Matrix_DC(ii,:));
    p_V=VP(:,ii); % punti del volume
    ed_V=p_V(loc_edg); % punti dei lati del volume
    for jj = 1:6
        ed=G1(:,VE(jj,ii))';
        check=(ed_V(:,1) == ed(1) & ed_V(:, 2) == ed(2));
        if sum(check)==0
           VE(jj,ii)=-VE(jj,ii);
           check2=(ed_V(:, 1) == ed(2) & ed_V(:, 2) == ed(1));
           VE_ind(jj,ii)=find(check2==1);
        elseif sum(check)==1
            VE_ind(jj,ii)=find(check==1);
        end
    end
end
toc
disp('... done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp('plot geo...')
fun_plot_geo2(F1,ind,Matrix_P0,xmin,xmax,ymin,ymax,zmin,zmax);
drawnow
disp('... done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp('tree cotree...')
[N,Q,ind] = fun_tree_cotree(ind,G1,Matrix_C,N);
disp('... done!')
disp('-------------------------------------------------------------------')
%%
disp('-------------------------------------------------------------------')
disp(['N.cotree (sys size) =',num2str(N.cotree)])
disp('-------------------------------------------------------------------')
%% Rhs-term External fields volu
disp('-------------------------------------------------------------------')
disp('rhs: Extrenal Electric Field...')
tic
[e_face] = fun_compute_ext_field_dual_edge_dot_tetra(Matrix_P0.',N.face,...
          F1,ind.face_shar,N.face_shar,ind.face_free,N.face_free,E_ext,7);
brhs=Q.'*e_face;
toc
disp('... done!')
disp('-------------------------------------------------------------------')
%%
rho_c=rrho_c*ones(N.volu,1);
w=2*pi*f;
norm_eps=0;
%%
bared=((Matrix_P0(G1(1,:),:)+Matrix_P0(G1(2,:),:))*0.5).';
end