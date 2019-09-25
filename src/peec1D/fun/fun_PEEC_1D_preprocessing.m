function [N,nSticks,PPgh,ll_h,ut_h,gauss_W,ne_cap_max,...
          Cap_Elem,NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,R,omega,Aobj,...
          Geapp,A_eapp,R_app,L_app,Cinv_app,A_napp,...
          beta,radius,b,Is,nNodes,lgt,radius_cap,bar_stick,barDoFs] ...
          = fun_PEEC_1D_preprocessing(Voltage_Source,...
          Current_Source,NN,G,radius,freq,np_gauss,rho,E_external,ret,E_ext)
%% Physical quantities
c_light = 299792458; %[m/s]
omega = 2*pi*freq; % angular frequency [Hz]
beta = omega/299792458; % propagation costant [1/m]
lambda = c_light/freq; % wave length [m]
%% Number of elements 
N.vol_sou=size([Voltage_Source.value],2); % number of appended barches
N.cur_sou=size([Current_Source.node],2); % number of injected currents
app=[[Voltage_Source.node_start],[Voltage_Source.node_end],[Current_Source.node]];
app=unique(app);
app=find(app<0);
N.node_app=length(app);% number of appended nodes
%  imposed currents are not supported in the matrix free code
%% plot geometry
disp('plotting geo...')
[~]=fun_plot_geo(NN,G); drawnow
disp('...done!')
%% arrangig data
disp('------------------------------------')
disp('Pre-processing of data...')
tic 
nSticks = size(G,2); % number of stick elements
nNodes = size(NN,2); % number of nodes
val_neg_g=G(1,:);
val_pos_g=G(2,:);
Matrix_G=sparse([1:nSticks],val_neg_g,-ones(nSticks,1),nSticks,nNodes);
Matrix_G=Matrix_G+sparse([1:nSticks],val_pos_g,ones(nSticks,1),nSticks,nNodes);
Aobj = Matrix_G'; %incidence matrix equivalent circuit
radius = radius*ones(1,nSticks);
% build capacitive elements
disp('...capacitive elements...')
ne_cap_max = max(full(sum(abs(Matrix_G),1)));
Cap_Elem = zeros(1+ne_cap_max,nNodes);
radius_cap=zeros(nNodes,1);
for k=1:nNodes
    idE = find(Matrix_G(:,k)); %edges connected to node k
    nE = length(idE); % number of edges connected to node k
    Cap_Elem(1,k) = nE; 
    Cap_Elem(2:nE+1,k) = idE;
    radius_cap(k)=mean(radius(idE)); % take a mean value for the radius related to the capacitive element
end
ne_cap_tot = sum(Cap_Elem(1,:));
toc
disp('...done!')
disp('------------------------------------')
%%
disp('------------------------------------')
disp('DETAILS:')
disp([num2str(freq) ' [Hz]'])
disp([num2str(nSticks) ' stick elements'])
disp([num2str(N.vol_sou) ' appended branches'])
disp([num2str(N.cur_sou) ' injected currents'])
disp(['N.B.: imposed currents are not supported'])
disp('------------------------------------')
%% pre-processing
disp('------------------------------------')
disp('Pre-processing...')
% L
disp('...pre-L...')
[gauss_P,gauss_W]=lgwt(np_gauss,-1,1); % gauss points and weights (local)
PPgh=zeros(3,np_gauss,nSticks);
ll_h=zeros(nSticks,1);
ut_h=zeros(3,nSticks);
for ii = 1:nSticks
    PPh=NN(1:3,G(1:2,ii)); 
    [PPgh(:,:,ii),ll_h(ii)] = Gauss_line_nvar(PPh,gauss_P,np_gauss); % gauss points global
    ut_h(:,ii)=(PPh(1:3,2)-PPh(1:3,1))/ll_h(ii);
end
% P 
disp('...pre-P...')
NN_xx1=zeros(3,ne_cap_max,nNodes);
NN_xx2=zeros(3,ne_cap_max,nNodes);
PPg_xx=zeros(3,np_gauss,ne_cap_max,nNodes);
ll_xx=zeros(ne_cap_max,nNodes);
ll_tot_xx=zeros(nNodes,1);
for hh=1:nNodes 
    ne_cap_hh=Cap_Elem(1,hh); % number of elements 
	idE_hh=Cap_Elem(2:ne_cap_max+1,hh);	
	glob_P=0.0;
   for jj = 1:ne_cap_hh
      NN_edge=NN(1:3,G(1:2,idE_hh(jj))); % take end points
      NN_xx1(1:3,jj,hh) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2)); % build capacitive element
      NN_xx2(1:3,jj,hh) = NN(1:3,hh); 
      [PPg_xx(1:3,1:np_gauss,jj,hh),ll_xx(jj,hh)] = Gauss_line_nvar([NN_xx1(:,jj,hh),NN_xx2(:,jj,hh)],gauss_P,np_gauss); %gauss points
      ll_tot_xx(hh)=ll_tot_xx(hh)+ll_xx(jj,hh); % length update
   end
end
disp('...done!')
disp('------------------------------------')
%% compute R (diagonal)
disp('------------------------------------')
disp('computing R (diagonal matrix)...')
tic
rho = rho*ones(1,nSticks);
lgt = compute_length(NN,G);
R = lgt.*(rho./(pi*radius.^2));
R = sparse(1:nSticks,1:nSticks,R,nSticks,nSticks);
toc
disp('...done!')
%% appended elements 
disp('------------------------------------')
disp('adding appended elements...')
tic
if N.vol_sou == 0 
    disp('no appended elements')
end
Us = zeros(nSticks,1);
Is = zeros(nNodes,1);
Is_napp = zeros(N.node_app,1);
%------------------ voltage generators -----------------
[A_eapp,A_napp,R_app,L_app,Cinv_app,Us_app] = set_VoltageSource(nNodes,N.vol_sou,N.node_app,Voltage_Source); 
% e_app: pre processing
disp('...pre-e_app...')
Geapp=zeros(2,N.vol_sou);
for ii = 1:N.vol_sou
    if ~isempty(find(A_eapp(:,ii)==-1, 1))
        Geapp(1,ii)=find(A_eapp(:,ii)==-1);
    end
    if ~isempty(find(A_eapp(:,ii)==1, 1))
        Geapp(2,ii)=find(A_eapp(:,ii)==1);
    end
end
%------------------ current generators -----------------
for k=1:N.cur_sou %conductors
    %current generator connected to a OBJ node
    if sign(Current_Source(k).node) == 1
        Is(Current_Source(k).node)=Current_Source(k).value;
    %current generator connected to a APPENDED node
    elseif sign(Current_Source(k).node) == -1
        Is_napp(abs(Current_Source(k).node))=Current_Source(k).value;
    end
end
toc
disp('...done!')
disp('------------------------------------')
%% external E
disp('------------------------------------')
disp('setting external electric field E0...')
if E_external == 1
    disp('...on...')
    tic
    Us = fun_compute_ext_field_dual_edge_dot_triang(NN,nSticks,G,E_ext,7);
    toc
else
    disp('...on...')
end
toc
disp('...done!')
disp('------------------------------------') 
%% system
disp('assembling rhs...')
tic
indIs=find(abs(Is)>0);
U_Is=zeros(nNodes,1);
if ret == 0
    for ii = 1:N.cur_sou
    U_Is=U_Is+P_stick_assemble_matlab_vec_MF(NN,G,radius,...
               np_gauss,ne_cap_max,Cap_Elem,...
               NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,...
               gauss_W,1:nNodes,indIs(ii))*Is(indIs(ii));
    end
else
    for ii = 1:N.cur_sou
    U_Is=U_Is+(P_stick_assemble_matlab_vec_MF(NN,G,radius,...
               np_gauss,ne_cap_max,Cap_Elem,...
               NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,...
               gauss_W,1:nNodes,indIs(ii))...
               .*exp(-1j*beta*dist_node_node_MF(NN,1:nNodes,indIs(ii))))*Is(indIs(ii));
    end
end
b = [ Us-Aobj.'*U_Is/(1j*omega); ... (e_obj)
      Us_app-A_eapp.'*U_Is/(1j*omega); ... (e_app)
      -Is_napp]; %(n_app)
toc
disp('...done!')
disp('------------------------------------')
%% barycenter related to DoFs
bar_stick = 0.5*(NN(:,G(1,:))+NN(:,G(2,:))).'; % barycenter of the sticks
barDoFs=999*ones(nSticks+N.vol_sou+N.node_app,3); % N.B. appended elements and nodes do not have a real geometric collocation, thus their barycenter is set to [999,999,999]
barDoFs(1:nSticks,1:3)=bar_stick;
end

