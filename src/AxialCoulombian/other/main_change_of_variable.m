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
Matrix_Q=blkdiag(speye(N.face_con),Matrix_G_m(:,phi_m_n0));
SYSTEM2= Matrix_Q.'*SYSTEM*Matrix_Q;
brhs2=Matrix_Q.'*brhs;
