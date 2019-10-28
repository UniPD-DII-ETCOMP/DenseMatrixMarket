function [J_r,J_i,J_norm_r,J_norm_i,bar] = fun_post_J(N_face,N_volu,N_volu_cond,Matrix_P0,VP,D1,X_I,ind_face_cond,ind_volu_cond)
XX_I=zeros(N_face,1);
XX_I(ind_face_cond,1)=X_I;
J_r = zeros(N_volu,3);
J_norm_r = zeros(N_volu,1);
J_i = zeros(N_volu,3);
bar=zeros(N_volu,3);
J_norm_i = zeros(N_volu,1);
for ii = 1:N_volu_cond
ind_volu=ind_volu_cond(ii);
pp = Matrix_P0(VP(:,ind_volu),:);
C = sum(pp)/4;
bar(ind_volu,:)=C;
for jj = 1:4
indf = abs(D1(jj,ind_volu));
[w] = fun_shape_tetra(pp,C,jj);
I_r = sign(D1(jj,ind_volu))*real(XX_I(indf));
J_r(ind_volu,:) = J_r(ind_volu,:)+w*I_r;
J_norm_r(ind_volu,1) = fun_my_norm(J_r(ind_volu,:));
I_i = sign(D1(jj,ind_volu))*imag(XX_I(indf));
J_i(ind_volu,:) = J_i(ind_volu,:)+w*I_i;
J_norm_i(ind_volu,1) = fun_my_norm(J_i(ind_volu,:));
end
end
end

function [w] = fun_shape_tetra(tetra,P,v)
J =   fun_my_vol_tetra(tetra)*6;%det([1 1 1 1;tetra']); %Jacobian
NP = size(P,1);
w = zeros(NP,3);
for ii = 1:NP
w(ii,:) = 2*(P(ii,:)-tetra(v,:))/J;
end
end

function [V] = fun_my_vol_tetra(tetra)
v1 =  tetra(2,:) - tetra(1,:);
v2 =  tetra(3,:) - tetra(1,:);
v3 =  tetra(4,:) - tetra(1,:);
V = det([v1;v2;v3])/6;
end

function [B] = fun_my_norm(A)
B = sqrt(A(:,1).^2 +A(:,2).^2 + A(:,3).^2);
end
