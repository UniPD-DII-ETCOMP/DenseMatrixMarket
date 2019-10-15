function [J_r,J_i,J_norm_r,J_norm_i,bar] = fun_post_J_Quad_face2(N_edge,N_face,Matrix_P0,F1,C1,X_I,Fac_Ed_loc,Aed,Led)
J_r = zeros(N_face,3);
J_norm_r = zeros(N_face,1);
J_i = zeros(N_face,3);
bar=zeros(N_face,3);
J_norm_i = zeros(N_face,1);
for ii = 1:N_face % 
    Face=Matrix_P0(F1(1:4,ii),:);
    Hexa(1:4,[1,3])=Face(1:4,[1,3]);
    Hexa(5:8,[1,3])=Face(1:4,[1,3]);
    Hexa(1:4,2)=-1;
    Hexa(5:8,2)=1;
C = sum(Hexa)/8; % 
bar(ii,:)=C;
[ww] = 2*fun_whitney_face_hexa_6_quad(Hexa,[0,-1,0].',1); % 
for jj = 1:4 % 
indf = abs(C1(jj,ii)); % 
ind_loc_face =  Fac_Ed_loc(jj,ii); % 
w=ww(1:3,1,ind_loc_face)*Led(indf)/Aed(indf); % 
I_r = sign(C1(jj,ii))*real(X_I(indf));
J_r(ii,:) = J_r(ii,:)+w.'*I_r;
J_norm_r(ii,1) = norm(J_r(ii,:));
I_i = sign(C1(jj,ii))*imag(X_I(indf));
J_i(ii,:) = J_i(ii,:)+w.'*I_i;
J_norm_i(ii,1) = norm(J_i(ii,:));
end
end
end
