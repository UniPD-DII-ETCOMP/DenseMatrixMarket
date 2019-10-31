function [G1,C1,D1,F1,Matrix_C,Matrix_G,Matrix_D] = fun_create_data(VP)
[G1,C1,D1,F1]=gcd_matlab(VP);
%% number of nodes edges faces volumes
nn = max(max(VP));
ne = size(G1,2);
nf = size(C1,2);
nv = size(D1,2);
%% G edges xnodes
Matrix_G = sparse(1:ne,G1(1,1:ne),-ones(ne,1),ne,nn);
Matrix_G = Matrix_G+sparse(1:ne,G1(2,1:ne),ones(ne,1),ne,nn);
%% C faces x edges
Cpos = (C1+abs(C1))/2;
Cneg = (C1-abs(C1))/2;
[r_pos_c,c_pos_c,val_pos_c] = find(Cpos);
[r_neg_c,c_neg_c,val_neg_c] = find(Cneg);
Matrix_C = sparse(c_pos_c,val_pos_c,ones(size(c_pos_c,1),1),nf,ne);
Matrix_C = Matrix_C+sparse(c_neg_c,abs(val_neg_c),-ones(size(c_neg_c,1),1),nf,ne);
%  Cnew = full(Cnew);
%% check 
disp('check if curl(grad()) = 0')
CxG=Matrix_C*Matrix_G;
ck=size(find(CxG),1);
if ck==0; disp('ok!'); else disp('no :('); end %disp(num2str(size(find(CxG),1)));
%% D volumes x faces
Dpos = (D1+abs(D1))/2;
Dneg = (D1-abs(D1))/2;
[r_pos_d,c_pos_d,val_pos_d] = find(Dpos);
[r_neg_d,c_neg_d,val_neg_d] = find(Dneg);
Matrix_D = sparse(c_pos_d,val_pos_d,ones(size(c_pos_d,1),1),nv,nf);
Matrix_D = Matrix_D+sparse(c_neg_d,abs(val_neg_d),-ones(size(c_neg_d,1),1),nv,nf);
%  Dnew = full(Dnew);
%% check 
disp('check if div(curl()) = 0')
DxC=Matrix_D*Matrix_C;
ck=size(find(DxC),1);
if ck==0; disp('ok!'); else disp('no :('); end %disp(num2str(size(find(DxC),1)));
%%
%% check
end
 
 