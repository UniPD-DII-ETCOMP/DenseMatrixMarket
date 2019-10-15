function Vext=fun_compute_ext_field_phi(Area,Matrix_P0,nFaces,F1,Field_ext,N_GAUSS)
%%
[xabsc,weig] = lgwt(N_GAUSS,-1, 1); %
NP=N_GAUSS^2; %
xi=zeros(NP,1);
eta=zeros(NP,1);
zet=zeros(NP,1);
WG_loc=zeros(NP,1);
hh=1;
for ii = 1:N_GAUSS % 
    for jj = 1:N_GAUSS
           xi(hh) =xabsc(ii);
           zet(hh)=xabsc(jj);
           WG_loc(hh)= weig(ii)*weig(jj); %
           hh=hh+1;
    end 
end
%%
Vext=zeros(nFaces,1);
detJ=zeros(NP,1);
PP=zeros(NP,3);
%%
for ii=1:nFaces
%%
    Face=Matrix_P0(F1(1:4,ii),:);
    Hexa(1:4,[1,3])=Face(1:4,[1,3]);
    Hexa(5:8,[1,3])=Face(1:4,[1,3]);
    Hexa(1:4,2)=-1;
    Hexa(5:8,2)=1;
%% 
    for jj = 1:NP % 
        [PP(jj,:)]=funTrilinear(xi(jj),0.0,zet(jj),Hexa);
    end
%% 
    for jj = 1:NP
        J(1:3,1)=(-(1.0d0-eta(jj))*(1.0d0+zet(jj))*Hexa(1,1:3)-(1.0d0-eta(jj))*(1.0d0-zet(jj))*Hexa(2,1:3) ...
                  +(1.0d0-eta(jj))*(1.0d0-zet(jj))*Hexa(3,1:3)+(1.0d0-eta(jj))*(1.0d0+zet(jj))*Hexa(4,1:3) ...
                  -(1.0d0+eta(jj))*(1.0d0+zet(jj))*Hexa(5,1:3)-(1.0d0+eta(jj))*(1.0d0-zet(jj))*Hexa(6,1:3) ...
                  +(1.0d0+eta(jj))*(1.0d0-zet(jj))*Hexa(7,1:3)+(1.0d0+eta(jj))*(1.0d0+zet(jj))*Hexa(8,1:3))*0.125d0;

        J(1:3,2)=(-(1.0d0- xi(jj))*(1.0d0+zet(jj))*Hexa(1,1:3)-(1.0d0- xi(jj))*(1.0d0-zet(jj))*Hexa(2,1:3) ...
                  -(1.0d0+ xi(jj))*(1.0d0-zet(jj))*Hexa(3,1:3)-(1.0d0+ xi(jj))*(1.0d0+zet(jj))*Hexa(4,1:3) ...
                  +(1.0d0- xi(jj))*(1.0d0+zet(jj))*Hexa(5,1:3)+(1.0d0- xi(jj))*(1.0d0-zet(jj))*Hexa(6,1:3) ...
                  +(1.0d0+ xi(jj))*(1.0d0-zet(jj))*Hexa(7,1:3)+(1.0d0+ xi(jj))*(1.0d0+zet(jj))*Hexa(8,1:3))*0.125d0;

        J(1:3,3)=(+(1.0d0- xi(jj))*(1.0d0-eta(jj))*Hexa(1,1:3)-(1.0d0- xi(jj))*(1.0d0-eta(jj))*Hexa(2,1:3) ...
                  -(1.0d0+ xi(jj))*(1.0d0-eta(jj))*Hexa(3,1:3)+(1.0d0+ xi(jj))*(1.0d0-eta(jj))*Hexa(4,1:3) ...
                  +(1.0d0- xi(jj))*(1.0d0+eta(jj))*Hexa(5,1:3)-(1.0d0- xi(jj))*(1.0d0+eta(jj))*Hexa(6,1:3) ...
                  -(1.0d0+ xi(jj))*(1.0d0+eta(jj))*Hexa(7,1:3)+(1.0d0+ xi(jj))*(1.0d0+eta(jj))*Hexa(8,1:3))*0.125d0;
       detJ(jj,1)=det(J); 
    end
%%

    res_1=0;
    for jj = 1:NP
        ext_phi=Field_ext(PP(jj,1),0,PP(jj,3));
        res_1=res_1+ext_phi(2)*PP(jj,1)*WG_loc(jj)*detJ(jj,1);
    end
    Vext(ii,1)=res_1*2*pi/Area(ii);
end
    
end
