function [Rdiag] = funRphi(N_face,rho,Area,N_GAUSS,Matrix_P0,F1)
%%
[xabsc,weig] = lgwt(N_GAUSS,-1, 1); % calcolo i punti e pesi di Gauss nella linea -[1 1]
NP=N_GAUSS^2; %numero di punti di Gauss nell'esaedro 
xi=zeros(NP,1);
eta=zeros(NP,1);
zet=zeros(NP,1);
WG_loc=zeros(NP,1);
hh=1;
for ii = 1:N_GAUSS % 
    for jj = 1:N_GAUSS
           xi(hh) =xabsc(ii);
           zet(hh)=xabsc(jj);
           WG_loc(hh)= weig(ii)*weig(jj); % peso di Gauss
           hh=hh+1;
    end 
end
%%
Rdiag=zeros(N_face,1);
Hexa=zeros(8,3);
PP=zeros(NP,3);
for ii = 1:N_face % ciclo sulle facce
    Face=Matrix_P0(F1(1:4,ii),:);
    Hexa(1:4,[1,3])=Face(1:4,[1,3]);
    Hexa(5:8,[1,3])=Face(1:4,[1,3]);
    Hexa(1:4,2)=-1;
    Hexa(5:8,2)=1;
%%
    for jj = 1:NP % punti di gauss in globale
        [PP(jj,:)]=funTrilinear(xi(jj),0.0,zet(jj),Hexa);
    end
%%
    Vol=0;
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
       detJ=det(J); 
       Vol =Vol+PP(jj,1)*WG_loc(jj)*detJ;
    end
    Vol=Vol*2*pi;
%%
    Rdiag(ii,1)=(rho(ii)/(Area(ii)^2))*Vol;
end
end

