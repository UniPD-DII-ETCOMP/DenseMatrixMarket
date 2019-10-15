function [PPP1,PPP2,NPP,NP1,NP2,WW] = funLphiphi3_acc_PrePro(N_face,N_GAUSS1,N_GAUSS2,Matrix_P0,F1)
%%
[xabsc1,weig1] = lgwt(N_GAUSS1,-1, 1); % calcolo i punti e pesi di Gauss nella linea -[1 1]
NP1=N_GAUSS1^2; %numero di punti di Gauss nell'esaedro 
xi1=zeros(NP1,1);
eta1=zeros(NP1,1);
zet1=zeros(NP1,1);
WG_loc1=zeros(NP1,1);
hh=1;
for ii = 1:N_GAUSS1 % 
    for jj = 1:N_GAUSS1
           xi1(hh) =xabsc1(ii);
           zet1(hh)=xabsc1(jj);
           WG_loc1(hh)= weig1(ii)*weig1(jj); % peso di Gauss
           hh=hh+1;
    end 
end
%%
[xabsc2,weig2] = lgwt(N_GAUSS2,-1, 1); % calcolo i punti e pesi di Gauss nella linea -[1 1]
NP2=N_GAUSS2^2; %numero di punti di Gauss nell'esaedro 
xi2=zeros(NP2,1);
eta2=zeros(NP2,1);
zet2=zeros(NP2,1);
WG_loc2=zeros(NP2,1);
hh=1;
for ii = 1:N_GAUSS2 % 
    for jj = 1:N_GAUSS2
           xi2(hh) =xabsc2(ii);
           zet2(hh)=xabsc2(jj);
           WG_loc2(hh)= weig2(ii)*weig2(jj); % peso di Gauss
           hh=hh+1;
    end 
end
%%
L=zeros(N_face,N_face);
Hexa1=zeros(8,3);
Hexa2=zeros(8,3);
PP1=zeros(NP1,3);
PP2=zeros(NP2,3);
J1=zeros(3,3);
J2=zeros(3,3);
%%
NPP=NP1*NP2;
WW=zeros(1,NPP);
hh=1;
for ii = 1:NP1
    for jj = 1:NP2
        WW(hh)=    WG_loc1(ii).*WG_loc2(jj); 
        hh=hh+1;
    end
end
%%
PPP1=zeros(NP1,3,N_face);
PPP2=zeros(NP2,3,N_face);
Lself_vec=zeros(1,N_face);
for ii = 1:N_face % ciclo sulle facce
    Face1=Matrix_P0(F1(:,ii),:);
    for kk = 1:NP1 % punti di gauss in globale
        [PPP1(kk,:,ii)]=funTrilinear2(xi1(kk),0.0,zet1(kk),Face1);
    end
   cent=sum(Face1)/4;
    dr=abs(2*(Face1(1,1)-cent(1)));
    dz=abs(2*(Face1(1,3)-cent(3)));
    rm=cent(1);
    [Lself] = L_self_Axial_rect_cross(rm,dr,dz);
    Lself_vec(ii)=Lself; 
end
for ii = 1:N_face % ciclo sulle facce
    Face1=Matrix_P0(F1(:,ii),:);
    for kk = 1:NP2 % punti di gauss in globale
        [PPP2(kk,:,ii)]=funTrilinear2(xi2(kk),0.0,zet2(kk),Face1);
    end    
end
end
