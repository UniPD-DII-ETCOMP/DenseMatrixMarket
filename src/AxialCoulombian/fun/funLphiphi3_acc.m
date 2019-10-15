function [L] = funLphiphi3_acc(N_face,N_GAUSS1,N_GAUSS2,Matrix_P0,F1)
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
%%
for ii = 1:N_face % ciclo sulle facce
    Face=Matrix_P0(F1(1:4,ii),:);
    PP1=PPP1(:,:,ii);
    for jj = 1:N_face
        PP2=PPP2(:,:,jj);
        hh=1;
        k2=zeros(1,NPP);
        Kg=zeros(1,NPP);
        ree=zeros(1,NPP);
        rqq=zeros(1,NPP);
        for qq = 1:NP1
            for ee = 1:NP2
                ree(hh) = PP2(ee,1); %source point r
                rqq(hh) = PP1(qq,1); %target point r
                h = PP2(ee,3)-PP1(qq,3);
                Kg2=(ree(hh)+rqq(hh))^2+h^2;
                k2(hh)=4*ree(hh)*rqq(hh)/Kg2;
                Kg(hh)=sqrt(Kg2);
                hh=hh+1;
            end
        end
        [e1,e2]=my_ellipke2(k2);
        G2Daxi=(4*ree./Kg).*((2-k2).*e1-2*e2)./(4*pi*k2);
        Lint=sum(rqq.*G2Daxi.*WW);        
        L(ii,jj)=Lint;
    end
end
L=(L+L.')/2;
L=L-diag(diag(L))+diag(Lself_vec/(2*pi/16));
L=L*2*pi/16;
end
%%
function [k,e] = my_ellipke2(m)
a0 = 1;
b0 = sqrt(1-m);
c0 = NaN;
s0 = m;
i1 = 0; mm = Inf;
while mm > 1e-12
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    i1 = i1 + 1;
    w1 = 2^i1*c1.^2;
    mm = max(w1(:));   
    % test for stagnation (may happen for TOL < machine precision)
%     if isequal(c0, c1)
%         error(message('MATLAB:ellipke:FailedConvergence'));
%     end   
    s0 = s0 + w1;  
    a0 = a1;  b0 = b1;  c0 = c1;
end
k = pi./(2*a1);
e = k.*(1-s0/2);
im = find(m==1);
if ~isempty(im)
    e(im) = ones(length(im),1);
    k(im) = inf;
end
end
