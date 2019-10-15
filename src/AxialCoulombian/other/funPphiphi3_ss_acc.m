function [P] = funPphiphi3_ss_acc(N_edge,N_GAUSS1,N_GAUSS2,Matrix_P0,G1,Aed)
%%
[xabsc1,weig1] = lgwt(N_GAUSS1,-1, 1); % calcolo i punti e pesi di Gauss nella linea -[1 1]
NP1=N_GAUSS1; %numero di punti di Gauss nell'esaedro 
xi1=zeros(NP1,1);
WG_loc1=zeros(NP1,1);
hh=1;
for ii = 1:N_GAUSS1 % 
   xi1(hh) =xabsc1(ii);
   WG_loc1(hh)= weig1(ii); % peso di Gauss
   hh=hh+1;
end
%%
[xabsc2,weig2] = lgwt(N_GAUSS2,-1, 1); % calcolo i punti e pesi di Gauss nella linea -[1 1]
NP2=N_GAUSS2; %numero di punti di Gauss nell'esaedro 
xi2=zeros(NP2,1);
WG_loc2=zeros(NP2,1);
hh=1;
for ii = 1:N_GAUSS2 % 
   xi2(hh) =xabsc2(ii);
   WG_loc2(hh)= weig2(ii); % peso di Gauss
   hh=hh+1;
end
%%
PPP1=zeros(NP1,3,N_edge);
lung=zeros(1,N_edge);
for ii = 1:N_edge
% lato target
    Ed1=Matrix_P0(G1(1:2,ii),1:3);
    Cent1=0.5*(Ed1(1,1:3)+Ed1(2,1:3));
    vec1=(Ed1(2,1:3)-Ed1(1,1:3))*0.5;
    lung1=norm(vec1)*2.0;
    lung(ii)=lung1;
% punti di gauss in globale
    for kk = 1:NP1 
        [PPP1(kk,:,ii)]=Cent1+vec1*xi1(kk);
    end
end
%%
PPP2=zeros(NP2,3,N_edge);
for ii = 1:N_edge
% lato target
    Ed1=Matrix_P0(G1(1:2,ii),1:3);
    Cent1=0.5*(Ed1(1,1:3)+Ed1(2,1:3));
    vec1=(Ed1(2,1:3)-Ed1(1,1:3))*0.5;
    lung1=norm(vec1)*2.0;
% punti di gauss in globale
    for kk = 1:NP2 
        [PPP2(kk,:,ii)]=Cent1+vec1*xi2(kk);
    end
end
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
P=zeros(N_edge,N_edge);
PP1=zeros(NP1,3);
PP2=zeros(NP2,3);
for ii = 1:N_edge % ciclo sui lati target
%% lato target
%     Ed1=Matrix_P0(G1(1:2,ii),1:3);
%     Cent1=0.5*(Ed1(1,1:3)+Ed1(2,1:3));
%     vec1=(Ed1(2,1:3)-Ed1(1,1:3))*0.5;
    lung1=lung(ii);%norm(vec1)*2.0;
%% punti di gauss in globale
%     for kk = 1:NP1 
%         [PP1(kk,1:3)]=Cent1+vec1*xi1(kk);
%     end
    PP1=PPP1(:,:,ii);
%%
    for jj = 1:N_edge % ciclo sui lati sorgente
    %% lato sou   
%         Ed2=Matrix_P0(G1(1:2,jj),1:3);
%         Cent2=0.5*(Ed2(1,1:3)+Ed2(2,1:3));
%         vec2=(Ed2(2,1:3)-Ed2(1,1:3))*0.5;
        lung2=lung(jj);%norm(vec2)*2.0;
    %% punti di gauss in globale
%         for kk = 1:NP2 
%             [PP2(kk,:)]=Cent2+vec2*xi2(kk);
%         end     
        PP2=PPP2(:,:,jj);
    %% 
        Pint=0.0d0;
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
%                 k=sqrt(k2);
%                 [e1]=my_ellipke(k2);
%                 G2Daxi=(ree)*(4*e1/Kg)/(4*pi);
%                 Pint=Pint+  rqq*G2Daxi*WG_loc1(qq)*WG_loc2(ee)*(lung1*0.5)*(lung2*0.5);
%                 Pint=Pint+  G2Daxi*WG_loc1(qq)*WG_loc2(ee)*(lung1*0.5)*(lung2*0.5);                
                hh=hh+1;
            end
        end
        [e1]=my_ellipke2(k2);
        G2Daxi=(ree).*(4*e1./Kg)/(4*pi);
        Pint=sum( G2Daxi.*WW*(lung1*0.5)*(lung2*0.5));
        Pint=Pint/(lung1*Aed(jj));
        P(ii,jj)=Pint;
%         P(jj,ii)=Pint;
    end
end
P=(P+P.')*0.5;
end
%%
function [k,e] = my_ellipke2(m)
a0 = 1;
b0 = sqrt(1-m);
c0 = NaN;
s0 = m;
i1 = 0; mm = Inf;
while mm > 1e-20
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