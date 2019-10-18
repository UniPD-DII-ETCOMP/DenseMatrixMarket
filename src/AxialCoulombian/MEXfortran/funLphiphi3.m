function [L] = funLphiphi3(N_face,Area,N_GAUSS1,N_GAUSS2,Matrix_P0,F1)
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
detJ1=zeros(NP1,1);
detJ2=zeros(NP2,1);
J1=zeros(3,3);
J2=zeros(3,3);
for ii = 1:N_face % ciclo sulle facce
% ii
    Face1=Matrix_P0(F1(1:4,ii),:);
%     Hexa1=zeros(8,3);
    Hexa1(1:4,[1,3])=Face1(1:4,[1,3]);
    Hexa1(5:8,[1,3])=Face1(1:4,[1,3]);
    Hexa1(1:4,2)=-1.0;
    Hexa1(5:8,2)=1.0;
%     PP1=zeros(NP1,3);
    for kk = 1:NP1 % punti di gauss in globale
        [PP1(kk,:)]=funTrilinear(xi1(kk),0.0,zet1(kk),Hexa1);
    end
%     detJ1=zeros(NP1,1);
%     J1=zeros(3,3);
    for kk = 1%:NP1
        J1(1:3,1)=(-(1.0d0-eta1(kk))*(1.0d0+zet1(kk))*Hexa1(1,1:3)-(1.0d0-eta1(kk))*(1.0d0-zet1(kk))*Hexa1(2,1:3) ...
                  +(1.0d0-eta1(kk))*(1.0d0-zet1(kk))*Hexa1(3,1:3)+(1.0d0-eta1(kk))*(1.0d0+zet1(kk))*Hexa1(4,1:3) ...
                  -(1.0d0+eta1(kk))*(1.0d0+zet1(kk))*Hexa1(5,1:3)-(1.0d0+eta1(kk))*(1.0d0-zet1(kk))*Hexa1(6,1:3) ...
                  +(1.0d0+eta1(kk))*(1.0d0-zet1(kk))*Hexa1(7,1:3)+(1.0d0+eta1(kk))*(1.0d0+zet1(kk))*Hexa1(8,1:3))*0.125d0;

        J1(1:3,2)=(-(1.0d0- xi1(kk))*(1.0d0+zet1(kk))*Hexa1(1,1:3)-(1.0d0- xi1(kk))*(1.0d0-zet1(kk))*Hexa1(2,1:3) ...
                  -(1.0d0+ xi1(kk))*(1.0d0-zet1(kk))*Hexa1(3,1:3)-(1.0d0+ xi1(kk))*(1.0d0+zet1(kk))*Hexa1(4,1:3) ...
                  +(1.0d0- xi1(kk))*(1.0d0+zet1(kk))*Hexa1(5,1:3)+(1.0d0- xi1(kk))*(1.0d0-zet1(kk))*Hexa1(6,1:3) ...
                  +(1.0d0+ xi1(kk))*(1.0d0-zet1(kk))*Hexa1(7,1:3)+(1.0d0+ xi1(kk))*(1.0d0+zet1(kk))*Hexa1(8,1:3))*0.125d0;

        J1(1:3,3)=(+(1.0d0- xi1(kk))*(1.0d0-eta1(kk))*Hexa1(1,1:3)-(1.0d0- xi1(kk))*(1.0d0-eta1(kk))*Hexa1(2,1:3) ...
                  -(1.0d0+ xi1(kk))*(1.0d0-eta1(kk))*Hexa1(3,1:3)+(1.0d0+ xi1(kk))*(1.0d0-eta1(kk))*Hexa1(4,1:3) ...
                  +(1.0d0- xi1(kk))*(1.0d0+eta1(kk))*Hexa1(5,1:3)-(1.0d0- xi1(kk))*(1.0d0+eta1(kk))*Hexa1(6,1:3) ...
                  -(1.0d0+ xi1(kk))*(1.0d0+eta1(kk))*Hexa1(7,1:3)+(1.0d0+ xi1(kk))*(1.0d0+eta1(kk))*Hexa1(8,1:3))*0.125d0;
       detJ1(:)=det(J1); 
    end
%%
cent=sum(Face1)/4;
dr=abs(2*(Face1(1,1)-cent(1)));
dz=abs(2*(Face1(1,3)-cent(3)));
rm=cent(1);
[Lself] = L_self_Axial_rect_cross(rm,dr,dz);
L(ii,ii)=Lself;
%% 
    for jj = ii+1:N_face
        Face2=Matrix_P0(F1(1:4,jj),:);
%         Hexa2=zeros(8,3);
        Hexa2(1:4,[1,3])=Face2(1:4,[1,3]);
        Hexa2(5:8,[1,3])=Face2(1:4,[1,3]);
        Hexa2(1:4,2)=-1.0;
        Hexa2(5:8,2)=1.0;
%         PP2=zeros(NP2,3);
        for kk = 1:NP2 % punti di gauss in globale
            [PP2(kk,:)]=funTrilinear(xi2(kk),0.0,zet2(kk),Hexa2);
        end
%         detJ2=zeros(NP2,1);
%         J2=zeros(3,3);
        for kk = 1%:NP2
            J2(1:3,1)=(-(1.0d0-eta2(kk))*(1.0d0+zet2(kk))*Hexa2(1,1:3)-(1.0d0-eta2(kk))*(1.0d0-zet2(kk))*Hexa2(2,1:3) ...
                      +(1.0d0-eta2(kk))*(1.0d0-zet2(kk))*Hexa2(3,1:3)+(1.0d0-eta2(kk))*(1.0d0+zet2(kk))*Hexa2(4,1:3) ...
                      -(1.0d0+eta2(kk))*(1.0d0+zet2(kk))*Hexa2(5,1:3)-(1.0d0+eta2(kk))*(1.0d0-zet2(kk))*Hexa2(6,1:3) ...
                      +(1.0d0+eta2(kk))*(1.0d0-zet2(kk))*Hexa2(7,1:3)+(1.0d0+eta2(kk))*(1.0d0+zet2(kk))*Hexa2(8,1:3))*0.125d0;

            J2(1:3,2)=(-(1.0d0- xi2(kk))*(1.0d0+zet2(kk))*Hexa2(1,1:3)-(1.0d0- xi2(kk))*(1.0d0-zet2(kk))*Hexa2(2,1:3) ...
                      -(1.0d0+ xi2(kk))*(1.0d0-zet2(kk))*Hexa2(3,1:3)-(1.0d0+ xi2(kk))*(1.0d0+zet2(kk))*Hexa2(4,1:3) ...
                      +(1.0d0- xi2(kk))*(1.0d0+zet2(kk))*Hexa2(5,1:3)+(1.0d0- xi2(kk))*(1.0d0-zet2(kk))*Hexa2(6,1:3) ...
                      +(1.0d0+ xi2(kk))*(1.0d0-zet2(kk))*Hexa2(7,1:3)+(1.0d0+ xi2(kk))*(1.0d0+zet2(kk))*Hexa2(8,1:3))*0.125d0;

            J2(1:3,3)=(+(1.0d0- xi2(kk))*(1.0d0-eta2(kk))*Hexa2(1,1:3)-(1.0d0- xi2(kk))*(1.0d0-eta2(kk))*Hexa2(2,1:3) ...
                      -(1.0d0+ xi2(kk))*(1.0d0-eta2(kk))*Hexa2(3,1:3)+(1.0d0+ xi2(kk))*(1.0d0-eta2(kk))*Hexa2(4,1:3) ...
                      +(1.0d0- xi2(kk))*(1.0d0+eta2(kk))*Hexa2(5,1:3)-(1.0d0- xi2(kk))*(1.0d0+eta2(kk))*Hexa2(6,1:3) ...
                      -(1.0d0+ xi2(kk))*(1.0d0+eta2(kk))*Hexa2(7,1:3)+(1.0d0+ xi2(kk))*(1.0d0+eta2(kk))*Hexa2(8,1:3))*0.125d0;
           detJ2(:)=det(J2); 
        end 
        % 
        Lint=0.0d0;
        for qq = 1:NP1
            for ee = 1:NP2
                ree = PP2(ee,1); %source point r
                rqq = PP1(qq,1); %target point r
                h = PP2(ee,3)-PP1(qq,3);
                Kg2=(ree+rqq)^2+h^2;
                k2=4*ree*rqq/Kg2;
                Kg=sqrt(Kg2);
                k=sqrt(k2);
                [e1,e2]=my_ellipke(k2);
                G2Daxi=(4*ree/Kg)*((2-k2)*e1-2*e2)/(4*pi*k2);
                Lint=Lint+  rqq*G2Daxi*WG_loc1(qq)*WG_loc2(ee)*detJ1(qq)*detJ2(ee);
            end
        end
        Lint=2*pi*Lint/(Area(ii)*Area(jj));
        L(ii,jj)=Lint;
        L(jj,ii)=Lint;
    end
end
end
%%
function [k,e] = my_ellipke(m,tol)
%ELLIPKE Complete elliptic integral.
%   [K,E] = ELLIPKE(M) returns the value of the complete elliptic
%   integrals of the first and second kinds, evaluated for each
%   element of M.  As currently implemented, M is limited to 0 <= M <= 1.
%   
%   [K,E] = ELLIPKE(M,TOL) computes the complete elliptic integrals to
%   the accuracy TOL instead of the default TOL = EPS(CLASS(M)).
%
%   Some definitions of the complete elliptic integrals use the modulus
%   k instead of the parameter M.  They are related by M = k^2.
%
%   Class support for input M:
%      float: double, single
%
%   See also ELLIPJ.

%   Modified to include the second kind by Bjorn Bonnevier
%   from the Alfven Laboratory, KTH, Stockholm, Sweden
%   Copyright 1984-2013 The MathWorks, Inc. 

%   ELLIPKE uses the method of the arithmetic-geometric mean
%   described in [1].

%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications", 1965, 17.6.

if nargin<1
  error(message('MATLAB:ellipke:NotEnoughInputs')); 
end

classin = superiorfloat(m);

if nargin<2, tol = eps(classin); end
if ~isreal(m) || ~isreal(tol)
    error(message('MATLAB:ellipke:ComplexInputs'))
end
if isempty(m), k = zeros(size(m),classin); e = k; return, end
if any(m(:) < 0) || any(m(:) > 1)
  error(message('MATLAB:ellipke:MOutOfRange'));
end
if ~isscalar(tol) || tol < 0 || ~isfinite(tol)
  error(message('MATLAB:ellipke:NegativeTolerance'));
end

a0 = 1;
b0 = sqrt(1-m);
c0 = NaN;
s0 = m;
i1 = 0; mm = Inf;
while mm > tol
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    i1 = i1 + 1;
    w1 = 2^i1*c1.^2;
    mm = max(w1(:));
    
    % test for stagnation (may happen for TOL < machine precision)
    if isequal(c0, c1)
        error(message('MATLAB:ellipke:FailedConvergence'));
    end
    
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

function [L] = L_self_Axial_rect_cross(rm,dr,dz)
%%
L = rm*(...
       log(8*rm/dr) + 1.0/12.0 - pi*dz/(3.0*dr) - 0.5*log(1+(dz^2)/(dr^2))...
       +(1.0/12.0)*(dr^2/dz^2)*log(1+dz^2/dr^2)+(1/12)*dz^2/dr^2*log(1+dr^2/dz^2)...
       +(2/3)*(dz/dr-dr/dz)*atan2(dz,dr)+dr^2/(96*rm^2)*((log(8*rm/dr)...
       -0.5*log(1+dz^2/dr^2))*(1+3*dz^2/dr^2)+3.45*dz^2/dr^2 ...
       +221/60 -1.6*pi*dz^3/dr^3 +3.2*dz^3/dr^3*atan2(dz,dr)-0.1...
       *dr^2/dz^2*log(1+dz^2/dr^2)+0.5*dz^4/dr^4*log(1+dr^2/dz^2)));

end
