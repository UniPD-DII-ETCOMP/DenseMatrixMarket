function [L] = funNme_acc(N_face_e,N_GAUSS1,N_point_m,Area_e,Matrix_P0_e,F1_e,Matrix_P0_m)
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
PPP1=zeros(NP1,3,N_face_e);
for ii = 1:N_face_e % ciclo sulle facce
    Face1=Matrix_P0_e(F1_e(:,ii),:);
    for kk = 1:NP1 % punti di gauss in globale
        [PPP1(kk,:,ii)]=funTrilinear2(xi1(kk),0.0,zet1(kk),Face1);
    end
end
%%
L=zeros(N_point_m,N_face_e);
for ii = 1:N_face_e % ciclo sulle facce
    PP1=PPP1(:,:,ii);
%% 
    for jj = 1:N_point_m
        Pm=Matrix_P0_m(jj,:).';
        ree = Pm(1); %source point r
        rqq = PP1(:,1); %target point r
        h = Pm(3)-PP1(:,3);
        Kg2=(ree(:)+rqq(:)).^2+h.^2;
        k2=4*ree(:)*rqq(:)./Kg2;
        Kg=sqrt(Kg2);
        [e1,e2]=my_ellipke2(k2);
        G2Daxi=(ree./Kg).*((2-k2).*e1-2*e2)./(k2);
        Lint=sum(rqq.*G2Daxi.*WG_loc1);
        Lint=Lint/2;
        L(jj,ii)=Lint;
    end
end
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