function [L] = funLphiphi3_acc_low_rank(PPP1,PPP2,NPP,NP1,NP2,WW,iif,jjf)
%%
L=zeros(length(iif),length(jjf));
%%
for iic = 1:length(iif) % ciclo sulle facce
    ii=iif(iic);
    PP1=PPP1(:,:,ii);
    for jjc = 1:length(jjf)
        jj=jjf(jjc);
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
        L(iic,jjc)=Lint;
    end
end
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
