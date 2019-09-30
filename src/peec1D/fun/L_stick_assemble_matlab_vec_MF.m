function [L]=L_stick_assemble_matlab_vec_MF(NN,G,radius,npg,Wint,PPghtot,ll_htot,ut_htot,gauss_W,hhout,kk)

L=zeros(length(hhout),length(kk));
% PPh=zeros(3,2);
% PPk=zeros(3,2);
% gauss_P=zeros(npg,1);
% gauss_W=zeros(npg,1);
% ut_h=zeros(3,1);
% ut_k=zeros(3,1);
% PPgh=zeros(3,npg);

%%
for hhh=1:length(hhout)
    hh=hhout(hhh);
    PPgh=PPghtot(:,:,hh);
    ll_h=ll_htot(hh);
    ut_h=ut_htot(:,hh);
    %mutual-inductance
%     kk=hh+1:nSticks;
        PPk1=NN(1:3,G(1,kk));
        PPk2=NN(1:3,G(2,kk));
        ll_k=ll_htot(kk);
        ut_k=ut_htot(:,kk);
        integ = zeros(length(kk),1);
        for ii=1:npg
            ri=fun_my_norm([PPgh(1,ii)-PPk1(1,:);...
                            PPgh(2,ii)-PPk1(2,:);...
                            PPgh(3,ii)-PPk1(3,:)]);
            rf=fun_my_norm([PPgh(1,ii)-PPk2(1,:);...
                            PPgh(2,ii)-PPk2(2,:);...
                            PPgh(3,ii)-PPk2(3,:)]);
            eps=ll_k./(ri+rf).';
            log_eps=log((1+eps)./(1-eps));
            integ=integ+gauss_W(ii)*log_eps.*fun_my_dot(ut_h,ut_k).';
        end 
        L(hhh,1:length(kk))=ll_h*integ.';
end 
L=1.0d-7*0.5*L;
%%
[C,IA,IB] = intersect(hhout,kk);
hhh=1;
if ~isempty(C)
    for hh=C(1):C(end)
    %self-inductance
        ll_h=ll_htot(hh);
        aa=log(ll_h/radius(hh)+sqrt((ll_h/radius(hh))^2+1));
        bb=sqrt(1+(radius(hh)/ll_h)^2);
        cc=radius(hh)/ll_h;
        L(IA(hhh),IB(hhh))=1.0d-7*2*ll_h*(aa-bb+cc+(Wint)*0.25);
        hhh=hhh+1;
    end
end
end 




