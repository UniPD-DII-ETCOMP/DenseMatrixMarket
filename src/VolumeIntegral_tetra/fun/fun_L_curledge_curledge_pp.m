function [Matrix_L] = fun_L_curledge_curledge_pp(EV,Neps,n_EV,X8e,...
    Y8e,Z8e,W8e,X11e,Y11e,Z11e,W11e,we,Nhh,Nkk,hhout,kkout)
Matrix_L=zeros(Nhh,Nkk);
for hhloc = 1:Nhh
    hh=hhout(hhloc);
    nv_hh=n_EV(hh); 
    for kkloc = 1:Nkk
        kk=kkout(kkloc);
        nv_kk=n_EV(kk);              
        for qq=1:nv_hh
            w_hh(1:3)=we(1:3,qq,hh);
            Xhq(1:8)=X8e(qq,1:8,hh);
            Yhq(1:8)=Y8e(qq,1:8,hh);
            Zhq(1:8)=Z8e(qq,1:8,hh);
            Whq(1:8)=W8e(qq,1:8,hh);        
            for ee=1:nv_kk
                w_kk(1:3)=we(1:3,ee,kk);
                Xke(1:11)=X11e(ee,1:11,kk);
                Yke(1:11)=Y11e(ee,1:11,kk);
                Zke(1:11)=Z11e(ee,1:11,kk);
                Wke(1:11)=W11e(ee,1:11,kk);            
                %
                fort_wh_wk=sign(EV(qq,hh))*sign(EV(ee,kk))*fun_my_dot(w_hh(1:3),w_kk(1:3)); % 
                L=0.0d0;
                for ii = 1:8 %
                   for jj = 1:11 % 
                       dist = fun_my_norm_eps([Xhq(ii),Yhq(ii),Zhq(ii)]-[Xke(jj),Yke(jj),Zke(jj)],Neps);
                       L=L+Whq(ii)*Wke(jj)/dist; 
                   end 
                end 
                L=L*fort_wh_wk;
                Matrix_L(hhloc,kkloc)=Matrix_L(hhloc,kkloc)+L;
            end 
        end        
    end
end
Matrix_L=Matrix_L*1e-7;

global nncoeff;
nncoeff=nncoeff+Nhh*Nkk;
end

