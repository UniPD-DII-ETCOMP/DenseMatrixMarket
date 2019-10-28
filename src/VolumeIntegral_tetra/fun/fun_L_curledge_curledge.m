function [Matrix_L] = fun_L_curledge_curledge(N_edge,nv_shar_max,EV,Matrix_P0,VP,EV_ind,Neps,n_EV)
Matrix_L=zeros(N_edge,N_edge);
vol_hh(1:nv_shar_max)=0.0d0;
vol_kk(1:nv_shar_max)=0.0d0;
w_hh(1:3,1:nv_shar_max)=0.0d0;
w_kk(1:3,1:nv_shar_max)=0.0d0;
Xhq(1:nv_shar_max,1:8)=0.0d0;
Yhq(1:nv_shar_max,1:8)=0.0d0;
Zhq(1:nv_shar_max,1:8)=0.0d0;
Whq(1:nv_shar_max,1:8)=0.0d0;
Xke(1:nv_shar_max,1:11)=0.0d0;
Yke(1:nv_shar_max,1:11)=0.0d0;
Zke(1:nv_shar_max,1:11)=0.0d0;
Wke(1:nv_shar_max,1:11)=0.0d0;
bar_hq(1:3,1:nv_shar_max)=0.0d0;
bar_ke(1:3,1:nv_shar_max)=0.0d0;
for hh = 1:N_edge
    nv_hh=n_EV(hh); % numero di volumi attaccati al lato hh
    idV_hh=abs(EV(1:nv_shar_max,hh)); % indice globale dei volumi attaccati al lato hh
    for qq=1:nv_hh %ciclo sui volumi comuni lato hh
        pvt(1:3,1:4)=Matrix_P0(1:3,VP(1:4,idV_hh(qq))); % punti del tetraefor
        [nnorm,area,vol_hh(qq)]= tetareavol(pvt); %trovo le grandezze del tetraedro
        w_hh(1:3,qq) = fun_whitney_tetra_rot_edge(nnorm,area,vol_hh(qq),EV_ind(qq,hh)); %trovo le whitney del lato hh nel volume qq
        [X8,Y8,Z8,W8] = fun_my_tetra_quad_fv_dim(pvt,vol_hh(qq));% punti e pesi di Gauss nel volume qq 
        Xhq(qq,1:8)=X8;
        Yhq(qq,1:8)=Y8;
        Zhq(qq,1:8)=Z8;
        Whq(qq,1:8)=W8;
        bar_hq(1:3,qq)=(pvt(1:3,1)+pvt(1:3,2)+pvt(1:3,3)+pvt(1:3,4))*0.25;% baricentro del tetraedro qq 
	end 
    for kk = 1:N_edge
        nv_kk=n_EV(kk); % numero di volumi attaccati al lato kk
        idV_kk=abs(EV(1:nv_shar_max,kk)); % indice globale dei volumi attaccati al lato kk
        for ee=1:nv_kk %ciclo sui volumi comuni lato kk
            pvt(1:3,1:4)=Matrix_P0(1:3,VP(1:4,idV_kk(ee))); % punti del tetraefor
            [nnorm,area,vol_kk(ee)]= tetareavol(pvt); %trovo le grandezze del tetraedro
            w_kk(1:3,ee) = fun_whitney_tetra_rot_edge(nnorm,area,vol_kk(ee),EV_ind(ee,kk)); %trovo le whitney del lato hh nel volume qq
            [X11,Y11,Z11,W11] =fun_tetra_my_quad_N11_D4_dim(pvt,vol_kk(ee));% punti e pesi di Gauss nel volume ee 
            Xke(ee,1:11)=X11;
            Yke(ee,1:11)=Y11;
            Zke(ee,1:11)=Z11;
            Wke(ee,1:11)=W11;
            bar_ke(1:3,ee)=(pvt(1:3,1)+pvt(1:3,2)+pvt(1:3,3)+pvt(1:3,4))*0.25;% baricentro del tetraedro ee        
        end              
%compute selected matrix entry
	for qq=1:nv_hh %ciclo su volumi lato hh
% 		Vqq=idV_hh(qq); % indice globale del volume qq
		for ee=1:nv_kk %ciclo volumi lato kk
%             Vee=idV_kk(ee); % indice globale del volume ee
	%------------- Quantità condivise -------------- %		
            fort_wh_wk=sign(EV(qq,hh))*sign(EV(ee,kk))*fun_my_dot(w_hh(1:3,qq),w_kk(1:3,ee)); % fort delle witney con segno giuto	
%             dist_qe = fun_my_norm(bar_hq(1:3,qq)-bar_ke(1:3,ee)); % distanza tra i baricentri dei tetraedr
	%------------- R -------------- %
%             if Vqq==Vee 
%                 R=vol_hh(qq)*fort_wh_wk;
%                 Matrix_RL(hh,kk)=Matrix_RL(hh,kk)+R*rho(Vqq);
% %                 SYS_re=SYS_re+R*rho_re(Vqq);
% %                 SYS_im=SYS_im+R*rho_im(Vqq);
%             end 	
	%------------- L -------------- %
            L=0.0d0;
            for ii = 1:8 % ciclo sui punti di Gauss del tetraedro qq 
               for jj = 1:11 % ciclo sui punti di Gauss del tetraedro ee
                   dist = fun_my_norm_eps([Xhq(qq,ii),Yhq(qq,ii),Zhq(qq,ii)]-[Xke(ee,jj),Yke(ee,jj),Zke(ee,jj)],Neps);
                   L=L+Whq(qq,ii)*Wke(ee,jj)/dist; 
               end 
            end 
            L=L*fort_wh_wk;
            Matrix_L(hh,kk)=Matrix_L(hh,kk)+L;
%              Matrix_RL(hh,kk)=Matrix_RL(hh,kk)-omega*L*sin(-beta*dist_qe)+1j*omega*L*cos(-beta*dist_qe);
%             SYS_re=SYS_re-omega*L*dsin(-beta*dist_qe);
%             SYS_im=SYS_im+omega*L*dcos(-beta*dist_qe);
		end 
    end        
    end
end
Matrix_L=Matrix_L*1e-7;
end

