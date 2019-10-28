function [we,X8e,Y8e,Z8e,W8e,X11e,Y11e,Z11e,W11e] = fun_L_curledge_curledge_PrePro(N_edge,nv_shar_max,EV,Matrix_P0,VP,EV_ind,n_EV)
vol_hh(1:nv_shar_max)=0.0d0;
we(1:3,1:nv_shar_max,1:N_edge)=0.0d0;
X8e(1:nv_shar_max,1:8,1:N_edge)=0.0d0;
Y8e(1:nv_shar_max,1:8,1:N_edge)=0.0d0;
Z8e(1:nv_shar_max,1:8,1:N_edge)=0.0d0;
W8e(1:nv_shar_max,1:8,1:N_edge)=0.0d0;
X11e(1:nv_shar_max,1:11,1:N_edge)=0.0d0;
Y11e(1:nv_shar_max,1:11,1:N_edge)=0.0d0;
Z11e(1:nv_shar_max,1:11,1:N_edge)=0.0d0;
W11e(1:nv_shar_max,1:11,1:N_edge)=0.0d0;
for hh = 1:N_edge
    nv_hh=n_EV(hh); % numero di volumi attaccati al lato hh
    idV_hh=abs(EV(1:nv_shar_max,hh)); % indice globale dei volumi attaccati al lato hh
    for qq=1:nv_hh %ciclo sui volumi comuni lato hh
        pvt(1:3,1:4)=Matrix_P0(1:3,VP(1:4,idV_hh(qq))); % punti del tetraefor
        [nnorm,area,vol_hh(qq)]= tetareavol(pvt); %trovo le grandezze del tetraedro
        we(1:3,qq,hh) = fun_whitney_tetra_rot_edge(nnorm,area,vol_hh(qq),EV_ind(qq,hh)); %trovo le whitney del lato hh nel volume qq
        [X8,Y8,Z8,W8] = fun_my_tetra_quad_fv_dim(pvt,vol_hh(qq));% punti e pesi di Gauss nel volume qq 
        X8e(qq,1:8,hh)=X8;
        Y8e(qq,1:8,hh)=Y8;
        Z8e(qq,1:8,hh)=Z8;
        W8e(qq,1:8,hh)=W8;
        %
        [X11,Y11,Z11,W11] =fun_tetra_my_quad_N11_D4_dim(pvt,vol_hh(qq));% punti e pesi di Gauss nel volume qq 
        X11e(qq,1:11,hh)=X11;
        Y11e(qq,1:11,hh)=Y11;
        Z11e(qq,1:11,hh)=Z11;
        W11e(qq,1:11,hh)=W11;
    end              
end
end

