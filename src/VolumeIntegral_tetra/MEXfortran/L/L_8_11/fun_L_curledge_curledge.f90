subroutine fun_L_curledge_curledge(N_volu,N_node,N_edge,nv_max,EV,Matrix_P0, & 
          VP,EV_ind,Neps,n_EV,N_thread,Nhh,Nkk,hhout,kkout,Matrix_L)
implicit none
integer*8 nv_max, EV(nv_max+1,N_edge), EV_ind(nv_max+1,N_edge)
integer*8 n_EV(N_edge), Nhh, Nkk, hhout(Nhh), kkout(Nkk)
real*8 Matrix_P0(3,N_node), Matrix_L(Nhh,Nkk)
integer*8 N_node, N_volu, N_thread,   VP(4,N_volu), N_edge
real*8 pvt(3,4), nnorm(3,4),area(4),vol
real*8 w1(3),w2(3)
integer*8 ii, jj, ind1, ind2, Neps
real*8 vol_hh(nv_max), vol_kk(nv_max)
real*8 w_hh(3,nv_max)
real*8 w_kk(3,nv_max)
real*8 Xhq(nv_max,8)
real*8 Yhq(nv_max,8)
real*8 Zhq(nv_max,8)
real*8 Whq(nv_max,8)
real*8 Xke(nv_max,11)
real*8 Yke(nv_max,11)
real*8 Zke(nv_max,11)
real*8 Wke(nv_max,11)
real*8 bar_hq(3,nv_max)
real*8 bar_ke(3,nv_max)
integer*8 nv_hh, idV_hh(nv_max), hh, kk
integer*8 qq,ee
real*8 X8(8), Y8(8), Z8(8), W8(8)
real*8 X11(11), Y11(11), Z11(11), W11(11)
integer*8 nv_kk, idV_kk(nv_max), kkloc, hhloc
real*8 fort_wh_wk, L, dist, vec(3)
Matrix_L(1:Nhh,1:Nkk)=0.0
vol_hh(1:nv_max)=0.0d0;
vol_kk(1:nv_max)=0.0d0;
w_hh(1:3,1:nv_max)=0.0d0;
w_kk(1:3,1:nv_max)=0.0d0;
Xhq(1:nv_max,1:8)=0.0d0;
Yhq(1:nv_max,1:8)=0.0d0;
Zhq(1:nv_max,1:8)=0.0d0;
Whq(1:nv_max,1:8)=0.0d0;
Xke(1:nv_max,1:11)=0.0d0;
Yke(1:nv_max,1:11)=0.0d0;
Zke(1:nv_max,1:11)=0.0d0;
Wke(1:nv_max,1:11)=0.0d0;
bar_hq(1:3,1:nv_max)=0.0d0;
bar_ke(1:3,1:nv_max)=0.0d0;
call omp_set_num_threads(N_thread)
!$OMP PARALLEL SHARED(N_edge,Matrix_P0,Matrix_L,n_EV,EV,nv_max,VP,EV_ind,Neps,Nhh,hhout,Nkk,kkout)
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(hh,nv_hh,idV_hh,qq,pvt,nnorm,area,vol_hh,w_hh,X8,Y8,Z8,W8,Xhq, & 
!$OMP                                Yhq,Zhq,Whq,bar_hq,kk,nv_kk,idV_kk,ee,vol_kk,w_kk, & 
!$OMP                                X11,Y11,Z11,W11,Xke,Yke,Zke,Wke,bar_ke,fort_wh_wk, & 
!$OMP                                L,ii,jj,dist,hhloc,kkloc)
do hhloc = 1,Nhh
	hh=hhout(hhloc)
    nv_hh=n_EV(hh); ! numero di volumi attaccati al lato hh
    idV_hh=abs(EV(1:nv_max,hh)); ! indice globale dei volumi attaccati al lato hh
    do qq=1,nv_hh !ciclo sui volumi comuni lato hh
        pvt(1:3,1:4)=Matrix_P0(1:3,VP(1:4,idV_hh(qq))); ! punti del tetraefor
		call tetareavol(pvt,nnorm,area,vol_hh(qq)); !trovo le grandezze del tetraedro
		call fun_whitney_tetra_rot_edge(nnorm,area,vol_hh(qq),EV_ind(qq,hh),w_hh(1:3,qq)); !trovo le whitney del lato hh nel volume qq
		call fun_my_tetra_quad_fv_dim(pvt,vol_hh(qq),X8,Y8,Z8,W8);! punti e pesi di Gauss nel volume qq 
		Xhq(qq,1:8)=X8;
        Yhq(qq,1:8)=Y8;
        Zhq(qq,1:8)=Z8;
        Whq(qq,1:8)=W8;
        bar_hq(1:3,qq)=(pvt(1:3,1)+pvt(1:3,2)+pvt(1:3,3)+pvt(1:3,4))*0.25;! baricentro del tetraedro qq 
	enddo 
    do kkloc = 1,Nkk
		kk=kkout(kkloc)
        nv_kk=n_EV(kk); ! numero di volumi attaccati al lato kk
        idV_kk=abs(EV(1:nv_max,kk)); ! indice globale dei volumi attaccati al lato kk
        do ee=1,nv_kk !ciclo sui volumi comuni lato kk
            pvt(1:3,1:4)=Matrix_P0(1:3,VP(1:4,idV_kk(ee))); ! punti del tetraefor
            call tetareavol(pvt,nnorm,area,vol_kk(ee)); !trovo le grandezze del tetraedro
            call fun_whitney_tetra_rot_edge(nnorm,area,vol_kk(ee),EV_ind(ee,kk),w_kk(1:3,ee)); !trovo le whitney del lato hh nel volume qq
            call fun_tetra_my_quad_N11_D4_dim(pvt,vol_kk(ee),X11,Y11,Z11,W11);! punti e pesi di Gauss nel volume ee 
            Xke(ee,1:11)=X11;
            Yke(ee,1:11)=Y11;
            Zke(ee,1:11)=Z11;
            Wke(ee,1:11)=W11;
            bar_ke(1:3,ee)=(pvt(1:3,1)+pvt(1:3,2)+pvt(1:3,3)+pvt(1:3,4))*0.25;! baricentro del tetraedro ee        
        enddo              
		do qq=1,nv_hh !ciclo su volumi lato hh
			do ee=1,nv_kk !ciclo volumi lato kk
				fort_wh_wk=sign(1.0d0,real(EV(qq,hh),8))*sign(1.0d0,real(EV(ee,kk),8))*dot_product(w_hh(1:3,qq),w_kk(1:3,ee)); ! fort delle witney con segno giuto	
				L=0.0d0;
				do ii = 1,8 ! ciclo sui punti di Gauss del tetraedro qq 
				   do jj = 1,11 ! ciclo sui punti di Gauss del tetraedro ee
					   call norm_eps([Xhq(qq,ii),Yhq(qq,ii),Zhq(qq,ii)]-[Xke(ee,jj),Yke(ee,jj),Zke(ee,jj)],Neps,dist);
					   L=L+Whq(qq,ii)*Wke(ee,jj)/dist; 
				   enddo 
				enddo 
				L=L*fort_wh_wk;
				Matrix_L(hhloc,kkloc)=Matrix_L(hhloc,kkloc)+L;
			enddo 
		enddo       
    enddo
enddo
!$OMP END DO NOWAIT 
!$OMP END PARALLEL 
Matrix_L=Matrix_L*1.0d-7;
end subroutine fun_L_curledge_curledge

subroutine norm_eps(a,eps,c) 
implicit none
real(kind=8), dimension(3), intent(in) :: a
real(kind=8)    :: c 
real(kind=8)    :: eps 
c=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3)+eps*eps)
end subroutine	norm_eps

subroutine tetareavol(v,norm,area,vol)
implicit none
double precision v(3,4),norm(3,4),area(4),vol
double precision r(3),s(3),t(3),vp(3),sp
integer*8 i
integer permut(4,4) ! face and nodes formin face
data permut / 1,2,3,4, 4,3,4,1, 2,1,1,3, 3,4,2,2 /
do i=1,4
r(:)=v(:,permut(i,3))-v(:,permut(i,2))
s(:)=v(:,permut(i,4))-v(:,permut(i,2))
vp(1)=r(2)*s(3)-r(3)*s(2)
vp(2)=-(r(1)*s(3)-r(3)*s(1))
vp(3)=r(1)*s(2)-r(2)*s(1)
area(i)=sqrt(vp(1)**2+vp(2)**2+vp(3)**2)
norm(:,i)=vp(:)/area(i)
enddo
area(:)=0.5d0*area(:)
t(:)=v(:,permut(4,2))-v(:,permut(4,1))
sp=dot_product(vp,t)
vol=abs(sp)/6.0d0
end subroutine tetareavol


subroutine fun_whitney_tetra_rot_edge(norm,area,vol,k,rotwi)
implicit none
integer*8 locedg(6,2)
data locedg /1,2,3,4,4,4,2,3,1,1,2,3/
double precision gradlambda(3,4),norm(3,4),area(4),vol,rotwi(3)
integer*8 k,nd1loc,nd2loc
call x2gradlambda(gradlambda,norm,area,vol)
!do k=1,6
nd1loc=locedg(k,1)
nd2loc=locedg(k,2)
call vecprod(gradlambda(1:3,nd1loc),gradlambda(1:3,nd2loc),rotwi)
rotwi=2.0d0*rotwi
!enddo
end subroutine fun_whitney_tetra_rot_edge

subroutine x2gradlambda(gradlambda,norm,area,vol)
implicit none
double precision gradlambda(3,4)
double precision norm(3,4),area(4),vol,vol3
integer*8 i
!call tetareavol(v,norm,area,vol)
vol3=3.0d0*vol
do i=1,4
gradlambda(1:3,i)=-norm(:,i)*area(i)/vol3
enddo
end subroutine x2gradlambda

subroutine vecprod(v1,v2,vp)
implicit none
double precision v1(3),v2(3),vp(3)
vp(1)=   v1(2)*v2(3)-v1(3)*v2(2)
vp(2)= -(v1(1)*v2(3)-v1(3)*v2(1))
vp(3)=   v1(1)*v2(2)-v1(2)*v2(1)
end subroutine 

subroutine fun_my_tetra_quad_fv_dim(tetra,V,X,Y,Z,W)
implicit none
real*8 tetra(3,4), V, X(8), Y(8), Z(8), W(8), F1(3), F2(3), F3(3), F4(3)
F1=(tetra(1:3,1)+tetra(1:3,2)+tetra(1:3,3))/3.0!sum(tetra(1:3,1:3))/3;
F2=(tetra(1:3,2)+tetra(1:3,3)+tetra(1:3,4))/3.0!sum(tetra(2:4,1:3))/3;
F3=(tetra(1:3,1)+tetra(1:3,3)+tetra(1:3,4))/3.0!sum(tetra([1,3,4],:))/3;
F4=(tetra(1:3,1)+tetra(1:3,2)+tetra(1:3,4))/3.0!sum(tetra([1,2,4],:))/3;
X = [tetra(1,1:4),F1(1),F2(1),F3(1),F4(1)]
Y = [tetra(2,1:4),F1(2),F2(2),F3(2),F4(2)]
Z = [tetra(3,1:4),F1(3),F2(3),F3(3),F4(3)]
W  = [V/40,V/40,V/40,V/40,9*V/40,9*V/40,9*V/40,9*V/40] ![V/40*ones(4,1);9*V/40*ones(4,1)]
end subroutine fun_my_tetra_quad_fv_dim

subroutine fun_tetra_my_quad_N11_D4_dim(tetra,Vol,Xg,Yg,Zg,WWg)
implicit none
real*8 tetra(3,4), Vol, Xg(11), Yg(11), Zg(11), WWg(11)
real*8 wg(11), p1(11), p2(11), p3(11), p4(11)
! 3D Gauss-Legendre weights for N = 11, D = 4;
wg=[-0.01315556,0.007622222,0.007622222,0.007622222,0.007622222,0.02488889,0.02488889,0.02488889,0.02488889,0.02488889,0.02488889]
p1=[0.2500000,0.7857143,0.07142857,0.07142857,0.07142857,0.1005964,0.1005964,0.1005964,0.3994034,0.3994034,0.3994034]
p2=[0.2500000,0.07142857,0.7857143,0.07142857,0.07142857,0.1005964,0.3994034,0.3994034,0.1005964,0.1005964,0.3994034]
p3=[0.2500000,0.07142857,0.07142857,0.7857143,0.07142857,0.3994034,0.1005964,0.3994034,0.1005964,0.3994034,0.1005964]
p4=[0.2500000,0.07142856,0.07142856,0.07142856,0.785714290000000,0.3994038,0.3994038,0.1005968,0.3994038,0.1005968,0.1005968]
!PG = p1*tetra(1,:)+p2*tetra(2,:)+p3*tetra(3,:)+p4*tetra(4,:)
Xg = p1*tetra(1,1)+p2*tetra(1,2)+p3*tetra(1,3)+p4*tetra(1,4)
Yg = p1*tetra(2,1)+p2*tetra(2,2)+p3*tetra(2,3)+p4*tetra(2,4)
Zg = p1*tetra(3,1)+p2*tetra(3,2)+p3*tetra(3,3)+p4*tetra(3,4)
WWG = wg*Vol*6;
end subroutine fun_tetra_my_quad_N11_D4_dim