subroutine fun_const_rot_edge_for90_sp2sp(N_edge,Matrix_P0,N_node,N_volu,N_thread,VE, &
VE_ind,VP,coef_mag,Rrcv,ne,elocsizout,nv_max,nnzRmaybe)
implicit none
!input e output
real*8 Matrix_P0(N_node,3), coef_mag(N_volu)
integer*8 N_node, N_volu, N_thread, VE(6,N_volu), VE_ind(6,N_volu), VP(4,N_volu), N_edge
!intirior variables
real*8 pv(4,3), pvt(3,4), norm(3,4),area(4),vol
real*8 w1(3),w2(3),Rloc
integer*8 ii, jj, hh, kk, ind1, ind2, nv_max
!------------------------------------ da cambiare in base ai punti di Gauss ------------------------------------------------------
real*8  Xgg(11), Ygg(11), Zgg(11), Wgg(11)
!---------------------------------------------hashtable
integer*8,allocatable,dimension(:) :: eht,elinks
integer*8,allocatable,dimension(:,:) :: g	
!real*8,allocatable,dimension(:,:) ::
integer*8 elocsizout, nnzRmaybe
real*8 Rrcv(3,elocsizout)	
real*8,allocatable,dimension(:) :: R_sparse
integer*8 esiz,ehtpsiz,elocsiz
integer*8 htsrcedge 
external htsrcedge
double precision log2
integer*8 ind, ied,ne
!real*8  X_ed_a(24,n_EV_max,N_edge), Y_ed_a(24,n_EV_max,N_edge), Z_ed_a(24,n_EV_max,N_edge), Wg_ed_a(24,n_EV_max,N_edge)
! --------------------------------------------------------------------------------------------------------------------------------
esiz=max(nnzRmaybe,10000)
ehtpsiz=int(log2(esiz))+1
elocsiz=2**ehtpsiz
allocate(eht(elocsiz),elinks(elocsiz))
allocate(g(2,elocsiz))	
allocate(R_sparse(elocsiz))
ne=0
eht(:)=0
elinks(:)=0
R_sparse(:) =0.0d0
g(:,:)=0
ind=0
ied=0
!----------------------------------------------------------------------------------------
!!call omp_set_num_threads(N_thread)
do ii=1,N_volu
   pv(1:4,1:3)=Matrix_P0(VP(:,ii),1:3) ! punti del tetraedo
   pvt(1:3,1)=pv(1,1:3)! traspongo il tetraedro
   pvt(1:3,2)=pv(2,1:3)
   pvt(1:3,3)=pv(3,1:3)
   pvt(1:3,4)=pv(4,1:3)
   call tetareavol(pvt,norm,area,vol) !trovo le grandezze del tetraedro
   !call fun_tetra_my_quad_N11_D4_2(wg,p1,p2,p3,p4,pv,vol,Xgg,Ygg,Zgg,Wgg) ! trovo i punti di Gauss del tetraedro
   do jj=1,6
      ind1=abs(VE(jj,ii))
      do hh = 1,6
	     ind2=abs(VE(hh,ii))
			call fun_whitney_tetra_rot_edge(norm,area,vol,pvt,VE_ind(jj,ii),w1)!trovo le whitney del lato hh esimo nel punto di gauss kk			
			call fun_whitney_tetra_rot_edge(norm,area,vol,pvt,VE_ind(hh,ii),w2)!trovo le whitney del lato hh esimo nel punto di gauss kk
			Rloc=coef_mag(ii)*vol*sign(1.0d0,real(VE(jj,ii),8))*sign(1.0d0,real(VE(hh,ii),8))*dot_product(w1,w2)
			!------------------------------------------------------------------------
			ind=htsrcedge([ind1,ind2],ehtpsiz,g,eht,elinks) ! cerco la coppia ind1,ind2
            ied=ind
			if(ind.eq.0)then ! controllo se ho già trovato la coppia, se non l'ho trovata ...
                ne=ne+1
                ied=ne
                call htinsedge (ied,[ind1,ind2],ehtpsiz,g,eht,elinks ) ! inserisco la coppia nella mia hash table
				R_sparse(ied)=R_sparse(ied)+Rloc ! aggungo il valore alla matrice R
			elseif(ind.ne.0)then
				R_sparse(ied)=R_sparse(ied)+Rloc	! aggungo il valore alla matrice R	 			
			end if 
			!------------------------------------------------------------------------
      enddo	  
   enddo
enddo 
!allocate(Rrcv(3,ne))
do hh = 1,ne
	Rrcv(1,hh)=real(g(1,hh),8)
	Rrcv(2,hh)=real(g(2,hh),8)
	ind=htsrcedge([g(1,hh),g(2,hh)],ehtpsiz,g,eht,elinks)
	if(ind.ne.0)then
		Rrcv(3,hh)=R_sparse(ind)
	endif
enddo
end subroutine fun_const_rot_edge_for90_sp2sp

subroutine fun_tetra_my_quad_N15_D5_2(wg,p1,p2,p3,p4,tetra,Vol,Xg,Yg,Zg,WWg)
implicit none
real*8 tetra(4,3), Vol, Xg(15), Yg(15), Zg(15), WWg(15)
real*8 wg(15), p1(15), p2(15), p3(15), p4(15)
Xg = p1*tetra(1,1)+p2*tetra(2,1)+p3*tetra(3,1)+p4*tetra(4,1)
Yg = p1*tetra(1,2)+p2*tetra(2,2)+p3*tetra(3,2)+p4*tetra(4,2)
Zg = p1*tetra(1,3)+p2*tetra(2,3)+p3*tetra(3,3)+p4*tetra(4,3)
WWG = wg*Vol*6;
end subroutine fun_tetra_my_quad_N15_D5_2

subroutine fun_tetra_my_quad_N24_D6_2(wg,p1,p2,p3,p4,tetra,Vol,Xg,Yg,Zg,WWg)
implicit none
real*8 tetra(4,3), Vol, Xg(24), Yg(24), Zg(24), WWg(24)
real*8 wg(24), p1(24), p2(24), p3(24), p4(24)
Xg = p1*tetra(1,1)+p2*tetra(2,1)+p3*tetra(3,1)+p4*tetra(4,1)
Yg = p1*tetra(1,2)+p2*tetra(2,2)+p3*tetra(3,2)+p4*tetra(4,2)
Zg = p1*tetra(1,3)+p2*tetra(2,3)+p3*tetra(3,3)+p4*tetra(4,3)
WWG = wg*Vol*6;
end subroutine fun_tetra_my_quad_N24_D6_2

subroutine fun_tetra_my_quad_N31_D7_2(wg,p1,p2,p3,p4,tetra,Vol,Xg,Yg,Zg,WWg)
implicit none
real*8 tetra(4,3), Vol, Xg(31), Yg(31), Zg(31), WWg(31)
real*8 wg(31), p1(31), p2(31), p3(31), p4(31)
Xg = p1*tetra(1,1)+p2*tetra(2,1)+p3*tetra(3,1)+p4*tetra(4,1)
Yg = p1*tetra(1,2)+p2*tetra(2,2)+p3*tetra(3,2)+p4*tetra(4,2)
Zg = p1*tetra(1,3)+p2*tetra(2,3)+p3*tetra(3,3)+p4*tetra(4,3)
WWG = wg*Vol*6;
end subroutine fun_tetra_my_quad_N31_D7_2

subroutine fun_tetra_my_quad_N45_D8_2(wg,p1,p2,p3,p4,tetra,Vol,Xg,Yg,Zg,WWg)
implicit none
real*8 tetra(4,3), Vol, Xg(45), Yg(45), Zg(45), WWg(45)
real*8 wg(45), p1(45), p2(45), p3(45), p4(45)
Xg = p1*tetra(1,1)+p2*tetra(2,1)+p3*tetra(3,1)+p4*tetra(4,1)
Yg = p1*tetra(1,2)+p2*tetra(2,2)+p3*tetra(3,2)+p4*tetra(4,2)
Zg = p1*tetra(1,3)+p2*tetra(2,3)+p3*tetra(3,3)+p4*tetra(4,3)
WWG = wg*Vol*6;
end subroutine fun_tetra_my_quad_N45_D8_2

subroutine fun_tetra_my_quad_N8_D4_2(wg,p1,p2,p3,p4,tetra,Vol,Xg,Yg,Zg,WWg)
implicit none
real*8 tetra(4,3), Vol, Xg(8), Yg(8), Zg(8), WWg(8)
real*8 wg(8), p1(8), p2(8), p3(8), p4(8)
Xg = p1*tetra(1,1)+p2*tetra(2,1)+p3*tetra(3,1)+p4*tetra(4,1)
Yg = p1*tetra(1,2)+p2*tetra(2,2)+p3*tetra(3,2)+p4*tetra(4,2)
Zg = p1*tetra(1,3)+p2*tetra(2,3)+p3*tetra(3,3)+p4*tetra(4,3)
WWG = wg*Vol*6;
end subroutine fun_tetra_my_quad_N8_D4_2


subroutine fun_tetra_my_quad_N11_D4_2(wg,p1,p2,p3,p4,tetra,Vol,Xg,Yg,Zg,WWg)
implicit none
real*8 tetra(4,3), Vol, Xg(11), Yg(11), Zg(11), WWg(11)
real*8 wg(11), p1(11), p2(11), p3(11), p4(11)
! 3D Gauss-Legendre weights for N = 11, D = 4;
! wg=[-0.01315556,0.007622222,0.007622222,0.007622222,0.007622222,0.02488889,0.02488889,0.02488889,0.02488889,0.02488889,0.02488889]
! p1=[0.2500000,0.7857143,0.07142857,0.07142857,0.07142857,0.1005964,0.1005964,0.1005964,0.3994034,0.3994034,0.3994034]
! p2=[0.2500000,0.07142857,0.7857143,0.07142857,0.07142857,0.1005964,0.3994034,0.3994034,0.1005964,0.1005964,0.3994034]
! p3=[0.2500000,0.07142857,0.07142857,0.7857143,0.07142857,0.3994034,0.1005964,0.3994034,0.1005964,0.3994034,0.1005964]
! p4=[0.2500000,0.07142856,0.07142856,0.07142856,0.785714290000000,0.3994038,0.3994038,0.1005968,0.3994038,0.1005968,0.1005968]
!PG = p1*tetra(1,:)+p2*tetra(2,:)+p3*tetra(3,:)+p4*tetra(4,:)
Xg = p1*tetra(1,1)+p2*tetra(2,1)+p3*tetra(3,1)+p4*tetra(4,1)
Yg = p1*tetra(1,2)+p2*tetra(2,2)+p3*tetra(3,2)+p4*tetra(4,2)
Zg = p1*tetra(1,3)+p2*tetra(2,3)+p3*tetra(3,3)+p4*tetra(4,3)
WWG = wg*Vol*6;
end subroutine fun_tetra_my_quad_N11_D4_2

subroutine fun_whitney_tetra_rot_edge_ori(v,rotwi)
implicit none
integer*8 locedg(6,2)
data locedg /1,2,3,4,4,4,2,3,1,1,2,3/
double precision v(3,4),gradlambda(3,4),norm(3,4),area(4),vol,rotwi(6,3)
integer*8 k,nd1loc,nd2loc
call x2gradlambda(v,gradlambda,norm,area,vol)
do k=1,6
nd1loc=locedg(k,1)
nd2loc=locedg(k,2)
call vecprod(gradlambda(:,nd1loc),gradlambda(:,nd2loc),rotwi(k,1:3))
rotwi(k,1:3)=2.0d0*rotwi(k,1:3)
enddo
end subroutine fun_whitney_tetra_rot_edge_ori

subroutine fun_whitney_tetra_rot_edge(norm,area,vol,v,k,rotwi)
implicit none
integer*8 locedg(6,2)
data locedg /1,2,3,4,4,4,2,3,1,1,2,3/
double precision v(3,4),gradlambda(3,4),norm(3,4),area(4),vol,rotwi(3)
integer*8 k,nd1loc,nd2loc
call x2gradlambda(v,gradlambda,norm,area,vol)
!do k=1,6
nd1loc=locedg(k,1)
nd2loc=locedg(k,2)
call vecprod(gradlambda(:,nd1loc),gradlambda(:,nd2loc),rotwi)
rotwi=2.0d0*rotwi
!enddo
end subroutine fun_whitney_tetra_rot_edge

subroutine vecprod(v1,v2,vp)
implicit none
double precision v1(3),v2(3),vp(3)
vp(1)=   v1(2)*v2(3)-v1(3)*v2(2)
vp(2)= -(v1(1)*v2(3)-v1(3)*v2(1))
vp(3)=   v1(1)*v2(2)-v1(2)*v2(1)
end subroutine   

subroutine fun_whitney_tetra_edge(norm,area,vol,v,lambda,k,wi)
implicit none
integer*8 locedg(6,2)
data locedg /1,2,3,4,4,4,2,3,1,1,2,3/
double precision v(3,4),gradlambda(3,4),norm(3,4),area(4),vol,lambda(4),wi(3)
integer*8 k,nd1loc,nd2loc
call x2gradlambda(v,gradlambda,norm,area,vol)
!do k=1,6
nd1loc=locedg(k,1)
nd2loc=locedg(k,2)
!wi(k,1:3)=lambda(nd1loc)*gradlambda(:,nd2loc)-lambda(nd2loc)*gradlambda(:,nd1loc)
wi(1:3)=lambda(nd1loc)*gradlambda(:,nd2loc)-lambda(nd2loc)*gradlambda(:,nd1loc)
!enddo
end subroutine fun_whitney_tetra_edge

subroutine x2gradlambda(v,gradlambda,norm,area,vol)
implicit none
double precision v(3,4),gradlambda(3,4)
double precision norm(3,4),area(4),vol,vol3
integer*8 i
!call tetareavol(v,norm,area,vol)
vol3=3.0d0*vol
do i=1,4
gradlambda(:,i)=-norm(:,i)*area(i)/vol3
enddo
end subroutine x2gradlambda

subroutine tetareavol(v,norm,area,vol)
implicit none
double precision v(3,4),norm(3,4),area(4),vol
double precision r(3),s(3),t(3),vp(3),sp
integer*8 i
integer*8 permut(4,4) ! face and nodes formin face
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


subroutine fun_my_vol_tetra(tetra,V)
implicit none
real*8 tetra(4,3), v1(3), v2(3), v3(3), V
v1 =  tetra(2,1:3) - tetra(1,1:3)
v2 =  tetra(3,1:3) - tetra(1,1:3)
v3 =  tetra(4,1:3) - tetra(1,1:3)
V = (v1(1)*(v2(2)*v3(3)-v2(3)*v3(2))-v1(2)*(v2(1)*v3(3)-v2(3)*v3(1))+v1(3)*(v2(1)*v3(2)-v2(2)*v3(1)))/6
end subroutine fun_my_vol_tetra

subroutine fun_my_tetra_quad__fv(tetra,V,X,Y,Z,W)
implicit none
real*8 tetra(4,3), V, X(8), Y(8), Z(8), W(8), F1(3), F2(3), F3(3), F4(3)
F1=(tetra(1,1:3)+tetra(2,1:3)+tetra(3,1:3))/3.0!sum(tetra(1:3,1:3))/3;
F2=(tetra(2,1:3)+tetra(3,1:3)+tetra(4,1:3))/3.0!sum(tetra(2:4,1:3))/3;
F3=(tetra(1,1:3)+tetra(3,1:3)+tetra(4,1:3))/3.0!sum(tetra([1,3,4],:))/3;
F4=(tetra(1,1:3)+tetra(2,1:3)+tetra(4,1:3))/3.0!sum(tetra([1,2,4],:))/3;
X = [tetra(1:4,1),F1(1),F2(1),F3(1),F4(1)]
Y = [tetra(1:4,2),F1(2),F2(2),F3(2),F4(2)]
Z = [tetra(1:4,3),F1(3),F2(3),F3(3),F4(3)]
W  = [V/40,V/40,V/40,V/40,9*V/40,9*V/40,9*V/40,9*V/40] ![V/40*ones(4,1);9*V/40*ones(4,1)]
end subroutine fun_my_tetra_quad__fv

subroutine fun_shape_tetra(tetra,P,v,NP,Vol,w) 
implicit none
real*8 tetra(4,3), Vol, P(NP,3), w(NP,3), J
integer*8 NP, ii, v
J =  Vol*6
w(1:NP,1:3) = 0.0 
do ii = 1,NP
w(ii,:) = 2*(P(ii,1:3)-tetra(v,1:3))/J
! w2(ii,:) = 2*(P(ii,:)-P0(2,:))/J
! w3(ii,:) = 2*(P(ii,:)-P0(3,:))/J
! w4(ii,:) = 2*(P(ii,:)-P0(4,:))/J
enddo
end subroutine fun_shape_tetra

subroutine fun_my_norm(a,c) 
implicit none
real(kind=8), dimension(3), intent(in) :: a
real(kind=8)    :: c
c=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
end subroutine fun_my_norm

subroutine fun_my_cross(a, b, c)
implicit none
real*8, DIMENSION(3) :: c
real*8, DIMENSION(3), INTENT(IN) :: a, b
c(1) = a(2) * b(3) - a(3) * b(2)
c(2) = a(3) * b(1) - a(1) * b(3)
c(3) = a(1) * b(2) - a(2) * b(1)
end subroutine fun_my_cross

!--------------------------------------------------------------------------------
subroutine htinsedge ( ind, ab, htedgpowersiz,edg, edght, edglinks )
  implicit none
  integer*8 ab(2)
  integer*8 ind,htedgpowersiz
  integer*8 k
  integer*8 leni
  integer*8 edg(2,*),edght(*),edglinks(*)
  integer*8 fhash
  external fhash
  leni=2
  k=fhash(ab, leni, htedgpowersiz, htedgpowersiz)
  edg(1:2,ind)=ab(1:2)
  edglinks(ind) = edght(k)
  edght(k) = ind
!  write(6,*) 'ins:',ind,ab,k,htedgpowersiz
end
!--------------------------------------------------------------------------------
function htsrcedge ( ab, htedgpowersiz,edg, edght, edglinks )
  implicit none
  integer*8 ab(2),htedgpowersiz
  integer*8 htsrcedge
  integer*8 ind
  integer*8 k
  integer*8 leni
  integer*8 edg(2,*),edght(*),edglinks(*)
  integer*8 fhash
  external fhash
  leni=2
  k=fhash(ab, leni, htedgpowersiz, htedgpowersiz)
!  write(6,*) 'src:',ab,k,htedgpowersiz
  ind = edght(k)
  do
    if ( ind == 0 ) then
      exit
    end if
    if ( edg(1,ind) == ab(1) .and. edg(2,ind) == ab(2) )then
            exit
    endif
    ind = edglinks(ind)
  end do
  htsrcedge = ind
end
!--------------------------------------------------------------------------------
integer*8 function fhash(arr,siz,seed,cut)
implicit none
integer*8 arr(*),siz,seed,cut
integer*8 i,shft
!FNV hash
fhash = 1966136261
do i=1,siz
   fhash=ieor(fhash*16777619,arr(i))
enddo
shft=ishft(1,cut)-1
fhash=iand(fhash,shft)+1
!write(6,*) 'fhash:',arr(1:siz),fhash,cut
end function
!--------------------------------------------------------------------------------
double precision function log2(n)
implicit none
integer*8 n
log2=log(dble(n))/log(dble(2))
end function
!--------------------------------------------------------------------------------
