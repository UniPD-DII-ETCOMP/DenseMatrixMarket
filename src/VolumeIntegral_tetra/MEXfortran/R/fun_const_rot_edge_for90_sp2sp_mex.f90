!---- never change ------

#include "fintrf.h"      

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    implicit none

! mwPointer mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
!---- end never change ------
	  
	  
! mwSize stuff for mexing
      mwSize mo,no,siz
	  
! mwPointer stuff for mexing
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
      mwPointer N_edge_pr,Matrix_P0_pr,N_node_pr,N_volu_pr,N_thread_pr,VE_pr,VE_ind_pr
	  mwPointer VP_pr,coef_mag_pr,Matrix_Gamma_pr, ne_pr, nnzRmaybe_pr
	  mwPointer nv_max_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN

!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
	  
! fortran subroutine arguments
      real*8,allocatable,dimension(:,:) :: Matrix_P0, Matrix_Gamma, VE_r, VE_ind_r, VP_r
      integer,allocatable,dimension(:,:) ::  VE, VE_ind, VP
      real*8,allocatable,dimension(:) :: coef_mag
     ! integer,allocatable,dimension(:) ::   
      integer*8 ne, N_node, N_thread, N_volu, N_edge
	  integer*8 esiz,ehtpsiz,elocsiz2,nv_max,nnzRmaybe
      real*8 N_node_r, N_thread_r, N_volu_r, N_edge_r, nv_max_r,nnzRmaybe_r
	  double precision log2
      character*80 msg
      logical debu
       
      debu = .true. ! .true. o .false. 
	  if(debu) open(unit=66,file='logRCurlEdSp2Sp.txt',status='unknown')


! check for proper number of arguments. 
      if (nrhs .ne. 11) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_const_rot_edge_for90_sp2sp:nInput', &
                                '11 input argument required.')
      elseif (nlhs .ne. 2) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_const_rot_edge_for90_sp2sp:nOutput', &
                                '2 output argument required.')
      endif
      if(debu) write(66,*) 'arguments checked'


!    Check to see inputs are numeric.
	  do ii = 1,11
        if (mxIsNumeric(prhs(ii)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:fun_const_rot_edge_for90_sp2sp:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
	  if(debu) write(66,*) 'check num'		
	  
	  
!     Check that input #11 is integer and fetch it
      m = mxGetM(prhs(11))
      n = mxGetN(prhs(11))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 11 must be scalar.')
      endif	  
      siz = m*n
      nnzRmaybe_pr = mxGetPr(prhs(11))
      call mxCopyPtrToReal8(nnzRmaybe_pr, nnzRmaybe_r, siz) ! da double precision a reale
      nnzRmaybe=int(nnzRmaybe_r) ! da reale a intero
	  if(debu) write(66,*) 'input 11' 
	  
!     Check that input #3 is integer and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 3 must be scalar.')
      endif	  
      siz = m*n
      N_node_pr = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(N_node_pr, N_node_r, siz) ! da double precision a reale
      N_node=int(N_node_r) ! da reale a intero
	  if(debu) write(66,*) 'input 3' 
      
	  
!     Check that input #2 is real matrix and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. N_node .or. n .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 2 must be N_nodex3')
      endif	  
      siz = m*n
      Matrix_P0_pr = mxGetPr(prhs(2))
	  allocate(Matrix_P0(N_node,3))
      call mxCopyPtrToReal8(Matrix_P0_pr, Matrix_P0, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 2'
	     	
	  
!     Check that input #5 is integer and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 5 must be scalar.')
      endif	  
      siz = m*n
      N_thread_pr = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(N_thread_pr, N_thread_r, siz) ! da double precision a reale
      N_thread=int(N_thread_r) ! da reale a intero
      if(debu) write(66,*) 'input 5' 	  	  

       
      
!     Check that input #4 is integer and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 4 must be scalar.')
      endif	  
      siz = m*n
      N_volu_pr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(N_volu_pr, N_volu_r, siz) ! da double precision a reale
      N_volu=int(N_volu_r) ! da reale a intero
      if(debu) write(66,*) 'input 4'       
      
      
!     Check that input #1 is integer and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be scalar.')
      endif	  
      siz = m*n
      N_edge_pr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(N_edge_pr, N_edge_r, siz) ! da double precision a reale
      N_edge=int(N_edge_r) ! da reale a intero
      if(debu) write(66,*) 'input 1'             
      
      
!     Check that input #6 is integer matrix and fetch it
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. 6 .or. n .ne. N_volu) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 6 must be 6 x N_volu')
      endif	  
      siz = m*n
      VE_pr = mxGetPr(prhs(6))
	  allocate(VE_r(6,N_volu))
      allocate(VE(6,N_volu))
      call mxCopyPtrToReal8(VE_pr, VE_r, siz) ! da double precision a reale
	  VE = int(VE_r)
	  if(debu) write(66,*) 'input 6'       
      
      
!     Check that input #8 is integer matrix and fetch it
      m = mxGetM(prhs(8))
      n = mxGetN(prhs(8))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. 4 .or. n .ne. N_volu) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 8 must be 4 x N_volu')
      endif	  
      siz = m*n
      VP_pr = mxGetPr(prhs(8))
	  allocate(VP_r(4,N_volu))
      allocate(VP(4,N_volu))
      call mxCopyPtrToReal8(VP_pr, VP_r, siz) ! da double precision a reale
	  VP = int(VP_r)
	  if(debu) write(66,*) 'input 8'        
      
      
!     Check that input #7 is integer matrix and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. 6 .or. n .ne. N_volu) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 7 must be 6 x N_volu')
      endif	  
      siz = m*n
      VE_ind_pr = mxGetPr(prhs(7))
	  allocate(VE_ind_r(6,N_volu))
      allocate(VE_ind(6,N_volu))
      call mxCopyPtrToReal8(VE_ind_pr, VE_ind_r, siz) ! da double precision a reale
	  VE_ind = int(VE_ind_r)
	  if(debu) write(66,*) 'input 7'       
        
		
!     Check that input #9 is real vector
      m = mxGetM(prhs(9))
      n = mxGetN(prhs(9))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. N_volu .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 9 must be N_volu x 1.')
      endif	  
      siz = m*n
	  if(debu) write(66,*) siz 
      coef_mag_pr = mxGetPr(prhs(9))
	  allocate(coef_mag(N_volu))
      call mxCopyPtrToReal8(coef_mag_pr, coef_mag, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 9'       
	  if(debu) write(66,*) size(coef_mag) 

!     Check that input #10 is integer and fetch it
      m = mxGetM(prhs(10))
      n = mxGetN(prhs(10))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be scalar.')
      endif	  
      siz = m*n
      nv_max_pr = mxGetPr(prhs(10))
      call mxCopyPtrToReal8(nv_max_pr, nv_max_r, siz) ! da double precision a reale
      nv_max=int(nv_max_r) ! da reale a intero
      if(debu) write(66,*) 'input 10' 



! call the computational subroutine.
	  if(debu) write(66,*) 'call!' 
	  ne=0
	  esiz=max(nnzRmaybe,10000)
	  ehtpsiz=int(log2(esiz))+1
	  elocsiz2=2**ehtpsiz
	  allocate(Matrix_Gamma(3,elocsiz2))
	  Matrix_Gamma(1:3,1:elocsiz2)=0.0d0 
      call fun_const_rot_edge_for90_sp2sp(N_edge,Matrix_P0,N_node,N_volu,N_thread, & 
      VE,VE_ind,VP,coef_mag,Matrix_Gamma,ne,elocsiz2,nv_max,nnzRmaybe)

       if(debu) write(66,*) 'matrix created'
      
      deallocate(Matrix_P0,VE,VE_r,VP,VP_r,VE_ind,VE_ind_r,coef_mag)
	  if(debu) write(66,*) 'deallocate'

! Create a matrix for the return argument 1
      mo=3
      no=elocsiz2
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
    if(debu) write(66,*) 'I am here1'
! Load the output 1 into a MATLAB array.
      Matrix_Gamma_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(Matrix_Gamma, Matrix_Gamma_pr, siz)
	  
	 
	  
! Create a matrix for the return argument 2
      mo=1
      no=1
	  ComplexFlag = 0
      plhs(2) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
    if(debu) write(66,*) 'I am here2'
! Load the output 1 into a MATLAB array.
      ne_pr = mxGetPr(plhs(2))
      siz=mo*no
      call mxCopyReal8ToPtr(real(ne,8), ne_pr, siz)	  
	  
	  
	  
	  

      if(debu) write(66,*) 'matrix converted for matlab'

      deallocate(Matrix_Gamma)
      
	  
      if(debu) write(66,*) 'closing lo file, bye bye'
      if(debu) close(66)
      return
      end
	  
	  !--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
