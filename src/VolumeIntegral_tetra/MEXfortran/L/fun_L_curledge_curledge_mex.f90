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
      mwPointer N_volu_pr,N_node_pr,N_edge_pr,nv_max_pr,EV_pr,Matrix_P0_pr 
      mwPointer VP_pr,EV_ind_pr,Neps_pr,n_EV_pr,N_thread_pr,Matrix_L_pr
	  mwPointer Nhh_pr,Nkk_pr,hhout_pr,kkout_pr,bared_pr,thr_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN

!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
	  
! fortran subroutine arguments
      real*8,allocatable,dimension(:,:) :: Matrix_P0, Matrix_L
	  real*8,allocatable,dimension(:,:) :: EV_r, EV_ind_r, VP_r, bared
      integer*8,allocatable,dimension(:,:) ::  EV, EV_ind, VP
      integer*8,allocatable,dimension(:) :: n_EV, hhout,kkout
	  real*8,allocatable,dimension(:) :: n_EV_r, hhout_r,kkout_r
     ! integer,allocatable,dimension(:) ::   
      integer*8 N_node, N_thread, N_volu, N_edge, nv_max, Nhh,Nkk
      real*8 N_node_r, N_thread_r, N_volu_r, N_edge_r, nv_max_r, Neps
	  real*8 Nhh_r,Nkk_r, thr
      character*80 msg
      logical debu
       
      debu = .true. ! .true. o .false. per attivare o disattivare il debug
	  !if(debu) open(unit=66,file='logL_curl_ed.txt',status='unknown')

! check for proper number of arguments. 
      if (nrhs .ne. 17) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:nInput', &
                                '17 input argument required.')
      elseif (nlhs .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:nOutput', &
                                '1 output argument required.')
      endif
      !if(debu) write(66,*) 'arguments checked'

!    Check to see inputs are numeric.
	  do ii = 1,17
        if (mxIsNumeric(prhs(ii)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
	  !if(debu) write(66,*) 'check num'		
	  
!     Check that input #1 is integer and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 1 must be scalar N_volu.')
      endif	  
      siz = m*n
      N_volu_pr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(N_volu_pr, N_volu_r, siz) ! da double precision a reale
      N_volu=int(N_volu_r,8) ! da reale a intero
	  !if(debu) write(66,*) 'input 1' 
      
!     Check that input #2 is integer and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 2 must be scalar N_node.')
      endif	  
      siz = m*n
      N_node_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(N_node_pr, N_node_r, siz) ! da double precision a reale
      N_node=int(N_node_r,8) ! da reale a intero
	  !if(debu) write(66,*) 'input 2' 	  
	  
!     Check that input #3 is integer and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 3 must be scalar N_edge.')
      endif	  
      siz = m*n
      N_edge_pr = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(N_edge_pr, N_edge_r, siz) ! da double precision a reale
      N_edge=int(N_edge_r,8) ! da reale a intero
	  !if(debu) write(66,*) 'input 3' 	

!     Check that input #4 is integer and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 4 must be scalar nv_max.')
      endif	  
      siz = m*n
      nv_max_pr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(nv_max_pr, nv_max_r, siz) ! da double precision a reale
      nv_max=int(nv_max_r,8) ! da reale a intero
	  !if(debu) write(66,*) 'input 4' 	

!     Check that input #5 is integer matrix and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
	  !if(debu) write(66,*) m
	  !if(debu) write(66,*) n
      if(m .ne. nv_max+1 .or. n .ne. N_edge) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 5 must be nv_max+1 x N_edge EV')
      endif	  
      siz = m*n
      EV_pr = mxGetPr(prhs(5))
	  allocate(EV_r(nv_max+1,N_edge))
      allocate(EV(nv_max+1,N_edge))
      call mxCopyPtrToReal8(EV_pr, EV_r, siz) ! da double precision a reale
	  EV = int(EV_r,8)
	  !if(debu) write(66,*) 'input 5'
	  deallocate(EV_r)
	  
!     Check that input #6 is real matrix and fetch it
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
	  !if(debu) write(66,*) m
	  !if(debu) write(66,*) n
      if(m .ne. 3 .or. n .ne. N_node) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 6 must be 3xN_node Matrix_P0')
      endif	  
      siz = m*n
      Matrix_P0_pr = mxGetPr(prhs(6))
	  allocate(Matrix_P0(3,N_node))
      call mxCopyPtrToReal8(Matrix_P0_pr, Matrix_P0, siz) ! da double precision a reale
	  !if(debu) write(66,*) 'input 6'
	     	
!     Check that input #7 is integer matrix and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
	  !if(debu) write(66,*) m
	  !if(debu) write(66,*) n
      if(m .ne. 4 .or. n .ne. N_volu) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 7 must be 4 x N_volu VP')
      endif	  
      siz = m*n
      VP_pr = mxGetPr(prhs(7))
	  allocate(VP_r(4,N_volu))
      allocate(VP(4,N_volu))
      call mxCopyPtrToReal8(VP_pr, VP_r, siz) ! da double precision a reale
	  VP = int(VP_r,8)
	  !if(debu) write(66,*) 'input 7'			
	  deallocate(VP_r)
			
!     Check that input #8 is integer matrix and fetch it
      m = mxGetM(prhs(8))
      n = mxGetN(prhs(8))
	  !if(debu) write(66,*) m
	  !if(debu) write(66,*) n
      if(m .ne. nv_max+1 .or. n .ne. N_edge) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 8 must be nv_max+1 x N_edge EV_ind')
      endif	  
      siz = m*n
      EV_ind_pr = mxGetPr(prhs(8))
	  allocate(EV_ind_r(nv_max+1,N_edge))
      allocate(EV_ind(nv_max+1,N_edge))
      call mxCopyPtrToReal8(EV_ind_pr, EV_ind_r, siz) ! da double precision a reale
	  EV_ind = int(EV_ind_r,8)
	  !if(debu) write(66,*) 'input 8'			
	  deallocate(EV_ind_r)	
			
!     Check that input #9 is integer and fetch it
      m = mxGetM(prhs(9))
      n = mxGetN(prhs(9))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 9 must be scalar Neps.')
      endif	  
      siz = m*n
      Neps_pr = mxGetPr(prhs(9))
      call mxCopyPtrToReal8(Neps_pr, Neps, siz) ! da double precision a reale
	  !if(debu) write(66,*) 'input 9' 			
			
			
!     Check that input #10 is integer matrix and fetch it
      m = mxGetM(prhs(10))
      n = mxGetN(prhs(10))
	  !if(debu) write(66,*) m
	  !if(debu) write(66,*) n
      if(m .ne. 1 .or. n .ne. N_edge) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 10 must be 1 x N_edge n_EV')
      endif	  
      siz = m*n
      n_EV_pr = mxGetPr(prhs(10))
	  allocate(n_EV(N_edge))
      allocate(n_EV_r(N_edge))
      call mxCopyPtrToReal8(n_EV_pr, n_EV_r, siz) ! da double precision a reale
	  n_EV = int(n_EV_r,8)
	  !if(debu) write(66,*) 'input 10'			
      deallocate(n_EV_r)			
		
!     Check that input #11 is integer and fetch it
      m = mxGetM(prhs(11))
      n = mxGetN(prhs(11))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 11 must be scalar N_thread.')
      endif	  
      siz = m*n
      N_thread_pr = mxGetPr(prhs(11))
      call mxCopyPtrToReal8(N_thread_pr, N_thread_r, siz) ! da double precision a reale
	  N_thread = int(N_thread_r,8)
	  !if(debu) write(66,*) 'input 11' 	

!     Check that input #12 is integer and fetch it
      m = mxGetM(prhs(12))
      n = mxGetN(prhs(12))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 12 must be scalar Nkk.')
      endif	  
      siz = m*n
      Nkk_pr = mxGetPr(prhs(12))
      call mxCopyPtrToReal8(Nkk_pr, Nkk_r, siz) ! da double precision a reale
      Nkk=int(Nkk_r,8) ! da reale a intero
	  !if(debu) write(66,*) 'input 12' 
	  
!     Check that input #13 is integer and fetch it
      m = mxGetM(prhs(13))
      n = mxGetN(prhs(13))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 13 must be scalar Nhh.')
      endif	  
      siz = m*n
      Nhh_pr = mxGetPr(prhs(13))
      call mxCopyPtrToReal8(Nhh_pr, Nhh_r, siz) ! da double precision a reale
      Nhh=int(Nhh_r,8) ! da reale a intero
	  !if(debu) write(66,*) 'input 13' 	  

!     Check that input #14 is integer matrix and fetch it
      m = mxGetM(prhs(14))
      n = mxGetN(prhs(14))
	  !if(debu) write(66,*) m
	  !if(debu) write(66,*) n
      if(m .ne. 1 .or. n .ne. Nhh) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 14 must be 1 x Nhh hhout')
      endif	  
      siz = m*n
      hhout_pr = mxGetPr(prhs(14))
	  allocate(hhout(Nhh))
      allocate(hhout_r(Nhh))
      call mxCopyPtrToReal8(hhout_pr, hhout_r, siz) ! da double precision a reale
	  hhout = int(hhout_r,8)
	  !if(debu) write(66,*) 'input 14'			
      deallocate(hhout_r)

!     Check that input #15 is integer matrix and fetch it
      m = mxGetM(prhs(15))
      n = mxGetN(prhs(15))
	  !if(debu) write(66,*) m
	  !if(debu) write(66,*) n
      if(m .ne. 1 .or. n .ne. Nkk) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 15 must be 1 x Nkk kkout')
      endif	  
      siz = m*n
      kkout_pr = mxGetPr(prhs(15))
	  allocate(kkout(Nkk))
      allocate(kkout_r(Nkk))
      call mxCopyPtrToReal8(kkout_pr, kkout_r, siz) ! da double precision a reale
	  kkout = int(kkout_r,8)
	  !if(debu) write(66,*) 'input 15'			
      deallocate(kkout_r)


!     Check that input #16 is real matrix and fetch it
      m = mxGetM(prhs(16))
      n = mxGetN(prhs(16))
	  !if(debu) write(66,*) m
	  !if(debu) write(66,*) n
      if(m .ne. 3 .or. n .ne. N_edge) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 16 must be 3xN_edge bared')
      endif	  
      siz = m*n
      bared_pr = mxGetPr(prhs(16))
	  allocate(bared(3,N_edge))
      call mxCopyPtrToReal8(bared_pr, bared, siz) ! da double precision a reale
	  !if(debu) write(66,*) 'input 16'

!     Check that input #17 is integer and fetch it
      m = mxGetM(prhs(17))
      n = mxGetN(prhs(17))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:fun_L_curledge_curledge:NonRowVector', &
                                'Input 17 must be scalar thr.')
      endif	  
      siz = m*n
      thr_pr = mxGetPr(prhs(17))
      call mxCopyPtrToReal8(thr_pr, thr, siz) ! da double precision a reale
	  !if(debu) write(66,*) 'input 17' 

! call the computational subroutine.
      allocate(Matrix_L(Nhh,Nkk))
	  Matrix_L(1:Nhh,1:Nkk)=0.0d0
	  !if(debu) write(66,*) 'call !' 
      call fun_L_curledge_curledge(N_volu,N_node,N_edge,nv_max,EV,Matrix_P0, & 
          VP,EV_ind,Neps,n_EV,N_thread,Nhh,Nkk,hhout,kkout,bared,thr,Matrix_L)
       !if(debu) write(66,*) 'matrix created'
      
      deallocate(EV,EV_ind,Matrix_P0,n_EV,hhout,kkout,bared)
	  !if(debu) write(66,*) 'deallocate'

! Create a matrix for the return argument 1
      mo=Nhh
      no=Nkk
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
    !if(debu) write(66,*) 'I am here1'
! Load the output 1 into a MATLAB array.
      Matrix_L_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(Matrix_L, Matrix_L_pr, siz)

      !if(debu) write(66,*) 'matrix convereted to matlab'

      deallocate(Matrix_L)
      
	  
      !if(debu) write(66,*) 'closing .txt'
      !if(debu) close(66)
      return
      end
