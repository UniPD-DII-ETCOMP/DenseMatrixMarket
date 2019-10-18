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
	  mwPointer N_node_pr,N_face_pr,Area_pr,N_GAUSS1_pr,N_GAUSS2_pr,Matrix_P0_pr,F1_pr,N_thread_pr,L_pr
	  mwPointer iiout_pr, jjout_pr, Nii_pr, Njj_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
! fortran subroutine arguments
	  real*8,allocatable,dimension(:,:) :: Matrix_P0, F1_r, L
      integer*8,allocatable,dimension(:,:) :: F1
	  real*8,allocatable,dimension(:) :: Area, iiout_r, jjout_r
	  integer*8,allocatable,dimension(:) :: iiout, jjout
!	  integer*8,allocatable,dimension(:) :: 
	  real*8 :: N_node_r, N_face_r, N_GAUSS1_r, N_GAUSS2_r, N_thread_r, Nii_r, Njj_r
      integer*8 :: N_node, N_face, N_GAUSS1, N_GAUSS2, N_thread, Nii, Njj
	  character*80 msg
      logical debu
       
      debu = .true. ! .true. o .false. per attivare o disattivare il debug
	  if(debu) open(unit=66,file='logL_axi.txt',status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
if (nrhs .ne. 12) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_auto:nInput', &
                                '12 input arguments required')
elseif (nlhs .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_auto:nOutput', &
                                '1 output argument required')
endif	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Check to see INPUTS are numeric.
	  do ii = 1,8
        if (mxIsNumeric(prhs(ii)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:L_axi_vv:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
	  if(debu) write(66,*) 'sono arrivato al check numerico'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #1 is INTEGER and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be scalar, N_node.')
      endif	  
      siz = m*n
      N_node_pr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(N_node_pr, N_node_r, siz) ! da double precision a reale
      N_node=int(N_node_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 1' 	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #2 is INTEGER and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 2 must be scalar, N_face.')
      endif	  
      siz = m*n
      N_face_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(N_face_pr, N_face_r, siz) ! da double precision a reale
      N_face=int(N_face_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 2' 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!     Check that input #3 is INTEGER vector and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. 1 .or. n .ne. N_face) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 3: Area: 1xN_face')
      endif	  
      siz = m*n
      Area_pr = mxGetPr(prhs(3))
      allocate(Area(N_face))
      call mxCopyPtrToReal8(Area_pr, Area, siz) ! da double precision a reale  
	  if(debu) write(66,*) 'input 3'	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #4 is INTEGER and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 4 must be scalar, N_GAUSS1.')
      endif	  
      siz = m*n
      N_GAUSS1_pr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(N_GAUSS1_pr, N_GAUSS1_r, siz) ! da double precision a reale
      N_GAUSS1=int(N_GAUSS1_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 4' 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #5 is INTEGER and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 5 must be scalar, N_GAUSS2.')
      endif	  
      siz = m*n
      N_GAUSS2_pr = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(N_GAUSS2_pr, N_GAUSS2_r, siz) ! da double precision a reale
      N_GAUSS2=int(N_GAUSS2_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 5' 	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	       
!	  Check that input #6 is REAL matrix and fetch it
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. N_node .or. n .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 6: Matrix_P0')
      endif	  
      siz = m*n
      Matrix_P0_pr = mxGetPr(prhs(6))
	  allocate(Matrix_P0(N_node,3))
      call mxCopyPtrToReal8(Matrix_P0_pr, Matrix_P0, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 6'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #7 is INTEGER matrix and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      if(m .ne. 4 .or. n .ne. N_face) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 7: F1')
      endif	  
      siz = m*n
      F1_pr = mxGetPr(prhs(7))
      allocate(F1_r(4,N_face))
      allocate(F1(4,N_face))
      call mxCopyPtrToReal8(F1_pr, F1_r, siz) ! da double precision a reale  
      F1 = int(F1_r,8)
      deallocate(F1_r)
      write(66,*) 'N_face',N_face
      !do ii=1,N_face
  	 !write(66,*) 'ii',ii, F1(1,ii),F1(2,ii),F1(3,ii),F1(4,ii)
      !enddo
      if(debu) write(66,*) 'input 7'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #8 is INTEGER and fetch it
      m = mxGetM(prhs(8))
      n = mxGetN(prhs(8))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 8 must be scalar, N_thread.')
      endif	  
      siz = m*n
      N_thread_pr = mxGetPr(prhs(8))
      call mxCopyPtrToReal8(N_thread_pr, N_thread_r, siz) ! da double precision a reale
      N_thread=int(N_thread_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 8' 	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #9 is INTEGER and fetch it
      m = mxGetM(prhs(9))
      n = mxGetN(prhs(9))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 9 must be scalar, Nii.')
      endif	  
      siz = m*n
      Nii_pr = mxGetPr(prhs(9))
      call mxCopyPtrToReal8(Nii_pr, Nii_r, siz) ! da double precision a reale
      Nii=int(Nii_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 9' 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #10 is INTEGER and fetch it
      m = mxGetM(prhs(10))
      n = mxGetN(prhs(10))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 10 must be scalar, Njj.')
      endif	  
      siz = m*n
      Njj_pr = mxGetPr(prhs(10))
      call mxCopyPtrToReal8(Njj_pr, Njj_r, siz) ! da double precision a reale
      Njj=int(Njj_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 10' 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!     Check that input #11 is INTEGER vector and fetch it
      m = mxGetM(prhs(11))
      n = mxGetN(prhs(11))
      if(m .ne. 1 .or. n .ne. Nii) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 11: iiout: 1xNii')
      endif	  
      siz = m*n
      iiout_pr = mxGetPr(prhs(11))
      allocate(iiout(Nii))
	  allocate(iiout_r(Nii))
      call mxCopyPtrToReal8(iiout_pr, iiout_r, siz) ! da double precision a reale  
	  iiout=int(iiout_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 11'	
	  deallocate(iiout_r)	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!     Check that input #12 is INTEGER vector and fetch it
      m = mxGetM(prhs(12))
      n = mxGetN(prhs(12))
      if(m .ne. 1 .or. n .ne. Njj) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 12: jjout: 1xNjj')
      endif	  
      siz = m*n
      jjout_pr = mxGetPr(prhs(12))
      allocate(jjout(Njj))
	  allocate(jjout_r(Njj))
      call mxCopyPtrToReal8(jjout_pr, jjout_r, siz) ! da double precision a reale  
	  jjout=int(jjout_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 12'	 
	  deallocate(jjout_r)	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      	  	  
! call the computational subroutine.
      allocate(L(Nii,Njj))
      call funLphiphi3(N_node,N_face,Area,N_GAUSS1,N_GAUSS2,Matrix_P0,F1,N_thread,Nii,Njj,iiout,jjout,L)
       if(debu) write(66,*) 'ho chiamato la subroutine e creato le matrici'
      deallocate(Area,Matrix_P0,F1,jjout,iiout)
	  if(debu) write(66,*) 'ho deallocato le cose inutili'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! Create a matrix for the return argument 1
      mo=Nii
      no=Njj
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! Load the output 1 into a MATLAB array.
      L_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(L, L_pr, siz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      if(debu) write(66,*) 'ho convertito la matrice per matlab'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      deallocate(L)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      if(debu) write(66,*) 'ho deallocato la matrice e ora chiudo il logLve.txt'
      if(debu) close(66)
      return
      end

