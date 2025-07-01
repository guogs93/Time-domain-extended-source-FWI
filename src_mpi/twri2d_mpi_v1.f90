!**************************************************************************************************
! 2D time domain extended-source full waveform inversion 
! author: Gaoshan guo at geoazur
! Reference: 
! [1] Ali Gholami, Hossein S. Aghamiry, and Stéphane Operto, (2022),
!     Extended-space full-waveform inversion in the time domain with the augmented Lagrangian method, GEOPHYSICS 87: R63-R77. 
! [2] Gaoshan Guo, Stéphane Operto, Ali Gholami, and Hossein S. Aghamiry, (2024),
!     Time-domain extended-source full-waveform inversion: Algorithm and practical workflow, GEOPHYSICS 89: R73-R94. 
! [3] Gaoshan Guo, Stéphane Operto, (2024)
!     Regional exploration of crustal structures from long-offset OBS data by time-domain extended-source FWI: 
!     a case study from eastern Nankai Trough, offshore Japan. Submitted to GRL.
!**************************************************************************************************

PROGRAM twri2d_mpi_v1
  USE mpi
  USE parameters
  USE kwsinc
  implicit none
  INCLUDE 'common.h'
  INCLUDE 'optim_type.h'
  INCLUDE 'smumps_struc.h'

  TYPE (optim_type)                 :: optim
  TYPE (grid_type)                  :: grid
  TYPE (inv_type)                   :: inv
  TYPE (smumps_struc)               :: id

  CHARACTER(LEN=4)                  :: FLAG

  integer                           :: infompi,i,i1,i2

  character(len=160)                :: name_true,name_init,name_acqui
  
  integer                           :: unit

  real(KIND=8)                      :: tbegin,tend
  integer                           :: niter,niter0
  integer                           :: ix1min,ix1max,ix2min,ix2max
  integer                           :: flaggrad
  integer, parameter                :: ibound=1
  real, parameter                   :: thresconv=1e-16
  logical, parameter                :: debug= .false.
  real, allocatable                 :: dm(:)

  ! parameters for ADMM
  logical                           :: is_admm
  real                              :: perc_admm,ro_admm, perc_ro
  integer                           :: reg
  real                              :: gradtlinnorm,gradnorm
  real, dimension(:,:,:), allocatable :: p, q
  real, dimension(:,:), allocatable :: grad_term
  real, dimension(:),   allocatable :: gradtlin

  integer:: nnzreg, nnzregpreco
  integer, allocatable :: reg_coo(:),i_reg_coo(:),j_reg_coo(:)

  ! parameters for hicks
  integer, allocatable :: src_tab(:,:),rec_tab(:,:,:)

  call MPI_INIT(infompi)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,infompi)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,infompi)

  tbegin=MPI_WTIME()

  id%COMM = MPI_COMM_WORLD
  id%SYM =  0                   ! Unsymmetric matrix
  id%PAR =  1                   ! Host working
  id%JOB = -1                   ! MUMPS initialization
  CALL smumps( id )

  ! ###############################################################################
  ! Read input parameters
  ! ###############################################################################

  call subreadpar(name_true,name_init,name_acqui,grid,inv,reg,perc_ro,perc_admm)

  if (reg.eq.1)then
     is_admm=.true.
  elseif(reg.eq.0)then
     is_admm=.false.
  endif

  if (mype.eq.0) then
     write(*,*) ''
     write(*,*) '=============================================================================='
     write(*,*) 'Run mode: ',grid%mode
     write(*,*) 'Initial velocity model file name (name_init)',name_init(1:LEN_TRIM(name_init))
     write(*,*) 'Acquisition file name (name_init)',name_acqui(1:LEN_TRIM(name_acqui))
     write(*,*) 'Space & time grid dimension and interval (n1 n2 h nt dt)',grid%n1,grid%n2,grid%h,grid%nt,grid%dt
     write(*,*) 'Penalty term', inv%mu
     write(*,*) 'Max. number of iterations',inv%nitermax
     write(*,*) 'Optimization algorithm ',inv%optimalgo
     write(*,*) 'perc = ',inv%perc
     write(*,*) 'x1min x1max x2min x2max for gradient computation = ',grid%x1min,grid%x1max,grid%x2min,grid%x2max
     write(*,*) '=============================================================================='
     write(*,*) ''
     open(NEWUNIT=unit,file='tfwi2d.log')
     write(unit,*) ''
     write(unit,*) '================================================================================'
     write(unit,*) 'Run mode: ',grid%mode
     write(unit,*) 'Initial velocity model file name (name_init)',name_init(1:LEN_TRIM(name_init))
     write(unit,*) 'Acquisition file name (name_init)',name_acqui(1:LEN_TRIM(name_acqui))
     write(unit,*) 'Space & time grid dimension and interval (n1 n2 h nt dt)',grid%n1,grid%n2,grid%h,grid%nt,grid%dt
     write(unit,*) 'Penalty term', inv%mu
     write(unit,*) 'Max. number of iterations',inv%nitermax
     write(unit,*) 'Optimization algorithm ',inv%optimalgo
     write(unit,*) 'perc = ',inv%perc
     write(unit,*) 'x1min x1max x2min x2max for gradient computation = ',grid%x1min,grid%x1max,grid%x2min,grid%x2max
     write(unit,*) '=============================================================================='
     write(unit,*) ''
     close(unit)
  end if

  ! ############################################################################
  ! SET CONSTANTS. CALL SUBHICKS TO COMPUTE coef_hicks

  call allocatekws2d
  call kwsinc2d(grid%ofs)

  ! ###############################################################################
  ! Set the optimization parameters by applying the mask defined by ix1min,ix1max,ix2min,ix2max
  ! ###############################################################################

  if (grid%x1max.lt.grid%x1min) stop 'x1max < x1min!!!!!!!!!'
  if (grid%x2max.lt.grid%x2min) stop 'x2max < x2min!!!!!!!!!'

  ix1min=int(grid%x1min/grid%h)+1
  ix1min=max(1,ix1min)
  ix1max=int(grid%x1max/grid%h)+1
  ix1max=min(ix1max,grid%n1)

  ix2min=int(grid%x2min/grid%h)+1
  ix2min=max(1,ix2min)
  ix2max=int(grid%x2max/grid%h)+1
  ix2max=min(ix2max,grid%n2)

  inv%nop=(ix1max-ix1min+1)*(ix2max-ix2min+1)
  grid%nn=grid%n1*grid%n2

  ! ###############################################################################
  ! dynamic memory allocation
  ! ###############################################################################
  
  allocate(grid%c0(grid%n1,grid%n2))
  allocate(grid%c(grid%n1,grid%n2))
  allocate(grid%ctrue(grid%n1,grid%n2))
  allocate(grid%slow(grid%n1,grid%n2))
  allocate(grid%slow0(grid%n1,grid%n2))
  allocate(grid%slowt(grid%n1,grid%n2))
  allocate(inv%m(inv%nop))
  allocate(inv%gradsum(inv%nop))                    !full gradient
  allocate(inv%gradprecosum(inv%nop))
  allocate(dm(inv%nop))

  ! ###############################################################################
  ! build Ricker wavelet
  ! ###############################################################################

  allocate(grid%s(grid%nt))
  
  ! ###############################################################################
  ! open and read input models
  ! ###############################################################################

  if (mype.eq.0) then
     open(NEWUNIT=unit,file=name_true,access='direct',recl=grid%n1*grid%n2*4)
     read(unit,rec=1) grid%ctrue
     close(unit)
     open(NEWUNIT=unit,file=name_init,access='direct',recl=grid%n1*grid%n2*4)
     read(unit,rec=1) grid%c0
     close(unit)
  end if

  call MPI_BCAST(grid%ctrue,grid%n1*grid%n2,MPI_REAL,0,MPI_COMM_WORLD,infompi)
  call MPI_BCAST(grid%c0,grid%n1*grid%n2,MPI_REAL,0,MPI_COMM_WORLD,infompi)

  grid%c(:,:)=grid%c0(:,:)

  grid%slow0(:,:)=1.0/grid%c0(:,:)**2
  grid%slowt(:,:)=1.0/grid%ctrue(:,:)**2
  grid%slow(:,:)=1.0/grid%c(:,:)**2

  ! ###############################################################################
  !  read acquisition
  ! ###############################################################################
  
  ! ################################################################
  ! READ ACQUISITION

  OPEN(unit=20,file=name_acqui)
  READ(20,*) grid%npatch
  allocate(grid%nsrc(grid%npatch))
  READ(20,*) grid%nsrc(1),grid%nrecmax
  allocate(grid%nrec(grid%nsrc(grid%npatch)))
  close(20)

  allocate(src_tab(grid%nsrc(1),5),rec_tab(grid%nsrc(1),grid%nrecmax,5))
  call subacqui(name_acqui,grid,src_tab,rec_tab)

  ! ###############################################################################
  !  compute first gradient
  ! ###############################################################################
 
  flaggrad=1
  niter0=0
  call subgradient(grid,inv,src_tab,rec_tab,ix1min,ix1max,ix2min,ix2max,flaggrad,niter0)

  if(mype.eq.0)then
    open(unit=40,file='cost.dat')
    write(40,*) niter0,inv%costsum,inv%costsum2
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,infompi)

  ! ###############################################################################
  ! gradient output
  ! ###############################################################################

  if (mype==0) then
  
    open(NEWUNIT=unit,file='gradient0.bin',access='direct',recl=inv%nop*4)
    write(unit,rec=1) inv%gradsum
    close(unit)

    open(NEWUNIT=unit,file='preco0.bin',access='direct',recl=inv%nop*4)
    write(unit,rec=1) inv%gradprecosum
    close(unit)
  end if
  
  if (grid%oio.eq.1) stop

  if (grid%mode==0) then
     go to 99
  end if

  ! ###############################################################################
  ! prepare ADMM parameters
  if (is_admm) then
     if (mype == 0) then
        perc_admm = 0.2

        ! Dual variable
        allocate(q(grid%n1, grid%n2, 2))
        q(:,:,:) = 0.

        ! Auxiliary variable
        allocate(p(grid%n1, grid%n2, 2))
        p(:,:,:) = 0.

        allocate(gradtlin(inv%nop))
        allocate(grad_term(grid%n1, grid%n2))
     end if
  end if

  ! ###############################################################################
  ! setting of the optimization toolbox
  ! ###############################################################################

  if (mype==0) then
     optim%niter_max = inv%nitermax
     optim%conv = thresconv
     optim%print_flag = 1
     optim%debug = .true.
  end if

  !***************************
  ! bound projection
  !***************************

  optim%bound= ibound
  allocate(optim%ub(inv%nop),optim%lb(inv%nop))

  optim%ub(:)=1.0/inv%minvel**2
  optim%lb(:)=1.0/inv%maxvel**2
  optim%threshold=1e-16

  if(mype==0)then

    call buildHess_reg(grid%n1,grid%n2,inv%nop,nnzreg)
    allocate(reg_coo(nnzreg),i_reg_coo(nnzreg), j_reg_coo(nnzreg))
    open(NEWUNIT=unit,file='mat_reg.bin',access='STREAM',form='UNFORMATTED')
    read(unit) (reg_coo(i),i=1,nnzreg)
    close(unit)
    open(NEWUNIT=unit,file='IRN_reg.bin',access='STREAM',form='UNFORMATTED')
    read(unit) (i_reg_coo(i),i=1,nnzreg)
    close(unit)
    open(NEWUNIT=unit,file='JCN_reg.bin',access='STREAM',form='UNFORMATTED')
    read(unit) (j_reg_coo(i),i=1,nnzreg)
    close(unit)

  endif

  ! ###############################################################################
  ! optimization loop
  ! ###############################################################################

  FLAG='INIT'       
  niter = 0

  if(mype==0)then

    id%ICNTL(20) = 0                ! Dense RHS
    id%ICNTL(21) = 0                ! Dense solution
    id%ICNTL(27) = 1                ! Blocking factor for multiplerhs
    id%ICNTL(7)  = 7                ! without the nested dissection step, MUMPS makes an automatic choice depending on the matrix characteristics
    id%ICNTL(14) = 80
    id%ICNTL(5) = 0
    id%ICNTL(18) = 0
    id%ICNTL(23) = 0
    id%ICNTL(35) = 0
    id%CNTL(7) = 1e-5
  endif

  do while (niter.le.inv%nitermax)

     ! ADMM
     if (is_admm) then
        call MPI_BARRIER(MPI_COMM_WORLD,infompi)
        if (mype == 0) then
           
           call subadmm_grad(grid, grid%slow, p, q, grad_term)
           call subc2m(grad_term, grid%n1, grid%n2, inv%nop,ix1min,ix1max,ix2min,ix2max, gradtlin)
           
           call subnorm(inv%gradsum,inv%nop,gradnorm) 
           call subnorm(gradtlin,inv%nop,gradtlinnorm)

           if(gradtlinnorm.eq.0)then
             ro_admm=0
           else
             ro_admm=perc_ro*gradnorm/gradtlinnorm
           endif

           write(120,*)'ro_admm',ro_admm
           inv%gradsum = inv%gradsum + ro_admm * gradtlin
                          
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,infompi)
     end if 

     !***********************************************
     ! update the model with preconditioned gradient
     !***********************************************

     if(mype==0)then

       call subperturbation(inv%nop,inv%waterlevel,inv%gradprecosum)

       call buildLHS(grid%n1,grid%n2,inv%nop,inv%gradprecosum,ro_admm,nnzreg,reg_coo,i_reg_coo,j_reg_coo,nnzregpreco)

       id%N=inv%nop
       id%NNZ=nnzregpreco
       id%NRHS=1
       id%LRHS=inv%nop
       allocate(id%IRN(id%NNZ),id%JCN(id%NNZ),id%A(id%NNZ),id%RHS(id%LRHS))

       open(NEWUNIT=unit,file='mat_regpreco.bin',access='STREAM',form='UNFORMATTED')
       read(unit) (id%A(i),i=1,id%NNZ)
       close(unit)
       open(NEWUNIT=unit,file='IRN_regpreco.bin',access='STREAM',form='UNFORMATTED')
       read(unit) (id%IRN(i),i=1,id%NNZ)
       close(unit)
       open(NEWUNIT=unit,file='JCN_regpreco.bin',access='STREAM',form='UNFORMATTED')
       read(unit) (id%JCN(i),i=1,id%NNZ)
       close(unit)
       
       id%RHS(:)=inv%gradsum(:)

     endif

      ! LU FACTORIZATION
      id%JOB = 4
      call smumps( id )

      ! SOLUTION STEP
      id%JOB = 3
      call smumps( id )

     if(mype==0)then
       dm(:)=id%RHS
       if(niter.eq.0)then
         open(NEWUNIT=unit,file='gradientpreco0.bin',access='direct',recl=inv%nop*4)
         write(unit,rec=1) dm
         close(unit)
       endif

       deallocate(id%IRN,id%JCN,id%A,id%RHS)

       if (niter==0) then
           WRITE(*,*) 'SCALE GRADIENT'
           call steplength(dm,inv%m,inv%nop,inv%perc,inv%alpha)
       end if

       inv%m(:)=inv%m(:)-inv%alpha*dm(:)
       do i=1,inv%nop
          if(inv%m(i).gt.optim%ub(i)) then
            inv%m(i)=optim%ub(i)
          endif
          if(inv%m(i).lt.optim%lb(i)) then
            inv%m(i)=optim%lb(i)
          endif
       enddo

     endif

     call MPI_BCAST(FLAG,4,MPI_CHARACTER,0,MPI_COMM_WORLD,infompi)
     call MPI_BCAST(inv%m,inv%nop,MPI_REAL,0,MPI_COMM_WORLD,infompi)

     ! obtain the new slowness

     call subm2c(inv%m,grid%n1,grid%n2,inv%nop,ix1min,ix1max,ix2min,ix2max,grid%slow0,grid%slow)
     
     ! transform the new slowness to velocity 

     grid%c(:,:)=1.0/sqrt(grid%slow(:,:))

     niter=niter+1

     if(mype==0)then
        open(NEWUNIT=unit,file='fwimodel.bin',access='direct',recl=grid%nn*4)
        write(unit,rec=niter) grid%c
        close(unit)
           
        open(NEWUNIT=unit,file='vfinal.bin',access='direct',recl=grid%nn*4)
        write(unit,rec=1) grid%c
        close(unit)
     endif

     flaggrad=1
     call subgradient(grid,inv,src_tab,rec_tab,ix1min,ix1max,ix2min,ix2max,flaggrad,niter) 

     if(mype.eq.0)then
        write(40,*)niter,inv%costsum,inv%costsum2
     endif

     !! ADMM STUFF
     if (is_admm) then
        call MPI_BARRIER(MPI_COMM_WORLD,infompi)
        if (mype == 0) then
           call subadmm_nste(grid,perc_admm, grid%slow, p, q)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,infompi)
     end if

  end do

  if(mype.eq.0)then
  close(40)
  endif

  deallocate(optim%ub,optim%lb)

  if (is_admm) then
     if(mype==0)then
       deallocate(q)
       deallocate(p)
       deallocate(gradtlin)
       deallocate(grad_term)
     endif
  endif

  if(mype==0)then
    deallocate(reg_coo,i_reg_coo, j_reg_coo)
  endif

99 continue
  
  deallocate(grid%nsrc)
  deallocate(grid%nrec)
  deallocate(grid%ctrue)
  deallocate(grid%c)
  deallocate(grid%c0)
  deallocate(grid%slow)
  deallocate(grid%slowt)
  deallocate(grid%slow0)
  deallocate(inv%gradsum)
  deallocate(inv%gradprecosum)
  deallocate(inv%m)
  deallocate(src_tab,rec_tab)
  deallocate(dm)

  ! ###########################################################
  ! Deallocate memory for kws coefficients
  call deallocatekws2d 

  ! CLOSE MUMPS SESSION
  id%JOB = -2
  CALL smumps( id )

  tend=MPI_WTIME()

  if(mype.eq.0) write(*,*) "ELAPSED TIME = ",tend-tbegin

  call MPI_FINALIZE(infompi)

end program twri2d_mpi_v1
