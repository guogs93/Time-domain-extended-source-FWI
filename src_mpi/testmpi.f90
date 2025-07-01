! ==========================================================================
!                            PROGRAM TFWI2D
! 2D full waveform inversion for wavespeed reconstruction.
! Velocity-stress acoustic wave equation is considered for forward problem.
! Explosve sources are used. Records are pressure.
! FWI is performed by means of local optimization. The misfit function is the 
! L2 norm of the difference bertween modeled and recorded seismic data.
! The gradient of the misfit function is computed with the adjoint-state method.
! The local optimization is implemented with the SEISCOPE optimization toolbox.
!
! Input parameters:
! name_init(char): initial velocity model (binary format, float samples)
! name_acqui(char): acquisition file name
! n1(int),n2(int),h(real),nt(int),dt(real): spatial grid dimensions (1 is the fast index). h: grid interval. 
! nt: number of time steps. dt: time step.
! type_src: pressure (0), vertical force (1)
! type_rec: hydrophone (0), vertical geophone (1)
! peakfreq(real): peak frequency of the Ricker wavelet (Hz).
! nitermax(int): maximum number of FWI iterations.
! optimalgo(int): optimization algorithm. 0: steepest descent; 1: lbfgs.
! perc(real): the gradient is scaled by perc x vel / grad_max where grad_max is the
! max. of the gradient and vel is the wavespeed in the initial model at the position
! where the max. of the gradient is reached.
! x1min (real) x1max(real) x2min(real) x2max(real): the medium is updated in the domain
! [x1min;x1max]x[x2min;x2max].

! Ref: L. Metivier and R. Brossier. The SEISCOPE Optimization Toolbox: A large-scale nonlinear
! optimization library based on reverse communication, Geophysics, 81(2), pages F11-F25, 2016.
! ==========================================================================

program testmpi
  use mpi
  implicit none

!  include 'common.h'
  integer:: nproc,mype

  integer:: infompi
  integer:: ista, iend, i, ierr
  real:: y_glob(10)
  real, allocatable :: y_loc(:)

  call MPI_INIT(infompi)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, infompi)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mype, infompi)
  call para_range(1, 10, nproc, mype, ista, iend)
  
  write(*,*)ista,iend,mype

  ! each process allocate its own local array
  allocate(y_loc(1:iend-ista+1))
  y_loc(:) = 0

  ! each process make is own stuff
  do i=ista, iend
    y_loc(i-ista+1)=i
  enddo
  !write(*,*)'y_loc:',y_loc

  ! each process communicate its part the master process (0)
  call mpi_gather(y_loc, iend-ista+1, MPI_REAL, y_glob, iend-ista+1, MPI_REAL,0, MPI_COMM_WORLD, infompi)

  if (mype == 0) then
    write(*,*)'y_glob:',y_glob
  end if

  ! or each process communicate to all other processes its part
  !call mpi_allgather(y_loc, iend-ista+1, MPI_REAL, y_glob, iend-ista+1,MPI_REAL, MPI_COMM_WORLD, infompi)

  !write(*,*) mype, 'y_glob :', y_glob

  call MPI_FINALIZE(infompi)

end program testmpi

subroutine para_range(n1, n2, nprocs, irank, ista, iend)
 iwork1 = (n2 -n1 +1) / nprocs
 iwork2 =  MOD(n2 - n1 + 1, nprocs)
 ista = irank*iwork1 + n1 + MIN(irank, iwork2)
 iend = ista + iwork1 - 1
 IF (iwork2 .GT. irank) iend = iend + 1
end subroutine para_range

