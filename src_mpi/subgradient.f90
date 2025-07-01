! =============================================================================================
! SUBROUTINE SUBGRADIENT
! Perform different action related to ES-FWI (parameter mode)
! === mode = 0: Perform forward modeling using as velocity model ctrue
! === mode = 1: ES-FWI based on MINRES
! ==============================================================================================
subroutine subgradient(grid,inv,src_tab,rec_tab,ix1min,ix1max,ix2min,ix2max,&
                      flaggrad,niter)
  use mpi
  use parameters
  use kwsinc
  implicit none
  include 'common.h'

  TYPE (grid_type)                  :: grid
  TYPE (inv_type)                   :: inv

  integer                           :: i,j,k,niter
  integer                           :: ntrace,ix1min,ix1max,ix2min,ix2max
  real, allocatable                 :: wavefieldr(:,:,:)
  real, allocatable                 :: wavefielde(:,:,:),wavefieldesc(:,:,:),awavefield(:,:,:)
  real, allocatable                 :: sumdeltadei(:,:,:),sumdeltadeo(:,:,:),wd(:,:,:),wb(:,:,:),source(:)
  real, allocatable                 :: gradtemp(:),gradprecotemp(:),precon(:),gradsumsmo(:)
  real, allocatable                 :: deltadet(:,:),deltadeadj(:,:)
  integer                           :: src_tab(grid%nsrc(1),5),rec_tab(grid%nsrc(1),grid%nrecmax,5)
  
  integer                           :: unit

  character(len=160)                :: namei

  integer                           :: istot,ipatch,isrc,irec0,bsize,kshot,niternum,niternumtemp
  integer                           :: ista,iend,ierr
  integer                           :: flaggrad

  real                              :: gradprecomax,gradpreconorm,gradnorm
  real                              :: ratio,ratiotemp,ratioaverage
  real                              :: costtemp,alpha,sum1,sum2,costtemp2,dto,sigmat,sigmar,fhig

  if(mype.eq.0) write(*,*) 'call subc2m',grid%n1,grid%n2,inv%nop,ix1min,ix1max,ix2min,ix2max
 
  call subc2m(grid%slow,grid%n1,grid%n2,inv%nop,ix1min,ix1max,ix2min,ix2max,inv%m)

  ! data samping for Nyquist theory

  grid%nto=nint(float(grid%nt-1)*grid%dt/grid%dto)+1

  allocate(wavefieldr(grid%n1,grid%n2,grid%nto))
  allocate(wavefieldesc(grid%n1,grid%n2,grid%nto))
  allocate(wavefielde(grid%n1,grid%n2,grid%nto))
  allocate(awavefield(grid%n1,grid%n2,grid%nto))
  allocate(inv%grad(inv%nop)) 
  allocate(inv%gradpreco(inv%nop))
  allocate(gradtemp(inv%nop))
  allocate(gradprecotemp(inv%nop))
  allocate(precon(inv%nop))
  allocate(source(grid%nt))
  allocate(gradsumsmo(inv%nop))

  inv%gradsum(:)=0.
  inv%costsum=0.
  inv%costsum2=0
  inv%gradprecosum(:)=0
  ratioaverage=0
  istot=0
  irec0=1
  kshot=0  
  niternumtemp=0
  niternum=0

  ! data samping for Nyquist theory
  
  if(mype==0) write(*,*)'data samping for Nyquist theory',grid%nto

  do ipatch=1,grid%npatch
     if (mype==0) write(*,*) "Mype Patch n ",mype,ipatch,grid%nsrc(ipatch),grid%nrec(ipatch)
     bsize=grid%nt*grid%nrec(ipatch)*4

     call para_range(1,grid%nsrc(ipatch),nproc,mype,ista,iend)
 
     if (ipatch.gt.1) kshot=kshot+grid%nsrc(ipatch-1)

     inv%grad(:)=0.
     inv%gradpreco(:)=0
     gradtemp(:)=0.
     gradprecotemp(:)=0
     inv%cost=0.
     inv%cost2=0
     costtemp=0.
     costtemp2=0

     do isrc=ista,iend

        istot=kshot+isrc

        allocate(sumdeltadei(grid%nt,grid%nrec(istot),grid%nsrc(ipatch)))
        allocate(sumdeltadeo(grid%nt,grid%nrec(istot),grid%nsrc(ipatch)))

        allocate(wd(grid%nt,grid%nrecmax,grid%nsrc(ipatch)))
        allocate(wb(grid%n1,grid%n2,grid%nsrc(ipatch)))

        if(grid%mode.ge.0)then
          open(NEWUNIT=unit,file='windowd.bin',access='stream',form='unformatted')
          read(unit) wd(:,:,:)
          close(unit)
          open(NEWUNIT=unit,file='windowb.bin',access='stream',form='unformatted')
          read(unit) wb(:,:,:)
          close(unit)
        endif

        sumdeltadei=0
        sumdeltadeo=0

        allocate(inv%dataobs(grid%nt,grid%nrec(istot)))
        allocate(inv%datar(grid%nt,grid%nrec(istot)))
        allocate(inv%deltadr(grid%nt,grid%nrec(istot)))
        allocate(inv%deltade(grid%nt,grid%nrec(istot)))
        allocate(inv%datae(grid%nt,grid%nrec(istot)))
        allocate(inv%dataesc(grid%nt,grid%nrec(istot)))
        allocate(deltadet(grid%nt,grid%nrec(istot)))
        allocate(deltadeadj(grid%nt,grid%nrec(istot)))

        write(*,*) istot, grid%nrec(istot), grid%nrecmax

        write(*,'(i4)', Advance = 'no') istot
        inv%dataobs(:,:)=0.
        inv%datar(:,:)=0.
        inv%deltadr(:,:)=0
        inv%deltade(:,:)=0
        inv%datae(:,:)=0
        inv%dataesc(:,:)=0
      
        ! read in wavelet
        call subname(istot,namei)

        open(NEWUNIT=unit,file='wavelet_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted',status='old')
        read(unit) source(:)
        close(unit)

        grid%s(:)=source(:)/grid%h**2

        ! =============================================================================
        ! MODE 0: forward simulation
        ! =============================================================================
        if (grid%mode==0) then

           call submodeling(istot,grid%ctrue,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,grid%nto,grid%dto,grid%s, &
                grid%type_src,grid%type_rec,src_tab(istot,:),rec_tab(istot,:,:),grid%nrecmax,grid%nrec(istot),inv%dataobs,wavefieldr)

           call subname(istot,namei)

           do j=1, grid%nrec(istot)
           do i=1, grid%nt
              inv%dataobs(i,j)=wd(i,j,isrc)*inv%dataobs(i,j)
           enddo
           enddo

           open(NEWUNIT=unit,file='dataobs_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
           write(unit) inv%dataobs
           close(unit)
           
        elseif (grid%mode==1) then

           !**********************************************************
           ! step 1: calculate reduced data residuals (delta d^r)
           !**********************************************************

           ! step 1.1: Forward -> reduced wavefield and data (u^r and d^r) 
           call submodeling(istot,grid%c,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,grid%nto,grid%dto,grid%s, &
                grid%type_src,grid%type_rec,src_tab(istot,:),rec_tab(istot,:,:),grid%nrecmax,grid%nrec(istot),inv%datar,wavefieldr)
                
           call subname(istot,namei)

           do j=1, grid%nrec(istot)
           do i=1, grid%nt
              inv%datar(i,j)=wd(i,j,isrc)*inv%datar(i,j)
           enddo
           enddo

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='datar_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) inv%datar
              close(unit)
           end if

           open(NEWUNIT=unit,file='dataobs_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
           read(unit) inv%dataobs
           close(unit)

           do j=1, grid%nrec(istot)
           do i=1, grid%nt
              inv%dataobs(i,j)=wd(i,j,isrc)*inv%dataobs(i,j)
           enddo
           enddo

           open(NEWUNIT=unit,file='dataobs_'//namei(1:LEN_TRIM(namei))//'.binw',access='stream',form='unformatted')
           write(unit) inv%dataobs
           close(unit)

            ! step 1.2: calculate reduced data residuals (time is reversed at here) (delta d^r=d^obs-d^r)

           call subresiduals(istot,inv%dataobs,inv%datar,inv%deltadr,grid%nt,grid%nrec(istot),inv%cost)

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='deltadr_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              do i=1, grid%nrec(istot)
              do j=1, grid%nt
                  write(unit) inv%deltadr(grid%nt-j+1,i)
              enddo
              enddo
              close(unit)
           end if

           !***********************************************************************************************************
           ! step 2: calculate extended data residuals using matching filter (delta d^e=(1\mu*S(m)S(m)^T+I) delta d^r)
           !***********************************************************************************************************
           
           do j=1, grid%nrec(istot)
           do i=1, grid%nt
              inv%deltadr(grid%nt-i+1,j)=wd(i,j,isrc)*inv%deltadr(grid%nt-i+1,j)
           enddo
           enddo

           call sublinearminres(istot,ipatch,grid,inv,rec_tab(istot,:,:),wd(:,:,isrc),wb(:,:,isrc),inv%deltade,niternum)

           do j=1, grid%nrec(istot)
           do i=1, grid%nt
              inv%deltade(grid%nt-i+1,j)=wd(i,j,isrc)*inv%deltade(grid%nt-i+1,j)
           enddo
           enddo



           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='deltade_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              do i=1, grid%nrec(istot)
              do j=1, grid%nt
                  write(unit) inv%deltade(grid%nt-j+1,i)
              enddo
              enddo
              close(unit)
           end if

           !**********************************************************************************
           ! step 3: adjoint -> calculate source extension  (deltab^e=1\mu*G^T delta d^e)
           !**********************************************************************************
           call submodeling_adjoint_esource(istot,grid%c,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,grid%nto,grid%dto, &
                inv%deltade,grid%type_rec,grid%nrec(istot),grid%nrecmax,rec_tab(istot,:,:),wb(:,:,isrc),inv%dataesc,wavefieldesc,inv%cost2)
          
           if (grid%oio.eq.1) then
              !open(NEWUNIT=unit,file='wavefield_esc'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              !write(unit) wavefieldesc(:,:,:)
              !close(unit)

              open(NEWUNIT=unit,file='dataesc_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) inv%dataesc(:,:)
              close(unit)
           end if
           !************************************************************
           ! summmation of reduced wave and scatter wave
           !************************************************************
           wavefielde(:,:,:)=wavefieldr(:,:,:) + wavefieldesc(:,:,:)
           inv%datae(:,:)=inv%datar(:,:) + inv%dataesc(:,:)
           
           do j=1, grid%nrec(istot)
           do i=1, grid%nt
              inv%datae(i,j)=wd(i,j,isrc)*inv%datae(i,j)
           enddo
           enddo

           if (grid%oio.eq.1)then
              open(NEWUNIT=unit,file='datae_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) inv%datae
              close(unit)
           end if

           !**************************************************
           ! step 5: adjoint -> extended adjoint wavfield
           !**************************************************

           if(niter.eq.0)then
             sumdeltadeo(:,:,isrc)=inv%deltade(:,:)
             deltadet(:,:)=inv%deltade(:,:)
           else
              ! read in the previous summation of extended data residuals
              open(NEWUNIT=unit,file='sumdeltade_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              do j=1, grid%nrec(istot)
              do i=1, grid%nt
                 read(unit)sumdeltadei(grid%nt-i+1,j,isrc)
              enddo
              enddo
              close(unit)

              ! calculate summation for save 
              sumdeltadeo(:,:,isrc)=sumdeltadei(:,:,isrc)+inv%deltade(:,:)

              ! calculate current data residuals
              deltadet(:,:)=inv%deltade(:,:)+sumdeltadeo(:,:,isrc)
           endif

           open(NEWUNIT=unit,file='sumdeltade_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
           do j=1, grid%nrec(istot)
           do i=1, grid%nt
           write(unit) sumdeltadeo(grid%nt-i+1,j,isrc)
           enddo
           enddo
           close(unit)

           do j=1, grid%nrec(istot)
           do i=1, grid%nt
              deltadeadj(grid%nt-i+1,j)=wd(i,j,isrc)*deltadet(grid%nt-i+1,j)
           enddo
           enddo

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='deltadew_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              do i=1, grid%nrec(istot)
              do j=1, grid%nt
                 write(unit) deltadeadj(grid%nt-j+1,i)
              enddo
              enddo
              close(unit)
           end if
           
           call submodeling_adjoint(istot,grid%c,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,grid%nto,grid%dto,&
                deltadeadj,grid%nrec(istot),grid%nrecmax,rec_tab(istot,:,:),wb(:,:,isrc),awavefield)

           !*************************************
           ! step 6: correlation -> gradient
           !*************************************
           call subcor(wavefielde,awavefield,wb(:,:,isrc),grid%n1,grid%n2,grid%nto,grid%dto,ix1min,ix1max,ix2min,ix2max,inv%nop,&
                inv%grad,inv%gradpreco)       

          end if
 
        deallocate(inv%dataobs)
        deallocate(inv%datar)
        deallocate(inv%deltadr)
        deallocate(inv%deltade)
        deallocate(inv%dataesc)
        deallocate(inv%datae)
        deallocate(deltadet)
        deallocate(deltadeadj)
 
     end do !end do is

     deallocate(sumdeltadei,sumdeltadeo,wd,wb)
          

     if (grid%mode.ge.1) then
        call MPI_REDUCE(inv%grad,gradtemp,inv%nop,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(inv%gradpreco,gradprecotemp,inv%nop,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(inv%cost,costtemp,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(inv%cost2,costtemp2,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if (mype==0) write(*,*) 'costtemp = ',costtemp
     end if

     !Summation over patchs

     if (mype.eq.0) then
        inv%gradsum(:)=inv%gradsum(:)+gradtemp(:)
        inv%gradprecosum(:)=inv%gradprecosum(:)+gradprecotemp(:)
        inv%costsum=inv%costsum+costtemp
        inv%costsum2=inv%costsum2+costtemp2
     end if

     irec0=irec0+grid%nrec(ipatch)
  end do
  !end do ip

  ! smooth the gradient
  
  if(mype==0)then
 
    allocate(grid%topo(grid%n2),grid%itopo(grid%n2))

    open(12,file='ftopo.ascii')
    do i=1, grid%n2
       read(12,*) grid%topo(i)
    end do
    close(12)
    do i=1, grid%n2
       grid%itopo(i)=nint(grid%topo(i)/grid%h)+2
    end do
    
    open(13,file='smooth.par')
    read(13,*)inv%tau1,inv%tau2
    read(13,*)inv%freq
    close(13)

    call smoothgauss2d(inv%gradsum,grid%c,inv%freq,inv%nop,grid%n1,grid%n2,grid%h,gradsumsmo,grid%itopo,inv%tau1,inv%tau2)

    inv%gradsum=gradsumsmo

    deallocate(grid%topo,grid%itopo)

  endif
  
  ! ============================================================================================
  !    Deallocation
  ! ============================================================================================

  deallocate(wavefieldr)
  deallocate(wavefielde,wavefieldesc)
  deallocate(awavefield)
  deallocate(inv%grad)
  deallocate(inv%gradpreco)
  deallocate(gradtemp)
  deallocate(gradprecotemp)
  deallocate(precon)
  deallocate(source)

end subroutine subgradient

! ============================================================================================
!             SUBROUTINE SUBNORM
! ============================================================================================

SUBROUTINE subnorm(x,n,x_norm)
  IMPLICIT NONE
  INTEGER :: i,n
  REAL :: x(n),x_norm
  x_norm=0.
  DO i=1,n
     x_norm=x_norm+ABS(x(i))
  END DO
  x_norm=x_norm/FLOAT(n)

END SUBROUTINE subnorm

! ============================================================================================
!             SUBROUTINE SUB_MAX
! ============================================================================================

SUBROUTINE submax(x,n,xmax)
  IMPLICIT NONE
  INTEGER :: n,i
  REAL :: x(n),xmax
  xmax=-1.
  DO i=1,n
     IF (ABS(x(i)).GT.xmax) xmax=ABS(x(i))
  END DO

END SUBROUTINE submax

! ============================================================================================
!             SUBROUTINE STEPLENGTH
! ============================================================================================

SUBROUTINE steplength(grad,x,n,perc,alpha)
  IMPLICIT NONE
  INTEGER :: n,i,imax
  REAL :: perc,gradmax,alpha
  REAL :: x(n),grad(n)

  gradmax=-1.
  imax=-1
  DO i=1,n
     IF (ABS(grad(i)).GT.gradmax) THEN
        gradmax=abs(grad(i))
        imax=i
     END IF
  END DO
  write(*,*) 'GRADMAX PERC = ',gradmax,perc
  alpha=x(imax)*perc/gradmax

END SUBROUTINE  steplength
