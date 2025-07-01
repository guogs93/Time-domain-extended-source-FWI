!***********************************************************************************
! linear minization residuals (MINRES) to solve Ax=b
! solve the normal system H_d delta d^e = delta d^*
! gaoshan guo at CNRS geoazur
!***********************************************************************************

subroutine sublinearminres(istot,ipatch,grid,inv,ir,wd,wb,deltade,itn)
  use kwsinc
  use minresModule
  implicit none
  type (grid_type)                  :: grid
  type (inv_type)                   :: inv
  integer                           :: i,j,k,niter,ipatch,istot,nn
  real                              :: deltade(grid%nt,grid%nrec(istot))
  real, allocatable                 :: deltade1d(:),deltadr1d(:)
  integer                           :: ir(grid%nrecmax,5)
  real                              :: wd(grid%nt,grid%nrecmax),wb(grid%n1,grid%n2)
  character(len=160)                :: namei

  ! Local variables
  logical             :: normal, precon, ifLS
  integer             :: n, nout
  real            :: shift, pertM

  logical   :: checkA 
  integer   :: itnlim, istop, itn, nprint
  real  :: Anorm, Acond, Arnorm, rnorm, rtol, r1norm, ynorm
  real  :: enorm, etol, wnorm, xnorm

  nn=grid%nt*grid%nrec(istot)
 
  allocate(deltade1d(nn),deltadr1d(nn))
 
  ! transform the observed data residuals into 1D
  call sub2d21d(inv%deltadr,grid%nrec(istot),grid%nt,nn,deltadr1d) 

  normal=.false.
  precon=.true.

  shift=0

  checkA=.false.
  itnlim = 100
  rtol   = 1.0e-5

  if(grid%oio.eq.1)then
    nout = 111
  else
    nout = -111
  endif

  call subname(istot,namei)

  if(nout.gt.0)then
    open(nout,file='MINRES_'//namei(1:LEN_TRIM(namei))//'.txt',status='unknown')
  endif

  call MINRES(nn, deltadr1d, shift, checkA, normal, &
           deltade1d, itnlim, nout, rtol,istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm, &
           istot,ipatch,grid,inv,ir,wd,wb)
 
  if(nout.gt.0)then
    close(nout)
  endif

  call sub1d22d(deltade1d,grid%nrec(istot),grid%nt,nn,deltade)

 deallocate(deltade1d,deltadr1d)

end subroutine

subroutine datanorm(datao,nt,nr,data_norm)
 implicit none
 integer:: nt,nr,i,j
 real:: datao(nt,nr)
 real:: data_norm

 data_norm=0

 do i=1, nt
 do j=1, nr
    data_norm=data_norm+datao(i,j)**2
 enddo
 enddo

 data_norm=sqrt(data_norm)

end subroutine

!***************************************
! The transform form 2d to 1d
!***************************************
subroutine sub2d21d(deltade,nr,nt,nn,deltade1d)
  implicit none
  integer :: nr,nn,nt,it,ir,k
  real :: deltade1d(nn),deltade(nt,nr)
  
  k=0
  do ir=1, nr
     do it=1, nt
        k=k+1
        deltade1d(k)=deltade(it,ir)
     end do
  end do

end subroutine sub2d21d

!******************************************
! The transform form 1d to 2d
!******************************************
subroutine sub1d22d(deltade1d,nr,nt,nn,deltade)
  implicit none
  integer :: nr,nn,nt,it,ir,k
  real :: deltade1d(nn),deltade(nt,nr)

  k=0
  do ir=1, nr
     do it=1, nt
        k=k+1
        deltade(it,ir)=deltade1d(k)
     end do
  end do

end subroutine sub1d22d

!*************************************
subroutine scalL2(n,x,y,scal_xy)

  implicit none

  !IN
  integer :: n
  real,dimension(n) :: x,y
  !IN/OUT
  real :: scal_xy
  !Local variables
  integer :: i

  scal_xy=0.
  do i=1,n
     scal_xy=scal_xy+x(i)*y(i)
  enddo

end subroutine scalL2

!**********************************
subroutine normL2(n,x,norm_x)

  implicit none

  !IN
  integer :: n
  real,dimension(n) :: x
  !IN/OUT
  real :: norm_x
  !Local variables
  integer :: i

  norm_x=0.
  do i=1,n
     norm_x=norm_x+x(i)**2
  enddo
  norm_x=sqrt(norm_x)

end subroutine normL2

!*********************************************
! The bandpass filter with Hanning branches
!*********************************************
subroutine filmartin(deltade1do,nr,nt,nn,dt,fhig,deltade1d)
  implicit none
  integer:: nr,nt,nn,i,j,k,nto,i2,nf,nord
  real:: deltade1do(nn),deltade1d(nn)
  real, allocatable:: deltadeo(:,:),deltade(:,:)
  real:: flow, fhig,llow,lhig 
  real:: dt, fe, fmax
  integer :: conv
  integer :: iflow, ifhig, illow, ilhig
  real, allocatable               :: trace(:)
  real, allocatable               :: han1(:),han2(:)
  real, allocatable               :: fmartin(:)
  real*8, allocatable             :: dmartin(:)
  real*8, allocatable             :: preal(:),pimag(:)
  real*8, allocatable             :: ybuf(:)
  
  integer, parameter              :: idnum=14
  integer inum(idnum)
  data inum /8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536/
  
  allocate(deltadeo(nt,nr),deltade(nt,nr))
  
  call sub1d22d(deltade1do,nr,nt,nn,deltadeo)  
  
  do i=1,idnum
     if (nt.le.inum(i)) then
         nto=inum(i)
         go to 100
     end if
  end do

100  nto=nto
  
  ! exponent of the Hanning taper
  nord=2
  
  ! bandwith (flow,fhigh) and length of the slopes (llow,lhigh) 
  flow=0
  llow=0
  lhig=2 
  
  !max frequency
  
  nf=nto
  fe=1./dt
  fmax=fe/2
  
  if (fhig.gt.fmax) then
     write(*,*) 'WARNING: fhig > fmax; fhig forced to fmax'
     fhig=fmax
     lhig=0.
  end if

  conv=float(nto)/fe
  iflow=int(flow*conv)+1
  ifhig=int(fhig*conv)+1
  illow=int(llow*conv)+1
  ilhig=int(lhig*conv)+1

  if (iflow-illow.lt.1.and.iflow.gt.0) then
     illow=iflow-1
  end if
        
  if (ifhig+ilhig.gt.nf) then
     ilhig=nf-ifhig
  end if

  !write(*,*) 'iflow ifhig illow ilhig = ',iflow,ifhig,illow,ilhig
  
  allocate(trace(nto))
  allocate(han1(illow+1))
  allocate(han2(ilhig+1))
  allocate(fmartin(nf+1))
  allocate(dmartin(nf+1))
  allocate(preal(nto))
  allocate(pimag(nto))
  allocate(ybuf(nto))
        
  trace(:)=0.
  han1(:)=0.
  han2(:)=0.

! Building filter

  fmartin(:)=0.
        
! Low-cut branch

  if (illow.gt.0) then
     call sphanning(han1,illow,nord)
     do j=1,illow
        fmartin(iflow-j+1)=han1(j)
     end do
  end if

! High-cut branch

  han2(:)=0.

  if (ilhig.gt.0) then
     call sphanning(han2,ilhig,nord)
     do j=1,ilhig
        fmartin(ifhig+j-1)=han2(j)
     end do
  end if
        
  do j=iflow,ifhig
     fmartin(j)=1.
  end do
        
  dmartin(:)=dble(fmartin(:)+0.0001)

  do i2=1, nr 
     pimag(:)=0.
     do j=1, nt
        trace(j)=deltadeo(j,i2)
     enddo
     preal(:)=dble(trace(:))
     call realft2s(nto,2,preal,pimag,ybuf)
     do j=1, nto
        preal(j)=preal(j)*dmartin(j)
        pimag(j)=pimag(j)*dmartin(j)
     end do
     call realft2s(nto,1,preal,pimag,ybuf)
     do j=1, nt
        deltade(j,i2)=sngl(preal(j))
     enddo
  end do
  
  ! transform the data residuals into 1D
  call sub2d21d(deltade,nr,nt,nn,deltade1d)
  
  deallocate(trace)
  deallocate(han1)
  deallocate(han2)
  deallocate(fmartin)
  deallocate(dmartin)
  deallocate(preal)
  deallocate(pimag)
  deallocate(ybuf)
  deallocate(deltadeo,deltade)
  
end subroutine

!**************************************
! hanning window
!**************************************
subroutine sphanning(han,nf,nord)
  implicit none
  integer:: j, nf, nord
  real:: han(nf+1), han1, han2, han3
  real,parameter:: pi=3.14159265
  
  do j=1, nf+1
     han1=2.*pi*float(j-1)/float(2*nf)
     han2=(cos(han1)+1.)/2.+0.000001
     han3=log(han2)*float(nord)
     han(j)=exp(han3)
  end do
  
end subroutine sphanning
