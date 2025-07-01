
subroutine filmartin2d(deltadeo,nr,nt,dt,flow,fhig,llow,lhig,deltade)
  implicit none
  integer:: nr,nt,nn,i,j,k,nto,i2,nf,nord
  real:: deltadeo(nt,nr),deltade(nt,nr)
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

do i=1,idnum
     if (nt.le.inum(i)) then
         nto=inum(i)
         go to 100
     end if
  end do

100  write(*,*) "nto = ",nto

  ! exponent of the Hanning taper
  nord=2


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

  write(*,*) 'iflow ifhig illow ilhig = ',iflow,ifhig,illow,ilhig

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

  deallocate(trace)
  deallocate(han1)
  deallocate(han2)
  deallocate(fmartin)
  deallocate(dmartin)
  deallocate(preal)
  deallocate(pimag)
  deallocate(ybuf)

end subroutine
