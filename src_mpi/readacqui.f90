! =====================================================================================================
! subroutine readacqui
! Read acquisition file
! format:
! xs1(1),xs2(1),t,unc,icode=0
! xr1(1,1),xs2(1,1),t,unc,icode=1
! xr1(2,1),xs2(2,1),t,unc,icode=1
! ...
! xr1(nr(1),1),xs2(nr(1),1),t,unc,icode=1
! xs1(2),xs2(2),t,unc,icode=0
! xr1(1,2),xs2(1,2),t,unc,icode=1
! xr1(2,2),xs2(2,2),t,unc,icode=1
! ...
! xr1(nr(2),2),xs2(nr(2),2),t,unc,icode=1
! ...
! ...
! =====================================================================================================

subroutine readacqui(name_acqui,dx,is1,is2,ir1,ir2,nr,ns,nrmax,tabs,ntrace)
  implicit none
  integer :: ns,nrmax,nr0,icode,ns0,ntrace
  real :: t,unc,x1,x2,dx
  character(len=160) name_acqui
  integer :: nr(ns),is1(ns),is2(ns),ir1(nrmax,ns),ir2(nrmax,ns),tabs(ns)
  integer, parameter :: unit=26

  open(unit,file=name_acqui)
  ns0=0
  ntrace=0
  tabs(1)=1
11 read(unit,*,end=91) x1,x2,t,unc,icode
  if (icode==0) then
     ns0=ns0+1
     nr(ns0)=0
     is1(ns0)=int(x1/dx)+1
     is2(ns0)=int(x2/dx)+1
     if (ns0.gt.1) tabs(ns0)=ntrace+1
  else
     ntrace=ntrace+1
     nr(ns0)=nr(ns0)+1
     ir1(nr(ns0),ns0)=int(x1/dx)+1
     ir2(nr(ns0),ns0)=int(x2/dx)+1
  end if
  go to 11
91 continue
  close(unit)
  if (ns0 .ne. ns) stop 'BUG: ns non equal to ns0'

end subroutine readacqui

! =====================================================================================================
! subroutine readacqui0 (just count the number of source and the maximum number of receiver per source)
! Read acquisition file
! format:
! xs1(1),xs2(1),t,unc,icode=0
! xr1(1,1),xs2(1,1),t,unc,icode=1
! xr1(2,1),xs2(2,1),t,unc,icode=1
! ...
! xr1(nr(1),1),xs2(nr(1),1),t,unc,icode=1
! xs1(2),xs2(2),t,unc,icode=0
! xr1(1,2),xs2(1,2),t,unc,icode=1
! xr1(2,2),xs2(2,2),t,unc,icode=1
! ...
! xr1(nr(2),2),xs2(nr(2),2),t,unc,icode=1
! ...
! ...
! =====================================================================================================

subroutine readacqui0(name_acqui,ns,nrmax)
  implicit none
  integer :: ns,nrmax,icode,nr
  real :: t,unc,x1,x2
  character(len=160) name_acqui
  integer, parameter :: unit=26

  open(unit,file=name_acqui)
  ns=0
  nrmax=-1
11 read(unit,*,end=91) x1,x2,t,unc,icode
  if (icode==0) then
     ns=ns+1
     nr=0
  else
     nr=nr+1
     nrmax=max(nr,nrmax)
  end if
  go to 11
91 continue
  close(unit)
  write(*,*) "Number of shots: ",ns
  write(*,*) "Maximum number of receivers / shot: ",nrmax

end subroutine readacqui0

! =================================================================================================
! subroutine sub_read_acqui
! Read fixed-spread (stationary recording) acquisition file
! The acquisition consists of several fixed-spread patches
! Each patch contains nsrc and nrec receivers
! the format of the acquisition file is as follow
! npatch
! nsrc(1) nrec(1)
! xs1(1), xs2(1)
! xs1(2), xs2(2)
! ...
! xs1(nsrc(1)), xs2(nrec(1))
! xr1(1), xr2(1)
! xr1(2), xr2(2)
! ...
! xr1(nsrc(1)), xr2(nrec(1))
! xs1(1+nsrc(1)), xs2(1+nsrc(1))
! xs1(2+nsrc(1)), xs2(2+nsrc(1))
! ...
! xs1(nsrc(2)+nsrc(1)), xs2(nrec(2)+nsrc(1))
! xr1(1+nsrc(1)), xr2(1+nsrc(1))
! xr1(2+nsrc(1)), xr2(2+nsrc(1))
! ...
! xr1(nrec(2)+nrec(1)), xr2(nrec(2)+nsrc(1))
! ....
! ....
! ======================================================================================================

SUBROUTINE sub_read_acqui(name_acqui,grid)
  USE parameters
  IMPLICIT NONE
  TYPE (grid_type) :: grid

  INTEGER :: i,j,k,l,unit
  REAL :: xs1,xs2,xr1,xr2,x1max,x2max
  CHARACTER(LEN=160) :: name_acqui

  OPEN(NEWUNIT=unit,file=name_acqui)

  READ(unit,*) grid%npatch
  write(*,*) grid%npatch
  
  grid%nsrc_tot=0
  grid%nrec_tot=0

  x1max=float(grid%n1-1)*grid%h
  x2max=float(grid%n2-1)*grid%h
  write(*,*) 'X1MAX X2MAX = ',x1max,x2max
! write(*,*)'I am here'
  ALLOCATE(grid%nsrc(grid%npatch))
  ALLOCATE(grid%nrec(grid%npatch))
!write(*,*)'I am here'
  grid%ntrace=0
  DO i=1,grid%npatch
     READ(unit,*) grid%nsrc(i),grid%nrec(i)
     ! write(*,*) grid%nsrc(i),grid%nrec(i)
     grid%nsrc_tot=grid%nsrc_tot+grid%nsrc(i)
     grid%nrec_tot=grid%nrec_tot+grid%nrec(i)
     grid%ntrace=grid%ntrace+grid%nsrc(i)*grid%nrec(i)
  END DO

  WRITE(*,*) ''
  WRITE(*,*) "==================================================="
  WRITE(*,*) "NUMBER OF TRACES: ",grid%ntrace
  WRITE(*,*) "==================================================="
  WRITE(*,*) ''

  ALLOCATE (grid%is(2,grid%nsrc_tot))
  ALLOCATE (grid%ir(2,grid%nrec_tot))

  k=0
  l=0
  DO i=1,grid%npatch
     DO j=1,grid%nsrc(i)
        k=k+1
        READ(unit,*) xs1,xs2
        ! write(*,*) xs1,xs2
        IF (xs1.LT.0..OR.xs1.GT.x1max) GOTO 94
        IF (xs2.LT.0..OR.xs2.GT.x2max) GOTO 95
        grid%is(1,k)=int(xs1/grid%h)+1+grid%npml
        grid%is(2,k)=int(xs2/grid%h)+1+grid%npml
     END DO
     DO j=1,grid%nrec(i)
        l=l+1
        READ(unit,*) xr1,xr2
        ! write(*,*) xr1,xr2
        IF (xr1.LT.0..OR.xr1.GT.x1max) GOTO 96
        IF (xr2.LT.0..OR.xr2.GT.x2max) then
           write(*,*) xr2,x2max
           GOTO 97
        end if
        grid%ir(1,l)=int(xr1/grid%h)+1+grid%npml
        grid%ir(2,l)=int(xr2/grid%h)+1+grid%npml
        !write(*,*) grid%ir(1,l),grid%ir(2,l)
     END DO
  END DO

  CLOSE(unit)

  RETURN

94 WRITE(*,*) 'SOURCE OUTSIDE THE DOMAIN (DIM 1)'
  STOP

95 WRITE(*,*) 'SOURCE OUTSIDE THE DOMAIN (DIM 2)'
  STOP

96 WRITE(*,*) 'RECEIVER OUTSIDE THE DOMAIN (DIM 1)'
  STOP

97 WRITE(*,*) 'RECEIVER OUTSIDE THE DOMAIN (DIM 2)'
  STOP

END SUBROUTINE sub_read_acqui
