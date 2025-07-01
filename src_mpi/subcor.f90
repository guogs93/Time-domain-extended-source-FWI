! ==============================================================================================
! subroutine subcor: perform correlation between incident and adjoint wavefields
! Here, wavefield is the time derivative of the incident wavefield in the frame of the first-order velocity stress equation
! awavefield is the adjoint wavefield
! ==============================================================================================

subroutine subcor(wavefield,awavefield,wb,n1,n2,nt,dt,ix1min,ix1max,ix2min,ix2max,nop,grad,gradpreco)
  implicit none
  integer               :: n1,n2,nt,ix1min,ix1max,ix2min,ix2max,nop
  integer               :: it,i1,i2,k
  real                  :: dt
  real                  :: wavefield(n1,n2,nt),awavefield(n1,n2,nt),grad(nop),gradpreco(nop),wb(n1,n2)
   
  do it=1, nt
     k=0
     do i2=ix2min,ix2max
        do i1=ix1min,ix1max
           k=k+1
           grad(k)=grad(k)+wavefield(i1,i2,it)*awavefield(i1,i2,nt-it+1)*wb(i1,i2)
           gradpreco(k)=gradpreco(k)+(wavefield(i1,i2,it)*wavefield(i1,i2,it)*wb(i1,i2))
        end do
     end do
  end do
   
end subroutine subcor

! ============================================================================================
!             SUBROUTINE DERIV
! ============================================================================================

subroutine deriv(wavefield,wavefieldd,n1,n2,nt,dt)
  implicit none
  integer :: n1, n2, nt, i1, i2, it
  real :: dt
  real :: wavefield(n1,n2,nt),wavefieldd(n1,n2,nt)
  
  wavefieldd(:,:,:)=0.
 
  do it=1,nt-1
     do i2=1,n2
        do i1=1,n1
           wavefieldd(i1,i2,it)=-(wavefield(i1,i2,it+1)-wavefield(i1,i2,it))/dt
        end do
     end do
  end do
 
end subroutine deriv
