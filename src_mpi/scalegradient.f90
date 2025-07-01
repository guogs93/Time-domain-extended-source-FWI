! ============================================================================
!  subroutine scalegradient: scale the gradient
! Find where (index imax) the max (gradmax) of the gradient is reached
! Defined alpha as alpha=m(imax) * perc/gradmax
! Scale the gradient and the misfit function  grad = grad * alpha
!                                             cost = cost * alpha
! ============================================================================

subroutine scalegradient(grad,cost,m,nop,perc,alpha)
  implicit none
  integer :: nx,nop,i,imax
  real :: perc,cost,gradmax,alpha
  real :: m(nop),grad(nop)

  gradmax=-1.
  imax=-1
  do i=1,nop
     if (abs(grad(i)).gt.gradmax) then
        gradmax=abs(grad(i))
        imax=i
     end if
  end do
  write(*,*) 'GRADMAX PERC = ',gradmax,perc
  alpha=m(imax)*perc/gradmax
  grad(:)=grad(:)*alpha
  cost=cost*alpha
  write(*,*) 'alpha = ',alpha
end subroutine scalegradient
