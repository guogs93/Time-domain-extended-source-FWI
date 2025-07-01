! ===============================================================================
! subroutine subm2c
! m(nop) ---> c(n1,n2)
! c(n1,n2): 2D velocity grid input in the code
! m(nop): vector loaded in the SEISCOPE tool box
! ===============================================================================
subroutine subm2c(m,n1,n2,nop,ix1min,ix1max,ix2min,ix2max,c0,c)
  implicit none
  integer :: n1,n2,nop,ix1min,ix1max,ix2min,ix2max,i1,i2,k
  real :: m(nop),c0(n1,n2),c(n1,n2)
  c(:,:)=c0(:,:)
  k=0
  do i2=ix2min,ix2max
     do i1=ix1min,ix1max
        k=k+1
        c(i1,i2)=m(k)
     end do
  end do
end subroutine subm2c

! ===============================================================================
! subroutine subc2m
! c(n1,n2)  ---> m(nop)
! c(n1,n2): 2D velocity grid input in the code
! m(nop): vector loaded in the SEISCOPE tool box
! ===============================================================================

subroutine subc2m(c,n1,n2,nop,ix1min,ix1max,ix2min,ix2max,m)
  implicit none
  integer :: n1,n2,nop,ix1min,ix1max,ix2min,ix2max,i1,i2,k
  real :: m(nop),c(n1,n2)
  k=0
  do i2=ix2min,ix2max
     do i1=ix1min,ix1max
        k=k+1
        m(k)=c(i1,i2)
     end do
  end do
end subroutine subc2m
