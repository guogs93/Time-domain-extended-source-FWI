!***************************************************************************
! compute the model perturbation with CG iteration Ax=b ( H*dm=g )
! gaoshan guo
!***************************************************************************

subroutine subperturbation(nop,waterlevel,preco1d)
 use parameters
 implicit none
 real :: preco1d(nop)
 real :: precomax,waterlevel
 integer :: nop
 
 call submax(preco1d,nop,precomax)

 preco1d(:)=(preco1d(:)+waterlevel*precomax)

end subroutine


