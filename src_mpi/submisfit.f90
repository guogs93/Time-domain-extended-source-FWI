! ======================================================================
! SUBROUTINE SUBMISFIT
! Compute the data residuals (res) and the l2 misfit function (cost)
! ======================================================================

subroutine submisfit(is,deltadr,wavebe,nt,nr,n1,n2,cost1,cost2)

  implicit none
  integer :: nt,nr,n1,n2
  integer :: i,j,k,is
  real :: cost1,cost2
  real :: deltadr(nt,nr),wavebe(n1,n2,nt)
  integer:: unit
  character(len=80) :: namei
  logical, parameter :: debug =.false.  

  cost1=0.
  
  do i=1,nr
  do j=1,nt
     cost1=cost1+0.5*deltadr(j,i)**2
  enddo
  enddo

  cost2=0

  do k=1, nt
  do j=1, n2
  do i=1, n1
     cost2=cost2+0.5*wavebe(i,j,k)**2
  enddo
  enddo
  enddo

  write(namei,FMT='(I8)') is
  namei = ADJUSTL(namei)

  if (debug) then
     open(NEWUNIT=unit,file='datares_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
     do i=1, nr
     do j=1, nt
        write(unit) deltadr(nt-j+1,i)
     enddo
     enddo
     close(unit)

     open(NEWUNIT=unit,file='srcres_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
     write(unit) wavebe
     close(unit)
  end if
  
end subroutine submisfit

