! transform an integer in a character string
subroutine subname(i,name)
  implicit none
  integer :: i
  character(len=160) :: name
  write(name,FMT='(I8)') i
  if (i.lt.10) then
     name = '0000'//ADJUSTL(name)
  else if (i.ge.10.and.i.lt.100) then
     name = '000'//ADJUSTL(name)
  else if (i.ge.100.and.i.lt.1000) then
     name = '00'//ADJUSTL(name)
  else if (i.ge.1000.and.i.lt.10000) then
     name = '0'//ADJUSTL(name)
  else if (i.ge.10000.and.i.lt.100000) then
     name = ADJUSTL(name)
  else
     write(*,*) 'i too big (subroutine subname)'
     stop
  end if
end subroutine subname
