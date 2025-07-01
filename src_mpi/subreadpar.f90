subroutine subreadpar(name_true,name_init,name_acqui,grid,inv,reg,perc_rho,perc_admm)
  USE parameters
  IMPLICIT NONE
          
  TYPE (grid_type)	    :: grid
  TYPE (inv_type)	    :: inv

  CHARACTER(LEN=160) :: name_true,name_init,name_acqui

  INTEGER, parameter :: unit=10

  integer:: niter_CG,style
  integer:: reg
  real:: perc_rho,perc_admm

  open(unit,file='tfwi2d.par')
  read(unit,*) grid%mode
  read(unit,*) name_true,name_init
  read(unit,*) name_acqui
  read(unit,*) grid%n1,grid%n2,grid%npml,grid%h,grid%nt,grid%dt,grid%dto
  read(unit,*) grid%type_src,grid%type_rec
  read(unit,*) grid%ofs
  read(unit,*) inv%mu
  read(unit,*) inv%nitermax
  read(unit,*) inv%optimalgo
  read(unit,*) inv%perc
  read(unit,*) inv%preco,inv%waterlevel
  read(unit,*) grid%x1min,grid%x1max,grid%x2min,grid%x2max
  read(unit,*) grid%oio
  read(unit,*) inv%maxvel,inv%minvel
  read(unit,*) reg
  read(unit,*) perc_rho,perc_admm
  read(unit,*) grid%flow, grid%fhig, grid%llow, grid%lhig
  close(unit)

end subroutine subreadpar

