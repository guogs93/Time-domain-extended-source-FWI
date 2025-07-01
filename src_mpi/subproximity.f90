  subroutine sub_proximity(grid, arg_prox, ro_admm, p)
    use parameters
    implicit none
    type (grid_type) :: grid
    real :: p(grid%n1,grid%n2,2)
    real :: arg_prox(grid%n1,grid%n2,2)
    real :: ro_admm
    integer :: admm_type

    p(:,:,:) = sign(1., arg_prox(:,:,:)) * max( abs(arg_prox(:,:,:)) - ro_admm, 0.) 
    
  end subroutine sub_proximity
