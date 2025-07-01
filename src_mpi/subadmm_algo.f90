  subroutine subadmm_grad(grid, m, p, q, grad_term)
    use parameters
    implicit none
    type (grid_type) :: grid
    real :: m(grid%n1,grid%n2)
    real :: p(grid%n1,grid%n2,2), q(grid%n1,grid%n2,2)
    real :: Am(grid%n1,grid%n2,2), arg_tmp(grid%n1,grid%n2,2)
    real :: grad_term(grid%n1,grid%n2)

    call compute_Am(grid, m, Am)
     
    arg_tmp = (p + q) - Am
    call compute_Am_adjoint(grid,arg_tmp, grad_term)
      
  end subroutine subadmm_grad
  
  subroutine subadmm_nste(grid,perc_admm, m, p, q)
    use parameters
    implicit none
    type (grid_type) :: grid
    real :: perc_admm, sigma_admm
    real :: m(grid%n1,grid%n2)
    real :: p(grid%n1,grid%n2,2), q(grid%n1,grid%n2,2)
    real :: Am(grid%n1,grid%n2,2), arg_prox(grid%n1,grid%n2,2)
    
    call compute_Am(grid, m, Am)
    arg_prox   = Am - q
    sigma_admm = perc_admm * maxval(abs(arg_prox))
    call sub_proximity(grid, arg_prox, sigma_admm, p)
    
    q = q - (Am - p)
    
  end subroutine subadmm_nste
  
  subroutine compute_Am(grid, m, Am)
    use parameters
    implicit none
    type (grid_type) :: grid                 
    real::  m(grid%n1,grid%n2)
    real:: Am(grid%n1,grid%n2,2)
    
    call subget_parameter_gradient(m, grid%n1, grid%n2, &
            grid%h, grid%h, Am)
            
  end subroutine compute_Am

  subroutine compute_Am_adjoint(grid, m, Am)
    use parameters
    implicit none
    type (grid_type) :: grid
    real :: m(grid%n1, grid%n2,2)
    real :: Am(grid%n1, grid%n2)
    
    call subget_parameter_adjoint_gradient(m, grid%n1, grid%n2, &
            grid%h, grid%h, Am)
            
  end subroutine compute_Am_adjoint
