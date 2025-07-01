! ===================================================================================
! SUBROUTINE SUBACQUI
! PROCESS ACQUISITION FILE
! FORMAT OF THE ACQUISITION FILE:
! xs1,xs2,icode
! xr1,xr2,icode
! ===================================================================================

  ! --------------------------------------------------------------
  ! SUBROUTINE SUBACQUI
  ! --------------------------------------------------------------
  SUBROUTINE subacqui(name_acqui,grid,src_tab,rec_tab)
    USE kwsinc
    USE parameters
    IMPLICIT NONE
    TYPE (grid_type)            :: grid
    INTEGER                     :: i,j,ix(2),ialpha(2),nalpha_frees,debug
    REAL                        :: xs(2), xr(2)
    REAL                        :: dalpha,xoversurf
    INTEGER                     :: nsrc,nrec  
    INTEGER                     :: src_tab(grid%nsrc(1),5),rec_tab(grid%nsrc(1),grid%nrecmax,5)
    character(len=26)           :: name_acqui
    integer                     :: unit_acqui=20
  
    debug = 1

    OPEN(unit_acqui,file=name_acqui)
    READ(unit_acqui,*) 
    READ(unit_acqui,*) !grid%nsrc(1),grid%nrec(1)
      
    DO i=1, grid%nsrc(1)
       READ(unit_acqui,*) 
       READ(unit_acqui,*) xs(1),xs(2)
       CALL subsourcekws(2,xs,grid%h,grid%ofs,ix,ialpha,nalpha_frees)
       src_tab(i,1) = ix(1) + grid%npml 
       src_tab(i,2) = ix(2) + grid%npml
       src_tab(i,3) = ialpha(1)
       src_tab(i,4) = ialpha(2)

       IF (grid%ofs.EQ.0) THEN
          src_tab(i,5) = 0
       ELSE
          src_tab(i,5) = nalpha_frees
       END IF
    
       READ(unit_acqui,*) grid%nrec(i)
       DO j=1, grid%nrec(i)
          READ(unit_acqui,*) xr(1),xr(2)
          CALL subsourcekws(2,xr,grid%h,grid%ofs,ix,ialpha,nalpha_frees)
          rec_tab(i,j,1) = ix(1) + grid%npml
          rec_tab(i,j,2) = ix(2) + grid%npml
          rec_tab(i,j,3) = ialpha(1)
          rec_tab(i,j,4) = ialpha(2)

          IF (grid%ofs.EQ.0) THEN
             rec_tab(i,j,5) = 0
          ELSE
             rec_tab(i,j,5) = nalpha_frees
          END IF
       END DO

    END DO

    CLOSE(unit_acqui)
    
  END SUBROUTINE subacqui
  
