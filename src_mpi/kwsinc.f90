! ======================================================================================================================================
! Build Kaiser-Windowed sinc (KSW) functions allowing for positioning of sources and extracting solution at arbitrary positions in a regular Cartesian grid
!
! Reference:
!@Article{Hicks_2002_ASR,
!  author  = {G. J. Hicks},
!  title   = {Arbitrary source and receiver positioning in finite-difference schemes using {K}aiser windowed sinc functions},
!  journal = {Geophysics},
!  year    = {2002},
!  volume  = {67},
!  pages   = {156--166},
!}
!
!
! subroutines subhicks1d, subhicks2d, subhicks3d: Build 1D/2D/3D Kaiser-windowed sinc function
! subroutine subsourcehicks: given the coordinates of the source, find the 7 attributes allowing to select the Hicks coefficients in the pre-computed table.
! Arrays coef_hicks_monopole(:,:,:,:,:,:,:), coef_hicks_dipole(:,:,:,:,:,:,:), coef_hicks_der2(:,:,:,:,:,:,:) have 7 entries.
! The first three are the index of the closest point in the grid, the next three are alpha (the shift with respect to the closest point), and the last
! store the number of grid points above free surface.
!
! ======================================================================================================================================
module kwsinc
  implicit none

  integer,parameter                 :: irad(1:3)=10
  integer                           :: nalpha(1:3)=100

  real, allocatable                 :: coef_hicks1D_monopole(:,:)
  real, allocatable                 :: coef_hicks1D_dipole(:,:)
  real, allocatable                 :: coef_hicks1D_der2(:,:)

  real, allocatable                 :: coef_hicks2D_monopole(:,:,:,:,:)
  real, allocatable                 :: coef_hicks2D_dipole(:,:,:,:,:)
  real, allocatable                 :: coef_hicks2D_der2(:,:,:,:,:)

  real, allocatable                 :: coef_hicks3D_monopole(:,:,:,:,:,:,:)
  real, allocatable                 :: coef_hicks3D_dipole(:,:,:,:,:,:,:)
  real, allocatable                 :: coef_hicks3D_der2(:,:,:,:,:,:,:)

  public BESSI0

  real, dimension(10)      :: bhicks0,bhicks1,bhicks2
  data bhicks0 /1.24,2.94,4.53,6.31,7.91,9.42,10.95,12.53,14.09,14.18/
  data bhicks1 /0.00,3.33,4.96,6.42,7.77,9.52,11.11,12.52,14.25,16.09/
  data bhicks2 /1.44,3.19,3.85,3.97,4.33,4.69,4.45,4.93,5.11,5.35/

contains

  SUBROUTINE allocatekws1d
    allocate(coef_hicks1D_monopole(-irad(1):irad(1), nalpha(1)))
    allocate(coef_hicks1D_dipole(-irad(1):irad(1), nalpha(1)))
    allocate(coef_hicks1D_der2(-irad(1):irad(1), nalpha(1)))
    coef_hicks1D_monopole(:,:) = 0.
    coef_hicks1D_dipole(:,:) = 0.
    coef_hicks1D_der2(:,:) = 0.
  END SUBROUTINE allocatekws1d

  SUBROUTINE deallocatekws1d
    deallocate(coef_hicks1D_monopole)
    deallocate(coef_hicks1D_dipole)
    deallocate(coef_hicks1D_der2)
  END SUBROUTINE deallocatekws1d

  SUBROUTINE allocatekws2d
    allocate(coef_hicks2D_monopole(-irad(1):irad(1), -irad(2):irad(2),nalpha(1), nalpha(2), 0:irad(1)))
    allocate(coef_hicks2D_dipole(-irad(1):irad(1), -irad(2):irad(2),nalpha(1), nalpha(2), 0:irad(1)))
    allocate(coef_hicks2D_der2(-irad(1):irad(1), -irad(2):irad(2),nalpha(1), nalpha(2), 0:irad(1)))
    coef_hicks2D_monopole(:,:,:,:,:) = 0.
    coef_hicks2D_dipole(:,:,:,:,:) = 0.
    coef_hicks2D_der2(:,:,:,:,:) = 0.
  END SUBROUTINE allocatekws2d

  SUBROUTINE deallocatekws2d
    deallocate(coef_hicks2D_monopole)
    deallocate(coef_hicks2D_dipole)
    deallocate(coef_hicks2D_der2)
  END SUBROUTINE deallocatekws2d

  SUBROUTINE allocatekws3d
    integer                           :: nalpha(1:3)=40
    allocate(coef_hicks3D_monopole(-irad(1):irad(1),-irad(2):irad(2),-irad(3):irad(3),nalpha(1),nalpha(2),nalpha(3),0:irad(1)))
    allocate(coef_hicks3D_dipole(-irad(1):irad(1),-irad(2):irad(2),-irad(3):irad(3),nalpha(1),nalpha(2),nalpha(3),0:irad(1)))
    allocate(coef_hicks3D_der2(-irad(1):irad(1),-irad(2):irad(2),-irad(3):irad(3),nalpha(1),nalpha(2),nalpha(3),0:irad(1)))
    coef_hicks3D_monopole(:,:,:,:,:,:,:) = 0.
    coef_hicks3D_dipole(:,:,:,:,:,:,:) = 0.
    coef_hicks3D_der2(:,:,:,:,:,:,:) = 0.
  END SUBROUTINE allocatekws3d

  SUBROUTINE deallocatekws3d
    deallocate(coef_hicks3D_monopole)
    deallocate(coef_hicks3D_dipole)
    deallocate(coef_hicks3D_der2)
  END SUBROUTINE deallocatekws3d

  SUBROUTINE kwsinc1d

    ! ==================================================================================================================================
    ! SUBROUTINE subhicks1d: 1D KSW coefficients
    ! ==================================================================================================================================

    real            :: bhicks_monopole,bhicks_dipole,bhicks_der2

    integer         :: i_alpha,nloop
    real            :: i0n_monopole,i0n_dipole,i0n_der2,i0d_monopole,i0d_dipole,i0d_der2
    real            :: xn,kaiser_monopole,kaiser_dipole,kaiser_der2

    real            :: alpha

    real, PARAMETER :: PI=3.14159265359
    real, PARAMETER :: HICKS_EPS_TEST = 1e-15,HICKS_EPS_KAISER = 1e-15

    ! Adaptation as a function of the dimension of propagation

    bhicks_monopole=bhicks0(irad(1))
    bhicks_dipole=bhicks1(irad(1))
    bhicks_der2=bhicks2(irad(1))

    ! compute Bessel func. for bb
    i0d_monopole = BESSI0(bhicks_monopole)  ! denominator, equation (6)
    i0d_dipole   = BESSI0(bhicks_dipole)
    i0d_der2     = BESSI0(bhicks_der2)

    ! Outer loop over the tabulated positions in the interval [0.5 ; -0.5[
    DO i_alpha = 1, nalpha(1)
       alpha = 0.5 - real(i_alpha - 1) / nalpha(1)

       ! Inner loop over the spatial support of the interpolating function
       DO nloop = -irad(1), irad(1)

          ! xn = algebraic distance between the source and the grid point of index n; x = n + alpha (see equation (3) in Hicks, Geophysics, 2002).
          ! nloop=0 corresponds to the nearest point to the source

          xn = real(nloop)  + alpha

          IF (ABS(xn) .le.  irad(1) ) THEN                  !This test prevents negative quantities in the square root of the argument of the Bessel function

             ! compute coef.
             IF (ABS(xn) > HICKS_EPS_TEST) THEN          ! xn ne 0 (the source is not on a grid point)

                i0n_monopole = BESSI0(bhicks_monopole * SQRT(1.-xn*xn/(irad(1)*irad(1))))         !Equation (6), numerator
                i0n_dipole   = BESSI0(bhicks_dipole   * SQRT(1.-xn*xn/(irad(1)*irad(1))))
                i0n_der2     = BESSI0(bhicks_der2     * SQRT(1.-xn*xn/(irad(1)*irad(1))))
                
                kaiser_monopole = i0n_monopole / i0d_monopole                                      !Equation 6
                kaiser_dipole   = i0n_dipole   / i0d_dipole
                kaiser_der2     = i0n_der2   / i0d_der2
                
                coef_hicks1D_monopole( nloop,  i_alpha )= kaiser_monopole * sin(PI*xn) / (PI*xn)            !Equations (4), (5)
                coef_hicks1D_dipole( nloop, i_alpha )  = kaiser_dipole * (cos(PI*xn)-(sin(PI*xn)/(PI*xn)))/(xn)   ! Equation (9)
                coef_hicks1D_der2( nloop, i_alpha ) = kaiser_der2 * ( -PI*sin(PI * xn)/(xn) &
                     - 2.*cos(PI * xn)/(xn**2) &
                     + 2.*sin(PI * xn)/(PI * xn**3) )

                ! check small values
                IF (ABS(coef_hicks1D_monopole(nloop,i_alpha)) <  HICKS_EPS_KAISER ) coef_hicks1D_monopole(nloop,i_alpha) = 0.
                IF (ABS(coef_hicks1D_dipole(nloop,i_alpha))   <  HICKS_EPS_KAISER ) coef_hicks1D_dipole(nloop,i_alpha)   = 0.
                IF (ABS(coef_hicks1D_der2(nloop,i_alpha))     <  HICKS_EPS_KAISER ) coef_hicks1D_der2(nloop,i_alpha)   = 0.

             ELSE                                                            ! xn = 0 (the source is on a grid point)

                i0n_monopole = BESSI0(bhicks_monopole)                       !Equation (6), numerator
                i0n_dipole   = BESSI0(bhicks_dipole)
                i0n_der2     = BESSI0(bhicks_der2)
                kaiser_monopole = 1.                                          !Equation 6
                kaiser_dipole   = i0n_dipole   / i0d_dipole
                kaiser_der2    = i0n_der2   / i0d_der2

                coef_hicks1D_monopole( nloop, i_alpha) = 1.
                coef_hicks1D_dipole( nloop, i_alpha)   = 0.
                coef_hicks1D_der2( nloop, i_alpha)   = -PI**2

             ENDIF

          end if

       ENDDO

    ENDDO !END DO ALPHA

  END SUBROUTINE kwsinc1d

  ! ==============================================================================================
  ! SUBROUTINE SUBHICKS2D
  ! Tabulate Kaiser-windowed sinc functions for source and receiver positioning in 2D grid
  !
  ! Reference:
  ! @Article{Hicks_2002_ASR,
  !  author  = {G. J. Hicks},
  !  title   = {Arbitrary source and receiver positioning in finite-difference schemes using {K}aiser windowed sinc functions},
  !  journal = {Geophysics},
  !  year    = {2002},
  !  volume  = {67},
  !  pages   = {156--166},
  !}
  !
  ! Author: S. Operto; email: operto@geoazur.unice.fr
  ! ==============================================================================================

  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !                                 SUBHICKS
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SUBROUTINE  kwsinc2d(ifreesurf)

    integer                 :: i_alpha,i1_alpha,i2_alpha,i1_sym
    integer                 :: nloop,nloop1,nloop2
    integer                 :: ifreesurf,ip

    logical, parameter      :: debug = .false.

    real, ALLOCATABLE       :: coef_hicks_temp_monopole(:,:,:,:)
    real, ALLOCATABLE       :: coef_hicks_temp_dipole(:,:,:,:)
    real, ALLOCATABLE       :: coef_hicks_temp_der2(:,:,:,:)

    write(*,*) 'NALPHA = ',nalpha(:)

    allocate(coef_hicks_temp_monopole(-irad(1):irad(1), -irad(2):irad(2),nalpha(1), nalpha(2)))
    allocate(coef_hicks_temp_dipole(-irad(1):irad(1), -irad(2):irad(2), nalpha(1), nalpha(2)))
    allocate(coef_hicks_temp_der2(-irad(1):irad(1), -irad(2):irad(2), nalpha(1), nalpha(2)))

    coef_hicks_temp_monopole(:,:,:,:) = 0.
    coef_hicks_temp_dipole(:,:,:,:) = 0.
    coef_hicks_temp_der2(:,:,:,:) = 0.

    ! ====================================================================================================================
    ! Tabulation of Hicks coefficients
    ! For alpha = 0, the source is on the grid node. alpha=0 for i_alpha=nalpha(1)/2+1
    ! alpha \in ]-0.5;0.5]
    ! alpha decrases from 0.5 to -0.5 from the left to the right (see the loop below)
    ! =====================================================================================================================

    !1D coefficients
    CALL allocatekws1d
    CALL kwsinc1d

    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !=============================================================
    ! compute the coef. in 2D without interaction of freee surface
    !==============================================================
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################

    ! Tensorial construction of the 3D interpolation function

    ! loop on the spatial samples
    DO i2_alpha = 1, nalpha(2)
       DO i1_alpha = 1, nalpha(1)

          ! loop on the support of the interpolation
          DO nloop2 = -irad(2) ,irad(2)
             DO nloop1 = -irad(1) ,irad(1)
                coef_hicks_temp_monopole(nloop1, nloop2, i1_alpha, i2_alpha) = &
                     coef_hicks1D_monopole(nloop2,i2_alpha) * coef_hicks1D_monopole(nloop1,i1_alpha)
                coef_hicks_temp_dipole(nloop1, nloop2, i1_alpha, i2_alpha) = &
                     coef_hicks1D_monopole(nloop2,i2_alpha) * coef_hicks1D_dipole(nloop1,i1_alpha)
                coef_hicks_temp_der2(nloop1, nloop2, i1_alpha, i2_alpha) = &
                     coef_hicks1D_monopole(nloop2,i2_alpha) * coef_hicks1D_der2(nloop1,i1_alpha)
             ENDDO
          ENDDO

       ENDDO
    ENDDO

    IF ( debug ) THEN
       WRITE(*,*) 'coef_hicks_temp_monopole(:,:,1,1)', coef_hicks_temp_monopole(:,:,1,1)
       WRITE(*,*) 'coef_hicks_temp_dipole(:,:,1,1)',   coef_hicks_temp_dipole(:,:,1,1)
       WRITE(*,*) 'coef_hicks_temp_der2(:,:,1,1)',     coef_hicks_temp_der2(:,:,1,1)
    ENDIF

    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    ! Dynamic allocation for the output array coef_hicks; the seven^th field is used to manage free surface interaction
    ! at most, there are irad - 1 non zero coefficients above free surface (the coefficient on the free surface is 0)
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################

    ! test if free surface

    ! No interaction with free surface; field 7 is 0
    !==========================================
    coef_hicks2D_monopole(:,:,:,:, 0) = coef_hicks_temp_monopole(:,:,:,:)
    coef_hicks2D_dipole(:,:,:,:, 0)  = coef_hicks_temp_dipole(:,:,:,:)
    coef_hicks2D_der2(:,:,:,:, 0)  = coef_hicks_temp_der2(:,:,:,:)

    !=========================================================================================
    ! Free-surface interaction; the field 7 is ip index, the number of points above the free surface
    !=========================================================================================

    IF ( iabs(ifreesurf) == 1 ) THEN

       ! loop over the possible number of nodes above free surface (nb interpolated points above free surface)
       ! ip: number of nodes above free surface

       DO ip = 1, irad(1)
          ! WRITE(*,*) 'ip = ',ip

          ! loop on the spatial samples
          DO i2_alpha = 1, nalpha(2)
             DO i1_alpha = 1, nalpha(1)
                !      alpha1 = -0.5 + real(i1_alpha) / real(nalpha(1))

                ! loop on the support of the interpolation
                DO nloop2 = -irad(2) ,irad(2)
                   DO nloop1 = -irad(1) ,irad(1)

                      ! Initialization
                      coef_hicks2D_monopole(nloop1, nloop2, i1_alpha, i2_alpha, ip) = &
                           coef_hicks_temp_monopole(nloop1, nloop2, i1_alpha, i2_alpha)
                      coef_hicks2D_dipole(nloop1, nloop2, i1_alpha, i2_alpha, ip)  = &
                           coef_hicks_temp_dipole(nloop1, nloop2, i1_alpha, i2_alpha)
                      coef_hicks2D_der2(nloop1, nloop2, i1_alpha, i2_alpha, ip)= &
                           coef_hicks_temp_der2(nloop1, nloop2, i1_alpha, i2_alpha)

                      ! Index of the node located at the free surface
                      i1_sym = -irad(1) + ip
                      coef_hicks2D_monopole(i1_sym,nloop2,i1_alpha,i2_alpha,ip)=0.   !set the source to 0 at the free surface

                      ! apply change of coef. below free surface
                      DO nloop = 1, ip

                         ! Mirroring of coefficients below free surface (anti-symmetry)
                         coef_hicks2D_monopole(i1_sym+nloop,nloop2, i1_alpha, i2_alpha,ip) = &
                              coef_hicks_temp_monopole(i1_sym+nloop,nloop2, i1_alpha, i2_alpha)  &
                              - coef_hicks_temp_monopole(i1_sym-nloop,nloop2,i1_alpha, i2_alpha)

                         coef_hicks2D_dipole(i1_sym+nloop,nloop2, i1_alpha, i2_alpha,ip) = &
                              coef_hicks_temp_dipole(i1_sym+nloop,nloop2, i1_alpha, i2_alpha) &
                              - coef_hicks_temp_dipole(i1_sym-nloop,nloop2, i1_alpha, i2_alpha)

                         coef_hicks2D_der2(i1_sym+nloop,nloop2, i1_alpha, i2_alpha,ip) =   &
                              coef_hicks_temp_der2(i1_sym+nloop,nloop2, i1_alpha, i2_alpha)  &
                              - coef_hicks_temp_der2(i1_sym-nloop,nloop2,i1_alpha, i2_alpha)

                      ENDDO

                   END DO
                END DO

             END DO
          END DO

       ENDDO   ! IP

       ! ========================

    ENDIF ! IF ( ifreesurf == 1) THEN

    ! Deallocate memory for 1D coefficients
    CALL deallocatekws1d

    ! deallocate local tables
    deallocate(coef_hicks_temp_monopole)
    deallocate(coef_hicks_temp_dipole)
    deallocate(coef_hicks_temp_der2)

  END SUBROUTINE kwsinc2d

  ! ==============================================================================================
  ! SUBROUTINE SUBHICKS3D
  ! Tabulate Kaiser-windowed sinc functions for source and receiver positioning in 3D grid
  !
  ! Reference:
  ! @Article{Hicks_2002_ASR,
  !  author  = {G. J. Hicks},
  !  title   = {Arbitrary source and receiver positioning in finite-difference schemes using {K}aiser windowed sinc functions},
  !  journal = {Geophysics},
  !  year    = {2002},
  !  volume  = {67},
  !  pages   = {156--166},
  !}
  !
  ! Author: S. Operto; email: operto@geoazur.unice.fr
  ! ==============================================================================================
  !
  SUBROUTINE kwsinc3d(ifreesurf)

    integer  :: i_alpha,i1_alpha,i2_alpha,i3_alpha,i1_sym
    integer  :: nloop,nloop1,nloop2,nloop3
    integer  :: ifreesurf,ip

    real       :: bhicks_monopole, bhicks_dipole, bhicks_der2
    real       :: i0d_monopole,i0d_dipole,i0d_der2

    logical, parameter  :: debug = .false.

    real, ALLOCATABLE  :: coef_hicks_temp_monopole(:,:,:,:,:,:)
    real, ALLOCATABLE  :: coef_hicks_temp_dipole(:,:,:,:,:,:)
    real, ALLOCATABLE  :: coef_hicks_temp_der2(:,:,:,:,:,:)

    allocate(coef_hicks_temp_monopole(                                     &
         -irad(1):irad(1), -irad(2):irad(2), -irad(3):irad(3), &
         nalpha(1), nalpha(2), nalpha(3)))
    allocate(coef_hicks_temp_dipole(                                       &
         -irad(1):irad(1), -irad(2):irad(2), -irad(3):irad(3), &
         nalpha(1), nalpha(2), nalpha(3)))
    allocate(coef_hicks_temp_der2(                                         &
         -irad(1):irad(1), -irad(2):irad(2), -irad(3):irad(3), &
         nalpha(1), nalpha(2), nalpha(3)))

    coef_hicks_temp_monopole(:,:,:,:,:,:) = 0.
    coef_hicks_temp_dipole(:,:,:,:,:,:) = 0.
    coef_hicks_temp_der2(:,:,:,:,:,:) = 0.

    ! Tabulation of Hicks coefficients
    ! For alpha = 0, the source is on the grid node. alpha=0 for i_alpha=nalpha(1)/2+1
    ! alpha \in ]-0.5;0.5]
    ! alpha decrases from 0.5 to -0.5 from the left to the right (see the loop below)

    !1D coefficients
    CALL allocatekws1d
    CALL kwsinc1d

    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !=============================================================
    ! compute the coef. in 3D without interaction of freee surface
    !==============================================================
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################

    ! Tensorial construction of the 3D interpolation function
    ! For 1D and 2D cases, note that nloop3 and/or nloop2 range between 0 and 0.

    WRITE(*,*) 'TENSORIAL CONSTRUCTION'

    ! loop on the spatial samples
    DO i3_alpha = 1, nalpha(3)
       DO i2_alpha = 1, nalpha(2)
          DO i1_alpha = 1, nalpha(1)

             ! loop on the support of the interpolation
             DO nloop3 = -irad(3) ,irad(3)
                DO nloop2 = -irad(2) ,irad(2)
                   DO nloop1 = -irad(1) ,irad(1)
                      coef_hicks_temp_monopole(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha) = &
                           coef_hicks1D_monopole(nloop3,i3_alpha) *  &
                           coef_hicks1D_monopole(nloop2,i2_alpha) * &
                           coef_hicks1D_monopole(nloop1,i1_alpha)
                      coef_hicks_temp_dipole(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha) = &
                           coef_hicks1D_monopole(nloop3,i3_alpha) * &
                           coef_hicks1D_monopole(nloop2,i2_alpha) * &
                           coef_hicks1D_dipole(nloop1,i1_alpha)
                      coef_hicks_temp_der2(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha) = &
                           coef_hicks1D_monopole(nloop3,i3_alpha) * &
                           coef_hicks1D_monopole(nloop2,i2_alpha) * &
                           coef_hicks1D_der2(nloop1,i1_alpha)
                   ENDDO
                ENDDO
             ENDDO

          ENDDO
       ENDDO
    ENDDO

    IF ( debug ) THEN
       WRITE(*,*) 'coef_hicks_temp_monopole(:,:,:,1,1,1)', coef_hicks_temp_monopole(:,:,:,1,1,1)
       WRITE(*,*) 'coef_hicks_temp_dipole(:,:,:,1,1,1)',   coef_hicks_temp_dipole(:,:,:,1,1,1)
       WRITE(*,*) 'coef_hicks_temp_der2(:,:,:,1,1,1)',     coef_hicks_temp_der2(:,:,:,1,1,1)
    ENDIF

    ! =============================================================================
    ! PROCESS FREE SURFACE
    ! =============================================================================

    ! No interaction with free surface; field 7 is 0
    !==========================================
    coef_hicks3D_monopole(:,:,:,:,:,:, 0) = coef_hicks_temp_monopole(:,:,:,:,:,:)
    coef_hicks3D_dipole(:,:,:,:,:,:, 0)  = coef_hicks_temp_dipole(:,:,:,:,:,:)
    coef_hicks3D_der2(:,:,:,:,:,:, 0)  = coef_hicks_temp_der2(:,:,:,:,:,:)

    !=========================================================================================
    ! Free-surface interaction; the field 7 is ip index, the number of points above the free surface
    !=========================================================================================

    IF ( iabs(ifreesurf) == 1 ) THEN

       ! loop over the possible number of nodes above free surface (nb interpolated points above free surface)
       ! ip: number of nodes above free surface

       DO ip = 1, irad(1)
          WRITE(*,*) 'ip = ',ip

          ! loop on the spatial samples
          DO i3_alpha = 1, nalpha(3)
             DO i2_alpha = 1, nalpha(2)
                DO i1_alpha = 1, nalpha(1)

                   ! loop on the support of the interpolation
                   DO nloop3 = -irad(3) ,irad(3)
                      DO nloop2 = -irad(2) ,irad(2)
                         DO nloop1 = -irad(1) ,irad(1)

                            ! Initialization
                            coef_hicks3D_monopole(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha, ip) = &
                                 coef_hicks_temp_monopole(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha)
                            coef_hicks3D_dipole(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha, ip)  = &
                                 coef_hicks_temp_dipole(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha)
                            coef_hicks3D_der2(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha, ip)  = &
                                 coef_hicks_temp_der2(nloop1, nloop2, nloop3, i1_alpha, i2_alpha, i3_alpha)

                            ! Index of the node located at the free surface
                            i1_sym = -irad(1) + ip
                            coef_hicks3D_monopole(i1_sym,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha,ip)=0.

                            ! apply change of coef. below free surface
                            DO nloop = 1, ip

                               ! Mirroring of coefficients below free surface (anti-symmetry)
                               coef_hicks3D_monopole(i1_sym+nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha,ip)=  &
                                    coef_hicks_temp_monopole(i1_sym+nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha)   &
                                    - coef_hicks_temp_monopole(i1_sym-nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha)
                               coef_hicks3D_dipole(i1_sym+nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha,ip)=          &
                                    coef_hicks_temp_dipole(i1_sym+nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha)     &
                                    - coef_hicks_temp_dipole(i1_sym-nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha)
                               coef_hicks3D_der2(i1_sym+nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha,ip)=            &
                                    coef_hicks_temp_der2(i1_sym+nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha)       &
                                    - coef_hicks_temp_der2(i1_sym-nloop,nloop2,nloop3,i1_alpha,i2_alpha,i3_alpha)

                            ENDDO

                         END DO
                      END DO
                   END DO

                END DO
             END DO
          END DO

       ENDDO   ! IP

    ENDIF ! IF ( ifreesurf == 1) THEN

    ! Deallocate memory for 1D coefficients
    call deallocatekws1d

    ! deallocate local tables
    deallocate(coef_hicks_temp_monopole )
    deallocate(coef_hicks_temp_dipole )
    deallocate(coef_hicks_temp_der2 )

  END SUBROUTINE kwsinc3d

  ! ==============================================================================================
  ! SUBROUTINE SUBSOURCEKWS
  ! Assign ix(ndim), ialpha(ndim) and nalpha_frees to a source in order to pick suitables values of KWS coefficients in pre-tabulated table
  ! ix(ndim): indexes of the closest grid point to the source
  ! ialpha(ndim): index of alpha in a array of dim nalpha
  ! nalpha_frees: number of grid points above the free surface
  ! Input: ndim, nalpha(ndim),xs(ndim),dz
  ! Output: ix(ndim), ialpha(ndim),nalpha_frees
  ! Author: S. Operto; email: operto@geoazur.unice.fr
  ! ==============================================================================================

  SUBROUTINE subsourcekws(ndim,xs,dz,ifs,ix,ialpha,nalpha_frees)
    integer, INTENT(IN)   :: ndim,ifs
    real, INTENT(IN)      :: xs(ndim),dz
    
    integer, INTENT(OUT)  :: ialpha(ndim),nalpha_frees
    integer, INTENT(OUT)  :: ix(ndim)
    
    integer               :: i
    real                  :: x0,alpha

    DO i=1,ndim
       ix(i)=NINT(xs(i)/dz)+1                          !index of the closest grid point
       x0=FLOAT(ix(i)-1)*dz                            !coordinates of the closest grid point
       alpha = ( x0 - xs(i) ) / dz                     !alpha
       ialpha(i)=NINT(( 0.5-alpha ) * nalpha(i) ) + 1  !index of alpha in a array od dim nalpha
       ialpha(i)=MIN(ialpha(i),nalpha(i))
    END DO

    nalpha_frees = irad(1)-NINT(xs(1)/dz)
    nalpha_frees = MAX(nalpha_frees,0)

    if (ifs==0) nalpha_frees = 0
    !   write(*,*) 'nalpha_frees = ',nalpha_frees

  END SUBROUTINE subsourcekws

  ! ==============================================================================================
  !            BESSEL FUNCTION
  ! ==============================================================================================

  FUNCTION BESSI0(X)

    !--------------------------------------- 
    ! Compute BESSEL function
    !---------------------------------------
    REAL BESSI0
    REAL X, AX
    REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
    DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
    DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,0.225319D-2,-0.157565D-2,0.916281D-2,&
         -0.2057706D-1,0.2635537D-1,-0.1647633D-1,0.392377D-2/

    IF (ABS(X).LT.3.75) THEN
       Y=(X/3.75)**2
       BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
    ELSE
       AX=ABS(X)
       Y=3.75/AX
       BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
    ENDIF

    RETURN 

  END FUNCTION BESSI0

end module kwsinc
