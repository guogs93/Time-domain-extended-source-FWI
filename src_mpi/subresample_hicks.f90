! =================================================================================
! PROGRAM HICKSINTERPOLATION
! Resampling 1D function with Kaiser-winded sinc function
! The code was originally writte, to resample at  éNyquist rate finite-sampled seismograms computed by finite-difference method
! A demo is provided to check the code
! namein nameout: input output files
! n1in n2: number of samples number of traces in the input file. Resampling is applied along dimension 1
! irad nalpha: Radius of the windowed Sinc function in number of points of the input signal
! nalpha: number of intervals sampling a sampling interval
!
! Reference
! @Article{Hicks_2002_ASR,
!  author  = {G. J. Hicks},
!  title   = {Arbitrary source and receiver positioning in finite-difference schemes using {K}aiser windowed sinc functions},
!  journal = {Geophysics},
!  year    = {2002},
!  volume  = {67},
!  pages   = {156--166},
!}
!
! https://www.geoazur.fr/svn/dav/wind_articles/TECHREPORT/TR11_Operto_HicksSource/TR011_2020_Operto_HicksSource_v1.pdf
!
! ========================================================================================

subroutine subresample(xin,xout,n1in,n2,n1out,d1in,d1out)
  implicit none

  real                       :: d1in,d1out,tmax
  integer                    :: n1in,n2,n1out
  integer                    :: irad,nalpha
  real                       :: xin(n1in,n2),xout(n1out,n2)

  real, allocatable          :: coef_hicks1D_monopole(:,:)

  real                       :: tout,t0,alpha,coeff,alpha_check,dalpha
  integer                    :: ialpha,iradmin,iradmax

  integer                    :: i,i1,i2,j,it0

  irad=4
  nalpha=10

  ! ==========================================================================
  
  dalpha=1./float(nalpha)

  allocate(coef_hicks1D_monopole(-irad:irad,nalpha))

  ! Precomputation of the (2*irad+1) windowed-sinc coefficients as a function of alpha
  call subhicks1d(coef_hicks1D_monopole,irad,nalpha)

  ! Loop over traces
  do i2=1,n2
     xout(:,i2)=0.
     do i1=1,n1out
        tout=float(i1-1)*d1out                                  !time of the sample where the function needs to be estimated
        it0=nint(tout/d1in)+1                                   !index of the closest sample in the original sampling
        t0=float(it0-1)*d1in                                    !time  of the closest sample
        alpha = ( t0 - tout ) / d1in                            !alpha
        ialpha=nint(( 0.5 - alpha ) * nalpha )+1                !index of alpha
        alpha_check = 0.5 - REAL(ialpha - 1) / nalpha
        if (ialpha.gt.nalpha) ialpha=nalpha

        iradmin=MAX(1-it0,   -irad)
        iradmax=MIN(n1in-it0, irad)

        ! internal loop over the points sampling the Kaiser windowed sinc
        
        DO j=iradmin,iradmax
           coeff = coef_hicks1D_monopole(j,ialpha)
           xout(i1,i2)=xout(i1,i2)+coeff*xin(it0+j,i2)
        END DO

     end do
  end do

  !write(*,*) "================================="
  !write(*,*) 'N1OUT = ',n1out
  !write(*,*) "================================="
  
  deallocate(coef_hicks1D_monopole)

end subroutine

! ====================================================================================================

subroutine subhicks1d(coef_hicks1D_monopole,irad,nalpha)
  implicit none

  real        :: alpha
  real        :: BESSI0
  integer     :: irad,nalpha
  integer     :: j
  integer     :: i_alpha,nloop
  real        :: i0n_monopole,i0d_monopole
  real        :: xn,kaiser_monopole
  real        :: bhicks_monopole
  real        :: coef_hicks1D_monopole(-irad:irad,nalpha)

  real, parameter   :: PI=3.14159265359
  real, parameter   :: HICKS_EPS_TEST=1e-15,HICKS_EPS_KAISER=1e-15

  REAL,DIMENSION(10)      :: bhicks0,bhicks1
  DATA bhicks0 /1.24,2.94,4.53,6.31,7.91,9.42,10.95,12.53,14.09,14.18/

  ! Hicks coefficient for Kaiser function
  bhicks_monopole=bhicks0(irad)

  !   write(*,*) 'BHICKS_MONOPOLE = ',bhicks_monopole

  ! compute Bessel func. for bb
  i0d_monopole = BESSI0(bhicks_monopole)  ! denominator, equation (6)

  ! Tabulation of Hicks coefficients
  ! For alpha = 0, the source is on the grid node. alpha=0 for i_alpha=nalpha(1)/2+1
  ! alpha \in ]-0.5;0.5]
  ! alpha decrases from 0.5 to -0.5 from the left to the right (see the loop below)

  ! Outer loop over the tabulated positions in the interval [0.5 ; -0.5[
  do i_alpha = 1, nalpha
     alpha = 0.5 - REAL(i_alpha - 1) / nalpha

     ! Inner loop over the spatial support of the interpolating function

     do nloop = -irad,irad

        ! xn = algebraic distance between the source and the grid point of index n; x = n + alpha (see equation (3) in Hicks, Geophysics, 2002).
        ! nloop=0 corresponds to the nearest point to the source

        xn = real(nloop)  + alpha

        if (abs(xn) .le.  irad ) then                  !This test prevents negative quantities in the square root of the argument of the Bessel function

           ! compute coef.
           if (abs(xn) > HICKS_EPS_TEST) then             ! xn ne 0 (the source is not on a grid point)

              i0n_monopole = BESSI0(bhicks_monopole * sqrt(1.-xn*xn/(irad*irad)))         		!Equation (6), numerator
              kaiser_monopole = i0n_monopole / i0d_monopole                                      	!Equation 6, denominator
              coef_hicks1D_monopole( nloop,  i_alpha )= kaiser_monopole * sin(PI*xn) / (PI*xn)            		!Equations (4), (5)

              ! check small values
              if (abs(coef_hicks1D_monopole(nloop,i_alpha)) <  HICKS_EPS_KAISER ) coef_hicks1D_monopole(nloop,i_alpha) = 0.

           else            ! xn = 0 (the source is on a grid point)

              i0n_monopole = BESSI0(bhicks_monopole)     				!Equation (6), numerator
              kaiser_monopole = 1.                                                     !Equation 6

              coef_hicks1D_monopole( nloop, i_alpha) = 1.
              if (i_alpha == nalpha/2) WRITE(*,*) "coeff1d ",coef_hicks1D_monopole(nloop,  i_alpha )


           end if

        end if

     end do

  end do !END DO ALPHA

end subroutine subhicks1d

! ---------------------------------------------------------------------------------------------------------------
! FUNCTION BESSI0 (FROM NUMERICAL RECIPES)
! ---------------------------------------------------------------------------------------------------------------

FUNCTION BESSI0(X)

  !---------------------------------------
  ! Compute BESSEL function
  !---------------------------------------

  REAL BESSI0, X, AX
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

