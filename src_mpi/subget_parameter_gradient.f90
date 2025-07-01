! ==========================================================================
! Compute the gradient of m
! ==========================================================================
subroutine subget_parameter_gradient(model, nx, ny, dx, dy, grad)
  implicit none
  integer :: nx, ny
  real :: dx, dy
  real :: model(nx,ny)
  real :: grad(nx,ny,2)
  integer :: ix, iy
  real    :: idx, idy
 
  !! Prepare
  idx = 1. / dx
  idy = 1. / dy
    
  grad=0
    
  !! Compute gradient
  do iy = 1, ny-1
     do ix = 1, nx-1
        grad(ix, iy, 1) = (model(ix+1, iy  ) - model(ix, iy))
        grad(ix, iy, 2) = (model(ix  , iy+1) - model(ix, iy))
     end do
  end do
   
end subroutine subget_parameter_gradient

!
! compute the adjoint gradient of m
!
subroutine subget_parameter_adjoint_gradient(grad, nx, ny, dx, dy, model)
  implicit none

  integer :: nx, ny
  real :: dx, dy
  real :: model(nx,ny)
  real :: grad(nx,ny,2)

  integer :: ix, iy
  real    :: idx, idy, d_dx, d_dy

  !! Prepare
  idx = 1. / dx
  idy = 1. / dy

  !! Compute adjoint of gradient (divergence grosso merdo)
  model(1,:) = 0.   ! check bnd condition
  model(:,1) = 0.   ! check bnd condition
  model(nx,:) = 0.
  model(:,ny) =0 

  do iy = 2, ny-1
     do ix = 2, nx-1
        d_dx          = (grad(ix, iy, 1) - grad(ix-1, iy  , 1))
        d_dy          = (grad(ix, iy, 2) - grad(ix  , iy-1, 2))
        model(ix, iy) = d_dx + d_dy
     end do
  end do
  
end subroutine subget_parameter_adjoint_gradient

SUBROUTINE subderivop(n1, n2, D1, iD1, jD1, nnz1, D2, iD2, jD2, nnz2, val)

  IMPLICIT NONE

  integer, INTENT(IN) :: n1, n2, val, nnz1, nnz2
  integer             :: i1, i2, k, l, nnz01, nnz02

  integer             :: iD1(nnz1),jD1(nnz1),iD2(nnz2),jD2(nnz2)
  real             :: D1(nnz1),D2(nnz2)

  !#############################################################################################################

  !val=1 ==> First order derivatives Dz, Dx
  !First order forward difference
  IF (val .EQ. 1 ) THEN

     nnz01=0
     DO i2=2,n2-1
        DO i1=2,n1-1

           nnz01=nnz01+1
           k=(i2-1)*n1 + i1

           D1(nnz01)=-1
           iD1(nnz01)=k
           jD1(nnz01)=k

           nnz01=nnz01+1
           l=(i2-1)*n1 + (i1+1)
           D1(nnz01)=1
           iD1(nnz01)=k
           jD1(nnz01)=l

        END DO

     END DO
     WRITE(*,*) "nnz01 and nnz1 ==",nnz01, nnz1
     IF ( nnz01 .NE. nnz1 ) THEN 
        WRITE(*,*) " ERROR nnz01 =/= nnz1 "
     END IF

     nnz02=0
     DO i2=2,n2-1
        DO i1=2,n1-1

           nnz02=nnz02+1
           k=(i2-1)*n1 + i1

           D2(nnz02)=-1
           iD2(nnz02)=k
           jD2(nnz02)=k

           nnz02=nnz02+1
           l=((i2+1)-1)*n1 + i1
           D2(nnz02)=1
           iD2(nnz02)=k
           jD2(nnz02)=l

        END DO
     END DO

     WRITE(*,*) "nnz02 and nnz2 ==",nnz02, nnz2
     IF ( nnz02 .NE. nnz2 ) THEN
        WRITE(*,*) " ERROR nnz02 =/= nnz2 "
     END IF

     !val=2 ==> Second order derivatives DDz, DDx
     !Second order forward difference
  ELSE

     nnz01=0
     DO i2=1,n2
        DO i1=1,n1-2

           nnz01=nnz01+1
           k=(i2-1)*n1 + i1

           D1(nnz01)=1
           iD1(nnz01)=k
           jD1(nnz01)=k

           nnz01=nnz01+1
           l=(i2-1)*n1 + (i1+1)
           D1(nnz01)=-2
           iD1(nnz01)=k
           jD1(nnz01)=l

           nnz01=nnz01+1
           l=(i2-1)*n1 + (i1+2)
           D1(nnz01)=1
           iD1(nnz01)=k
           jD1(nnz01)=l

        END DO

     END DO
     WRITE(*,*) "nnz01 and nnz1 ==",nnz01, nnz1
     IF ( nnz01 .NE. nnz1 ) THEN 
        WRITE(*,*) " ERROR nnz01 =/= nnz1 "
     END IF

     nnz02=0
     DO i2=1,n2-1
        DO i1=1,n1

           nnz02=nnz02+1
           k=(i2-1)*n1 + i1

           D2(nnz02)=1
           iD2(nnz02)=k
           jD2(nnz02)=k

           nnz02=nnz02+1
           l=((i2+1)-1)*n1 + i1
           D2(nnz02)=-2
           iD2(nnz02)=k
           jD2(nnz02)=l

           nnz02=nnz02+1
           l=((i2+2)-1)*n1 + i1
           D2(nnz02)=1
           iD2(nnz02)=k
           jD2(nnz02)=l

        END DO
     END DO
     WRITE(*,*) "nnz02 and nnz2 ==",nnz02, nnz2
     IF ( nnz02 .NE. nnz2 ) THEN 
        WRITE(*,*) " ERROR nnz02 =/= nnz2 "
     END IF

  END IF

END SUBROUTINE subderivop


SUBROUTINE subderivopt(n1, n2, D1, iD1, jD1, nnz1, D2, iD2, jD2, nnz2, val)

  IMPLICIT NONE

  integer, INTENT(IN) :: n1, n2, val, nnz1, nnz2
  integer             :: i1, i2, k, l, nnz01, nnz02

  integer             :: iD1(nnz1),jD1(nnz1),iD2(nnz2),jD2(nnz2)
  real             :: D1(nnz1),D2(nnz2)

  !#############################################################################################################

  !val=1 ==> First order derivatives Dz, Dx
  !First order forward difference
  IF (val .EQ. 1 ) THEN

     nnz01=0
     DO i2=2,n2
        DO i1=2,n1
           nnz01=nnz01+1
           k=(i2-1)*n1 + (i1)
           l=(i2-1)*n1 + (i1-1)
           
           D1(nnz01)=-1
           iD1(nnz01)=k
           jD1(nnz01)=l

           nnz01=nnz01+1
           
           D1(nnz01)=1
           iD1(nnz01)=k
           jD1(nnz01)=k
        END DO

     END DO
     WRITE(*,*) "nnz01 and nnz1 ==",nnz01, nnz1
     IF ( nnz01 .NE. nnz1 ) THEN 
        WRITE(*,*) " ERROR nnz01 =/= nnz1 "
     END IF

     nnz02=0
     DO i2=2,n2
        DO i1=2,n1
           nnz02=nnz02+1
           k=(i2-1)*n1 + i1 
           l=(i2-2)*n1 + i1
           D2(nnz02)=1
           iD2(nnz02)=k
           jD2(nnz02)=k
           
           nnz02=nnz02+1
           
           D2(nnz02)=-1
           iD2(nnz02)=k
           jD2(nnz02)=l

        END DO
     END DO

     WRITE(*,*) "nnz02 and nnz2 ==",nnz02, nnz2
     IF ( nnz02 .NE. nnz2 ) THEN
        WRITE(*,*) " ERROR nnz02 =/= nnz2 "
     END IF

     !val=2 ==> Second order derivatives DDz, DDx
     !Second order forward difference
  ELSE

     nnz01=0
     DO i2=1,n2
        DO i1=1,n1-2

           nnz01=nnz01+1
           k=(i2-1)*n1 + i1

           D1(nnz01)=1
           iD1(nnz01)=k
           jD1(nnz01)=k

           nnz01=nnz01+1
           l=(i2-1)*n1 + (i1+1)
           D1(nnz01)=-2
           iD1(nnz01)=k
           jD1(nnz01)=l

           nnz01=nnz01+1
           l=(i2-1)*n1 + (i1+2)
           D1(nnz01)=1
           iD1(nnz01)=k
           jD1(nnz01)=l

        END DO

     END DO
     WRITE(*,*) "nnz01 and nnz1 ==",nnz01, nnz1
     IF ( nnz01 .NE. nnz1 ) THEN 
        WRITE(*,*) " ERROR nnz01 =/= nnz1 "
     END IF

     nnz02=0
     DO i2=1,n2-1
        DO i1=1,n1

           nnz02=nnz02+1
           k=(i2-1)*n1 + i1

           D2(nnz02)=1
           iD2(nnz02)=k
           jD2(nnz02)=k

           nnz02=nnz02+1
           l=((i2+1)-1)*n1 + i1
           D2(nnz02)=-2
           iD2(nnz02)=k
           jD2(nnz02)=l

           nnz02=nnz02+1
           l=((i2+2)-1)*n1 + i1
           D2(nnz02)=1
           iD2(nnz02)=k
           jD2(nnz02)=l

        END DO
     END DO
     WRITE(*,*) "nnz02 and nnz2 ==",nnz02, nnz2
     IF ( nnz02 .NE. nnz2 ) THEN 
        WRITE(*,*) " ERROR nnz02 =/= nnz2 "
     END IF

  END IF

END SUBROUTINE subderivopt

subroutine buildLHS(n1,n2,nn,preco,ro,nnzreg,reg_coo,i_reg_coo,j_reg_coo,nnzregpreco)
 implicit none
 integer:: n1, n2, nn, nnzreg, nnzregpreco, unit, i, job, ierr, nnzp
 real:: ro
 real:: preco(nn),reg_coo(nnzreg)
 integer:: i_reg_coo(nnzreg),j_reg_coo(nnzreg)

 real,   dimension(:),allocatable:: reg_csr, preco_csr, regpreco_csr,regpreco_coo
 integer,dimension(:),allocatable:: i_reg_csr, j_reg_csr, ipreco, jpreco, i_preco_csr, j_preco_csr
 integer,dimension(:),allocatable:: i_regpreco_csr, j_regpreco_csr, i_regpreco_coo, j_regpreco_coo
 integer,dimension(:),allocatable:: iw , ndegr

 allocate(reg_csr(nnzreg),i_reg_csr(nn+1),j_reg_csr(nnzreg))

 !convert reg into CSR format ( COO ==> CSR)
 call coocsr2 (nn, nnzreg, reg_coo,i_reg_coo,j_reg_coo, reg_csr,j_reg_csr,i_reg_csr)

 reg_csr=ro*reg_csr

 allocate(ipreco(nn),jpreco(nn))

 do i=1, nn
     ipreco(i)=i
     jpreco(i)=i
 enddo

 allocate(preco_csr(nn),i_preco_csr(nn+1),j_preco_csr(nn))
 !convert preco into CSR format ( COO ==> CSR)
 call coocsr2 (nn, nn, preco,ipreco,jpreco, preco_csr,j_preco_csr,i_preco_csr)
 
 allocate(ndegr(nn))
 allocate(iw(nn))
 ! nnzregpreco with regpreco = preco + reg
 call aplbdg(nn,nn,j_preco_csr,i_preco_csr,nn,j_reg_csr,i_reg_csr,nnzreg,ndegr,nnzregpreco,iw)
 deallocate(ndegr)
 deallocate(iw)
 write(*,*) 'nnz(preco + reg) = ',nnzregpreco

 allocate(regpreco_csr(nnzregpreco),j_regpreco_csr(nnzregpreco),i_regpreco_csr(nn+1))
 allocate(iw(nn))
 job=1
 call aplb(nn,nn,job,preco_csr,j_preco_csr,i_preco_csr,nn,reg_csr,j_reg_csr,i_reg_csr,nnzreg,&
      regpreco_csr,j_regpreco_csr,i_regpreco_csr,nnzregpreco,iw,ierr)           
 deallocate(iw)
 write(*,*) 'preco + reg, OK !!!'

 !convert preco_csr to COO format ( COO --> CSR)
 allocate(regpreco_coo(nnzregpreco),i_regpreco_coo(nnzregpreco),j_regpreco_coo(nnzregpreco))

 job=3
 call csrcoo (nn,job,nnzregpreco,regpreco_csr,j_regpreco_csr,i_regpreco_csr,nnzp,&
       regpreco_coo,i_regpreco_coo,j_regpreco_coo,ierr)

 if (ierr .ne. 0) print *,'ERROR while converting preco_csr, CSR ==> COO '
 if (nnzp .ne. nnzregpreco) print *,'ERROR while converting preco_csr, nnzp =/=nnzPreco '
 write(*,*) 'convert preco_csr en format COO, OK !!!'

 open(NEWUNIT=unit,file='mat_regpreco.bin',access='STREAM',form='UNFORMATTED',action='WRITE')
 write(unit) (regpreco_coo(i),i=1,nnzregpreco)
 close(unit)
 open(NEWUNIT=unit,file='IRN_regpreco.bin',access='STREAM',form='UNFORMATTED',action='WRITE')
 write(unit) (i_regpreco_coo(i),i=1,nnzregpreco)
 close(unit)
 open(NEWUNIT=unit,file='JCN_regpreco.bin',access='STREAM',form='UNFORMATTED',action='WRITE')
 write(unit) (j_regpreco_coo(i),i=1,nnzregpreco)
 close(unit)

end subroutine


subroutine buildHess_reg(n1,n2,nn,nnzR1)
 implicit none
 integer:: n1, n2
 integer:: Hessreg(nn,nn)
 integer:: nn, i, j, val=1, ierr, job, unit
 integer:: nnzDh, nnzDv, nnzDDh, nnzDDv, nnzR1, nnzR12

 integer,dimension(:),allocatable :: iDh,jDh,iDv,jDv,iDht,jDht,iDvt,jDvt
 real,dimension(:),allocatable :: Dh, Dv, Dht, Dvt

 integer,dimension(:),allocatable :: jDv_csr,iDv_csr,jDvt_csr,iDvt_csr
 integer,dimension(:),allocatable :: jDh_csr,iDh_csr,jDht_csr,iDht_csr
 integer,dimension(:),allocatable :: jDDh_csr,iDDh_csr,jDDv_csr,iDDv_csr 
 real,dimension(:),allocatable :: Dv_csr,Dh_csr,Dvt_csr,Dht_csr,DDh_csr,DDv_csr
 
 integer,dimension(:),allocatable :: iDDh_coo,jDDv_coo,iDDv_coo,iR1,jR1,i_reg_coo,j_reg_coo 
 real,dimension(:),allocatable :: R1,reg_coo 

 integer,dimension(:),allocatable :: iw, ndegr

 nnzDh = 2*(n2-2)*(n1-2)
 nnzDv = 2*(n2-2)*(n1-2)

 allocate(Dh(nnzDh),iDh(nnzDh),jDh(nnzDh)) 
 allocate(Dv(nnzDv),iDv(nnzDv),jDv(nnzDv))

 allocate(Dht(nnzDh),iDht(nnzDh),jDht(nnzDh))   
 allocate(Dvt(nnzDv),iDvt(nnzDv),jDvt(nnzDv))

 call subderivop(n1, n2, Dv, iDv, jDv, nnzDv, Dh, iDh, jDh, nnzDh, val) 
 !call subderivopt(n1, n2, Dvt, iDvt, jDvt, nnzDv, Dht, iDht, jDht, nnzDh, val)

 allocate(Dh_csr(nnzDh),iDh_csr(nn+1),jDh_csr(nnzDh))
 allocate(Dv_csr(nnzDv),iDv_csr(nn+1),jDv_csr(nnzDv))

 allocate(Dht_csr(nnzDh),iDht_csr(nn+1),jDht_csr(nnzDh))
 allocate(Dvt_csr(nnzDv),iDvt_csr(nn+1),jDvt_csr(nnzDv))
 
 call coocsr2(nn, nnzDv, Dv, iDv, jDv,  Dv_csr, jDv_csr, iDv_csr)
 call coocsr2(nn, nnzDv, Dv, jDv, iDv, Dvt_csr,jDvt_csr,iDvt_csr)
 call coocsr2(nn, nnzDh, Dh, iDh, jDh,  Dh_csr, jDh_csr, iDh_csr)
 call coocsr2(nn, nnzDh, Dh ,jDh,iDh, Dht_csr,jDht_csr,iDht_csr)

 !Product Dvt*Dv in CSR format
 allocate(ndegr(nn))
 allocate(iw(nn))
 call amubdg_test(nn,nn,nn,jDvt_csr,iDvt_csr,nnzDv,jDv_csr,iDv_csr,nnzDv,ndegr,nnzDDv,iw)
 deallocate(ndegr)
 deallocate(iw)
 write(*,*) " nnzDv, nnzDDv ==",nnzDv, nnzDDv

 allocate(DDv_csr(nnzDDv),jDDv_csr(nnzDDv),iDDv_csr(nn+1))
 allocate(iw(nn))
 job = 1
 call amubtest (nn,nn,nn,job,Dvt_csr,jDvt_csr,iDvt_csr,nnzDv,Dv_csr,jDv_csr,&
      iDv_csr,nnzDv,DDv_csr,jDDv_csr,iDDv_csr,nnzDDv,iw,ierr)
 if (ierr .ne. 0) print *,'Problem Dv^t*Dv, ierr = ',ierr
 write(*,*) ' Product DDv = Dv^t * Dv, OK !!!'
 deallocate(iw)

 !Product Dht*Dh in CSR format
 allocate(ndegr(nn))
 allocate(iw(nn))
 call amubdg_test(nn,nn,nn,jDht_csr,iDht_csr,nnzDh,jDh_csr,iDh_csr,nnzDh,ndegr,nnzDDh,iw)
 deallocate(ndegr)
 deallocate(iw)
 write(*,*) " nnzDh, nnzDDh ==",nnzDh, nnzDDh

 allocate(DDh_csr(nnzDDh),jDDh_csr(nnzDDh),iDDh_csr(nn+1))
 allocate(iw(nn))
 job = 1
 call amubtest (nn,nn,nn,job,Dht_csr,jDht_csr,iDht_csr,nnzDh,Dh_csr,jDh_csr,&
       iDh_csr,nnzDh,DDh_csr,jDDh_csr,iDDh_csr,nnzDDh,iw,ierr)
 if (ierr .ne. 0) print *,'Probleme Dh^t*Dh,  ierr = ',ierr
 write(*,*) ' Product DDh = Dh^t * Dh, OK !!!'
 deallocate(iw)

 ! nnzR1 with R1 =  DDv + DDh
 allocate(ndegr(nn))
 allocate(iw(nn))
 call aplbdg(nn,nn,jDDv_csr,iDDv_csr,nnzDDv,jDDh_csr,iDDh_csr,nnzDDh,ndegr,nnzR1,iw)
 deallocate(ndegr)
 deallocate(iw)
 write(*,*) 'nnz(DDv + DDh) = ', nnzR1

 allocate(R1(nnzR1),jR1(nnzR1),iR1(nn+1))
 allocate(iw(nn))
 job = 1
 call aplb (nn,nn,job,DDv_csr,jDDv_csr,iDDv_csr,nnzDDv,DDh_csr,jDDh_csr,iDDh_csr,nnzDDh, &
       R1,jR1,iR1,nnzR1,iw,ierr)
 deallocate(iw)
 write(*,*) ' DDv + DDh , OK !!!'

 !convert preco_csr to COO format ( COO --> CSR)
 job = 3
 allocate(reg_coo(nnzR1),i_reg_coo(nnzR1),j_reg_coo(nnzR1))
 call csrcoo (nn,job,nnzR1,R1,jR1,iR1,nnzR12,reg_coo,i_reg_coo,j_reg_coo,ierr)
 if (ierr .ne. 0) print *,'ERROR while converting preco_csr, CSR ==> COO '
 if (nnzR12 .ne. nnzR1) print *,'ERROR while converting preco_csr, nnzp =/=nnzPreco '
 write(*,*) 'Convert reg_csr to format COO, OK !!!'
  
 open(NEWUNIT=unit,file='mat_reg.bin',access='STREAM',form='UNFORMATTED',action='WRITE')
 write(unit) (reg_coo(i),i=1,nnzR1)
 close(unit)
 open(NEWUNIT=unit,file='IRN_reg.bin',access='STREAM',form='UNFORMATTED',action='WRITE')
 write(unit) (i_reg_coo(i),i=1,nnzR1)
 close(unit)
 open(NEWUNIT=unit,file='JCN_reg.bin',access='STREAM',form='UNFORMATTED',action='WRITE')
 write(unit) (j_reg_coo(i),i=1,nnzR1)
 close(unit)

 deallocate(Dh_csr,iDh_csr,jDh_csr)
 deallocate(Dv_csr,iDv_csr,jDv_csr)

 deallocate(Dht_csr,iDht_csr,jDht_csr)
 deallocate(Dvt_csr,iDvt_csr,jDvt_csr)
 
 deallocate(DDh_csr,jDDh_csr,iDDh_csr)
 deallocate(R1,jR1,iR1)

 deallocate(reg_coo,i_reg_coo,j_reg_coo)

 deallocate(Dh,iDh,jDh)   
 deallocate(Dv,iDv,jDv)
 deallocate(Dht,iDht,jDht)
 deallocate(Dvt,iDvt,jDvt)

end subroutine
