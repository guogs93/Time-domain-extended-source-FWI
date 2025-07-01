SUBROUTINE subderivop(n1, n2, D1, iD1, jD1, nnz1, D2, iD2, jD2, nnz2, val)
 
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: n1, n2, val, nnz1, nnz2
    INTEGER             :: i1, i2, k, l, nnz01, nnz02

    INTEGER,ALLOCATABLE :: iD1(nnz1),jD1(nnz1),iD2(nnz2),jD2(nnz2)
    integer,ALLOCATABLE :: D1(nnz1),D2(nnz2)
    
    !#############################################################################################################
    
    !val=1 ==> First order derivatives Dz, Dx
    !First order forward difference
    IF (val .EQ. 1 ) THEN
        
        nnz01=0
        DO i2=1:n2
            DO i1=1:n1-1
            
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
        IF ( nnz01 .NE. nnz1 ) THEN WRITE(*,*) " ERROR nnz01 =/= nnz1 "
        
        nnz02=0
        DO i2=1:n2-1
            DO i1=1:n1
            
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
        IF ( nnz02 .NE. nnz2 ) THEN WRITE(*,*) " ERROR nnz02 =/= nnz2 "
        
    !val=2 ==> Second order derivatives DDz, DDx
    !Second order forward difference
    ELSE
    
        nnz01=0
        DO i2=1:n2
            DO i1=1:n1-2
            
            nnz01=nnz01+1
            k=(i2-1)*n1 + i1
            
            D1(nnz01)=1 !cmplx(+1.0,0.0)
            iD1(nnz01)=k
            jD1(nnz01)=k
            
            nnz01=nnz01+1
            l=(i2-1)*n1 + (i1+1)
            D1(nnz01)=-2 !cmplx(-2.0,0.0)
            iD1(nnz01)=k
            jD1(nnz01)=l
            
            nnz01=nnz01+1
            l=(i2-1)*n1 + (i1+2)
            D1(nnz01)=1 !cmplx(+1.0,0.0)
            iD1(nnz01)=k
            jD1(nnz01)=l
            
            END DO
            
        END DO
        WRITE(*,*) "nnz01 and nnz1 ==",nnz01, nnz1
        IF ( nnz01 .NE. nnz1 ) THEN WRITE(*,*) " ERROR nnz01 =/= nnz1 "
        
        nnz02=0
        DO i2=1:n2-1
            DO i1=1:n1
            
            nnz02=nnz02+1
            k=(i2-1)*n1 + i1
            
            D2(nnz02)=1 !cmplx(+1.0,0.0)
            iD2(nnz02)=k
            jD2(nnz02)=k
            
            nnz02=nnz02+1
            l=((i2+1)-1)*n1 + i1
            D2(nnz02)=-2 !cmplx(-2.0,0.0)
            iD2(nnz02)=k
            jD2(nnz02)=l
            
            nnz02=nnz02+1
            l=((i2+2)-1)*n1 + i1
            D2(nnz02)=1 !cmplx(+1.0,0.0)
            iD2(nnz02)=k
            jD2(nnz02)=l
            
            END DO
        END DO
        WRITE(*,*) "nnz02 and nnz2 ==",nnz02, nnz2
        IF ( nnz02 .NE. nnz2 ) THEN WRITE(*,*) " ERROR nnz02 =/= nnz2 "
    
    END IF

END SUBROUTINE subderivop


SUBROUTINE build2ndsubpb_TVreg_LHS(grid,inv)
 
    IMPLICIT NONE
 
    INCLUDE 'mpif.h'
    INCLUDE 'common.h'
    INCLUDE 'ffwi2d.h'
    
    TYPE (grid_type)                :: grid
    TYPE (inv_type)                 :: inv
    
    INTEGER             :: unit=10, ierr
    CHARACTER(len=160)   :: nameout
    
    INTEGER             :: i, j, i1, i2, k, l, n1, n2, nn, val=1, job
    INTEGER,ALLOCATABLE :: iw(:) , ndegr(:)
    
    INTEGER,ALLOCATABLE :: iDh(:),jDh(:),iDv(:),jDv(:),iDh_csr(:),jDh_csr(:),iDv_csr(:),jDv_csr(:)
    INTEGER,ALLOCATABLE :: iDht_csr(:),jDht_csr(:),iDvt_csr(:),jDvt_csr(:),iId_gamm0(:),jId_gamm0(:)
    INTEGER,ALLOCATABLE :: iDDh_csr(:),jDDh_csr(:),iDDv_csr(:),jDDv_csr(:),iId_gamm0_csr(:), jId_gamm0_csr(:)
    integer,ALLOCATABLE :: Dh(:),Dv(:), Dh_csr(:), Dv_csr(:), Dht_csr(:),Dvt_csr(:), DDh_csr(:), DDv_csr(:) , Id_gamm0(:), Id_gamm0_csr(:)
    integer             :: R2(:), R1(:), R1_tmp(:)
    INTEGER,ALLOCATABLE :: iR1(:), jR1(:), iR1_tmp(:), jR1_tmp(:)
    INTEGER             :: nnzDv, nnzDh, nnzDDh, nnzDDv, nnzR1_tmp, nnzR1_tmp,nnzPreco, nnzp
    integer,ALLOCATABLE :: gradp_csr(:), preco_csr(:), preco_reg_coo(:)
    INTEGER,ALLOCATABLE :: igradp(:), jgradp(:), igradp_csr(:), jgradp_csr(:),ipreco_csr(:), jpreco_csr(:)
    INTEGER,ALLOCATABLE :: ipreco_reg_coo(:), jpreco_reg_coo(:)

    nnzDh = 2(grid%n2 - 1)*grid%n1
    nnzDv = 2(grid%n2)*(grid%n1 - 1)

    ALLOCATE(Dh(nnzDh))
    ALLOCATE(iDh(nnzDh))  
    ALLOCATE(jDh(nnzDh)) 
    ALLOCATE(Dv(nnzDv))
    ALLOCATE(iDv(nnzDv))  
    ALLOCATE(jDv(nnzDv))

    var=1

    Dh(:)  = 0
    iDh(:) = 0
    jDh(:) = 0
    Dv(:)  = 0
    iDv(:) = 0
    jDv(:) = 0

    !Compute the first order derivative matrics Dv et Dh
    CALL subderivop(grid%n1, grid%n2, Dv, iDv, jDv, nnzDv, Dh, iDh, jDh, nnzDh, val)

    !convert Dv into CSR format ( COO ==> CSR)
    CAll coocsr2 (grid%nn,nnzDv,Dv,iDv,jDv, Dv_csr,jDv_csr,iDv_csr)

    !convert Dvt into CSR format ( COO ==> CSR)
    CAll coocsr2 (grid%nn, nnzDv, Dv,jDv,iDv, Dvt_csr,jDvt_csr,iDvt_csr)
    
    !convert Dh into CSR format ( COO ==> CSR)
    CAll coocsr2 (grid%nn, nnzDh, Dh,iDh,jDh, Dh_csr,jDh_csr,iDh_csr)
    
    !convert Dht into CSR format ( COO ==> CSR)
    CAll coocsr2 (grid%nn, nnzDh, Dh,jDh,iDh, Dht_csr,jDht_csr,iDht_csr)

    

end
