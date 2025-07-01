!****************************************
! global parameter
!****************************************

module parameters
! ==============================================================================
! PARAMETERS RELATED TO THE FORWARD PROBLEM (grid, acquisition)
! ==============================================================================

    TYPE grid_type
    SEQUENCE
    INTEGER                     :: mode
    INTEGER                     :: n1,n2,npml,nt,nn,nto
    REAL                        :: h,dt,t0,dto
    REAL                        :: x1min,x1max,x2min,x2max
    REAL,ALLOCATABLE            :: ctrue(:,:),c(:,:),c0(:,:),s(:)
    REAL,ALLOCATABLE            :: slow(:,:),slow0(:,:),slowt(:,:)
    INTEGER                     :: type_src,type_rec
    INTEGER                     :: ofs
    INTEGER                     :: oio
    REAL                        :: flow, fhig, llow, lhig
    INTEGER                     :: npatch
    INTEGER                     :: ntrace,nsrc_tot,nrec_tot,nrecmax
    INTEGER, ALLOCATABLE        :: nsrc(:),nrec(:)
    REAL, ALLOCATABLE           :: xs(:,:),xr(:,:),topo(:)
    INTEGER, ALLOCATABLE        :: is(:,:),ir(:,:),itopo(:)

    END TYPE grid_type

! ==============================================================================
! PARAMETERS RELATED TO THE INVERSE PROBLEM
! ==============================================================================

    TYPE inv_type

    SEQUENCE

    REAL, ALLOCATABLE           :: dataobs(:,:),datar(:,:),deltadr(:,:),deltade(:,:),dataesc(:,:)
    REAL, ALLOCATABLE           :: datae(:,:)
    INTEGER                     :: nitermax,optimalgo,nop,preco
    REAL                        :: perc,cost,costsum,waterlevel,alpha,mu,cost2,costsum2,tau1,tau2,freq,maxvel,minvel
    REAL, ALLOCATABLE           :: grad(:),gradsum(:),gradpreco(:),gradprecosum(:),m(:),precon(:)

    END TYPE inv_type

end module parameters

