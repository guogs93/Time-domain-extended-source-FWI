! ==========================================================================
! LOOP PARTIONING FOR MPI (FROM IBM MPI DOCUMENTTATION)
! ==========================================================================
! n1: the lowest value of the iteration variable (IN)
! n2: the highest value of the iteration variable (IN)
! nprocs: the number of processes (IN)
! irank: the rank for which you want to know the range of iterations (IN) 
! ista: the lowest value of the iteration variable that processes irank executes (OUT)
! iend: the highest value of the iteration variable the processes irank executes (OUT)

 SUBROUTINE para_range(n1, n2, nprocs, irank, ista, iend)
 iwork1 = (n2 -n1 +1) / nprocs
 iwork2 =  MOD(n2 - n1 + 1, nprocs)
 ista = irank*iwork1 + n1 + MIN(irank, iwork2)
 iend = ista + iwork1 - 1
 IF (iwork2 .GT. irank) iend = iend + 1
 END SUBROUTINE para_range
