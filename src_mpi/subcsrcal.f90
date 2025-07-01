 subroutine  coocsr2 (ni,nn,a,ir,jc,ao,jao,iao)
  implicit none
  integer , intent(in) :: ni , nn
  integer              :: i , j , k0 , k , iad
  real , dimension(nn) , intent(in) :: a
  real , dimension(nn), intent(out) :: ao
  integer , dimension(nn) , intent(in) :: ir,jc
  integer , dimension(nn) , intent(out) :: jao
  integer , dimension(ni+1) , intent(out) ::iao
  real :: x

  ao(:) = 0
  jao(:) = 0
  iao(:) = 0

  do 1 k=1,ni+1
       iao(k) = 0
1      continue
       ! determine row-lengths.
         do 2 k=1, nn
              iao(ir(k)) = iao(ir(k))+1
2             continue
              ! starting position of each row..
              k = 1
              do 3 j=1, ni+1
                   k0 = iao(j)
                   iao(j) = k
                   k = k+k0
3                  continue
              ! go through the structure  once more. Fill in output matrix.
                   do 4 k=1, nn
                        i = ir(k)
                        j = jc(k)
                        x = a(k)
                        iad = iao(i)
                        ao(iad) =  x
                        jao(iad) = j
                        iao(i) = iad+1
4                       continue
                        ! shift back iao
                        do 5 j=ni,1,-1
                             iao(j+1) = iao(j)
5                            continue
                             iao(1) = 1
                             return

 end subroutine coocsr2

subroutine amubdg_test (ni,nj,njb,ja,ia,nnza,jb,ib,nnzb,ndegr,nn,iw)
  ! Détermination du nombre de coefficients non nuls de la matrice C = A*B
  ! renvoie nnz_C et ndegr (le nombre de coefficients non nuls par ligne)
  !subroutine amubdg (nrow,ncol,ncolb,ja,ia,jb,ib,ndegr,nnz,iw)

  integer, intent(in)  :: ni, nj, njb,nnza,nnzb
  integer, intent(out) :: nn
  integer, intent(in)  :: ja(nnza), ia(ni+1), jb(nnzb), ib(nj+1)
  integer, intent(out) :: ndegr(ni) 
  integer, intent(out) :: iw(njb)
  integer              :: j, ii, jr, jc, k, ldg, last

  do 1 k=1, njb
     iw(k) = 0
1    continue

     do 2 k=1, ni
        ndegr(k) = 0
2       continue

        !     method used: Transp(A) * A = sum [over i=1, nrow]  a(i)^T a(i)
        !     where a(i) = i-th row of  A. We must be careful not to add  the
        !     elements already accounted for.

        do 7 ii=1,ni

           !     for each row of A
           ldg = 0

           !    end-of-linked list
           last = -1
           do 6 j = ia(ii),ia(ii+1)-1

              !     row number to be added:
              jr = ja(j)
              do 5 k=ib(jr),ib(jr+1)-1
                 jc = jb(k)
                 if (iw(jc) .eq. 0) then

                    !     add one element to the linked list
                    ldg = ldg + 1
                    iw(jc) = last
                    last = jc
                 endif
5                continue
6                continue
                 ndegr(ii) = ldg

                 !     reset iw to zero
                 do 61 k=1,ldg
                    j = iw(last)
                    iw(last) = 0
                    last = j
61                  continue

7                   continue
                    nn = 0
                    do 8 ii=1, ni
                       nn = nn+ndegr(ii)
8                      continue
                       return

 end  subroutine amubdg_test

 subroutine amubtest(ni,nib,nj,job,a,ja,ia,nnza,b,jb,ib,nnzb,c,jc,ic,nn,iw,ierr)
   ! Effectue le produit des matrices creuses A et B en
   ! format CSR
   ! renvoie C = A*B en format CSR, on a C, jC et iC

   integer, intent(in) :: nj, ni, nib, nnza, nnzb, nn ,job
   integer, intent(out) :: ierr
   integer              :: l, j,ii,ka, jj, kb, jcol, jpos,k
   real, intent(in) :: a(nnza) , b(nnzb)
   real, intent(out) :: c(nn)
   integer, intent(in) :: ja(nnza),jb(nnzb),ia(ni+1),ib(nib+1)
   integer, intent(out) :: jc(nn), ic(ni+1) , iw(nj)
   real :: scal
   logical :: values

   c(:) = 0
   jc(:) = 0
   ic(:) = 0

   values = (job .ne. 0)
   l = 0
   ic(1) = 1
   ierr = 0
   !     initialize array iw.
   do 1 j=1, nj
        iw(j) = 0
1       continue

        do 500 ii=1, ni
               !     row i
               do 200 ka=ia(ii), ia(ii+1)-1
                      if (values) scal = a(ka)
                         jj   = ja(ka)
                         do 100 kb=ib(jj),ib(jj+1)-1
                                jcol = jb(kb)
                                jpos = iw(jcol)
                                if (jpos .eq. 0) then
                                   l = l+1
                                   if (l .gt. nn) then
                                      ierr = ii
                                      return
                                   endif
                                   jc(l) = jcol
                                   iw(jcol)= l
                                   if (values) c(l)  = scal*b(kb)
                                else
                                   if (values) c(jpos) = c(jpos) + scal*b(kb)
                                endif
100                             continue
200                             continue
                                do 201 k=ic(ii), l
                                      iw(jc(k)) = 0
201                                   continue
                                      ic(ii+1) = l+1
500                                   continue
                                      return

 end  subroutine amubtest

subroutine aplb (ni,nj,job,a,ja,ia,nnza,b,jb,ib,nnzb,c,jc,ic,m,iw,ierr)
   ! Effectue la somme de deux matrices creuses C = A + B
! renvoie C en format CSR

   integer, intent(in) :: nj, ni,nnza, nnzb
   integer, intent(in) :: m
   integer , intent(in)  :: job
   integer , intent(out) :: ierr
   integer              :: ll, j, k , ii, ka, kb, jcol, jpos
   real , intent(in) :: a(nnza), b(nnzb)
   real , intent(out) :: c(m)
   integer , intent(in) :: ja(nnza),jb(nnzb),ia(ni+1),ib(ni+1)
   integer , intent(out) :: jc(m), ic(ni+1) , iw(nj)
   logical :: values


   values = (job .ne. 0) 
   ierr = 0
   ll = 0
   ic(1) = 1 
   do 1 j=1, nj
        iw(j) = 0
1       continue

        do 500 ii=1, ni
           !     row i 
           do 200 ka=ia(ii), ia(ii+1)-1 
                  ll = ll+1
                  jcol    = ja(ka)
                  if (ll .gt. m) goto 999
                     jc(ll) = jcol 
                     if (values) c(ll)  = a(ka) 
                        iw(jcol)= ll
200                     continue     
                        do 300 kb=ib(ii),ib(ii+1)-1
                               jcol = jb(kb)
                               jpos = iw(jcol)
                               if (jpos .eq. 0) then
                                  ll = ll+1
                                  if (ll .gt. m) goto 999
                                     jc(ll) = jcol
                                     if (values) c(ll)  = b(kb)
                                        iw(jcol)= ll
                                     else
                                        if (values) c(jpos) = c(jpos) + b(kb)
                                     endif
300                                  continue
                                     do 301 k=ic(ii), ll
                                            iw(jc(k)) = 0
301                                         continue
                                            ic(ii+1) = ll+1
500                                         continue
                                            return
999                                         ierr = ii
                                            return

    end subroutine aplb

subroutine aplbdg (ni,nj,ja,ia,nnza,jb,ib,nnzb,ndegr,m,iw) 
   ! Détermine le nombre de coefficients non nuls de C = A + B
   ! renvoie m = nnz_C et ndegr (le nombre de coefficients non nuls par ligne)

   integer, intent(in) :: nj, ni, nnza, nnzb
   integer, intent(out) :: m
   integer              :: k , ii, ldg , last , j, jr , jc
   integer , intent(in) :: ja(nnza),jb(nnzb),ia(ni+1),ib(ni+1)
   integer , intent(out) :: ndegr(ni) , iw(nj)


   do 1 k=1, nj 
        iw(k) = 0 
1       continue

        do 2 k=1, ni
             ndegr(k) = 0 
2            continue

             do 7 ii=1,ni
                  ldg = 0 

                  !    end-of-linked list

                last = -1 

                !     row of A

            do 5 j = ia(ii),ia(ii+1)-1 
                 jr = ja(j) 
                 !     add element to the linked list 

                 ldg = ldg + 1
                 iw(jr) = last 
                 last = jr
5                continue

                 !     row of B

                 do 6 j=ib(ii),ib(ii+1)-1
                      jc = jb(j)
                      if (iw(jc) .eq. 0) then 
                         !     add one element to the linked list 

                          ldg = ldg + 1
                          iw(jc) = last 
                          last = jc
                      endif
6                     continue
                      !     done with row ii. 
                      ndegr(ii) = ldg
                      !     reset iw to zero

                      do 61 k=1,ldg 
                            j = iw(last) 
                            iw(last) = 0
                            last = j
61                          continue

7                           continue

                            m = 0
                            do 8 ii=1, ni 
                                 m = m+ndegr(ii) 
8                                continue
                                 return

     end subroutine aplbdg


SUBROUTINE csrcoo (ni,job,nzmax,a,ja,ia,nnz_a,ao,ir,jc,ierr)
   ! Convertit une matrice stockée en format CSR en format COO (CSR --> COO)
   ! Renvoie A, jc (indices de colonne ), ir ( indices de ligne)

   integer, intent(in) :: ni , nzmax
   integer , intent(in)  :: job
   integer , intent(out) :: ierr
   integer              :: i ,k , k1, k2 
   integer, intent(out)  :: nnz_a
   real , intent(in) :: a(nzmax) 
   real , intent(out) :: ao(nzmax)
   integer , intent(in) :: ja(nzmax), ia(ni+1) 
   integer , intent(out) :: ir(nzmax), jc(nzmax)

     ierr = 0
     nnz_a = ia(ni+1)-1
     if (nnz_a .gt. nzmax) then
        ierr = 1
        return
     endif

     goto (3,2,1) job
1         do 10 k=1,nnz_a
                ao(k) = a(k)
10              continue
2         do 11 k=1,nnz_a
                jc(k) = ja(k)
11              continue

                                            !     copy backward to allow for
                                            !     in-place processing. 

3               do 13 i=ni,1,-1
                      k1 = ia(i+1)-1
                      k2 = ia(i)
                      do 12 k=k1,k2,-1
                            ir(k) = i
12                          continue
13                          continue
                      return

     end SUBROUTINE csrcoo
