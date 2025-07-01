! ----------------------------------------------------------------------------------------------
! nz: nbre de points d'entree de la fonction
! npi: nbre de coefficients de b-spline
! deltat: position of the spline knots
! x: values of the function
! nz: number of samples of x and deltat
! bcoef, tcoef, npi: output vectors and parameter

 subroutine subresample(x,x1,n1,n2,n1out,d1in,d1out)

 IMPLICIT NONE

 integer :: ntmax0
 parameter (ntmax0=500000)

 INTEGER :: i1,i2,i3,i,j,k,n1,n2,n1out,n2out,npi,mflag

 REAL :: d1in,d1out

 REAL :: x(n1,n2),x1(n1out,n2)
 !REAL :: t(ntmax0),xval(ntmax0),deltat(ntmax0),bcoeft(ntmax0),tcoeft(ntmax0),bcoefz(ntmax0),tcoefz(ntmax0)
 real :: t(n1out),xval(n1out),deltat(n1),bcoefz(n1+2),tcoefz(n1+6)

 k=4

! ----------------------------------------------------------------------
! RESAMPLING ALONG DIRECTION 1

! set vector t

 do j=1,n1out
 t(j)=float(j-1)*d1out
 end do

 do j=1,n1
 deltat(j)=float(j-1)*d1in
 end do

 do i2=1,n2

 call binit(deltat,x(:,i2),bcoefz,tcoefz,n1,npi)

! Interpolate the value at each time sample

 do j=1,n1out
 call bsplval(tcoefz,bcoefz,npi,k,t(j),xval(j),mflag)
 end do

 do j=1,n1out
 x1(j,i2)=xval(j)
 end do

 enddo

 return
 end
! ======================================================================
          subroutine binit(xpas,prof,bcoef,t,np,npi)
! ************************************************************

! subroutine pour calculer les matrices t et bcoef contenant les noeuds
! et les coefficients des bsplines

! variables d entree
! xpas , prof,np
! variables de sortie: bcoef,t,npi

! dimensionnement
! la matrice t doit avoir une dimension > np+6
! la matrice bcoef doit avoir une dimension > np+2
! le nombre de points np doit etre superieur a 4

! bibliotheque des variables
! np  : nbre de points d entree de la fonction
! xpas: matrice de dimension > np contenant les points d echantillonage
! prof: matrice de dimension > np contenant les valeurs de la fonction
!        aux abscisses xpas(i)
! npi : nbre de coefficients de bsplines
! t   : matrice contenant les noeuds des bsplines
! bcoef : matrice contenant les poids des bsplines

! ************************************************************

      dimension prof(*),bcoef(*),t(*),xpas(*)

      fonc(x)=(x-x2)*(x-x3)*(x-x4)*f1/((x1-x2)*(x1-x3)*(x1-x4))               &
      +(x-x1)*(x-x3)*(x-x4)*f2/((x2-x1)*(x2-x3)*(x2-x4))                      &    
      +(x-x1)*(x-x2)*(x-x4)*f3/((x3-x1)*(x3-x2)*(x3-x4))                      &
      +(x-x1)*(x-x2)*(x-x3)*f4/((x4-x1)*(x4-x2)*(x4-x3))                    

      bcoef(1)=prof(1)
      xminii=xpas(1)
      t(1)=xminii
      t(2)=xminii
      t(3)=xminii
      t(4)=xminii
      ncond=4
      ncond2=2
      npm1=np-1
      npm2=np-2
      npp1=np+1
      npp2=np+2
      do 210 j=2,npm1
      ncond=ncond+1
      t(ncond)=xpas(j)
210   continue
      pt=xpas(np)
      ncond=ncond+1
      t(ncond)=pt

      y1=xminii
      y2=xminii
      y3=t(5)
      x1=xpas(1)
      x2=xpas(2)
      x3=xpas(3)
      x4=xpas(4)
      f1=prof(1)
      f2=prof(2)
      f3=prof(3)
      f4=prof(4)
      xi=(y1+y2+y3)/3.
      bcoef(2)=fonc(xi)
      do 217 j=2,npm2
      ncond2=ncond2+1
      x1=t(ncond2+1)
      x2=t(ncond2+2)
      x3=t(ncond2+3)
      x4=t(ncond2+4)
      f1=prof(ncond2-2)
      f2=prof(ncond2-1)
      f3=prof(ncond2)
      f4=prof(ncond2+1)
      xi=(x1+x2+x3)/3.
      bcoef(ncond2)=fonc(xi)
217   continue
      ncond2=ncond2+1
      xi=(x4+x2+x3)/3.
      bcoef(ncond2)=fonc(xi)
      do 216 j=1,3
      ncond=ncond+1
      t(ncond)=pt
216   continue
      ncond2=ncond2+1
      y1=t(ncond2+1)
      y2=t(ncond2+2)
      y3=t(ncond2+3)
      xi=(y1+y2+y3)/3.
      bcoef(ncond2)=fonc(xi)
      ncond2=ncond2+1
      bcoef(ncond2)=prof(np)
      npi=ncond2
      return
      end

! ***************************************************************

      subroutine bsplval(t,bcoef,n,k,x,bvalue,mflag)
      common/bspl/ajc,dl,dr,xi,km1

! cf de boor :a practical guide to splines ,page 144
! calculates value at x of the  spline from b-repr
! the spline is taken to be continuous from the right
! calcule de plus les coefficients necessaires pour le prog bspldps.f77
! (ajc,dl,dr,km1)
!  ,programme qui calcule les derivees premieres et secondes
! necessite de mettre le common/bspl/ dans le programme
! appelant les deux subroutines  si on desire calculer
! les derivees premiere et seconde ensuite

! ***********input *************
! t,bcoef,n,k...forms the b-representation of the spline f to be
! evaluated.specifically
! t....knot sequence,of length n+k,assumed non decreasing.
! bcoef...b_coeeficient sequence , of length n.
! n....length of bcoef and dimension of spline(k,t),
!           assumed positive.
! k.....order of the spline .
!
! x....the point at which to evaluate .
! jderiv....integer giving the order of the derivative to be evaluated
!       assumed to be zero
!
! ***********output************
! bvalue.....the value of f at x .
!
! *****************************************************************
      parameter (kmax=5)
      real bcoef(*),t(*),x,aj(kmax),dl(kmax),dr(kmax),ajc(kmax)
      bvalue=0.
      xi=x
      if(0.ge.k)go to 99

! find i / 1 .le. i .lt. n+k and t(i).lt.t(i+1) and t(i).le.x.lt.t(i+1)
! if no such i can be found ,x lies outside the support of the spline f
! and bvalue =0.
! (the assymetry in this choice of i makes f right continuous)
      call interv(t,n+k,x,i,mflag)
      if(mflag.ne.0)go to 99
! ***if k=1 bvalue=bcoef(i)
      km1=k-1
      if(km1.gt.0)go to 1
      bvalue=bcoef(i)
      go to 99

! ***store the k b_spline coeeficients relevant for the knot interval
!   (t(i),t(i+1))in aj(1),....,aj(k) and compute dl(j)=x-t(i+1-j),
!   dr(j)=t(i+j)-x ,j=1,...,k-1 . set any of the aj not obtainable
!    from input to zero. set any t.s not obtainable equal to t(1) or
!   to t(n+k) appropriately.
1     jcmin=1
      imk=i-k
      if(imk.ge.0)go to 8
      jcmin=1-imk
      do 5 j=1,i
5     dl(j)=x-t(i+1-j)
      do 6 j=i,km1
      aj(k-j)=0.
      ajc(k-j)=0.
6     dl(j)=dl(i)
      go to 10
8     do 9 j=1,km1
9     dl(j)=x-t(i+1-j)

10    jcmax=k
      nmi=n-i
      if(nmi.ge.0)go to 18
      jcmax=k+nmi
      do 15 j=1,jcmax
15    dr(j)=t(i+j)-x
      do 16 j=jcmax,km1
      aj(j+1)=0.
      ajc(j+1)=0.
16    dr(j)=dr(jcmax)
      go to 20
18    do 19 j=1,km1
19    dr(j)=t(i+j)-x

20    do 21 jc=jcmin,jcmax
      aj(jc)=bcoef(imk+jc)
21    ajc(jc)=aj(jc)

! ***compute value at x in (t(i),t(i+1))of the function
!   given its relevant b-spline coeffs in aj(1),...,aj(k).
30    if(0.eq.km1)go to 39
      do 33 j=1,km1
      kmj=k-j
      ilo=kmj
      do 33 jj=1,kmj
      aj(jj)=(aj(jj+1)*dl(ilo)+aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
33    ilo=ilo-1
39    bvalue=aj(1)

99    return
      end
! ***************************************************************
        subroutine filter(fil,ntnew1,fmin,fmax,f01,f02,dt)

!   if01 <---> ifmin <-------------> ifmax <---> if02

        implicit none
	integer nn,ntnew1
        real fil(ntnew1)
        real fs,dt,fnyq,df
        real a1,b1,a2,b2
        real f
        real fmin,fmax,f01,f02
        integer i

        nn=ntnew1-1
        fs=1./dt
        fnyq=0.5*fs
        df=fnyq/float(nn-1)

        if (fmin.eq.f01) then
        a1=-9999.
        else
        a1=1./(fmin-f01)
        b1=-a1*f01
        end if

        if (fmax.eq.f02) then
        a2=9999.
        else
        a2=1./(fmax-f02)
        b2=-a2*f02
        end if

        do i=1,ntnew1
        f=float(i-1)*df

        if (f.lt.f01) then
        fil(i)=0.
        elseif (f.ge.f01.and.f.le.fmin) then
                if (a1.eq.-9999.) then
                fil(i)=1.
                else
                fil(i)=a1*f+b1
                end if
        elseif (f.gt.fmin.and.f.lt.fmax) then
        fil(i)=1.
        elseif (f.ge.fmax.and.f.le.f02) then
                if (a2.eq.9999.) then
                fil(i)=1.
                else
                fil(i)=a2*f+b2
                end if
        elseif (f.gt.f02) then
        fil(i)=0.
        end if

        end do

        return
        end
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine interv(xt,lxt,x,left,mflag)
! computes   left=max(i,, 1.le.i.le.lxt .and. xt(i).le.x )

! ************input ****
! xt...assumed non decreasing
! lxt...number of terms in the sequence xt.
! x...the point whose location with respect to the sequence xt is
!     to be determined.

! ************output*******
! left,mflag.....both integers,whose value is

!  1     -1        if       x.lt.xt(1)
!  i      0        if     xt(i).le.x.lt.xt(i+1)
! lxt    1        if     xt(lxt).le.x

      real xt(lxt)
      data ilo/1/
!     save ilo
      ihi=ilo+1
      if(ihi.lt.lxt)go to 20
      if(x.ge.xt(lxt))go to 110
      if(lxt.le.1)go to 90
      ilo=lxt-1
      ihi=lxt

20    if(x.ge.xt(ihi))go to 40
      if(x.ge.xt(ilo))go to 100

!     **** now x.lt.xt(ilo).decrease ilo to capture x .
      istep=1
31    ihi=ilo
      ilo=ihi-istep
      if(ilo.le.1)go to 35
      if(x.ge.xt(ilo))go to 50
      istep=istep*2
      go to 31
35    ilo=1
      if(x.lt.xt(1))go to 90
      go to 50
!     ***8now x.ge.xt(ihi). increase ihi to capture x .
40    istep=1
41    ilo=ihi
      ihi=ilo+istep
      if(ihi.ge.lxt)go to 45
      if(x.lt.xt(ihi))go to 50
      istep=istep*2
      go to 41
45    if(x.ge.xt(lxt))go to 110
      ihi=lxt

!     ****now xt(ilo).le.x.lt.xt(ihi) . narrow the interval
50    middle=(ilo+ihi)/2
      if(middle.eq.ilo)go to 100
!   note. it is assumed that middle=ilo in case ihi=ilo+1
      if(x.lt.xt(middle))go to 53
      ilo=middle
      go to 50
53    ihi=middle
      go to 50
!       *****set output and return
90    mflag=-1
      left=1
      return
100   mflag=0
      left=ilo
      return
110   mflag=1
      left=lxt
      return
      end

