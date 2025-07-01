
!####################################################################
! FOURIER TRANSFORM SUBROUTINES (SEE NUMERICAL RECIPES)

! subroutine fourn.f

  subroutine fourn(data,nn,ndim,isign)
  integer isign,ndim
  integer nn(ndim)
  real*4 data(*)
  real*4 tempi,tempr
  real*8 theta,wi,wpi,wpr,wr,wtemp
  integer i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot

 ntot=1
 do idim=1,ndim
 ntot=ntot*nn(idim)
 end do

 nprev=1
! boucle sur les dimensions
 do idim=1,ndim
 n=nn(idim)
 nrem=ntot/(n*nprev)
 ip1=2*nprev
 ip2=ip1*n
 ip3=ip2*nrem
 i2rev=1
     do i2=1,ip2,ip1
     if(i2.lt.i2rev)then
     do i1=i2,i2+ip1-2,2
     do i3=i1,ip3,ip2
     i3rev=i2rev+i3-i2
     tempr=data(i3)
     tempi=data(i3+1)
     data(i3)=data(i3rev)
     data(i3+1)=data(i3rev+1)
     data(i3rev)=tempr
     data(i3rev+1)=tempi
     end do
     end do
     endif
     ibit=ip2/2
1    if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
     i2rev=i2rev-ibit
     ibit=ibit/2
     goto 1
     endif
     i2rev=i2rev+ibit
     end do
     ifp1=ip1
2    if(ifp1.lt.ip2)then
     ifp2=2*ifp1
     theta=isign*6.28318530717959d0/(ifp2/ip1)
     wpr=-2.d0*sin(0.5d0*theta)**2
     wpi=sin(theta)
     wr=1.d0
     wi=0.d0
     do i3=1,ifp1,ip1
     do i1=i3,i3+ip1-2,2
     do i2=i1,ip3,ifp2
     k1=i2
     k2=k1+ifp1
     tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
     tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
     data(k2)=data(k1)-tempr
     data(k2+1)=data(k1+1)-tempi
     data(k1)=data(k1)+tempr
     data(k1+1)=data(k1+1)+tempi
     end do
     end do
     wtemp=wr
     wr=wr*wpr-wi*wpi+wr
     wi=wi*wpr+wtemp*wpi+wi
     end do
     ifp1=ifp2
     goto 2
     endif
     nprev=n*nprev
     end do

     end subroutine fourn

!####################################################################

!-- subroutine four1.f

   SUBROUTINE FOUR1(DATA,NN,ISIGN)
   REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
   DIMENSION DATA(*)
   N=2*NN
   J=1
   DO 11 I=1,N,2
   IF(J.GT.I)THEN
   TEMPR=DATA(J)
   TEMPI=DATA(J+1)
   DATA(J)=DATA(I)
   DATA(J+1)=DATA(I+1)
   DATA(I)=TEMPR
   DATA(I+1)=TEMPI
   ENDIF
   M=N/2
1               IF ((M.GE.2).AND.(J.GT.M)) THEN
   J=J-M
   M=M/2
   GO TO 1
   ENDIF
   J=J+M
11      CONTINUE
   MMAX=2
2       IF (N.GT.MMAX) THEN
   ISTEP=2*MMAX
   THETA=6.28318530717959D0/(ISIGN*MMAX)
   WPR=-2.D0*DSIN(0.5D0*THETA)**2
   WPI=DSIN(THETA)
   WR=1.D0
   WI=0.D0
   DO 13 M=1,MMAX,2
   DO 12 I=M,N,ISTEP
   J=I+MMAX
   TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
   TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
   DATA(J)=DATA(I)-TEMPR
   DATA(J+1)=DATA(I+1)-TEMPI
   DATA(I)=DATA(I)+TEMPR
   DATA(I+1)=DATA(I+1)+TEMPI
12  CONTINUE
   WTEMP=WR
   WR=WR*WPR-WI*WPI+WR
   WI=WI*WPR+WTEMP*WPI+WI
13              CONTINUE
    MMAX=ISTEP
    GO TO 2
    ENDIF
    RETURN
    END SUBROUTINE four1

!####################################################################

!-- subroutine realft.f

  subroutine realft(data,n,isign)
  real*8 wr,wi,wpr,wpi,wtemp,theta
  dimension data(*)
  theta=3.141592653589793d0/dble(n)
  c1=0.5
  if (isign.eq.1) then
  c2=-0.5
  call four1(data,n,+1)
  else
  c2=0.5
  theta=-theta
  endif
  wpr=-2.0d0*dsin(0.5d0*theta)**2
  wpi=dsin(theta)
  wr=1.0d0+wpr
  wi=wpi
  n2p3=2*n+3

  do i=2,n/2
  i1=2*i-1
  i2=i1+1
  i3=n2p3-i2
  i4=i3+1
  wrs=sngl(wr)
  wis=sngl(wi)
  h1r=c1*(data(i1)+data(i3))
  h1i=c1*(data(i2)-data(i4))
  h2r=-c2*(data(i2)+data(i4))
  h2i=c2*(data(i1)-data(i3))
  data(i1)=h1r+wrs*h2r-wis*h2i
  data(i2)=h1i+wrs*h2i+wis*h2r
  data(i3)=h1r-wrs*h2r+wis*h2i
  data(i4)=-h1i+wrs*h2i+wis*h2r
  wtemp=wr
  wr=wr*wpr-wi*wpi+wr
  wi=wi*wpr+wtemp*wpi+wi
  end do

 if (isign.eq.1) then
 h1r=data(1)
 data(1)=h1r+data(2)
 data(2)=h1r-data(2)
 else
 h1r=data(1)
 data(1)=c1*(h1r+data(2))
 data(2)=c1*(h1r-data(2))
 call four1(data,n,-1)
 endif

 return
 end subroutine realft

!####################################################################

	! -----------------------------------------------------------------------------------------
! -- subroutine realft.f	
! ----------------------------------------------------------------------------------------
	SUBROUTINE REALFT2S(N1,ISIGN1,DATA,DATAI,DATABUF)
      	REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      	REAL*8 DATA(N1),DATAI(N1),DATABUF(N1)
	N=N1/2
	IF (ISIGN1.EQ.2) THEN
	ISIGN=1
	ELSE
	ISIGN=-1
	END IF
      	THETA=3.141592653589793D0/DBLE(N)
      	C1=0.5
	IF (ISIGN.EQ.-1) THEN
	DO I=1,N
	DATABUF(I)=DATA(I)
	END DO
	K=0
	DO I=1,2*N,2
	K=K+1
	DATA(I)=DATABUF(K)
	DATA(I+1)=DATAI(K)
	END DO
	END IF
      	IF (ISIGN.EQ.1) THEN
        	C2=-0.5
        	CALL FOUR2(DATA,N,+1)
      	ELSE
        	C2=0.5
        	THETA=-THETA
      	ENDIF
      	WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      	WPI=DSIN(THETA)
      	WR=1.0D0+WPR
     	 WI=WPI
      	N2P3=2*N+3
      	DO 11 I=2,N/2
        	I1=2*I-1
       	 I2=I1+1
        	I3=N2P3-I2
       	 I4=I3+1
       	 WRS=SNGL(WR)
        	WIS=SNGL(WI)
        	H1R=C1*(DATA(I1)+DATA(I3))
        	H1I=C1*(DATA(I2)-DATA(I4))
        	H2R=-C2*(DATA(I2)+DATA(I4))
        	H2I=C2*(DATA(I1)-DATA(I3))
       	 DATA(I1)=H1R+WRS*H2R-WIS*H2I
        	DATA(I2)=H1I+WRS*H2I+WIS*H2R
        	DATA(I3)=H1R-WRS*H2R+WIS*H2I
        	DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        	WTEMP=WR
        	WR=WR*WPR-WI*WPI+WR
        	WI=WI*WPR+WTEMP*WPI+WI
11    	CONTINUE
      	IF (ISIGN.EQ.1) THEN
        	H1R=DATA(1)
        	DATA(1)=H1R+DATA(2)
        	DATA(2)=H1R-DATA(2)
      	ELSE
        	H1R=DATA(1)
        	DATA(1)=C1*(H1R+DATA(2))
        	DATA(2)=C1*(H1R-DATA(2))
        	CALL FOUR2(DATA,N,-1)
      	ENDIF
	IF (ISIGN.EQ.1) THEN
		K=0
        	DO I=1,2*N,2
		K=K+1
		DATAI(K)=DATA(I+1)
		DATA(K)=DATA(I)
		END DO
      	ENDIF
	IF (ISIGN.EQ.-1) THEN
        	DO I=1,2*N
		DATA(I)=DATA(I)/FLOAT(N)
		END DO
      	ENDIF
      	RETURN
      	END
!
! ------------------------------------------------------------------------------------
!-- calcul de TF et TFI pour un tableau de complexes
!-- Cf. numerical recipes p.394	
! ------------------------------------------------------------------------------------	
		subroutine four2(data,nn,isign)
		real*8 wr,wi,wpr,wpi,wtemp,theta
		real*8 tempr,tempi
		real*8 data(2*nn)
		n=2*nn
		j=1
		do 11 i=1,n,2
		if (j.gt.i) then
		tempr=data(j)
		tempi=data(j+1)
		data(j)=data(i)
		data(j+1)=data(i+1)
		data(i)=tempr
		data(i+1)=tempi
		end if
		m=n/2
1		if ((m.ge.2).and.(j.gt.m)) then
		j=j-m
		m=m/2
		goto 1
		end if
		j=j+m
11		continue
		mmax=2
2		if (n.gt.mmax) then
		istep=2*mmax
		theta=6.28318530717959d0/(isign*mmax)
		wpr=-2.d0*dsin(0.5d0*theta)**2
		wpi=dsin(theta)
		wr=1.d0
		wi=0.d0
		do 13 m=1,mmax,2
		do 12 i=m,n,istep
		j=i+mmax
		tempr=wr*data(j)-wi*data(j+1)
		tempi=wr*data(j+1)+wi*data(j)
		data(j)=data(i)-tempr
		data(j+1)=data(i+1)-tempi
		data(i)=data(i)+tempr
		data(i+1)=data(i+1)+tempi
12		continue
		wtemp=wr
		wr=wr*wpr-wi*wpi+wr
		wi=wi*wpr+wtemp*wpi+wi
13		continue
		mmax=istep
		goto 2
		end if
		return
		end

