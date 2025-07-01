subroutine smoothgauss2d(x1d,vel0,freq,nop,n1,n2,h,x21d,itopo,frac1,frac2)
 implicit none

  real:: x1d(nop), x21d(nop), vel0(n1,n2)

  real, allocatable:: x(:,:), x2(:,:)
  !------------------- arrays to be smoothed 
  REAL,ALLOCATABLE ::x1(:,:)
  !------------------- arrays for topography manipulation
  INTEGER:: itopo(n2)
  !------------------- arrays for gaussian distribution
  REAL :: beta2,beta1,vel,frac1,frac2,freq

  INTEGER :: n1,n2,itypsmo,i1,i2,il2,il1,k2,k1,j1,j2,jj1,jj2,it,nop,k
  REAL :: h,tau1,tau2,xl2,xl1,tau12,tau22,d,betatot,tmp,betatot1

  ALLOCATE (x1(n1,n2),x2(n1,n2),x(n1,n2))

  k=0
  do i2=1, n2
     do i1=1, n1
        k=k+1
        x(i1,i2)=x1d(k)
     end do
  end do 
 
  !vel0(:,:)=1./sqrt(x(:,:))

  x2(:,:)=x(:,:)

  !==================================   1D loop
  do i2=1,n2
  !write(*,*) 'i2 / n2 = ',i2,' / ',n2
     do i1=itopo(i2),n1
        vel=vel0(i1,i2)
        tau2=vel*frac2/freq
        tau22=tau2*tau2
        xl2=3.*tau2
        il2=int(xl2/h)+1
        x1(i1,i2)=0.
        betatot=0.
        jj2=0
        do j2=-il2,il2
           jj2=jj2+1
           k2=j2+i2
           if(k2 >= 1 .and. k2 <= n2) then
              d=float(j2)*h
              beta2=exp(-d**2/tau22)
              x1(i1,i2)=x1(i1,i2)+beta2*x(i1,k2)
              betatot=betatot+beta2
           endif
        end do
        x1(i1,i2)=x1(i1,i2)/betatot   ! get the weighted new value
     end do
  end do

  !================================= The other loop

  do i2=1,n2
     do i1=itopo(i2),n1
        vel=vel0(i1,i2)
        tau1=vel*frac1/freq
        tau12=tau1*tau1
        xl1=3.*tau1
        il1=int(xl1/h)+1
        x2(i1,i2)=0.
        betatot=0.
        jj1=0
        do j1=-il1,il1
           jj1=jj1+1
           k1=j1+i1
           if(k1 >= itopo(i2) .and. k1 <= n1) then ! only values below the topography
              d=float(j1)*h
              beta1=exp(-d**2/tau12)
              x2(i1,i2)=x2(i1,i2)+beta1*x1(k1,i2)
              betatot=betatot+beta1
           endif
        end do
        x2(i1,i2)=x2(i1,i2)/betatot ! get the average value
     end do
  end do

  do i2=1, n2
  do i1=1, n1
     if(i1.le.itopo(i2))then
        x2(i1,i2)=0
     endif
  enddo
  enddo

  k=0
  do i2=1, n2
     do i1=1, n1
        k=k+1
        x21d(k)=x2(i1,i2)
     end do
  end do

  DEALLOCATE (x1,x,x2)

end subroutine
