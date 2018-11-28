module Integrals

      use Constants
      use Variables
      use Vectors

contains

real(dp) function Ana_in_in(m,A1,A2,k1,k2,a)

      implicit none
      double precision, external:: s13adf, ei
      integer, intent(in) :: m
      integer :: i, ifail
      real(dp) :: x, integral,integral_err, f, f2, interval
      real(dp), intent(in) :: A1,A2,k1,k2,a

      Ana_in_in =  &
      (A1*A2/(4*pi*k1*k2))*(a*(k1-k2)*s13adf((k1-k2)*a,ifail)+a*(k1+k2)*s13adf((k1+k2)*a,ifail) - &
              2*sin(k1*a)*sin(k2*a))/(2*a)
      return

end function Ana_in_in

!Overlap integral on simgle dot, numerical integration
real(dp) function OverlapNum(A1,A2,B1,B2,kin1,kin2,kout1,kout2,r)

      implicit none
      integer :: i, m
      real(dp) :: x, fin, fout, interval
      real(dp), intent(in) :: A1,A2,B1,B2,kin1,kin2,kout1,kout2,r

      i=0
     
      m=1000

      interval= r/m

      do i=0,m-1
      x=0.5*interval+i*interval
      fin = fin + A1*sin(kin1*x)*A2*sin(kin2*x)
      enddo

      do i=m,2*m
      x=0.5*interval+i*interval
      fout = fout + B1*B2*exp(-kout1*x)*exp(-kout2*x)
      enddo

      OverlapNum=abs(fin+fout)*interval

end function OverlapNum

!Overlap integral on simgle dot, analytical integration
real(dp) function OverlapAna(A1,A2,B1,B2,kin1,kin2,kout1,kout2,r)

      implicit none
      real(dp), intent(in) :: A1,A2,B1,B2,kin1,kin2,kout1,kout2,r
      
      OverlapAna=(A1*A2*(kin2*sin(kin1*r)*cos(kin2*r)-kin1*cos(kin1*r)*sin(kin2*r))/(kin1**2-kin2**2)  &
                 + 1.0*B1*B2*exp(-1*(kout1+kout2)*r)/(kout1+kout2))

end function OverlapAna

real(dp) function O_in_out_dimer(m,A1,B1,kin1,kout1,A2,B2,kin2,kout2,a,b,l)

      implicit none
      integer, intent(in) :: m
      integer :: i
      real(dp) :: x, y, integral,integral_err, f1, f2, f3, interval
      real(dp), intent(in) :: A1,B1,kin1,kout1,A2,B2,kin2,kout2,a, b, l

      open(40,file='Overlap-dimer.dat')

      f1= 0.0
      f2= 0.0
      f3= 0.0
      i=0
      interval= (a+b+l)/m


      do i=0,int((a*m)/(a+b+l))-1
      x=0.5*interval+i*interval
      y=a+b+l-x
      f1 = f1 + A1*sin(kin1*x)*B2*exp(-1*kout2*y)/(2*y*x)
      enddo
      do i=int((a*m)/(a+b+l)),int(((a+l)*m)/(a+b+l))
      x=0.5*interval+i*interval
      y=a+b+l-x
      f2 = f2 + B1*exp(-1*kout1*x)*B2*exp(-1*kout2*y)/(2*x*y)
      enddo
      do i=int(((a+l)*m)/(a+b+l))+1,m-1
      x=0.5*interval+i*interval
      y=a+b+l-x
      f3 = f3 + B1*exp(-1*kout1*x)*A2*sin(kin2*y)/(2*x*y)
      enddo

      write(40,*) f1*interval, f2*interval, f3*interval, (f1+f2+f3)*interval

      O_in_out_dimer=abs(f1+f2+f3)*interval

end function O_in_out_dimer

real(dp) function O_in_out_dimer_cart(m,A1,B1,kin1,kout1,A2,B2,kin2,kout2,a,b,l)

      implicit none
      integer :: m, x1,x2,x3
      integer :: xx,yy,zz,mm
      real(dp) :: x,y,z,integral,integral_err, f1, f2, f3, interval, f
      real(dp), intent(in) :: A1,B1,kin1,kout1,A2,B2,kin2,kout2,a, b, l
      real(dp),dimension(3) :: r1, r2 

      open(40,file='Overlap-dimer.dat')

      m=100
       
      interval= (a+b+l)/m
      f=0.0
      do xx=0,(m-1)
      x=0.5*interval+xx*interval
      f1 = 0.0 
      f2 = 0.0 
      f3 = 0.0 
      do yy=-1*int(b/interval),int(b/interval)
      y=0.5*interval+yy*interval
      do zz=-1*int(b/interval),int(b/interval)
      z=0.5*interval+zz*interval
      r1(1)=x
      r1(2)=y
      r1(3)=z
      r2(1)=a+b+l-x
      r2(2)=y
      r2(3)=z 
        
      if ((norm2(r1) .le. a) .and. (norm2(r2) .gt. b)) then
      f1 = f1 + A1*sin(kin1*norm2(r1))*B2*exp(-1*kout2*norm2(r2))/(4*pi*norm2(r1)*norm2(r2))
      else if ((norm2(r1) .gt. a) .and. (norm2(r2) .le. b)) then
      f2 = f2 + B1*exp(-1*kout1*norm2(r1))*A2*sin(kin2*norm2(r2))/(4*pi*norm2(r1)*norm2(r2)) 
      else if ((norm2(r1) .gt. a) .and. (norm2(r2) .gt. b)) then
      f3 = f3 + B1*exp(-1*kout1*norm2(r1))*B2*exp(-1*kout2*norm2(r2))/(4*pi*norm2(r1)*norm2(r2))
      endif
      enddo
      enddo
      f = f + f1 + f2 + f3
!      write (40,*) x, f1*interval**3, f2*interval**3, f3*interval**3, (f1+f2+f3)*interval**3
      enddo

      O_in_out_dimer_cart=abs(f)*interval**3

end function O_in_out_dimer_cart

!Computation of dimer transition dipole moments integrating in a grid
real(dp) function TransDip_dimer_cart(A1,B1,kin1,kout1,A2,B2,kin2,kout2,a,b,l)

      implicit none
      integer :: x1,x2,x3
      integer :: m, xx,yy,zz,mm
      real(dp) :: x,y,z, f1, f2, f3,  interval, f
      real(dp), intent(in) :: A1,B1,kin1,kout1,A2,B2,kin2,kout2,a, b, l
      real(dp),dimension(3) :: r1, r2 

      !write(6,*) a, b, l, A1,B1,kin1,kout1,A2,B2,kin2,kout2

      m=100
       
      interval= (a+b+l)/m
      f=0.0
      do xx=0,(m-1)
      x=0.5*interval+xx*interval
      f1 = 0.0 
      f2 = 0.0 
      f3 = 0.0 
      do yy=-1*int(b/interval),int(b/interval)
      y=0.5*interval+yy*interval
      do zz=-1*int(b/interval),int(b/interval)
      z=0.5*interval+zz*interval
      r1(1)=x
      r1(2)=y
      r1(3)=z
      r2(1)=a+b+l-x
      r2(2)=y
      r2(3)=z 
        
      if ((norm2(r1) .le. a) .and. (norm2(r2) .gt. b)) then
      f1 = f1 + A1*sin(kin1*norm2(r1))*B2*exp(-1*kout2*norm2(r2))/(4*pi*norm2(r2))
      else if ((norm2(r1) .gt. a) .and. (norm2(r2) .le. b)) then
      f2 = f2 + B1*exp(-1*kout1*norm2(r1))*A2*sin(kin2*norm2(r2))/(4*pi*norm2(r2)) 
      else if ((norm2(r1) .gt. a) .and. (norm2(r2) .gt. b)) then
      f3 = f3 + B1*exp(-1*kout1*norm2(r1))*B2*exp(-1*kout2*norm2(r2))/(4*pi*norm2(r2))
      endif
      enddo
      enddo
      f = f + f1 + f2 + f3
      enddo
      TransDip_dimer_cart=elec*abs(f)*interval**3

end function TransDip_dimer_cart

!Computation of dimer transition dipole moments integrating using Monte-Carlo
real(dp) function TransDip_dimer_MC(A1,B1,kin1,kout1,A2,B2,kin2,kout2,a,b,l)

      implicit none
      integer :: m, i, i1, i2, i3
      real(dp) :: f1, f2, f3, interval, f, maxrad, vola, volb, volout, side
      real(dp), intent(in) :: A1,B1,kin1,kout1,A2,B2,kin2,kout2,a, b, l
      real(dp), allocatable:: r1(:,:) , r2(:,:) , r1Anorm(:), r1Bnorm(:), RAB(:), r1Nnorm(:)

      
      !open(40,file='Rtest.dat')

      m=500000

      allocate(RAB(3))
      allocate(r1(m,3))
      allocate(r1Nnorm(m))
      allocate(r1Anorm(m))
      allocate(r1Bnorm(m))

      RAB(1)=a+b+l
      RAB(2)=0
      RAB(3)=0

      i1=0
      i2=0
      i3=0
      f1 = 0.0
      f2 = 0.0
      f3 = 0.0
      maxrad=max(a,b)

      vola=(4.0/3)*pi*a**3
      volb=(4.0/3)*pi*b**3
      side=0.5d-9
      volout=abs((2*maxrad+2*side)**2*(2*a+l+2*b+2*side)-(vola+volb))

      call random_seed()
      call random_number(r1(:,1))
      call random_number(r1(:,2))
      call random_number(r1(:,3))

r1Anorm(:)=sqrt(((2*a+2*b+l+2*side)*r1(:,1)-(a+side))**2+((2*maxrad+2*side)*r1(:,2)-(maxrad+side))**2+ &
               ((2*maxrad+2*side)*r1(:,3)-(maxrad+side))**2)
r1Bnorm(:)=sqrt(((2*a+2*b+l+2*side)*r1(:,1)-(RAB(1)+a+side))**2+((2*maxrad+2*side)*r1(:,2)-(maxrad+side))**2+&
                ((2*maxrad+2*side)*r1(:,3)-(maxrad+side))**2)
r1Nnorm(:)=sqrt(((2*a+2*b+l+2*side)*r1(:,1)-(l/2+2*a+side))**2+((2*maxrad+2*side)*r1(:,2)-(maxrad+side))**2+&
                ((2*maxrad+2*side)*r1(:,2)-(maxrad+side))**2)

      do i=1,m
      if ((r1Anorm(i) .le. a) .and. (r1Bnorm(i) .gt. b)) then
      f1 = f1 + r1Anorm(i)*exp(-1.0*kout2*RAB(1))*sin(kin1*r1Anorm(i))*exp(-1.0*kout2*r1Bnorm(i))/&
                (4.0*pi*r1Anorm(i)*(r1Bnorm(i)+RAB(1)))
      i1 = i1 + 1
      else if ((r1Anorm(i) .gt. a) .and. (r1Bnorm(i) .le. b)) then
      f2 = f2 + r1Anorm(i)*exp(-1.0*kout1*r1Anorm(i))*sin(kin2*(r1Bnorm(i)+RAB(1)))/&
                (4.0*pi*r1Anorm(i)*(r1Bnorm(i)+RAB(1)))
      i2 = i2 + 1
      else if ((r1Anorm(i) .gt. a) .and. (r1Bnorm(i) .gt. b)) then
      f3 = f3 + r1Anorm(i)*exp(-1.0*kout2*RAB(1))*exp(-1.0*kout1*r1Anorm(i))*exp(-1.0*kout2*r1Bnorm(i))/&
                (4.0*pi*r1Anorm(i)*(r1Bnorm(i)+RAB(1)))
      i3 = i3 + 1
      endif
      enddo


write(6,*) abs(elec*(A1*B2*f1*vola/i1)), abs(elec*(B1*A2*f2*volb/i2)), abs(elec*(B1*B2*f3*volout/i3))

      TransDip_dimer_MC=elec*(abs((A1*B2*f1*vola/i1))+abs((B1*A2*f2*volb/i2))+abs((B1*B2*f3*volout/i3)))

end function TransDip_dimer_MC

!Computation of dimer transition dipole moments integrating using Monte-Carlo
real(dp) function TransDip_dimer_MC_off(A1,B1,kin1,kout1,A2,B2,kin2,kout2,a,b,l)

      implicit none
      integer :: m, i, i1, i2, i3
      real(dp) :: fx1, fx2, fx3, fy1, fy2, fy3,fz1, fz2, fz3,interval, f, maxrad, vola, volb, volout
      real(dp), intent(in) :: A1,B1,kin1,kout1,A2,B2,kin2,kout2,a, b, l
      real(dp), allocatable:: r1(:,:) , r2(:,:) , r1A(:,:), r1B(:,:), RAB(:), r1N(:,:)

      
      open(41,file='Rtest.dat')

      m=1000000

      allocate(RAB(3))
      allocate(r1(m,3))
      allocate(r1N(m,3))
      allocate(r1A(m,3))
      allocate(r1B(m,3))

      RAB(1)=a+b+l
      RAB(2)=0
      RAB(3)=0

      i1=0
      i2=0
      i3=0
      fx1 = 0.0
      fx2 = 0.0
      fx3 = 0.0
      fy1 = 0.0
      fy2 = 0.0
      fy3 = 0.0
      fz1 = 0.0
      fz2 = 0.0
      fz3 = 0.0
      maxrad=max(a,b)

      vola=(4.0/3)*pi*a**3
      volb=(4.0/3)*pi*b**3
      volout=abs(((2*maxrad+2*side)**2*(2*a+l+2*b+2*side)-(vola+volb)))

      call random_seed()
      call random_number(r1(:,1))
      call random_number(r1(:,2))
      call random_number(r1(:,3))

!bottom back corner
!r1N(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))
!r1N(:,2)=(r1(:,2)*(2*maxrad+2*side))
!r1N(:,3)=(r1(:,3)*(2*maxrad+2*side))
!r1A(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+a)
!r1A(:,2)=(r1(:,2)*(2*maxrad+2*side))-(side+maxrad)
!r1A(:,3)=(r1(:,3)*(2*maxrad+2*side))-(side+maxrad)
!r1B(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+2*a+l+b)
!r1B(:,2)=(r1(:,2)*(2*maxrad+2*side))-(side+maxrad)
!r1B(:,3)=(r1(:,3)*(2*maxrad+2*side))-(side+maxrad)

!Point between QD surfaces
r1N(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+2*a+l/2)
r1N(:,2)=(r1(:,2)*(2*maxrad+2*side))-(maxrad+side)
r1N(:,3)=(r1(:,3)*(2*maxrad+2*side))-(maxrad+side)
r1A(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+a)
r1A(:,2)=(r1(:,2)*(2*maxrad+2*side)-(maxrad+side))
r1A(:,3)=(r1(:,3)*(2*maxrad+2*side)-(maxrad+side))
r1B(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+2*a+l+b)
r1B(:,2)=(r1(:,2)*(2*maxrad+2*side)-(maxrad+side))
r1B(:,3)=(r1(:,3)*(2*maxrad+2*side)-(maxrad+side))

!center dot A
!r1N(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+a)
!r1N(:,2)=(r1(:,2)*(2*maxrad+2*side))-(maxrad+side)
!r1N(:,3)=(r1(:,3)*(2*maxrad+2*side))-(maxrad+side)
!r1A(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+a)
!r1A(:,2)=(r1(:,2)*(2*maxrad+2*side)-(maxrad+side))
!r1A(:,3)=(r1(:,3)*(2*maxrad+2*side)-(maxrad+side))
!r1B(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+2*a+l+b)
!r1B(:,2)=(r1(:,2)*(2*maxrad+2*side)-(maxrad+side))
!r1B(:,3)=(r1(:,3)*(2*maxrad+2*side)-(maxrad+side))

      do i=1,m
      if ((norm2(r1A(i,:)) .le. a) .and. (norm2(r1B(i,:)) .gt. b)) then
      fx1 = fx1 + r1N(i,1)*sin(kin1*norm2(r1A(i,:)))*exp(-1.0*kout2*norm2(r1B(i,:)))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      fy1 = fy1 + r1N(i,2)*sin(kin1*norm2(r1A(i,:)))*exp(-1.0*kout2*norm2(r1B(i,:)))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      fz1 = fz1 + r1N(i,3)*sin(kin1*norm2(r1A(i,:)))*exp(-1.0*kout2*norm2(r1B(i,:)))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      i1 = i1 + 1
      else if ((norm2(r1A(i,:)) .gt. a) .and. (norm2(r1B(i,:)) .le. b)) then
      fx2 = fx2 + r1N(i,1)*exp(-1.0*kout1*norm2(r1A(i,:)))*sin(kin2*(norm2(r1B(i,:))))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      fy2 = fy2 + r1N(i,2)*exp(-1.0*kout1*norm2(r1A(i,:)))*sin(kin2*(norm2(r1B(i,:))))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      fz2 = fz2 + r1N(i,3)*exp(-1.0*kout1*norm2(r1A(i,:)))*sin(kin2*(norm2(r1B(i,:))))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      i2 = i2 + 1
      else if ((norm2(r1A(i,:)) .gt. a) .and. (norm2(r1B(i,:)) .gt. b)) then
      fx3 = fx3 + r1N(i,1)*exp(-1.0*kout1*norm2(r1A(i,:)))*exp(-1.0*kout2*norm2(r1B(i,:)))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      fy3 = fy3 + r1N(i,2)*exp(-1.0*kout1*norm2(r1A(i,:)))*exp(-1.0*kout2*norm2(r1B(i,:)))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      fz3 = fz3 + r1N(i,3)*exp(-1.0*kout1*norm2(r1A(i,:)))*exp(-1.0*kout2*norm2(r1B(i,:)))/&
                (norm2(r1A(i,:))*norm2(r1B(i,:)))
      i3 = i3 + 1
      endif
      enddo

TransDip_dimer_MC_off=sqrt(((elec/(4*pi))*(((A1*B2*fx1*vola/i1)+(B1*A2*fx2*volb/i2)+(B1*B2*fx3*volout/i3))))**2+ &
                           ((elec/(4*pi))*(((A1*B2*fy1*vola/i1)+(B1*A2*fy2*volb/i2)+(B1*B2*fy3*volout/i3))))**2+ &
                           ((elec/(4*pi))*(((A1*B2*fz1*vola/i1)+(B1*A2*fz2*volb/i2)+(B1*B2*fz3*volout/i3))))**2)


end function TransDip_dimer_MC_off

!Normalization of WF analytical
real(dp) function Norm_Ana(r,A,B,kin,kout)

      implicit none
      real(dp) :: r,A,B,kin,kout

      Norm_Ana=-A**2*(sin(2*kin*r)-2*kin*r)/(4*kin) + B**2*exp(-2*kout*r)/(2*kout)

end function Norm_Ana

!Normalization of WF in radial coordinates
real(dp) function Norm_Num(r,A,B,kin,kout)

      implicit none
      integer :: i, m
      real(dp) :: x, fin, fout, interval
      real(dp), intent(in) :: A,B,kin,kout,r

      fin= 0.0
      fout= 0.0
      m=2000
      i=0
      x=0

      interval= r/m

      do i=0,m-1
      x=0.5*interval+i*interval
      fin = fin + (A*sin(kin*x))**2
      enddo

      do i=m,4*m
      x=0.5*interval+i*interval
      fout = fout + (B*exp(-1*kout*x))**2
      enddo

      Norm_Num=abs(fin+fout)*interval

end function Norm_Num

!Normalization of WF in radial coordinates
real(dp) function Norm_MC_off(A1,B1,kin,kout,a,b,l)

      implicit none
      integer :: i, m, i1, i2, i3
      real(dp) :: fx1, fx2, fx3, fy1, fy2, fy3,fz1, fz2, fz3,interval, f, maxrad, vola, volb, volout
      real(dp), allocatable:: r1(:,:) , r2(:,:) , r1A(:,:), r1B(:,:), RAB(:), r1N(:)
      real(dp) :: x, fin, fout
      real(dp), intent(in) :: A1,B1,kin,kout,a,b,l

      open(41,file='Rtest.dat')

      m=10000

      allocate(RAB(3))
      allocate(r1(m,3))
      allocate(r1N(3))
      allocate(r1A(m,3))
      allocate(r1B(m,3))

      RAB(1)=a+b+l
      RAB(2)=0
      RAB(3)=0

      i1=0
      i2=0
      i3=0
      fx1 = 0.0
      fx2 = 0.0
      fx3 = 0.0
      fy1 = 0.0
      fy2 = 0.0
      fy3 = 0.0
      fz1 = 0.0
      fz2 = 0.0
      fz3 = 0.0
      maxrad=max(a,b)

      vola=(4.0/3)*pi*a**3
      volb=(4.0/3)*pi*b**3
      volout=abs(((2*maxrad+2*side)**2*(2*a+l+2*b+2*side)-(vola+volb)))

      call random_seed()
      call random_number(r1(:,1))
      call random_number(r1(:,2))
      call random_number(r1(:,3))

r1A(:,1)=(r1(:,1)*(2*a+2*b+l+2*side))-(side+a)
r1A(:,2)=(r1(:,2)*(2*maxrad+2*side)-(maxrad+side))
r1A(:,3)=(r1(:,3)*(2*maxrad+2*side)-(maxrad+side))
r1N(1)=((maxrad+side))
r1N(2)=((maxrad+side))
r1N(3)=((maxrad+side))


      do i=1,m
      if ((norm2(r1A(i,:)) .le. a)) then
      !if ((norm2(r1A(i,:)) .le. a) .and. (norm2(r1B(i,:)) .gt. b)) then
      fx1 = fx1 + (sin(kin*norm2(r1A(i,:))))**2/(norm2(r1A(i,:)))**2
      !fx1 = fx1 + (sin(kin*norm2(r1A(i,:)))*cos(kin*norm2(r1N(:)))+cos(kin*norm2(r1A(i,:)))*sin(kin*norm2(r1N(:))))**2*&
      !            (norm2(r1A(i,:))+norm2(r1N(:)))**2
      !fx1 = fx1 + (sin(kin*(r1A(i,1)))*cos(kin*(r1N(1)))+cos(kin*(r1A(i,1)))*sin(kin*(r1N(1))))**2*&
      !            ((r1A(i,1))*(r1N(1)))**2
      !fy1 = fy1 + (sin(kin*(r1A(i,2)))*cos(kin*(r1N(2)))+cos(kin*(r1A(i,2)))*sin(kin*(r1N(2))))**2*&
      !            ((r1A(i,2))*(r1N(2)))**2
      !fz1 = fz1 + (sin(kin*(r1A(i,3)))*cos(kin*(r1N(3)))+cos(kin*(r1A(i,3)))*sin(kin*(r1N(3))))**2*&
      !            ((r1A(i,3))*(r1N(3)))**2
      !write(41,'(6f10.4)') r1N(i,:)*1e9, r1N(i,:)*1e9, r1N(i,:)*1e9, r1A(i,:)*1e9, r1A(i,:)*1e9, r1A(i,:)*1e9
      i1 = i1 + 1
      else if ((norm2(r1A(i,:)) .gt. a)) then
!write(41,*) R1A(i,1), R1A(i,2), R1A(i,3)
      !else if ((norm2(r1A(i,:)) .gt. a) .and. (norm2(r1B(i,:)) .le. b)) then
      fx2 = fx2 + (exp(-1.0*kout*norm2(r1A(i,:))))**2/(norm2(r1A(i,:)))**2
      !fx2 = fx2 + (exp(-1.0*kout*norm2(r1A(i,:)))*exp(-1.0*kout*norm2(r1N(:))))**2/(norm2(r1A(i,:))+norm2(r1N(:)))**2
      !fx2 = fx2 + (exp(-1.0*kout*(r1A(i,1)))*exp(-1.0*kout*(r1N(1))))**2/((r1A(i,1))*(r1N(1)))**2
      !fy2 = fy2 + (exp(-1.0*kout*(r1A(i,2)))*exp(-1.0*kout*(r1N(2))))**2/((r1A(i,2))*(r1N(2)))**2 
      !fz2 = fz2 + (exp(-1.0*kout*(r1A(i,3)))*exp(-1.0*kout*(r1N(3))))**2/((r1A(i,3))*(r1N(3)))**2
      i2 = i2 + 1
      endif
      enddo

write(6,*) i1, fx1, i2, fx2

      Norm_MC_off=(A1*A1*(fx1*vola/i1)+B1*B1*(fx2*volout/i2))/(4*pi)
      !Norm_MC_off=(A1*A1*(fx1*vola/i1+fy1*vola/i1+fz1*vola/i1)+B1*B1*(fx2*volout/i2+fy2*volout/i2+fz2*volout/i2))/(4*pi)

end function Norm_MC_off

!Normalization of WF in cartesian coordinates
real(dp) function Norm_cart(AB1,AB2,k1,k2,a)

      implicit none
      double precision, external:: s13adf, ei
      integer :: m
      integer :: i, xx1,yy1,zz1
      real(dp) :: f1, f2, interval
      real(dp) :: AB1,AB2,k1,k2,a,x1,y1,z1, b
      real(dp), dimension(3) :: r1

      m=100
      f1= 0.0
      f2= 0.0
      i=0

      interval= 2*a/m

      do xx1=0,m
      x1=0.5*interval+xx1*interval
      do yy1=0,m
      y1=0.5*interval+yy1*interval
      do zz1=0,m
      z1=0.5*interval+zz1*interval
      r1(1)=x1
      r1(2)=y1
      r1(3)=z1


      if (norm2(r1) .le. a) then
      f1 = f1 + (AB1*sin(k1*norm2(r1))/(sqrt(4*pi)*norm2(r1)))**2
      else if (norm2(r1) .gt. a) then
      f2 = f2 + (AB2*exp(-1*k2*norm2(r1))/(sqrt(4*pi)*norm2(r1)))**2
      endif
      
      enddo
      enddo
      enddo

      Norm_cart=8*(f1+f2)*interval**3

end function Norm_cart

!Normalization of WF in cartesian coordinates
real(dp) function Norm_cart_off(AB1,AB2,k1,k2,a)

      implicit none
      double precision, external:: s13adf, ei
      integer :: m
      integer :: i, xx1,yy1,zz1
      real(dp) :: f1, f2, interval
      real(dp) :: AB1,AB2,k1,k2,a,x1,y1,z1, b
      real(dp), dimension(3) :: r1

      m=500
      f1= 0.0
      f2= 0.0
      i=0

      interval= 4*a/m

      do xx1=0,m
      x1=0.5*interval+xx1*interval
      do yy1=0,m
      y1=0.5*interval+yy1*interval
      do zz1=0,m
      z1=0.5*interval+zz1*interval
      r1(1)=2*a-x1
      r1(2)=2*a-y1
      r1(3)=2*a-z1


      if (norm2(r1) .le. a) then
      f1 = f1 + (AB1*sin(k1*norm2(r1))/(sqrt(4*pi)*norm2(r1)))**2
      else if (norm2(r1) .gt. a) then
      f2 = f2 + (AB2*exp(-1*k2*norm2(r1))/(sqrt(4*pi)*norm2(r1)))**2
      endif
      
      enddo
      enddo
      enddo

      Norm_cart_off=(f1+f2)*interval**3

end function Norm_cart_off

!Normalization of WF in cartesian coordinates using Monte Carlo method 
real(dp) function Norm_cart_Rdm(AB1,AB2,k1,k2,a,b)
      implicit none
      integer(kind=8) :: m, i, iin, iout, j
      real(dp) :: x, y, z, f, fin, fout, f1in, f2in, f1out, f2out, AB1,AB2,k1,k2,a,b , n, maxv
      real(dp), dimension(3) :: r
      integer(kind=8), dimension(18) :: list
      
      f= 0.0
      f1in= 0.0
      f2in= 0.0
      f1out= 0.0
      f2out= 0.0
      i=0
      x=0
      iin =0
      iout =0
      !write(6,*) huge(m)
      !m=1000000000
      list = (/10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000, 100000000, &
             200000000, 500000000, 1000000000, 2000000000, 5000000000 /)

      do j = 1,size(list)

      do i=1,list(j)
      call random_number(x)
      call random_number(y)
      call random_number(z)
      r(1)=2*b*x-b
      r(2)=2*b*y-b
      r(3)=2*b*z-b

      if (norm2(r) .le. a) then
      fin = (AB1*sin(k1*norm2(r))/(sqrt(4*pi)*norm2(r)))**2
      iin = iin +1
      f1in = f1in + fin
      f2in = f2in + fin**2
      else if (norm2(r) .gt. a) then
      fout = (AB2*exp(-1*k2*norm2(r))/(sqrt(4*pi)*norm2(r)))**2
      iout = iout + 1
      f1out = f1out + fout
      f2out = f2out + fout**2
      endif

      enddo

!write(6,*) (f/iin)*(4*pi/3)*a**3 , (f1/iout)*((2*b)**3-(4*pi/3)*a**3)  , (f/iin)*(4*pi/3)*a**3 + (f1/iout)*((2*b)**3-(4*pi/3)*a**3) 
write(6,*)list(j),(f1in/iin)*(4*pi/3)*a**3,sqrt(((f2in/iin)-(f1in/iin)**2)/iin)*(4*pi/3)*a**3, &
                  (f1out/iout)*((2*b)**3-(4*pi/3)*a**3),sqrt(((f2out/iout)-(f1out/iout)**2)/iout)*((2*b)**3-(4*pi/3)*a**3),&
                  (f1in/iin)*(4*pi/3)*a**3 + (f1out/iout)*((2*b)**3-(4*pi/3)*a**3),&
                  sqrt((sqrt(((f2in/iin)-(f1in/iin)**2)/iin)*(4*pi/3)*a**3)**2+&
                  (sqrt(((f2out/iout)-(f1out/iout)**2)/iout)*((2*b)**3-(4*pi/3)*a**3))**2)

enddo


      Norm_cart_Rdm=4*(f/m)*(b)**3

end function Norm_cart_Rdm

!Normalization of WF in cartesian coordinates using Monte Carlo method 
real(dp) function Norm_cart_MC(m,AB1,AB2,k1,k2,a,b)
      implicit none
      integer :: m, i, i1, i2, iin, iout
      real(dp) :: x, y, z, f, f1, f2, AB1,AB2,k1,k2,a,b , n, maxv
      real(dp), dimension(3) :: r

      f= 0.0
      f1= 0.0
      f2= 0.0
      i=0
      x=0
      m=1000000
      maxv=1.78e26 !1.7711381564138995e26
      i1=0 
      do i=1,m
      call random_number(x)
      call random_number(y)
      call random_number(z)
      call random_number(n)
      r(1)=x*b
      r(2)=y*b
      r(3)=z*b
      n=n*maxv

      if (norm2(r) .le. a) then
      f = (AB1*sin(k1*norm2(r))/(sqrt(4*pi)*norm2(r)))**2
      !f1 = f1 + f
      !f2 = f2 + f**2
      iin = iin + 1
        if (n .le. f) then
        i1=i1+1
        endif
      else if (norm2(r) .gt. a) then
      f = (AB2*exp(-1*k2*norm2(r))/(sqrt(4*pi)*norm2(r)))**2
      !f1 = f1 + f
      !f2 = f2 + f**2
      iout = iout + 1 
        if (n .le. f) then
        i2=i2+1
        endif
      endif
 
      enddo

write(6,*) (1.0*i1/iin)*(4*pi/24)*maxv*a**3 , (1.0*i2/iout)*(b**3-(4*pi/24)*a**3)*maxv  ,&
                                          (1.0*i1/iin)*(4*pi/24)*a**3*maxv + (1.0*i2/iout)*(b**3-(4*pi/24)*a**3)*maxv 
!write(6,*) (f1*i1/m)*maxv*b**3 !, (sqrt((f2/m)-(f1/m))/(m-1))*b**3*maxv   !sqrt((f2-f1**2)/m)*maxv*(2*b)**3

      Norm_cart_MC=(f/m)*(b)**3

end function Norm_cart_MC


!Analytical term of overlap between hout and eout
real(dp) function Ana_out_out(m,B1,B2,k1,k2,a)

      implicit none
      double precision, external:: ei
      integer :: m, i, ifail
      real(dp) :: B1,B2,k1,k2,a, x, integral,integral_err, f, f2

      Ana_out_out = (B1*B2/(4*pi*k1*k2))*(-1*(k1+k2)*ei(-1*(k1+k2)*a)+exp(-1*a*(k1+k2))/a)
      return

end function Ana_out_out

!Transition dipole inside dot volume, spherical symmetry
real(dp) function TransDip_Num(A1,A2,B1,B2,kin1,kin2,kout1,kout2,r)

      implicit none
      integer :: i, m
      real(dp) :: x, fin, fout, interval
      real(dp), intent(in) :: A1,A2,B1,B2,kin1,kin2,kout1,kout2,r

      i=0
      m=1000

      interval= r/m

      do i=0,m-1
      x=0.5*interval+i*interval
      fin = fin + elec*x*A1*sin(kin1*x)*A2*sin(kin2*x)
      enddo

      do i=m,2*m
      x=0.5*interval+i*interval
      fout = fout + elec*x*B1*exp(-1.0*kout1*x)*B2*exp(-1.0*kout2*x)
      enddo

      TransDip_Num=abs(fin+fout)*interval

end function TransDip_Num

!Transition dipole one dot, analytic, spherical symmetry
real(dp) function TransDip_Ana(A1,A2,B1,B2,kin1,kin2,kout1,kout2,r)

      implicit none
      real(dp) :: A1,A2,B1,B2,kin1,kin2,kout1,kout2,r,rmin,rmax

      rmin=0.0
      rmax=2*r

      TransDip_Ana=-1.0*elec*(((A1*A2*0.5*((r*sin(r*(kin1-kin2))/(kin1-kin2))-(r*sin(r*(kin1+kin2))/(kin1+kin2))+&
                       (cos(r*(kin1-kin2))/(kin1-kin2)**2)-(cos(r*(kin1+kin2))/(kin1+kin2)**2))) - &
                       A1*A2*0.5*((cos(rmin*(kin1-kin2))/(kin1-kin2)**2)-(cos(rmin*(kin1+kin2))/(kin1+kin2)**2))) + &
                       (B1*B2*((-1.0*(exp(-1.0*rmax*(kout1+kout2))*(kout1*rmax+kout2*rmax+1)/(kout1+kout2)**2))+&
                       ((exp(-1.0*r*(kout1+kout2))*(kout1*r+kout2*r+1)/(kout1+kout2)**2)))))

end function TransDip_Ana

real(dp) function TransDip_EMA(Egap,r)

      implicit none
      real(dp) :: Egap, r

      TransDip_EMA=sqrt(((elec**2*hbar**2)/(6*Egap**2*m0))* &
       ((m0/me)-1)*((V0*(V0+(0.42*elec)))/(V0+(2/3)*(0.42*elec))))

end function TransDip_EMA

real(dp) function TransDip_cart(A1,B1,kin1,kout1,A2,B2,kin2,kout2,a)

      implicit none      
      integer :: m, i1,i2,i3,i4, i
      integer :: xx1,yy1,zz1,xx2,yy2,zz2,mm, oo,o
      real(dp) :: x1,y1,z1,x2,y2,z2,integral,integral_err, f1, f2, f3, f4, interval, f
      real(dp), intent(in) :: A1,B1,kin1,kout1,A2,B2,kin2,kout2,a
      real(dp),dimension(3) :: r1, r2, r, RAB

      open(40,file='Dinin.dat')
      open(41,file='Dinout.dat')
      open(42,file='Doutin.dat')
      open(43,file='Doutout.dat')

      m=10

      interval= a/m
      f=0.0
      f1 = 0.0
      f2 = 0.0
      f3 = 0.0
      f4 = 0.0

      do xx1=-1*nint(a/interval),nint(a/interval)
      x1=0.5*interval+xx1*interval
      do yy1=-1*nint(a/interval),nint(a/interval)
      y1=0.5*interval+yy1*interval
      do zz1=-1*nint(a/interval),nint(a/interval)
      z1=0.5*interval+zz1*interval
      r1(1)=x1
      r1(2)=y1
      r1(3)=z1

         do xx2=-1*nint(a/interval),nint(a/interval)
         x2=0.5*interval+xx2*interval
         do yy2=-1*nint(a/interval),nint(a/interval)
         y2=0.5*interval+yy2*interval
         do zz2=-1*nint(a/interval),nint(a/interval)
         z2=0.5*interval+zz2*interval
         r2(1)=x2
         r2(2)=y2
         r2(3)=z2

         r = r2 - r1

      if (norm2(r) .lt. 1e-12 ) then
      exit
      else if ((norm2(r1) .le. a) .and. (norm2(r2) .le. a)) then
         f1 = f1 + A1*A2*norm2(r)*(sin(kin1*norm2(r1))*sin(kin2*norm2(r2)))/(norm2(r1)*norm2(r2))
      else if ((norm2(r1) .le. a) .and. (norm2(r2) .gt. a)) then
         f2 = f2 + A1*B2*norm2(r)*(sin(kin1*norm2(r1))*exp(-1*kout2*norm2(r2)))/(norm2(r1)*norm2(r2))
      else if ((norm2(r1) .gt. a) .and. (norm2(r2) .le. a)) then
         f3 = f3 + B1*A2*norm2(r)*(exp(-1*kout1*norm2(r1))*sin(kin2*norm2(r2)))/(norm2(r1)*norm2(r2))
      else if ((norm2(r1) .gt. a) .and. (norm2(r2) .gt. a)) then
         f4 = f4 + B1*B2*norm2(r)*(exp(-1*kout1*norm2(r1))*exp(-1*kout2*norm2(r2)))/(norm2(r1)*norm2(r2))
      endif
         enddo
         enddo
         enddo
     enddo
     enddo
     enddo

      TransDip_cart=(elec/(4*pi))*abs(f1+f2+f3+f4)*interval**6

end function TransDip_cart

!!Exciton coupling, direct coulomb integral on one dot, inside volume, h1-e/h2-e
real(dp) function D12_in_in(oo,A1,k1,A2,k2,A3,k3,A4,k4,r)

      implicit none
      double precision, external:: s13adf, ei
      integer :: oo, l, k,  i1, i2,m
      real(dp) :: x, y, Integii, inte, inth, ende, endh, angl, A1, A2, A3, A4, k1, k2, k3, k4,&
      eps1,eps2,interval,r,epsR
      real(dp), dimension(3) :: r1, r2

      oo=100

      m=100

      interval=r/m
      
      i1=0
      i2=0
      Integii= 0.0
      x=0
      y=0

      do i1=0,m-1
      x=0.5*interval + i1*interval
      do i2=0,m-1
      y=0.5*interval + i2*interval
      epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-x/rhoe)+exp(-y/rhoh))/2))
      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
      !write(6,*) eps1, eps2, epsR
      do k=1,oo
      r1 = vector(x)
      r2 = vector(y)
      Integii= Integii + (sin(k1*x)*sin(k2*y)*sin(k3*x)*sin(k4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
      enddo
      enddo
      enddo
    
      D12_in_in = (A1*A2*A3*A4*elec*Integii*interval**2)/(4*pi*oo*eps0)

end function D12_in_in

!Exciton coupling, direct coulomb integral on one dot, inside volume, h1-e/h2-e
real(dp) function D12_out_in(oo,B1,k1,A2,k2,B3,k3,A4,k4,r)

      implicit none
      double precision, external:: s13adf, ei
      integer :: oo, l, k,  i1, i2,m
      real(dp) :: x, y, Integii, inte, inth, ende, endh, angl, B1, A2, B3, A4, k1, k2, k3, k4,&
      eps1,eps2,interval, epsR, r
      real(dp), dimension(3) :: r1, r2

      oo=100

      m=100

      interval=r/m
      
      i1=0
      i2=0
      Integii= 0.0
      x=0
      y=0

      do i1=m,2*m
      x=0.5*interval + i1*interval
      do i2=0,m-1
      y=0.5*interval + i2*interval
      epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-x/rhoe)+exp(-y/rhoh))/2))
      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
      !write(6,*) r, x, y, eps1, eps2, epsR
      do k=1,oo
      r1 = vector(x)
      r2 = vector(y)
      Integii= Integii + (exp(-1.0*k1*x)*sin(k2*y)*exp(-1.0*k3*x)*sin(k4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
      enddo
      enddo
      enddo
    
      D12_out_in = (B1*A2*B3*A4*elec*Integii*interval**2)/(4*pi*oo*eps0)

end function D12_out_in

!Exciton coupling, direct coulomb integral on one dot, inside volume, h1-e/h2-e
real(dp) function D12_in_out(oo,A1,k1,B2,k2,A3,k3,B4,k4,r)

      implicit none
      double precision, external:: s13adf, ei
      integer :: oo, l, k,  i1, i2,m
      real(dp) :: x, y,  Integii, inte, inth, ende, endh, angl, A1, B2, A3, B4, k1, k2, k3, k4,&
      eps1,eps2,interval,r,epsR
      real(dp), dimension(3) :: r1, r2

      oo=100

      m=100

      interval=r/m
      
      i1=0
      i2=0
      Integii= 0.0
      x=0
      y=0

      do i1=0,m-1
      x=0.5*interval + i1*interval
      do i2=m,2*m
      y=0.5*interval + i2*interval
      epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-x/rhoe)+exp(-y/rhoh))/2))
      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
      !write(6,*) r, x, y, eps1, eps2, epsR
      do k=1,oo
      r1 = vector(x)
      r2 = vector(y)
      Integii= Integii + (sin(k1*x)*exp(-1.0*k2*y)*sin(k3*x)*exp(-1.0*k4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
      enddo
      enddo
      enddo
    
      D12_in_out = (A1*B2*A3*B4*elec*Integii*interval**2)/(4*pi*oo*eps0)

end function D12_in_out

!Exciton coupling, direct coulomb integral on one dot, inside volume, h1-e/h2-e
real(dp) function D12_out_out(oo,B1,k1,B2,k2,B3,k3,B4,k4,r)

      implicit none
      double precision, external:: s13adf, ei
      integer :: oo, l, k,  i1, i2,m
      real(dp) :: x, y, Integii, inte, inth, ende, endh, angl, B1, B2, B3, B4, k1, k2, k3, k4,&
      eps1,eps2,interval,r,epsR
      real(dp), dimension(3) :: r1, r2

      oo=100

      m=100

      interval=r/m
      
      i1=0
      i2=0
      Integii= 0.0
      x=0
      y=0
 
      do i1=m,2*m
      x=0.5*interval + i1*interval
      do i2=m,2*m
      y=0.5*interval + i2*interval
      epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-x/rhoe)+exp(-y/rhoh))/2))
      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
      do k=1,oo
      r1 = vector(x)
      r2 = vector(y)
      Integii= Integii + (exp(-1.0*k1*x)*exp(-1.0*k2*y)*exp(-1.0*k3*x)*exp(-1.0*k4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
      enddo
      enddo
      enddo
    
      D12_out_out = (B1*B2*B3*B4*elec*Integii*interval**2)/(4*pi*oo*eps0)

end function D12_out_out

!!Exciton coupling, exchange coulomb integral on one dot, inside volume, h1-e/h2-e
!real(dp) function D12ex_in_out(oo,A1,k1,A2,k2,B3,k3,B4,k4,r)
!
!      implicit none
!      double precision, external:: s13adf, ei
!      integer :: oo, l, k,  i1, i2,m
!      real(dp) :: x, y, pi, Integii, inte, inth, ende, endh, angl, A1, A2, B3, B4, k1, k2, k3, k4,&
!      eps,eps0,elec,epsout,epsinfa,eps1,eps2,slope,interval,r,h,hbar,omegaLO,m0,me,mh,rhoe,rhoh, epsR
!      real(dp), dimension(3) :: r1, r2
!
!      pi=  4.D0*DATAN(1.D0)
!      eps=    9.56
!      eps0=   8.854187817620000e-12
!      elec=   1.60217662e-19
!      epsout= 2.5 
!      slope=50e9
!      h=      6.626070040e-34
!      hbar=   h/(2*pi)
!      omegaLO= 5.99585e12
!      m0=  9.10938356e-31
!      me= 0.13*m0
!      mh= 0.82*m0 
!      rhoe  = 1.0/sqrt((2*me*omegaLO)/hbar)
!      rhoh  = 1.0/sqrt((2*mh*omegaLO)/hbar)
!
!      oo=100
!
!      m=100
!
!      interval=r/m
!      
!      i1=0
!      i2=0
!      Integii= 0.0
!      x=0
!      y=0
!
!      epsinfa = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*r))**1.2)
!      !write(6,*) r, epsinfa
!      
!      do i1=0,m-1
!      x=0.5*interval + i1*interval
!      do i2=m,2*m
!      y=0.5*interval + i2*interval
!      epsR = 1.0/((1.0/epsinfa)-((1.0/epsinfa)-(1.0/(epsinfa+3.5)))*(1-(exp(-x/rhoe)+exp(-y/rhoh))/2))
!      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      !write(6,*) r, x, y, eps1, eps2, epsR
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integii= Integii + (sin(k1*x)*sin(k2*x)*exp(-1.0*k3*y)*exp(-1.0*k4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!    
!      D12ex_in_out = (A1*A2*B3*B4*elec*Integii*interval**2)/(4*pi*oo*eps0)
!
!end function D12ex_in_out
!
!!Exciton coupling, exchange coulomb integral on one dot, inside volume, h1-e/h2-e
!real(dp) function D12ex_out_in(oo,B1,k1,B2,k2,A3,k3,A4,k4,r)
!
!      implicit none
!      double precision, external:: s13adf, ei
!      integer :: oo, l, k,  i1, i2,m
!      real(dp) :: x, y, pi, Integii, inte, inth, ende, endh, angl, B1, B2, A3, A4, k1, k2, k3, k4,&
!      eps,eps0,elec,epsout,epsinfa,eps1,eps2,slope,interval,r,h,hbar,omegaLO,m0,me,mh,rhoe,rhoh, epsR
!      real(dp), dimension(3) :: r1, r2
!
!      pi=  4.D0*DATAN(1.D0)
!      eps=    9.56
!      eps0=   8.854187817620000e-12
!      elec=   1.60217662e-19
!      epsout= 2.5 
!      slope=50e9
!      h=      6.626070040e-34
!      hbar=   h/(2*pi)
!      omegaLO= 5.99585e12
!      m0=  9.10938356e-31
!      me= 0.13*m0
!      mh= 0.82*m0 
!      rhoe  = 1.0/sqrt((2*me*omegaLO)/hbar)
!      rhoh  = 1.0/sqrt((2*mh*omegaLO)/hbar)
!
!      oo=100
!
!      m=100
!
!      interval=r/m
!      
!      i1=0
!      i2=0
!      Integii= 0.0
!      x=0
!      y=0
!
!      epsinfa = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*r))**1.2)
!      !write(6,*) r, epsinfa
!      
!      do i1=m,2*m
!      x=0.5*interval + i1*interval
!      do i2=0,m-1
!      y=0.5*interval + i2*interval
!      epsR = 1.0/((1.0/epsinfa)-((1.0/epsinfa)-(1.0/(epsinfa+3.5)))*(1-(exp(-x/rhoe)+exp(-y/rhoh))/2))
!      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      !write(6,*) r, x, y, eps1, eps2, epsR
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integii= Integii + (exp(-1.0*k1*x)*exp(-1.0*k2*x)*sin(k3*y)*sin(k4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!    
!      D12ex_out_in = (B1*B2*A3*A4*elec*Integii*interval**2)/(4*pi*oo*eps0)
!
!end function D12ex_out_in
!
!!Exciton coupling, exchange coulomb integral on one dot, inside volume, h1-e/h2-e
!real(dp) function D12ex_out_out(oo,B1,k1,B2,k2,B3,k3,B4,k4,r)
!
!      implicit none
!      double precision, external:: s13adf, ei
!      integer :: oo, l, k,  i1, i2,m
!      real(dp) :: x, y, pi, Integii, inte, inth, ende, endh, angl, B1, B2, B3, B4, k1, k2, k3, k4,&
!      eps,eps0,elec,epsout,epsinfa,eps1,eps2,slope,interval,r, epsR
!      real(dp), dimension(3) :: r1, r2
!
!      oo=100
!
!      m=100
!
!      interval=r/m
!      
!      i1=0
!      i2=0
!      Integii= 0.0
!      x=0
!      y=0
!
!      epsinfa = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*r))**1.2)
!      !write(6,*) r, epsinfa
!      
!      do i1=m,2*m
!      x=0.5*interval + i1*interval
!      do i2=m,2*m
!      y=0.5*interval + i2*interval
!      epsR = 1.0/((1.0/epsinfa)-((1.0/epsinfa)-(1.0/(epsinfa+3.5)))*(1-(exp(-x/rhoe)+exp(-y/rhoh))/2))
!      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      !write(6,*) r, x, y, eps1, eps2, epsR
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integii= Integii + (exp(-1.0*k1*x)*exp(-1.0*k2*x)*exp(-1.0*k3*y)*exp(-1.0*k4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!    
!      D12ex_out_out = (B1*B2*B3*B4*elec*Integii*interval**2)/(4*pi*oo*eps0)
!
!end function D12ex_out_out

!Exciton coupling exchange same dot 8 loops
real(dp) function D12ex_8loops(oo,A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4,r)

      implicit none
      integer :: oo, l, k,  i1, i2,m
      real(dp) :: x, y, Integii,Integio,Integoi,Integoo, inte, inth, ende, endh, angl,&
      eps1,eps2,interval,r,epsR, A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4
      real(dp), dimension(3) :: r1, r2

      oo=200
      m=100
      interval=r/m
      i1=0
      i2=0
      Integii= 0.0
      Integio= 0.0
      Integoi= 0.0
      Integoo= 0.0
      x=0
      y=0

      do i1=0,m-1
      x=0.5*interval + i1*interval
      do i2=0,m-1
      y=0.5*interval + i2*interval
      do k=1,oo
      r1 = vector(x)
      r2 = vector(y)
      epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-norm2(r1-r2)/rhoe)+exp(-norm2(r1-r2)/rhoh))/2))
      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
      Integii= Integii + (sin(kin1*x)*sin(kin2*x)*sin(kin3*y)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
      enddo
      enddo
      enddo

      do i1=0,m-1
      x=0.5*interval + i1*interval
      eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      do i2=m,2*m
      y=0.5*interval + i2*interval
      eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
      do k=1,oo
      r1 = vector(x)
      r2 = vector(y)
      Integio= Integio + (sin(kin1*x)*sin(kin2*x)*exp(-1.0*kout3*y)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
      enddo
      enddo
      enddo

      do i1=m,2*m
      x=0.5*interval + i1*interval
      eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      do i2=0,m-1
      y=0.5*interval + i2*interval
      eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
      do k=1,oo
      r1 = vector(x)
      r2 = vector(y)
      Integoi= Integoi + (exp(-1.0*kout1*x)*exp(-1.0*kout2*x)*sin(kin3*y)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
      enddo
      enddo
      enddo

      do i1=m,2*m
      x=0.5*interval + i1*interval
      eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      do i2=m,2*m
      y=0.5*interval + i2*interval
      eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
      do k=1,oo
      r1 = vector(x)
      r2 = vector(y)
      Integoo= Integoo + (exp(-1.0*kout1*x)*exp(-1.0*kout2*x)*exp(-1.0*kout3*y)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
      enddo
      enddo
      enddo

      D12ex_8loops = (elec*interval**2)/(4*pi*oo*eps0)*( &
              (A1*A2*A3*A4*Integii) + &
              (A1*A2*B3*B4*Integio) + &
              (B1*B2*A3*A4*Integoi) + &
              (B1*B2*B3*B4*Integoo))

end function D12ex_8loops
!
!!Exciton coupling same dot direct 8 loops, slow
!real(dp) function D12dir(oo,A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4,r)
!
!      implicit none
!      integer :: oo, l, k,  i1, i2,m
!      real(dp) :: x, y, Integii,Integio,Integoi,Integoo, inte, inth, ende, endh, angl,&
!      eps1,eps2,interval,r,epsR, A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4, start, finish
!      real(dp), dimension(3) :: r1, r2
!
!      oo=200
!      m=100
!      interval=r/m
!      i1=0
!      i2=0
!      Integii= 0.0
!      Integio= 0.0
!      Integoi= 0.0
!      Integoo= 0.0
!      x=0
!      y=0
!
!call cpu_time(start)
!
!      do i1=0,m-1
!      x=0.5*interval + i1*interval
!      do i2=0,m-1
!      y=0.5*interval + i2*interval
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-norm2(r1-r2)/rhoe)+exp(-norm2(r1-r2)/rhoh))/2))
!      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      Integii= Integii + (sin(kin1*x)*sin(kin2*y)*sin(kin3*x)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!
!      do i1=0,m-1
!      x=0.5*interval + i1*interval
!      eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      do i2=m,2*m
!      y=0.5*interval + i2*interval
!      eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integio= Integio + (sin(kin1*x)*exp(-kout2*y)*sin(kin3*x)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!
!      do i1=m,2*m
!      x=0.5*interval + i1*interval
!      eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      do i2=0,m-1
!      y=0.5*interval + i2*interval
!      eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integoi= Integoi + (exp(-1.0*kout1*x)*sin(kin2*y)*exp(-1.0*kout3*x)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!
!      do i1=m,2*m
!      x=0.5*interval + i1*interval
!      eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      do i2=m,2*m
!      y=0.5*interval + i2*interval
!      eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integoo= Integoo + (exp(-1.0*kout1*x)*exp(-1.0*kout2*y)*exp(-1.0*kout3*x)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!
!write(6,*) "in_in",  (A1*A2*A3*A4*elec*Integii*interval**2)/(4*pi*oo*eps0)
!write(6,*) "in_out",  (A1*B2*A3*B4*elec*Integio*interval**2)/(4*pi*oo*eps0)
!write(6,*) "out_in",  (B1*A2*B3*A4*elec*Integoi*interval**2)/(4*pi*oo*eps0)
!write(6,*) "out_out",  (B1*B2*B3*B4*elec*Integoo*interval**2)/(4*pi*oo*eps0)
!
!      D12dir = (elec*interval**2)/(4*pi*oo*eps0)*( &
!              (A1*A2*A3*A4*Integii) + &
!              (A1*B2*A3*B4*Integio) + &
!              (B1*A2*B3*A4*Integoi) + &
!              (B1*B2*B3*B4*Integoo))
!
!call cpu_time(finish)
!
!write(6,*) "Time", finish-start
!
!end function D12dir
!
!!Exciton coupling same dot direct 6 loops
!real(dp) function D12dir2(oo,A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4,r)
!
!      implicit none
!      integer :: oo, l, k,  i1, i2,m
!      real(dp) :: x, y, Integii,Integio,Integoi,Integoo, inte, inth, ende, endh, angl,&
!      eps1,eps2,interval,r,epsR, A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4, start, finish
!      real(dp), dimension(3) :: r1, r2
!
!      oo=200
!      m=100
!      interval=r/m
!      i1=0
!      i2=0
!      Integii= 0.0
!      Integio= 0.0
!      Integoi= 0.0
!      Integoo= 0.0
!      x=0
!      y=0
!
!call cpu_time(start)
!
!      do i1=0,m-1
!      x=0.5*interval + i1*interval
!      eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      do i2=0,m-1
!      y=0.5*interval + i2*interval
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-norm2(r1-r2)/rhoe)+exp(-norm2(r1-r2)/rhoh))/2))
!      eps1 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      eps2 = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      Integii= Integii + (sin(kin1*x)*sin(kin2*y)*sin(kin3*x)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      do i2=m,2*m
!      y=0.5*interval + i2*interval
!      eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integio= Integio + (sin(kin1*x)*exp(-kout2*y)*sin(kin3*x)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!
!      do i1=m,2*m
!      x=0.5*interval + i1*interval
!      eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
!      do i2=0,m-1
!      y=0.5*interval + i2*interval
!      eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integoi= Integoi + (exp(-1.0*kout1*x)*sin(kin2*y)*exp(-1.0*kout3*x)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      do i2=m,2*m
!      y=0.5*interval + i2*interval
!      do k=1,oo
!      r1 = vector(x)
!      r2 = vector(y)
!      Integoo= Integoo + (exp(-1.0*kout1*x)*exp(-1.0*kout2*y)*exp(-1.0*kout3*x)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
!      enddo
!      enddo
!      enddo
!
!!write(6,*) "in_in",  (A1*A2*A3*A4*elec*Integii*interval**2)/(4*pi*oo*eps0)
!!write(6,*) "in_out",  (A1*B2*A3*B4*elec*Integio*interval**2)/(4*pi*oo*eps0)
!!write(6,*) "out_in",  (B1*A2*B3*A4*elec*Integoi*interval**2)/(4*pi*oo*eps0)
!!write(6,*) "out_out",  (B1*B2*B3*B4*elec*Integoo*interval**2)/(4*pi*oo*eps0)
!
!      D12dir2 = (elec*interval**2)/(4*pi*oo*eps0)*( &
!              (A1*A2*A3*A4*Integii) + &
!              (A1*B2*A3*B4*Integio) + &
!              (B1*A2*B3*A4*Integoi) + &
!              (B1*B2*B3*B4*Integoo))
!
!call cpu_time(finish)
!
!write(6,*) "Time2", finish-start
!
!end function D12dir2


!Exciton coupling same dot, direct, 2 loops and ifs, fast
real(dp) function DXXdir(A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4,r)

      implicit none
      integer :: oo, i1, i2, m, k
      real(dp) :: x, y, Integii,Integio,Integoi,Integoo, &
      eps1,eps2,eps1R,eps2R,interval,r,epsR, A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4, start, finish
      real(dp), dimension(3) :: r1, r2

      oo=200
      m=100
      interval=r/m
      i1=0
      i2=0
      Integii= 0.0
      Integio= 0.0
      Integoi= 0.0
      Integoo= 0.0
      x=0
      y=0

!call cpu_time(start)

      do i1=0,2*m
      x=0.5*interval + i1*interval
            eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      do i2=0,2*m
      y=0.5*interval + i2*interval
            eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
         if ( (i1 .lt. m) .and. (i2 .lt. m) ) then
            do k=1,oo
            r1 = vector(x)
            r2 = vector(y)
           epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-norm2(r1-r2)/rhoe)+exp(-norm2(r1-r2)/rhoh))/2))
            eps1R = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
            eps2R = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
            Integii= Integii + (sin(kin1*x)*sin(kin2*y)*sin(kin3*x)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1R*eps2R))
            enddo
         elseif ( (i1 .lt. m) .and. (i2 .ge. m) ) then
            do k=1,oo
            r1 = vector(x)
            r2 = vector(y)
            Integio= Integio + (sin(kin1*x)*exp(-kout2*y)*sin(kin3*x)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
            enddo
         elseif ( (i1 .ge. m) .and. (i2 .lt. m) ) then
            do k=1,oo
            r1 = vector(x)
            r2 = vector(y)
            Integoi= Integoi + (exp(-1.0*kout1*x)*sin(kin2*y)*exp(-1.0*kout3*x)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
            enddo
         elseif ( (i1 .ge. m) .and. (i2 .ge. m) ) then
            do k=1,oo
            r1 = vector(x)
            r2 = vector(y)
      Integoo= Integoo + (exp(-1.0*kout1*x)*exp(-1.0*kout2*y)*exp(-1.0*kout3*x)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
            enddo
         endif
      enddo
      enddo

      DXXdir = (elec*interval**2)/(4*pi*oo*eps0)*( &
              (A1*A2*A3*A4*Integii) + &
              (A1*B2*A3*B4*Integio) + &
              (B1*A2*B3*A4*Integoi) + &
              (B1*B2*B3*B4*Integoo))

!call cpu_time(finish)

!write(6,*) "Time", finish-start

end function DXXdir

!Exciton coupling same dot, exchange, 2 loops and ifs, fast
real(dp) function DXXex(A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4,r)

      implicit none
      integer :: oo, i1, i2, m, k
      real(dp) :: x, y, Integii,Integio,Integoi,Integoo, &
      eps1,eps2,eps1R,eps2R,interval,r,epsR, A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4, start, finish
      real(dp), dimension(3) :: r1, r2

      oo=200
      m=100
      interval=r/m
      i1=0
      i2=0
      Integii= 0.0
      Integio= 0.0
      Integoi= 0.0
      Integoo= 0.0
      x=0
      y=0

!call cpu_time(start)

      do i1=0,2*m
      x=0.5*interval + i1*interval
            eps1 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(x-r)))/pi
      do i2=0,2*m
      y=0.5*interval + i2*interval
            eps2 = epsout + (epsin(n)-epsout)*((pi/2) - atan(slope*(y-r)))/pi
         if ( (i1 .lt. m) .and. (i2 .lt. m) ) then
            do k=1,oo
            r1 = vector(x)
            r2 = vector(y)
           epsR = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-norm2(r1-r2)/rhoe)+exp(-norm2(r1-r2)/rhoh))/2))
            eps1R = epsout + (epsR-epsout)*((pi/2) - atan(slope*(x-r)))/pi
            eps2R = epsout + (epsR-epsout)*((pi/2) - atan(slope*(y-r)))/pi
            Integii= Integii + (sin(kin1*x)*sin(kin2*x)*sin(kin3*y)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1R*eps2R))
            enddo
         elseif ( (i1 .lt. m) .and. (i2 .ge. m) ) then
            do k=1,oo
            r1 = vector(x)
            r2 = vector(y)
            Integio= Integio + (sin(kin1*x)*sin(kin2*x)*exp(-1.0*kout3*y)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
            enddo
         elseif ( (i1 .ge. m) .and. (i2 .lt. m) ) then
            do k=1,oo
            r1 = vector(x)
            r2 = vector(y)
            Integoi= Integoi + (exp(-1.0*kout1*x)*exp(-1.0*kout2*x)*sin(kin3*y)*sin(kin4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
            enddo
         elseif ( (i1 .ge. m) .and. (i2 .ge. m) ) then
            do k=1,oo
            r1 = vector(x)
            r2 = vector(y)
      Integoo= Integoo + (exp(-1.0*kout1*x)*exp(-1.0*kout2*x)*exp(-1.0*kout3*y)*exp(-1.0*kout4*y))/(norm2(r1-r2)*sqrt(eps1*eps2))
            enddo
         endif
      enddo
      enddo

      DXXex = (elec*interval**2)/(4*pi*oo*eps0)*( &
              (A1*A2*A3*A4*Integii) + &
              (A1*A2*B3*B4*Integio) + &
              (B1*B2*A3*A4*Integoi) + &
              (B1*B2*B3*B4*Integoo))

!call cpu_time(finish)

!write(6,*) "Time", finish-start

end function DXXex


end module Integrals
