module Integrals_Dmax

      use Constants
      use Variables

contains

real*8 function D_MCMCoff(A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4,a,b,l,flag)

      implicit none

      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
      character*2 :: flag
      character*4 :: inA_inA, inA_inB, inA_out, inB_inA, inB_inB, inB_out, out_inA, out_inB, out_out
      integer :: m, ii, i, i1,i2,i3,i4, j, j1,j2,j3,j4,j5,j6,j7,j8,j9, num
      real(dp) :: f1, f2, f3, f4, f5, f6, f7, f8, f9, interval, f, &
                f1mc, f2mc, f3mc, f1mcmc, f2mcmc, f3mcmc, add, maxrad, vola, volb, volouta, voloutb, volout, start, finish, &
                epsR, epsinf, epsinfa, epsinfb, eps1, eps2, side, a, b, l
      real(dp) :: A1,B1,kin1,kout1,A2,B2,kin2,kout2,A3,B3,kin3,kout3,A4,B4,kin4,kout4
      real(dp), allocatable:: r1(:,:) , r2(:,:) , ra(:), rb(:), rc(:), rd(:), &
                              r1Anorm(:), r2Anorm(:), rnorm(:), RAB(:), r1Bnorm(:), r2Bnorm(:)

      include '../inputs/GetVar.f90'

      m=100000

      allocate(RAB(3))
      allocate(r1(m,3))
      allocate(r1Anorm(m))
      allocate(r1Bnorm(m))
      allocate(r2(m,3))
      allocate(r2Bnorm(m))
      allocate(r2Anorm(m))
      allocate(rnorm(m))
      allocate(ra(m))
      allocate(rb(m))
      allocate(rc(m))
      allocate(rd(m))

      RAB(1)=a+b+l
      RAB(2)=0
      RAB(3)=0

      i=0
      j=0
      maxrad=max(a,b)
      i1 = 0
      i2 = 0
      i3 = 0
      f1mc=0.0
      f2mc=0.0
      f3mc=0.0
      f1mcmc=0.0
      f2mcmc=0.0
      f3mcmc=0.0
      vola=(4.0/3)*PI*a**3
      volb=(4.0/3)*PI*b**3
      volout=abs(((2*maxrad+2*side)**2*(2*a+l+2*b+2*side))-(vola+volb))

     call random_seed()
     call random_number(r1)

r1Anorm(:)=sqrt(((2*a+2*b+l+2*side)*r1(:,1)-(a+side))**2+((2*maxrad+2*side)*r1(:,2)-(maxrad+side))**2+ &
               ((2*maxrad+2*side)*r1(:,3)-(maxrad+side))**2)
r1Bnorm(:)=sqrt(((2*a+2*b+l+2*side)*r1(:,1)-(RAB(1)+a+side))**2+((2*maxrad+2*side)*r1(:,2)-(maxrad+side))**2+&
                ((2*maxrad+2*side)*r1(:,3)-(maxrad+side))**2)

epsinfa = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*a))**1.2)
epsinfb = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*b))**1.2)

     do i=1,m

     call random_number(r2)

     r2Anorm(:)=sqrt(((2*a+2*b+l+2*side)*r2(:,1)-(a+side))**2+((2*maxrad+2*side)*r2(:,2)-(maxrad+side))**2+&
                   ((2*maxrad+2*side)*r2(:,3)-(maxrad+side))**2)
     r2Bnorm(:)=sqrt(((2*a+2*b+l+2*side)*r2(:,1)-(RAB(1)+a+side))**2+((2*maxrad+2*side)*r2(:,2)-(maxrad+side))**2+&
                   ((2*maxrad+2*side)*r2(:,3)-(maxrad+side))**2)
     rnorm(:)=sqrt(((2*a+2*b+l+2*side)*(r2(:,1)-r1(i,1)))**2+((2*maxrad+2*side)*(r2(:,2)-r1(i,2)))**2+&
                  ((2*maxrad+2*side)*(r2(:,3)-r1(i,3)))**2)

      j=0
      j1=0
      j2=0
      j3=0
      j4=0
      j5=0
      j6=0
      j7=0
      j8=0
      j9=0
      f1=0.0
      f2=0.0
      f3=0.0
      f4=0.0
      f5=0.0
      f6=0.0
      f7=0.0
      f8=0.0
      f9=0.0

      if ((r1Anorm(i) .le. a) .and. (r1Bnorm(i) .gt. b)) then !inA
               i1=i1+1
         do j=1,m

               include '../inputs/GetVec.f90'

               if (rnorm(j) .lt. 1e-12 ) then
               exit
               else if ((r2Bnorm(j) .le. b)) then !inA-inB
                  eps1 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r1Anorm(i)-a)))/PI
                  eps2 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r2Bnorm(j)-b)))/PI
                  f1 = f1 + wf(inA_inB,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j1=j1+1
              else if ((r2Bnorm(j) .gt. b) .and. (r2Anorm(j) .lt. a)) then !inA-inA
                  epsR = 1.0/((1.0/epsinfa)-((1.0/epsinfa)-(1.0/(epsinfa+3.5)))*(1-(exp(-rnorm(j)/rhoe)+exp(-rnorm(j)/rhoh))/2))
                  eps1 = epsout + (epsR-epsout)*((PI/2) - atan(slope*(r1Anorm(i)-a)))/PI
                  eps2 = epsout + (epsR-epsout)*((PI/2) - atan(slope*(r2Anorm(j)-a)))/PI
                  f2 = f2 + wf(inA_inA,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j2=j2+1
              else if ((r2Anorm(j) .gt. a) .and. (r2Bnorm(j) .gt. b)) then !inA-out
                  eps1 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r1Anorm(i)-a)))/PI
                  if ((r2Anorm(j)-a) .lt. (r2Bnorm(j)-b)) then
                     eps2 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r2Anorm(j)-a)))/PI 
                     else
                     eps2 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r2Bnorm(j)-b)))/PI
                  endif
                  f3 = f3 + wf(inA_out,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j3=j3+1
              endif
         enddo

         f1mc = f1mc + (f1/j1)*nor(inA_inB,A1,B1,A2,B2,A3,B3,A4,B4)*volb +&
                       (f2/j2)*nor(inA_inA,A1,B1,A2,B2,A3,B3,A4,B4)*vola +&
                       (f3/j3)*nor(inA_out,A1,B1,A2,B2,A3,B3,A4,B4)*volout

      else if ((r1Anorm(i) .gt. a) .and. (r1Bnorm(i) .lt. b)) then !inB
               i2=i2+1
         do j=1,m

               include '../inputs/GetVec.f90'

               if (rnorm(j) .lt. 1e-12 ) then
               exit
               else if ((r2Bnorm(j) .le. b)) then !inB-inB
                  epsR = 1.0/((1.0/epsinfb)-((1.0/epsinfb)-(1.0/(epsinfb+3.5)))*(1-(exp(-rnorm(j)/rhoe)+exp(-rnorm(j)/rhoh))/2))
                  eps1 = epsout + (epsR-epsout)*((PI/2) - atan(slope*(r1Bnorm(i)-b)))/PI
                  eps2 = epsout + (epsR-epsout)*((PI/2) - atan(slope*(r2Bnorm(j)-b)))/PI
                  f4 = f4 + wf(inB_inB,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j4=j4+1
              else if ((r2Bnorm(j) .gt. b) .and. (r2Anorm(j) .lt. a)) then !inB-inA
                  eps1 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r1Bnorm(i)-b)))/PI
                  eps2 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r2Anorm(j)-a)))/PI
                  f5 = f5 + wf(inB_inA,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j5=j5+1
              else if ((r2Anorm(j) .gt. a) .and. (r2Bnorm(j) .gt. b)) then !inB-out
                  eps1 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r1Bnorm(i))-b))/PI
                  if ((r2Anorm(j)-a) .lt. (r2Bnorm(j)-b)) then
                     eps2 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r2Anorm(j)-a)))/PI 
                     else
                     eps2 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r2Bnorm(j)-b)))/PI
                  endif
                  f6 = f6 + wf(inB_out,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j6=j6+1
              endif

          enddo

           f2mc = f2mc + (f4/j4)*nor(inB_inB,A1,B1,A2,B2,A3,B3,A4,B4)*volb +&
                         (f5/j5)*nor(inB_inA,A1,B1,A2,B2,A3,B3,A4,B4)*vola +&
                         (f6/j6)*nor(inB_out,A1,B1,A2,B2,A3,B3,A4,B4)*volout

      else if ((r1Anorm(i) .gt. a) .and. (r1Bnorm(i) .gt. b)) then !out
               i3=i3+1
         do j=1,m

               include '../inputs/GetVec.f90'

               if (rnorm(j) .lt. 1e-12 ) then
               exit
               else if ((r2Bnorm(j) .le. b)) then !out-inB
                  if ((r1Anorm(i)-a) .lt. (r1Bnorm(i)-b)) then
                     eps1 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r1Anorm(i)-a)))/PI 
                     else
                     eps1 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r1Bnorm(i)-b)))/PI
                  endif
                  eps2 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r2Bnorm(j)-b)))/PI
                  f7 = f7 + wf(out_inB,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j7=j7+1
              else if ((r2Bnorm(j) .gt. b) .and. (r2Anorm(j) .lt. a)) then !out-inA
                  if ((r1Anorm(i)-a) .lt. (r1Bnorm(i)-b)) then
                     eps1 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r1Anorm(i)-a)))/PI 
                     else
                     eps1 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r1Bnorm(i)-b)))/PI
                  endif
                  eps2 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(min(r2Anorm(j),r2Bnorm(j))-b)))/PI
                  f8 = f8 + wf(out_inA,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j8=j8+1
              else if ((r2Anorm(j) .gt. a) .and. (r2Bnorm(j) .gt. b)) then !out-out
                  if ((r1Anorm(i)-a) .lt. (r1Bnorm(i)-b)) then
                     eps1 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r1Anorm(i)-a)))/PI 
                     else
                     eps1 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r1Bnorm(i)-b)))/PI
                  endif
                  if ((r2Anorm(j)-a) .lt. (r2Bnorm(j)-b)) then
                     eps2 = epsout + (epsinfa-epsout)*((PI/2) - atan(slope*(r2Anorm(j)-a)))/PI 
                     else
                     eps2 = epsout + (epsinfb-epsout)*((PI/2) - atan(slope*(r2Bnorm(j)-b)))/PI
                  endif
                  f9 = f9 + wf(out_out,kin1,kout1,ra(i),kin2,kout2,rb(j),kin3,kout3,rc(i),kin4,kout4,rd(j)) /&
                            (rnorm(j)*sqrt(eps1*eps2))
                  j9=j9+1
              endif

          enddo

           f3mc = f3mc + (f7/j7)*nor(out_inB,A1,B1,A2,B2,A3,B3,A4,B4)*volb +&
                         (f8/j8)*nor(out_inA,A1,B1,A2,B2,A3,B3,A4,B4)*vola +&
                         (f9/j9)*nor(out_out,A1,B1,A2,B2,A3,B3,A4,B4)*volout

!write(6,*) f7, f8, f9


      endif


     enddo

           f1mcmc = (f1mc/i1)*vola
           f2mcmc = (f2mc/i2)*volb
           f3mcmc = (f3mc/i3)*volout

!call cpu_time(finish)

!write(6,*) finish-start

deallocate(r2Anorm,rnorm,r1Anorm,r2Bnorm, r1Bnorm)

!enddo

      D_MCMCoff=(elec/(64*PI**3*eps0))*(f1mcmc+f2mcmc+f3mcmc)

end function D_MCMCoff

real*8 function wf(integral,ki1,ko1,ra,ki2,ko2,rb,ki3,ko3,rc,ki4,ko4,rd)
      implicit none
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
      character*4 :: integral
      real(dp) :: ki1,ko1,ra,ki2,ko2,rb,ki3,ko3,rc,ki4,ko4,rd
! 4in
if (integral .eq. 'iiii') then
wf=sin(ki1*ra)*sin(ki2*rb)*sin(ki3*rc)*sin(ki4*rd)/(ra*rb*rc*rd)
! 3in
else if (integral .eq. 'iiio') then
wf=sin(ki1*ra)*sin(ki2*rb)*sin(ki3*rc)*exp(-1.0*ko4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'iioi') then
wf=sin(ki1*ra)*sin(ki2*rb)*exp(-1.0*ko3*rc)*sin(ki4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'ioii') then
wf=sin(ki1*ra)*exp(-1.0*ko2*rb)*sin(ki3*rc)*sin(ki4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'oiii') then
wf=exp(-1.0*ko1*ra)*sin(ki2*rb)*sin(ki3*rc)*sin(ki4*rd)/(ra*rb*rc*rd)
! 2in
else if (integral .eq. 'iioo') then
wf=sin(ki1*ra)*sin(ki2*rb)*exp(-1.0*ko3*rc)*exp(-1.0*ko4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'ooii') then
wf=exp(-1.0*ko1*ra)*exp(-1.0*ko2*rb)*sin(ki3*rc)*sin(ki4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'iooi') then
wf=sin(ki1*ra)*exp(-1.0*ko2*rb)*exp(-1.0*ko3*rc)*sin(ki4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'oiio') then
wf=exp(-1.0*ko1*ra)*sin(ki2*rb)*sin(ki3*rc)*exp(-1.0*ko4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'oioi') then
wf=exp(-1.0*ko1*ra)*sin(ki2*rb)*exp(-1.0*ko3*rc)*sin(ki4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'ioio') then
wf=sin(ki1*ra)*exp(-1.0*ko2*rb)*sin(ki3*rc)*exp(-1.0*ko4*rd)/(ra*rb*rc*rd)
! 1in
else if (integral .eq. 'iooo') then
wf=sin(ki1*ra)*exp(-1.0*ko2*rb)*exp(-1.0*ko3*rc)*exp(-1.0*ko4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'oioo') then
wf=exp(-1.0*ko1*ra)*sin(ki2*rb)*exp(-1.0*ko3*rc)*exp(-1.0*ko4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'ooio') then
wf=exp(-1.0*ko1*ra)*exp(-1.0*ko2*rb)*sin(ki3*rc)*exp(-1.0*ko4*rd)/(ra*rb*rc*rd)
else if (integral .eq. 'oooi') then
wf=exp(-1.0*ko1*ra)*exp(-1.0*ko2*rb)*exp(-1.0*ko3*rc)*sin(ki4*rd)/(ra*rb*rc*rd)
! 4in
else if (integral .eq. 'oooo') then
wf=exp(-1.0*ko1*ra)*exp(-1.0*ko2*rb)*exp(-1.0*ko3*rc)*exp(-1.0*ko4*rd)/(ra*rb*rc*rd)
endif
end function wf

real*8 function nor(integral,A1,B1,A2,B2,A3,B3,A4,B4)
      implicit none
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
      character*4 :: integral
      real(dp) :: A1,B1,A2,B2,A3,B3,A4,B4
! 4in
if (integral .eq. 'iiii') then
nor=A1*A2*A3*A4
! 3in
else if (integral .eq. 'iiio') then
nor=A1*A2*A3*B4
else if (integral .eq. 'iioi') then
nor=A1*A2*B3*A4
else if (integral .eq. 'ioii') then
nor=A1*B2*A3*A4
else if (integral .eq. 'oiii') then
nor=B1*A2*A3*A4
! 2in
else if (integral .eq. 'iioo') then
nor=A1*A2*B3*B4
else if (integral .eq. 'ooii') then
nor=B1*B2*A3*A4
else if (integral .eq. 'iooi') then
nor=A1*B2*B3*A4
else if (integral .eq. 'oiio') then
nor=B1*A2*A3*B4
else if (integral .eq. 'oioi') then
nor=B1*A2*B3*A4
else if (integral .eq. 'ioio') then
nor=A1*B2*A3*B4
! 1in
else if (integral .eq. 'iooo') then
nor=A1*B2*B3*B4
else if (integral .eq. 'oioo') then
nor=B1*A2*B3*B4
else if (integral .eq. 'ooio') then
nor=B1*B2*A3*B4
else if (integral .eq. 'oooi') then
nor=B1*B2*B3*A4
! 4in
else if (integral .eq. 'oooo') then
nor=B1*B2*B3*B4
endif
end function nor

end module Integrals_Dmax
