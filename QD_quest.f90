include 'specfun.f90'

program ModelQD_one

use Constants
use Variables
use Integrals
use Vectors
use Output
use Make_Ham

implicit none

double precision, external:: s13adf, ei, eone, nag_bessel_j0

character*64 :: arg1, arg2, filename, argA

integer :: je,jh,k,nsteps,r,ifail, r1, r2

real(dp) :: Ef,le,re,lh,rh,delta, start, finish, mu, A, epsinf, aR1, aR2
real(dp) :: Rine, Route, Rinh1, Routh1, Rinh2, Routh2, RadProbe, RadProbh1, RadProbh2, r0
real(dp),allocatable :: Ae(:), Ah1(:), Ah2(:), Be(:), Bh1(:), Bh2(:)
real(dp),allocatable :: I1eh1(:), I1eh2(:), I2eh1(:), I2eh2(:), I3eh1(:), I3eh2(:), kine(:), kinh1(:), kinh2(:)
real(dp),allocatable :: koute(:), kouth1(:), kouth2(:),diffe(:), diffh(:), E(:)

call getVariables

delta=  0.00001d-18
Ef=     1.28174d-18
nsteps= int(Ef/delta)
ifail=  1

allocate(E(nsteps))
allocate(diffe(nsteps))
allocate(diffh(nsteps))
allocate(minEe(nsteps,rmax+1))
allocate(minEh(nsteps,rmax+1))
allocate(Eeh1(rmax+1)) 
allocate(Eeh2(rmax+1)) 
allocate(Ae(rmax+1)) 
allocate(Ah1(rmax+1)) 
allocate(Ah2(rmax+1)) 
allocate(Be(rmax+1)) 
allocate(Bh1(rmax+1)) 
allocate(Bh2(rmax+1)) 
allocate(Cb_eh1(rmax+1)) 
allocate(Cb_eh2(rmax+1))
allocate(I1eh1(rmax+1)) 
allocate(I1eh2(rmax+1)) 
allocate(I2eh1(rmax+1)) 
allocate(I2eh2(rmax+1)) 
allocate(I3eh1(rmax+1))
allocate(I3eh2(rmax+1)) 
allocate(kine(rmax+1)) 
allocate(kinh1(rmax+1)) 
allocate(kinh2(rmax+1)) 
allocate(koute(rmax+1)) 
allocate(kouth1(rmax+1)) 
allocate(kouth2(rmax+1))
allocate(Norm_Ana_e(rmax+1))
allocate(Norm_Ana_h1(rmax+1))
allocate(Norm_Ana_h2(rmax+1))
allocate(OverlapAna_h1e(rmax+1))
allocate(OverlapAna_h2e(rmax+1))
allocate(Cb_Num_eh1(rmax+1))
allocate(Cb_Num_eh1_eh2(rmax+1))
allocate(Cb_Num_eh2(rmax+1))
allocate(TransDip_Ana_h1e(rmax+1))
allocate(TransDip_Ana_h2e(rmax+1))
allocate(Oscillator_Ana_h1e(rmax+1))
allocate(Oscillator_Ana_h2e(rmax+1))
allocate(ExctCoef_h1e(rmax+1))
allocate(ExctCoef_h2e(rmax+1))
!allocate((rmax+1))
!allocate((rmax+1))
!allocate((rmax+1))
!allocate((rmax+1))
!allocate((rmax+1))

open(11,file='Eh1e.dat')
open(12,file='Eh2e.dat')
open(13,file='wavefunctionA.dat')
open(14,file='wavefunctionB.dat')
open(15,file='radialdisA.dat')
open(16,file='radialdisB.dat')
open(17,file='TransDip-local.dat')
open(18,file='TransDip-CT.dat')

k=1

!Computation of energies
n=0

do n = rmin,rmax

i=0
r=0

je=1
jh=1

do i=1,nsteps
E(i)=delta*i
diffe(i) = abs(sqrt(2*me*E(i))/hbar * aR(n) * 1/tan(sqrt(2*me*E(i))/hbar * aR(n)) - 1 + (me/m0) + (me*aR(n))/(hbar) &
           * sqrt((2/m0)*(V0h(n)-E(i))))
if ((diffe(0) .eq. 0.000) .and. (diffh(0) .eq. 0.00)) then 
        diffe(0)=diffe(i)
        diffh(0)=diffe(i)
endif

if (diffe(i) .le. diffe(i-1)) then
        minEe(je,n) = E(i)
elseif ( (diffe(i) .ge. diffe(i-1)) .and. (E(i-1) .eq. minEe(je,n)) ) then
        je=je+1
endif

diffh(i) = abs(sqrt(2*mh*E(i))/hbar * aR(n) * 1/tan(sqrt(2*mh*E(i))/hbar * aR(n)) - 1 + (mh/m0) + (mh*aR(n))/(hbar) &
           * sqrt((2/m0)*(V0h(n)-E(i))))
if (diffh(i) .le. diffh(i-1)) then
        minEh(jh,n) = E(i)
elseif ( (diffh(i) .ge. diffh(i-1)) .and. (E(i-1) .eq. minEh(jh,n)) ) then
        jh=jh+1
endif
enddo

!wave vectors in and out
kine(n)=sqrt(2*me*minEe(1,n))/hbar
koute(n)=sqrt(2*m0*(V0h(n)-minEe(1,n)))/hbar

kinh1(n)=sqrt(2*mh*minEh(1,n))/hbar
kouth1(n)=sqrt(2*m0*(V0h(n)-minEh(1,n)))/hbar

kinh2(n)=sqrt(2*mh*minEh(2,n))/hbar
kouth2(n)=sqrt(2*m0*(V0h(n)-minEh(2,n)))/hbar

!normalization factors
Ae(n)=1/sqrt(aR(n)/2-sin(2*kine(n)*aR(n))/(4*kine(n))+sin(kine(n)*aR(n))**2/(2*koute(n)))
Be(n)=Ae(n)*sin(kine(n)*aR(n))*exp(koute(n)*aR(n)) 

Ah1(n)=1/sqrt(aR(n)/2-sin(2*kinh1(n)*aR(n))/(4*kinh1(n))+sin(kinh1(n)*aR(n))**2/(2*kouth1(n)))
Bh1(n)=Ah1(n)*sin(kinh1(n)*aR(n))*exp(kouth1(n)*aR(n)) 

Ah2(n)=1/sqrt(aR(n)/2-sin(2*kinh2(n)*aR(n))/(4*kinh2(n))+sin(kinh2(n)*aR(n))**2/(2*kouth2(n)))
Bh2(n)=Ah2(n)*sin(kinh2(n)*aR(n))*exp(kouth2(n)*aR(n))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Comlomb correction (analytical) between electron and holes

I1eh1(n)=(Ae(n)**2*Ah1(n)**2*aR(n)/4)*(1-sin(2*kine(n)*aR(n))/(2*kine(n)*aR(n))-s13adf(2*kine(n)*aR(n),ifail)/(2*kine(n)*aR(n))+&
         (s13adf(2*kine(n)*aR(n)-2*kinh1(n)*aR(n),ifail)+s13adf(2*kine(n)*aR(n)+2*kinh1(n)*aR(n),ifail))/(4*kine(n)*aR(n))) + &
        (Ae(n)**2*Ah1(n)**2*aR(n)/4)*(1-sin(2*kinh1(n)*aR(n))/(2*kinh1(n)*aR(n))-s13adf(2*kinh1(n)*aR(n),ifail)/(2*kinh1(n)*aR(n))+&
         (s13adf(2*kinh1(n)*aR(n)-2*kine(n)*aR(n),ifail)+s13adf(2*kinh1(n)*aR(n)+2*kine(n)*aR(n),ifail))/(4*kinh1(n)*aR(n))) 

I2eh1(n)=(Ae(n)**2*Bh1(n)**2*aR(n)/2)*(1-sin(2*kine(n)*aR(n))/(2*kine(n)*aR(n)))*eone(2*kouth1(n)*aR(n))+&
         (Ah1(n)**2*Be(n)**2*aR(n)/2)*(1-sin(2*kinh1(n)*aR(n))/(2*kinh1(n)*aR(n)))*eone(2*koute(n)*aR(n))

I3eh1(n)=(Be(n)**2*Bh1(n)**2*aR(n)/(2*koute(n)*aR(n))*(exp(-2*koute(n)*aR(n))*& 
      eone(2*kouth1(n)*aR(n))-eone(2*koute(n)*aR(n)+2*kouth1(n)*aR(n))))+&
      (Be(n)**2*Bh1(n)**2*aR(n)/(2*kouth1(n)*aR(n))*(exp(-2*kouth1(n)*aR(n))*&
      eone(2*koute(n)*aR(n))-eone(2*kouth1(n)*aR(n)+2*koute(n)*aR(n))))

I1eh2(n)=(Ae(n)**2*Ah2(n)**2*aR(n)/4)*(1-sin(2*kine(n)*aR(n))/(2*kine(n)*aR(n))-s13adf(2*kine(n)*aR(n),ifail)/(2*kine(n)*aR(n))+&
         (s13adf(2*kine(n)*aR(n)-2*kinh2(n)*aR(n),ifail)+s13adf(2*kine(n)*aR(n)+2*kinh2(n)*aR(n),ifail))/(4*kine(n)*aR(n))) + &
        (Ae(n)**2*Ah2(n)**2*aR(n)/4)*(1-sin(2*kinh2(n)*aR(n))/(2*kinh2(n)*aR(n))-s13adf(2*kinh2(n)*aR(n),ifail)/(2*kinh1(n)*aR(n))+&
         (s13adf(2*kinh2(n)*aR(n)-2*kine(n)*aR(n),ifail)+s13adf(2*kinh2(n)*aR(n)+2*kine(n)*aR(n),ifail))/(4*kinh2(n)*aR(n))) 

I2eh2(n)=(Ae(n)**2*Bh2(n)**2*aR(n)/2)*(1-sin(2*kine(n)*aR(n))/(2*kine(n)*aR(n)))*eone(2*kouth2(n)*aR(n))+&
         (Ah2(n)**2*Be(n)**2*aR(n)/2)*(1-sin(2*kinh2(n)*aR(n))/(2*kinh2(n)*aR(n)))*eone(2*koute(n)*aR(n))

I3eh2(n)=(Be(n)**2*Bh2(n)**2*aR(n)/(2*koute(n)*aR(n))*(exp(-2*koute(n)*aR(n))*&
       eone(2*kouth2(n)*aR(n))-eone(2*koute(n)*aR(n)+2*kouth2(n)*aR(n))))+&
         (Be(n)**2*Bh2(n)**2*aR(n)/(2*kouth2(n)*aR(n))*(exp(-2*kouth2(n)*aR(n))*&
       eone(2*koute(n)*aR(n))-eone(2*kouth2(n)*aR(n)+2*koute(n)*aR(n))))


Cb_eh1(n)=(elec**2/(4*pi*epsin(n)*eps0))*(I1eh1(n)+I2eh1(n)+I3eh1(n))

Cb_eh2(n)=(elec**2/(4*pi*epsin(n)*eps0))*(I1eh2(n)+I2eh2(n)+I3eh2(n))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Eeh1(n) = (minEe(1,n)+minEh(1,n))+V0-Cb_eh1(n) 
Eeh2(n) = (minEe(1,n)+minEh(2,n))+V0-Cb_eh2(n)

!write(6,*)  aR(n), "Egap e-h1", Eeh1(n)/elec, "including correction", Cb_eh1(n)/elec
!write(6,*) "Egap e-h2", Eeh2(n)/elec, "including correction", Cb_eh2(n)/elec 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Dipole moment 
if ( o_DipS == 'y' ) then
TransDip_Ana_h1e(n) = abs(TransDip_Ana(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))
!TransDip_Num_h1e(n) = abs(TransDip_Num(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))/Cm_to_D 
!TransDip_EMA_h1e(n) = abs(TransDip_EMA(Eeh1(n),aR(n)))/Cm_to_D
TransDip_Ana_h2e(n) = abs(TransDip_Ana(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))
!TransDip_Num_h2e(n) = abs(TransDip_Num(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))/Cm_to_D
!TransDip_EMA_h2e(n) = abs(TransDip_EMA(Eeh2(n),aR(n)))/Cm_to_D
endif

if ( ( vers .eq. 'singl' ) .or. ( vers .eq. 'dimer') ) then 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Wave function QDA and radial probability
!do r=1,rsteps
!
!r0=r*stepr0
!
!Rine=Ae(n)*sin(kine(n)*r0)/(sqrt(4*pi)*r0)
!Route=Be(n)*exp(-1*koute(n)*r0)/(sqrt(4*pi)*r0)
!
!Rinh1=Ah1(n)*sin(kinh1(n)*r0)/(sqrt(4*pi)*r0)
!Routh1=Bh1(n)*exp(-1*kouth1(n)*r0)/(sqrt(4*pi)*r0)
!
!Rinh2=Ah2(n)*sin(kinh2(n)*r0)/(sqrt(4*pi)*r0)
!Routh2=Bh2(n)*exp(-1*kouth2(n)*r0)/(sqrt(4*pi)*r0)
!
!if (r0 .le. aR(n)) then
!        RadProbe=Rine
!        RadProbh1=Rinh1
!        RadProbh2=Rinh2
!else if (r0 .gt. aR(n)) then
!        RadProbe=Route
!        RadProbh1=Routh1
!        RadProbh2=Routh2
!endif
!
!write(17,*) r0, RadProbe, RadProbh1, RadProbh2
!write(15,*) r0, r0**2*RadProbe**2, r0**2*RadProbh1**2, r0**2*RadProbh2**2
!
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Normalization integrals
if ( o_Norm == 'y' ) then
Norm_Ana_e(n) = Norm_Ana(aR(n),Ae(n),Be(n),kine(n),koute(n))
!Norm_Num_e(n) = Norm_Num(aR(n),Ae(n),Be(n),kine(n),koute(n))
Norm_Ana_h1(n) = Norm_Ana(aR(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n))
!Norm_Num_e(n) = Norm_Num(aR(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n))
Norm_Ana_h2(n) = Norm_Ana(aR(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n))
!Norm_Num_e(n) = Norm_Num(aR(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n))
endif

!Overlap integrals
if ( o_Over == 'y' ) then
OverlapAna_h1e(n) = abs(OverlapAna(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))
!OverlapNum_h1e(n) = abs(OverlapNum(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))
OverlapAna_h2e(n) = abs(OverlapAna(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))
!OverlapNum_h2e(n) = abs(OverlapNum(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))
endif

!Coulomb correction
if ( o_Coul == 'y' ) then
Cb_Num_eh1(n) = DXXdir(Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
                                       Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),aR(n)) + &
                                DXXex(Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
                                      Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),aR(n))
Cb_Num_eh1_eh2(n) =  DXXdir(Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
                                       Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n)) + &
                                DXXex(Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
                                      Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n))
Cb_Num_eh2(n) = DXXdir(Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),&
                                       Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n)) + &
                                DXXex(Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),&
                                      Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n))
endif

!Oscilator strength
if ( o_Osci == 'y' ) then
Oscillator_Ana_h1e(n) = ((2*me)/(3*elec**2*hbar**2))*Eeh1(n)*&
                                       TransDip_Ana(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n))**2
Oscillator_Ana_h2e(n) =  ((2*me)/(3*elec**2*hbar**2))*Eeh2(n)*&
                                       TransDip_Ana(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n))**2
endif

!molar extinction coefficient
if ( o_Exti == 'y' ) then
A=1.25d5*elec**2
mu=(me*mh)/(me+mh)
ExctCoef_h1e(n) = (A*mu*aR(n)**2)/(sqrt(2*pi)*Eeh1(n)*hbar**2*0.1)
ExctCoef_h2e(n) = (A*mu*aR(n)**2)/(sqrt(2*pi)*Eeh2(n)*hbar**2*0.1)
endif

endif

enddo

!Dipole moment h1-e dimer
!write(6,*) aR(n), linker, TransDip_dimer_MC(Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
!                                           Ae(n),Be(n),kine(n),koute(n),aR(n),aR(n),linker)

!write(6,*) aR(n), linker, TransDip_dimer_cart(m,Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
!                                              Ae(n),Be(n),kine(n),koute(n),aR(n),aR(n),linker)
!                 TransDip_dimer_cart(m,Ah2(n),Bh2(n),kinh2(n),kouth2(n),&
!                                              Ae(n),Be(n),kine(n),koute(n),aR(n),aR(n),linker)
!                    TransDip_dimer_cart(m,Ae(n),Be(n),kine(n),koute(n),&
!                                              Ah1(n),Bh1(n),kinh1(n),kouth1(n),aR(n),aR(n),linker),&
!                    TransDip_dimer_cart(m,Ae(n),Be(n),kine(n),koute(n),&
!                                              Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n),aR(n),linker)


!if ( o_DipD == 'y' ) then
!do r1=rmin,rmax
!do r1=rmin,rmax
!
!!do r2=rmin,rmax
!
!do r2=rmin,rmax
!
!!linker = r2*rsteps
!
!write(24,*) r1, r2, aR(r1), aR(r2), &
! TransDip_dimer_MC_off(Ah1(r1),Bh1(r1),kinh1(r1),kouth1(r1),Ae(r2),Be(r2),kine(r2),koute(r2),aR(r1), aR(r2),linker), &
! TransDip_dimer_MC_off(Ah2(r1),Bh2(r1),kinh2(r1),kouth2(r1),Ae(r2),Be(r2),kine(r2),koute(r2),aR(r1), aR(r2),linker)
!
!enddo
!
!write(23,*) '     '
!
!enddo

!endif

!enddo

!D12 Coulomb integral

!write(6,*) aR(n), D12_in_in(oo,Ah1(n),kinh1(n),Ae(n),kine(n),Ah2(n),kinh2(n),Ae(n),kine(n),aR(n)) + &
!                  D12_out_in(oo,Bh1(n),kouth1(n),Ae(n),kine(n),Bh2(n),kouth2(n),Ae(n),kine(n),aR(n)) + &
!                  D12_in_out(oo,Ah1(n),kinh1(n),Be(n),koute(n),Ah2(n),kinh2(n),Be(n),koute(n),aR(n)) + &
!                  D12_out_out(oo,Bh1(n),kouth1(n),Be(n),koute(n),Bh2(n),kouth2(n),Be(n),koute(n),aR(n)), & 
!                  D12ex_in_in(oo,Ah1(n),kinh1(n),Ae(n),kine(n),Ah2(n),kinh2(n),Ae(n),kine(n),aR(n)) + & 
!                  D12ex_out_in(oo,Bh1(n),kouth1(n),Be(n),koute(n),Ah2(n),kinh2(n),Ae(n),kine(n),aR(n)) + &
!                  D12ex_in_out(oo,Ah1(n),kinh1(n),Ae(n),kine(n),Bh2(n),kouth2(n),Be(n),koute(n),aR(n)) + &
!                  D12ex_out_out(oo,Bh1(n),kouth1(n),Be(n),koute(n),Bh2(n),kouth2(n),Be(n),koute(n),aR(n))
!write(50,*) aR(n), &!D12_in_in(oo,Ah1(n),kinh1(n),Ae(n),kine(n),Ah1(n),kinh1(n),Ae(n),kine(n),aR(n)) + &

!write(6,*) aR(n), D12dir(oo,Ah1(n),Bh1(n),kinh1(n),kouth1(n),Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),Ae(n),Be(n),kine(n),koute(n),aR(n)), &
!               D12ex(oo,Ah1(n),Bh1(n),kinh1(n),kouth1(n),Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),Ae(n),Be(n),kine(n),koute(n),aR(n)), &
!               D12ex_8loops(oo,Ah1(n),Bh1(n),kinh1(n),kouth1(n),Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),Ae(n),Be(n),kine(n),koute(n),aR(n))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if ( dyn .eq. 'y' ) then 

open(21,file='Ham0.dat')
open(22,file='Pulse.dat')
open(23,file='Hamt.dat')
open(24,file='Realc.dat')
open(25,file='Popc.dat')
open(26,file='Norm.dat')

Ham      = 0.0
TransHam = 0.0

if ( aA .eq. aB ) then
call make_Ham_ho
TransHam(0,1) = TransDip_Ana_h1e(1)
TransHam(0,2) = TransDip_Ana_h2e(2)
TransHam(0,3) = TransHam(0,1)
TransHam(0,4) = TransHam(0,2)
TransHam(0,5) = TransDip_Fit_h1e_ho(aR(1),link)*1d-33
TransHam(0,6) = TransDip_Fit_h2e_ho(aR(2),link)*1d-33
TransHam(0,7) = TransHam(0,5)
TransHam(0,8) = TransHam(0,6)
do i=1,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

elseif ( aA .ne. aB ) then
call make_Ham_he
TransHam(0,1) = TransDip_Ana_h1e(1)
TransHam(0,2) = TransDip_Ana_h2e(1)
TransHam(0,3) = TransDip_Ana_h1e(2)
TransHam(0,4) = TransDip_Ana_h2e(2)
TransHam(0,5) = TransDip_Fit_h1e_he(aR(1),aR(2))*1d-33
TransHam(0,6) = TransDip_Fit_h2e_he(aR(1),aR(2))*1d-33
TransHam(0,7) = TransDip_Fit_h1e_he(aR(1),aR(2))*1d-33
TransHam(0,8) = TransDip_Fit_h2e_he(aR(1),aR(2))*1d-33
do i=1,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

endif

!do i=0,nstates-1
!write(21,'(9es14.6e2)') (real(Ham(i,j)), j=0,nstates-1)
!enddo
!
!do i=0,nstates-1
!write(6,'(9es14.6e2)') (real(TransHam(i,j)), j=0,nstates-1)
!enddo

!!!!!INITIAL POPULATIONS
c0(0) = 1.0
c0(1) = 0.0
c0(2) = 0.0
c0(3) = 0.0
c0(4) = 0.0
c0(5) = 0.0
c0(6) = 0.0
c0(7) = 0.0
c0(8) = 0.0

xHam = dcmplx(Ham,0.0)
xHamt(:,:,0) = xHam(:,:)
xTransHam = dcmplx(TransHam,0.0)
xc0 = dcmplx(c0,0.0)
xc = 0.0
xc(:,0) = xc0(:)
xt01 = dcmplx(t01,0.0)
xt02 = dcmplx(t02,0.0)
xt03 = dcmplx(t03,0.0)

do t=0,ntime
 
xtime = dcmplx(t*timestep,0.0)

write(23,*) real(xtime)

do i=0,nstates-1
   do j=0,nstates-1
xHamt(i,j,t)  = xHam(i,j) - xTransHam(i,j) * xEd * cos(xomega*(xtime-xt01)+xphase) * exp(-1.0*(xtime-xt01)**2/(2*(xwidth**2))) &
                          - xTransHam(i,j) * xEd * cos(xomega*(xtime-xt01)+xphase) * exp(-1.0*(xtime-xt02)**2/(2*(xwidth**2))) &
                          - xTransHam(i,j) * xEd * cos(xomega*(xtime-xt01)+xphase) * exp(-1.0*(xtime-xt03)**2/(2*(xwidth**2)))
enddo

write(23,'(9es14.6e2)') (real(xHamt(i,k,t)), k=0,nstates-1)

enddo

write(23,*) 

enddo

do t=0,ntime

xtime = dcmplx(t*timestep,0.0)

write(22,*) real(xtime) , real(xEd * cos(xomega*(xtime-xt01)+xphase) * exp(-1.0*(xtime-xt01)**2/(2*(xwidth**2)))) , &
                          real(xEd * cos(xomega*(xtime-xt02)+xphase) * exp(-1.0*(xtime-xt02)**2/(2*(xwidth**2)))) , &
                          real(xEd * cos(xomega*(xtime-xt03)+xphase) * exp(-1.0*(xtime-xt03)**2/(2*(xwidth**2)))) 

k1 =dcmplx(0.0,0.0)
k2 =dcmplx(0.0,0.0)
k3 =dcmplx(0.0,0.0)
k4 =dcmplx(0.0,0.0)

do i=0,nstates-1
do j=0,nstates-1
k1(i) = k1(i) + (-1.0)*(im/xhbar)*&
(xHam(i,j) - xTransHam(i,j) * xEd * cos(xomega*(xtime-xt01)+xphase) * exp(-1.0*(xtime-xt01)**2/(2*(xwidth**2))) - &
             xTransHam(i,j) * xEd * cos(xomega*(xtime-xt02)+xphase) * exp(-1.0*(xtime-xt02)**2/(2*(xwidth**2))) - &
             xTransHam(i,j) * xEd * cos(xomega*(xtime-xt03)+xphase) * exp(-1.0*(xtime-xt03)**2/(2*(xwidth**2))))*&
xc(j,t)
enddo
enddo

do i=0,nstates-1
do j=0,nstates-1
k2(i) = k2(i) + (-1.0)*(im/xhbar)*&
(xHam(i,j) - xTransHam(i,j) * xEd * cos(xomega*((xtime-xt01)+xh/2)+xphase) * exp(-1.0*((xtime+xh/2)-xt01)**2/(2*(xwidth)**2)) - &
             xTransHam(i,j) * xEd * cos(xomega*((xtime-xt02)+xh/2)+xphase) * exp(-1.0*((xtime+xh/2)-xt02)**2/(2*(xwidth)**2)) - &
             xTransHam(i,j) * xEd * cos(xomega*((xtime-xt03)+xh/2)+xphase) * exp(-1.0*((xtime+xh/2)-xt03)**2/(2*(xwidth)**2)))*&
(xc(j,t)+(xh/2)*k1(j))
enddo
enddo

do i=0,nstates-1
do j=0,nstates-1
k3(i) = k3(i) + (-1.0)*(im/xhbar)*&
(xHam(i,j) - xTransHam(i,j) * xEd * cos(xomega*((xtime-xt01)+xh/2)+xphase) * exp(-1.0*((xtime+xh/2)-xt01)**2/(2*(xwidth)**2)) - &
             xTransHam(i,j) * xEd * cos(xomega*((xtime-xt02)+xh/2)+xphase) * exp(-1.0*((xtime+xh/2)-xt02)**2/(2*(xwidth)**2)) - &
             xTransHam(i,j) * xEd * cos(xomega*((xtime-xt03)+xh/2)+xphase) * exp(-1.0*((xtime+xh/2)-xt03)**2/(2*(xwidth)**2)))*&
(xc(j,t)+(xh/2)*k2(j))
enddo
enddo

do i=0,nstates-1
do j=0,nstates-1
k4(i) = k4(i) + (-1.0)*(im/xhbar)*&
(xHam(i,j) - xTransHam(i,j) * xEd * cos(xomega*((xtime-xt01)+xh)+xphase) * exp(-1.0*((xtime+xh)-xt01)**2/(2*(xwidth)**2)) - &
             xTransHam(i,j) * xEd * cos(xomega*((xtime-xt02)+xh)+xphase) * exp(-1.0*((xtime+xh)-xt02)**2/(2*(xwidth)**2)) -&
             xTransHam(i,j) * xEd * cos(xomega*((xtime-xt03)+xh)+xphase) * exp(-1.0*((xtime+xh)-xt03)**2/(2*(xwidth)**2)))*&
(xc(j,t)+xh*k3(j))
enddo
enddo

do i=0,nstates-1
xc(i,t+1) = xc(i,t)+(xh/6)*(k1(i)+2*k2(i)+2*k3(i)+k4(i))
enddo

!!!NORM
write(26,*) real(xtime), &
real(xc(1,t))**2+aimag(xc(1,t))**2+real(xc(2,t))**2+aimag(xc(2,t))**2+real(xc(3,t))**2+aimag(xc(3,t))**2+&
real(xc(4,t))**2+aimag(xc(4,t))**2+real(xc(5,t))**2+aimag(xc(5,t))**2+real(xc(6,t))**2+aimag(xc(6,t))**2+&
real(xc(7,t))**2+aimag(xc(7,t))**2+real(xc(8,t))**2+aimag(xc(8,t))**2+real(xc(0,t))**2+aimag(xc(0,t))**2

!!!!POPULATIONS
write(25,*) real(xtime), real(xc(0,t))**2+aimag(xc(0,t))**2, real(xc(1,t))**2+aimag(xc(1,t))**2,&
                         real(xc(2,t))**2+aimag(xc(2,t))**2, real(xc(3,t))**2+aimag(xc(3,t))**2,&
                         real(xc(4,t))**2+aimag(xc(4,t))**2, real(xc(5,t))**2+aimag(xc(5,t))**2,&
                         real(xc(6,t))**2+aimag(xc(6,t))**2, real(xc(7,t))**2+aimag(xc(7,t))**2,&
                         real(xc(8,t))**2+aimag(xc(8,t))**2

enddo

endif

if ( vers .eq. 'singl') then
call makeOutputSingle
elseif ( vers .eq. 'dimer') then
call makeOutputDimer
elseif ( vers .eq. 'range') then
call makeOutputRange
elseif ( vers .eq. 'rdmho') then
call makeOutputRdmho
elseif ( vers .eq. 'rdmhe') then
call makeOutputRdmhe
endif

deallocate(E,diffe,diffh,minEe,minEh)

end
