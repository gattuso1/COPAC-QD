include 'specfun.f90'

program ModelQD_one

use Constants
use Variables
use Integrals
use Vectors

implicit none

double precision, external:: s13adf, ei, eone, nag_bessel_j0

character*64 :: arg1, arg2, filename, argA

character*2 :: flag

integer :: je,jh,k,nsteps,r,rsteps,ifail, r2

real(dp) :: maxr0,stepr0,Ef,le,re,lh,rh,delta, start, finish, mu, A, epsinf, aR1, aR2, &
        Rine, Route, Rinh1, Routh1, Rinh2, Routh2, RadProbe, RadProbh1, RadProbh2, r0

real(dp),allocatable :: Eeh1(:), Eeh2(:), Ae(:), Ah1(:), Ah2(:), Be(:), Bh1(:), Bh2(:), Cb_eh1(:), Cb_eh2(:), &
                        I1eh1(:), I1eh2(:), I2eh1(:), I2eh2(:), I3eh1(:), I3eh2(:), kine(:), kinh1(:), kinh2(:), &
                        koute(:), kouth1(:), kouth2(:)

real(dp),allocatable :: minEe(:,:),minEh(:,:),diffe(:), diffh(:), E(:), Ham(:,:)

call getVariables

maxr0=  3d-9
stepr0= 0.05d-9
delta=  0.00001d-18
Ef=     1.28174d-18
nsteps= int(Ef/delta)
rsteps= int(maxr0/stepr0)+1
ifail=  1

write(6,*) rsteps

allocate(E(nsteps))
allocate(diffe(nsteps))
allocate(diffh(nsteps))
allocate(minEe(nsteps,ndots))
allocate(minEh(nsteps,ndots))
allocate(Eeh1(ndots)) 
allocate(Eeh2(ndots)) 
allocate(Ae(ndots)) 
allocate(Ah1(ndots)) 
allocate(Ah2(ndots)) 
allocate(Be(ndots)) 
allocate(Bh1(ndots)) 
allocate(Bh2(ndots)) 
allocate(Cb_eh1(ndots)) 
allocate(Cb_eh2(ndots))
allocate(I1eh1(ndots)) 
allocate(I1eh2(ndots)) 
allocate(I2eh1(ndots)) 
allocate(I2eh2(ndots)) 
allocate(I3eh1(ndots))
allocate(I3eh2(ndots)) 
allocate(kine(ndots)) 
allocate(kinh1(ndots)) 
allocate(kinh2(ndots)) 
allocate(koute(ndots)) 
allocate(kouth1(ndots)) 
allocate(kouth2(ndots))
allocate(Ham(9,9))

open(10,file='Ee.dat')
open(11,file='E1se-1sh.dat')
open(12,file='E1se-2sh.dat')
open(13,file='E1se-1sh-Cb.dat')
open(14,file='E1se-2sh-Cb.dat')
open(15,file='Renv-QDA.dat')
open(16,file='Renv-QDB.dat')
open(17,file='wavefunctionA.dat')
open(18,file='wavefunctionB.dat')
open(19,file='V.dat')
open(20,file='Elevels.dat')
open(21,file='TransDiph1e.dat')
open(22,file='TransDiph2e.dat')
open(23,file='TransDip-dimer-eh1.dat')
open(24,file='TransDip-dimer-eh2.dat')

k=1

!Computation of energies

do n = 1,ndots
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

write(6,*) "Egap e-h1", Eeh1(n)/elec, "including correction", Cb_eh1(n)/elec
write(6,*) "Egap e-h2", Eeh2(n)/elec, "including correction", Cb_eh2(n)/elec 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Wave function QDA and radial probability

do r=1,rsteps

r0=r*stepr0

Rine=Ae(n)*sin(kine(n)*r0)/(sqrt(4*pi)*r0)
Route=Be(n)*exp(-1*koute(n)*r0)/(sqrt(4*pi)*r0)

Rinh1=Ah1(n)*sin(kinh1(n)*r0)/(sqrt(4*pi)*r0)
Routh1=Bh1(n)*exp(-1*kouth1(n)*r0)/(sqrt(4*pi)*r0)

Rinh2=Ah2(n)*sin(kinh2(n)*r0)/(sqrt(4*pi)*r0)
Routh2=Bh2(n)*exp(-1*kouth2(n)*r0)/(sqrt(4*pi)*r0)

if (r0 .le. aR(n)) then
        RadProbe=Rine
        RadProbh1=Rinh1
        RadProbh2=Rinh2
else if (r0 .gt. aR(n)) then
        RadProbe=Route
        RadProbh1=Routh1
        RadProbh2=Routh2
endif

write(17,*) r0, RadProbe, RadProbh1, RadProbh2
write(15,*) r0, r0**2*RadProbe**2, r0**2*RadProbh1**2, r0**2*RadProbh2**2

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Normalization integrals
if ( o_Norm == 'y' ) then
write(6,*) "Norm e  ana", Norm_Ana(aR(n),Ae(n),Be(n),kine(n),koute(n))
write(6,*) "Norm e  num", Norm_Num(aR(n),Ae(n),Be(n),kine(n),koute(n))
write(6,*) "Norm h1 ana", Norm_Ana(aR(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n))
write(6,*) "Norm h1 num", Norm_Num(aR(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n))
write(6,*) "Norm h2 ana", Norm_Ana(aR(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n))
write(6,*) "Norm h2 num", Norm_Num(aR(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n))
endif

!Overlap integrals
if ( o_Over == 'y' ) then
write(6,*) "Overlap e-h1 ana", abs(OverlapAna(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))
write(6,*) "Overlap e-h1 num", abs(OverlapNum(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))
write(6,*) "Overlap e-h2 ana", abs(OverlapAna(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))
write(6,*) "Overlap e-h2 num", abs(OverlapNum(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))
endif

!Coulomb correction
if ( o_Coul == 'y' ) then
write(6,*) "Cb dir+ex eh1-eh1", DXXdir(Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
                                       Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),aR(n)) + &
                                DXXex(Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
                                       Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),aR(n))
write(6,*) "Cb dir+ex eh1-eh2", DXXdir(Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
                                       Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n)) + &
                                DXXex(Ae(n),Be(n),kine(n),koute(n),Ah1(n),Bh1(n),kinh1(n),kouth1(n),&
                                      Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n))
write(6,*) "Cb dir+ex eh2-eh2", DXXdir(Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),&
                                       Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n)) + &
                                DXXex(Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),&
                                      Ae(n),Be(n),kine(n),koute(n),Ah2(n),Bh2(n),kinh2(n),kouth2(n),aR(n))
endif

!Dipole moment 
if ( o_DipS == 'y' ) then
write(6,*) "Dipole e-h1 ana",  abs(TransDip_Ana(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))/Cm_to_D
write(6,*) "Dipole e-h1 num",  abs(TransDip_Num(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n)))/Cm_to_D 
write(6,*) "Dipole e-h1 EMA",  abs(TransDip_EMA(Eeh1(n),aR(n)))/Cm_to_D
write(6,*) "Dipole e-h2 ana",  abs(TransDip_Ana(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))/Cm_to_D
write(6,*) "Dipole e-h2 num",  abs(TransDip_Num(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n)))/Cm_to_D
write(6,*) "Dipole e-h1 EMA",  abs(TransDip_EMA(Eeh2(n),aR(n)))/Cm_to_D
endif

!Oscilator strength
if ( o_Osci == 'y' ) then
write(6,*) "Oscillator strength e-h1", ((2*me)/(3*elec**2*hbar**2))*Eeh1(n)*&
                                       TransDip_Ana(Ae(n),Ah1(n),Be(n),Bh1(n),kine(n),kinh1(n),koute(n),kouth1(n),aR(n))**2
write(6,*) "Oscillator strength e-h2", ((2*me)/(3*elec**2*hbar**2))*Eeh2(n)*&
                                       TransDip_Ana(Ae(n),Ah2(n),Be(n),Bh2(n),kine(n),kinh2(n),koute(n),kouth2(n),aR(n))**2
endif

!molar extinction coefficient
if ( o_Exti == 'y' ) then
A=1.25e5*elec**2
mu=(me*mh)/(me+mh)
write(6,*) "Exct. Coef. e-h1", (A*mu*aR(n)**2)/(sqrt(2*pi)*Eeh1(n)*hbar**2*0.1)
write(6,*) "Exct. Coef. e-h2", (A*mu*aR(n)**2)/(sqrt(2*pi)*Eeh2(n)*hbar**2*0.1)
endif

enddo

!Dielectric constant size dependent
!epsin(n)f = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*aR(n)))**1.2)
!epsin(n)f = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*1.35e-9))**1.2)
!disteh = (36.0/35.0)*aR(n)
!disteh = (36.0/35.0)*1.35e-9
!epsin(n)frR= 1.0/((1.0/epsin(n)f)-((1.0/epsin(n)f)-(1.0/(epsin(n)f+3.5)))*(1-(exp(-disteh/rhoe)+exp(-disteh/rhoh))/2))

!write(6,*) aR(n), epsin(n)f , 1.0/((1.0/epsin(n)f)-((1.0/epsin(n)f)-(1.0/(epsin(n)f+3.5)))*(1-(exp(-disteh/rhoe)+exp(-disteh/rhoh))/2)), &
!                    (exp(-disteh/rhoe)+exp(-disteh/rhoh))/2


!write(6,*) aR(n), epsout + (epsin(n)f-epsout)*((pi/2) - atan(slope*(aR(n)-1.35e-9)))/pi , &
!                  epsout + (epsin(n)frR-epsout)*((pi/2) - atan(slope*(aR(n)-1.35e-9)))/pi

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

!enddo
!write(6,*) "         "


!enddo

!do r=1,rsteps
!
!aR1=r*stepr0
!
!do r2=1,rsteps
!
!aR2=r2*stepr0

write(6,*) epsin(1), epsin(2), aR(1), aR(2)

!if ( o_DipD == 'y' ) then
write(6,*) "e-h1", aR(1), aR(2), &
                3.33564d30*TransDip_dimer_MC(Ah1(1),Bh1(1),kinh1(1),kouth1(1),Ae(2),Be(2),kine(2),koute(2),aR(1),aR(2),linker)
write(6,*) "h1-e", aR(2), aR(1), &
               3.33564d30*TransDip_dimer_MC(Ae(2),Be(2),kine(2),koute(2),Ah1(1),Bh1(1),kinh1(1),kouth1(1),aR(2),aR(1),linker)
!write(24,*) aR1, aR2,  TransDip_dimer_MC(Ah2(1),Bh2(1),kinh2(1),kouth2(1),Ae(2),Be(2),kine(2),koute(2),aR1,aR2,linker)
!write(6,*) 'e-h1',  TransDip_dimer_MC(Ae(2),Be(2),kine(2),koute(2),Ah1(1),Bh1(1),kinh1(1),kouth1(1),aR(2),aR(1),linker)
!write(6,*) 'h2-e',  TransDip_dimer_MC(Ah2(1),Bh2(1),kinh2(1),kouth2(1),Ae(2),Be(2),kine(2),koute(2),aR(1),aR(2),linker)
!write(6,*) 'e-h2',  TransDip_dimer_MC(Ae(2),Be(2),kine(2),koute(2),Ah2(1),Bh2(1),kinh2(1),kouth2(1),aR(2),aR(1),linker)

!enddo
!
!enddo

                  ! TransDip_dimer_MC(Ah1(n),Bh1(n),kinh1(n),kouth1(n),Ae(n),Be(n),kine(n),koute(n),aR(n),aB,linker)
!write(6,*) linker, TransDip_dimer_MC(Ah2(n),Bh2(n),kinh2(n),kouth2(n),Ae(n),Be(n),kine(n),koute(n),aB,aR(n),linker)
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

!  OA= Ana_in_in(m,Ah1(n),Ae(n),kinh1(n),kine(n),aR(n)) + & 
!      Ana_out_out(m,Bh1(n),Be(n),kouth1(n),koute(n),aR(n)) 
!
!  OB= Ana_in_in(m,Ah2(n),Ae(n),kinh2(n),kine(n),aR(n)) + & 
!      Ana_out_out(m,Bh2(n),Be(n),kouth2(n),koute(n),aR(n)) 

!write(6,*) Ana_in_in(m,Ah1(n),Ae(n),kinh1(n),kine(n),aR(n)), & 
!      Ana_out_out(m,Bh1(n),Be(n),kouth1(n),koute(n),aR(n)) 

!write(6,*) Ana_in_in(m,Ah2(n),Ae(n),kinh2(n),kine(n),aR(n)), & 
!      Ana_out_out(m,Bh2(n),Be(n),kouth2(n),koute(n),aR(n)) 

!if ( o_DCou == 'y' ) then
!include 'fetchIntegrals.f90'
!
!do i = 1,nstates
!write(14,"(9f6.3)") (Ham(i,j), j=0,nterms-1)
!enddo
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(E,diffe,diffh,minEe,minEh)

END
