include 'specfun.f90'

program ModelQD_one

use Constants
use Variables
use Integrals
use Vectors

implicit none

double precision, external:: s13adf, ei, eone, nag_bessel_j0

character*64:: arg1, arg2, filename, argA

integer:: i,je,jh,k,nsteps,r,rsteps,ifail

real(dp):: maxr0,stepr0,Ef,le,re,lh,rh,delta, & 
        start, finish, mu, A, epsinf, &
        Rine, Route, Rinh1, Routh1, Rinh2, Routh2, RadProbe, RadProbh1, RadProbh2, &
        r0, Ae, Ah1, Ah2, Be, Bh1, Bh2, Cb_eh1, Cb_eh2,&
        Eeh1, Eeh2, &
        I1eh1, I1eh2, I2eh1, I2eh2, I3eh1, I3eh2, kine, kinh1, kinh2, koute, kouth1, kouth2

real(dp),allocatable:: minEe(:),minEh(:),diffe(:), diffh(:), E(:)

call getVariables

maxr0=  4d-9
stepr0= 0.05d-9
delta=  0.00001d-18
Ef=     1.28174d-18
nsteps= int(Ef/delta)
rsteps= int(maxr0/stepr0)+1
ifail=  1

allocate(E(nsteps))
allocate(diffe(nsteps))
allocate(diffh(nsteps))
allocate(minEe(nsteps))
allocate(minEh(nsteps))

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
open(23,file='TransDip-dimer.dat')

k=1

!Computation of energies

!diffe(:)=0
!diffh(:)=0
!E(:)=0
!minEe(:)=0 
!minEh(:)=0 
i=0
r=0

je=1
jh=1

do i=1,nsteps
E(i)=delta*i
diffe(i) = abs(sqrt(2*me*E(i))/hbar * aA * 1/tan(sqrt(2*me*E(i))/hbar * aA) - 1 + (me/m0) + (me*aA)/(hbar) &
           * sqrt((2/m0)*(V0e-E(i))))
if ((diffe(0) .eq. 0.000) .and. (diffh(0) .eq. 0.00)) then 
        diffe(0)=diffe(i)
        diffh(0)=diffe(i)
endif


if (diffe(i) .le. diffe(i-1)) then
        minEe(je) = E(i)
elseif ( (diffe(i) .ge. diffe(i-1)) .and. (E(i-1) .eq. minEe(je)) ) then
        je=je+1
endif

diffh(i) = abs(sqrt(2*mh*E(i))/hbar * aA * 1/tan(sqrt(2*mh*E(i))/hbar * aA) - 1 + (mh/m0) + (mh*aA)/(hbar) &
           * sqrt((2/m0)*(V0h-E(i))))
if (diffh(i) .le. diffh(i-1)) then
        minEh(jh) = E(i)
elseif ( (diffh(i) .ge. diffh(i-1)) .and. (E(i-1) .eq. minEh(jh)) ) then
        jh=jh+1
endif
enddo

!wave vectors in and out
kine=sqrt(2*me*minEe(1))/hbar
koute=sqrt(2*m0*(V0e-minEe(1)))/hbar

kinh1=sqrt(2*mh*minEh(1))/hbar
kouth1=sqrt(2*m0*(V0h-minEh(1)))/hbar

kinh2=sqrt(2*mh*minEh(2))/hbar
kouth2=sqrt(2*m0*(V0h-minEh(2)))/hbar

!normalization factors
Ae=1/sqrt(aA/2-sin(2*kine*aA)/(4*kine)+sin(kine*aA)**2/(2*koute))
Be=Ae*sin(kine*aA)*exp(koute*aA) 

Ah1=1/sqrt(aA/2-sin(2*kinh1*aA)/(4*kinh1)+sin(kinh1*aA)**2/(2*kouth1))
Bh1=Ah1*sin(kinh1*aA)*exp(kouth1*aA) 

Ah2=1/sqrt(aA/2-sin(2*kinh2*aA)/(4*kinh2)+sin(kinh2*aA)**2/(2*kouth2))
Bh2=Ah2*sin(kinh2*aA)*exp(kouth2*aA)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Comlomb correction (analytical) between electron and holes

I1eh1=(Ae**2*Ah1**2*aA/4)*(1-sin(2*kine*aA)/(2*kine*aA)-s13adf(2*kine*aA,ifail)/(2*kine*aA)+&
         (s13adf(2*kine*aA-2*kinh1*aA,ifail)+s13adf(2*kine*aA+2*kinh1*aA,ifail))/(4*kine*aA)) + &
         (Ae**2*Ah1**2*aA/4)*(1-sin(2*kinh1*aA)/(2*kinh1*aA)-s13adf(2*kinh1*aA,ifail)/(2*kinh1*aA)+&
         (s13adf(2*kinh1*aA-2*kine*aA,ifail)+s13adf(2*kinh1*aA+2*kine*aA,ifail))/(4*kinh1*aA)) 

I2eh1=(Ae**2*Bh1**2*aA/2)*(1-sin(2*kine*aA)/(2*kine*aA))*eone(2*kouth1*aA)+&
         (Ah1**2*Be**2*aA/2)*(1-sin(2*kinh1*aA)/(2*kinh1*aA))*eone(2*koute*aA)

I3eh1=(Be**2*Bh1**2*aA/(2*koute*aA)*(exp(-2*koute*aA)*eone(2*kouth1*aA)-eone(2*koute*aA+2*kouth1*aA)))+&
         (Be**2*Bh1**2*aA/(2*kouth1*aA)*(exp(-2*kouth1*aA)*eone(2*koute*aA)-eone(2*kouth1*aA+2*koute*aA)))

I1eh2=(Ae**2*Ah2**2*aA/4)*(1-sin(2*kine*aA)/(2*kine*aA)-s13adf(2*kine*aA,ifail)/(2*kine*aA)+&
         (s13adf(2*kine*aA-2*kinh2*aA,ifail)+s13adf(2*kine*aA+2*kinh2*aA,ifail))/(4*kine*aA)) + &
         (Ae**2*Ah2**2*aA/4)*(1-sin(2*kinh2*aA)/(2*kinh2*aA)-s13adf(2*kinh2*aA,ifail)/(2*kinh1*aA)+&
         (s13adf(2*kinh2*aA-2*kine*aA,ifail)+s13adf(2*kinh2*aA+2*kine*aA,ifail))/(4*kinh2*aA)) 

I2eh2=(Ae**2*Bh2**2*aA/2)*(1-sin(2*kine*aA)/(2*kine*aA))*eone(2*kouth2*aA)+&
         (Ah2**2*Be**2*aA/2)*(1-sin(2*kinh2*aA)/(2*kinh2*aA))*eone(2*koute*aA)

I3eh2=(Be**2*Bh2**2*aA/(2*koute*aA)*(exp(-2*koute*aA)*eone(2*kouth2*aA)-eone(2*koute*aA+2*kouth2*aA)))+&
         (Be**2*Bh2**2*aA/(2*kouth2*aA)*(exp(-2*kouth2*aA)*eone(2*koute*aA)-eone(2*kouth2*aA+2*koute*aA)))


Cb_eh1=(elec**2/(4*pi*epsin*eps0))*(I1eh1+I2eh1+I3eh1)

Cb_eh2=(elec**2/(4*pi*eps*eps0))*(I1eh2+I2eh2+I3eh2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Eeh1 = (minEe(1)+minEh(1))+V0-Cb_eh1 
Eeh2 = (minEe(1)+minEh(2))+V0-Cb_eh2

write(6,*) "Egap e-h1", Eeh1/elec, "including correction", Cb_eh1/elec
write(6,*) "Egap e-h2", Eeh2/elec, "including correction", Cb_eh2/elec 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Wave function QDA and radial probability

do r=1,rsteps

r0=r*stepr0

Rine=Ae*sin(kine*r0)/(sqrt(4*pi)*r0)
Route=Be*exp(-1*koute*r0)/(sqrt(4*pi)*r0)

Rinh1=Ah1*sin(kinh1*r0)/(sqrt(4*pi)*r0)
Routh1=Bh1*exp(-1*kouth1*r0)/(sqrt(4*pi)*r0)

Rinh2=Ah2*sin(kinh2*r0)/(sqrt(4*pi)*r0)
Routh2=Bh2*exp(-1*kouth2*r0)/(sqrt(4*pi)*r0)

if (r0 .le. aA) then
        RadProbe=Rine
        RadProbh1=Rinh1
        RadProbh2=Rinh2
else if (r0 .gt. aA) then
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
write(6,*) "Norm e  ana", Norm_Ana(aA,Ae,Be,kine,koute)
write(6,*) "Norm e  num", Norm_Num(aA,Ae,Be,kine,koute)
write(6,*) "Norm h1 ana", Norm_Ana(aA,Ah1,Bh1,kinh1,kouth1)
write(6,*) "Norm h1 num", Norm_Num(aA,Ah1,Bh1,kinh1,kouth1)
write(6,*) "Norm h2 ana", Norm_Ana(aA,Ah2,Bh2,kinh2,kouth2)
write(6,*) "Norm h2 num", Norm_Num(aA,Ah2,Bh2,kinh2,kouth2)
endif

!Overlap integrals
if ( o_Over == 'y' ) then
write(6,*) "Overlap e-h1 ana", abs(OverlapAna(Ae,Ah1,Be,Bh1,kine,kinh1,koute,kouth1,aA))
write(6,*) "Overlap e-h1 num", abs(OverlapNum(Ae,Ah1,Be,Bh1,kine,kinh1,koute,kouth1,aA))
write(6,*) "Overlap e-h2 ana", abs(OverlapAna(Ae,Ah2,Be,Bh2,kine,kinh2,koute,kouth2,aA))
write(6,*) "Overlap e-h2 num", abs(OverlapNum(Ae,Ah2,Be,Bh2,kine,kinh2,koute,kouth2,aA))
endif

!Coulomb correction
if ( o_Coul == 'y' ) then
write(6,*) "Cb dir+ex eh1-eh1", DXXdir(Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,aA) + &
                                DXXex(Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,aA)
write(6,*) "Cb dir+ex eh1-eh2", DXXdir(Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah2,Bh2,kinh2,kouth2,aA) + &
                                DXXex(Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah2,Bh2,kinh2,kouth2,aA)
write(6,*) "Cb dir+ex eh2-eh2", DXXdir(Ae,Be,kine,koute,Ah2,Bh2,kinh2,kouth2,Ae,Be,kine,koute,Ah2,Bh2,kinh2,kouth2,aA) + &
                                DXXex(Ae,Be,kine,koute,Ah2,Bh2,kinh2,kouth2,Ae,Be,kine,koute,Ah2,Bh2,kinh2,kouth2,aA)
endif

!Dipole moment 
if ( o_DipS == 'y' ) then
write(6,*) "Dipole e-h1 ana",  abs(TransDip_Ana(Ae,Ah1,Be,Bh1,kine,kinh1,koute,kouth1,aA))/Cm_to_D
write(6,*) "Dipole e-h1 num",  abs(TransDip_Num(Ae,Ah1,Be,Bh1,kine,kinh1,koute,kouth1,aA))/Cm_to_D 
write(6,*) "Dipole e-h1 EMA",  abs(TransDip_EMA(Eeh1,aA))/Cm_to_D
write(6,*) "Dipole e-h2 ana",  abs(TransDip_Ana(Ae,Ah2,Be,Bh2,kine,kinh2,koute,kouth2,aA))/Cm_to_D
write(6,*) "Dipole e-h2 num",  abs(TransDip_Num(Ae,Ah2,Be,Bh2,kine,kinh2,koute,kouth2,aA))/Cm_to_D
write(6,*) "Dipole e-h1 EMA",  abs(TransDip_EMA(Eeh2,aA))/Cm_to_D
endif

!Oscilator strength
if ( o_Osci == 'y' ) then
write(6,*) "Oscillator strength e-h1", ((2*me)/(3*elec**2*hbar**2))*Eeh1*&
                                       TransDip_Ana(Ae,Ah1,Be,Bh1,kine,kinh1,koute,kouth1,aA)**2
write(6,*) "Oscillator strength e-h2", ((2*me)/(3*elec**2*hbar**2))*Eeh2*&
                                       TransDip_Ana(Ae,Ah2,Be,Bh2,kine,kinh2,koute,kouth2,aA)**2
endif

!molar extinction coefficient
if ( o_Exti == 'y' ) then
A=1.25e5*elec**2
mu=(me*mh)/(me+mh)
write(6,*) "Exct. Coef. e-h1", (A*mu*aA**2)/(sqrt(2*pi)*Eeh1*hbar**2*0.1)
write(6,*) "Exct. Coef. e-h2", (A*mu*aA**2)/(sqrt(2*pi)*Eeh2*hbar**2*0.1)
endif

!Dielectric constant size dependent
!epsinf = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*aA))**1.2)
!epsinf = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*1.35e-9))**1.2)
!disteh = (36.0/35.0)*aA
!disteh = (36.0/35.0)*1.35e-9
!epsinfrR= 1.0/((1.0/epsinf)-((1.0/epsinf)-(1.0/(epsinf+3.5)))*(1-(exp(-disteh/rhoe)+exp(-disteh/rhoh))/2))

!write(6,*) aA, epsinf , 1.0/((1.0/epsinf)-((1.0/epsinf)-(1.0/(epsinf+3.5)))*(1-(exp(-disteh/rhoe)+exp(-disteh/rhoh))/2)), &
!                    (exp(-disteh/rhoe)+exp(-disteh/rhoh))/2


!write(6,*) aA, epsout + (epsinf-epsout)*((pi/2) - atan(slope*(aA-1.35e-9)))/pi , &
!                  epsout + (epsinfrR-epsout)*((pi/2) - atan(slope*(aA-1.35e-9)))/pi

!Dipole moment h1-e dimer
!write(6,*) aA, linker, TransDip_dimer_MC(Ah1,Bh1,kinh1,kouth1,&
!                                           Ae,Be,kine,koute,aA,aA,linker)

!write(6,*) aA, linker, TransDip_dimer_cart(m,Ah1,Bh1,kinh1,kouth1,&
!                                              Ae,Be,kine,koute,aA,aA,linker)
!                 TransDip_dimer_cart(m,Ah2,Bh2,kinh2,kouth2,&
!                                              Ae,Be,kine,koute,aA,aA,linker)
!                    TransDip_dimer_cart(m,Ae,Be,kine,koute,&
!                                              Ah1,Bh1,kinh1,kouth1,aA,aA,linker),&
!                    TransDip_dimer_cart(m,Ae,Be,kine,koute,&
!                                              Ah2,Bh2,kinh2,kouth2,aA,aA,linker)

!enddo
!write(6,*) "         "


!enddo

!if ( o_DipD == 'y' ) then
!write(6,*)  TransDip_dimer_MC(Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,aA,aB,linker), &
!                 TransDip_dimer_MC(Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,aB,aA,linker)
                  ! TransDip_dimer_MC(Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,aA,aB,linker)
!write(6,*) linker, TransDip_dimer_MC(Ah2,Bh2,kinh2,kouth2,Ae,Be,kine,koute,aB,aA,linker)
!endif

!enddo

!D12 Coulomb integral

!write(6,*) aA, D12_in_in(oo,Ah1,kinh1,Ae,kine,Ah2,kinh2,Ae,kine,aA) + &
!                  D12_out_in(oo,Bh1,kouth1,Ae,kine,Bh2,kouth2,Ae,kine,aA) + &
!                  D12_in_out(oo,Ah1,kinh1,Be,koute,Ah2,kinh2,Be,koute,aA) + &
!                  D12_out_out(oo,Bh1,kouth1,Be,koute,Bh2,kouth2,Be,koute,aA), & 
!                  D12ex_in_in(oo,Ah1,kinh1,Ae,kine,Ah2,kinh2,Ae,kine,aA) + & 
!                  D12ex_out_in(oo,Bh1,kouth1,Be,koute,Ah2,kinh2,Ae,kine,aA) + &
!                  D12ex_in_out(oo,Ah1,kinh1,Ae,kine,Bh2,kouth2,Be,koute,aA) + &
!                  D12ex_out_out(oo,Bh1,kouth1,Be,koute,Bh2,kouth2,Be,koute,aA)
!write(50,*) aA, &!D12_in_in(oo,Ah1,kinh1,Ae,kine,Ah1,kinh1,Ae,kine,aA) + &

!write(6,*) aA, D12dir(oo,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,aA), &
!               D12ex(oo,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,aA), &
!               D12ex_8loops(oo,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,aA)

!  OA= Ana_in_in(m,Ah1,Ae,kinh1,kine,aA) + & 
!      Ana_out_out(m,Bh1,Be,kouth1,koute,aA) 
!
!  OB= Ana_in_in(m,Ah2,Ae,kinh2,kine,aA) + & 
!      Ana_out_out(m,Bh2,Be,kouth2,koute,aA) 

!write(6,*) Ana_in_in(m,Ah1,Ae,kinh1,kine,aA), & 
!      Ana_out_out(m,Bh1,Be,kouth1,koute,aA) 

!write(6,*) Ana_in_in(m,Ah2,Ae,kinh2,kine,aA), & 
!      Ana_out_out(m,Bh2,Be,kouth2,koute,aA) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(E,diffe,diffh,minEe,minEh)

END
