include 'specfun.f90'

program ModelQD_one

use Constants
use Variables
use Integrals
use Vectors

implicit none

double precision, external:: s13adf, ei, eone, nag_bessel_j0

character*64:: arg1, arg2, filename, argA

integer:: i,je,jh,k,l,nsteps,r,rsteps,ifail,raA,raB, m, o, i1, i2, oo, iii, iio, ioi, ioo, rx

real(dp):: maxr0,stepr0,Ef,le,re,lh,rh,delta,x,y,z,D,Dev,lambda, & 
        aex, J , E2, integral, truc, truc2, truc3, a1, a2, a3, OA, OB, dcouplA, &
        InteA, Inth1A, Inth2A, InteB, Inth1B, Inth2B, Overlapeh1A, dcouplB, epsinfrR, &
        C_in_in, C_in_out, C_out_out, prod,  inie, ende, inih, endh, inte, inth, LePo, ra, angl, &
        rand1, rand2, start, finish, mu, A, epsinf, disteh, &
        Rine, Route, Rinh1, Routh1, Rinh2, Routh2, RadProbe, RadProbh1, RadProbh2, &
        r0, Ae, Ah1, Ah2, Be, Bh1, Bh2, &
        Eeh1, Eeh2, alphae, alphah1, alphah2, betae, betah1, betah2, &
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
o=      1000000

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

diffe(:)=0
diffh(:)=0
E(:)=0
minEe(:)=0 
minEh(:)=0 
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

!write(10,*) E(i), sqrt(2*me*E(i))/hbar * aA * 1/tan(sqrt(2*me*E(i))/hbar * aA), - 1 + (me/m0) + (me*aA)/(hbar) &
!           * sqrt((2/m0)*(V0e-E(i)))

if (diffe(i) .le. diffe(i-1)) then
        minEe(je) = E(i)
elseif ( (diffe(i) .ge. diffe(i-1)) .and. (E(i-1) .eq. minEe(je)) ) then
        je=je+1
endif

!diffh(i) = abs(sqrt(2*mh*E(i))/hbar * aA * 1/tan(sqrt(2*mh*E(i))/hbar * aA) - 1 + (mh/m0) + (mh*aA)/(hbar) &
!           * sqrt((2/m0)*(V0h-E(i))))
!
!if (diffh(i) .le. diffh(i-1)) then
!        minEh(jh) = E(i)
!elseif ( (diffh(i) .ge. diffh(i-1)) .and. (E(i-1) .eq. minEh(jh)) ) then
!        jh=jh+1
!endif

enddo

!write(6,*) V0h, m0, minEe(1),minEe(2), hbar

!wave vectors in and out

kine=sqrt(2*me*minEe(1))/hbar
koute=sqrt(2*m0*(V0e-minEe(1)))/hbar

kinh1=sqrt(2*mh*minEh(1))/hbar
kouth1=sqrt(2*m0*(V0h-minEh(1)))/hbar

kinh2=sqrt(2*mh*minEh(2))/hbar
kouth2=sqrt(2*m0*(V0h-minEh(2)))/hbar

alphae=kine*aA
alphah1=kinh1*aA
alphah2=kinh2*aA

betae=koute*aA
betah1=kouth1*aA
betah2=kouth2*aA

!normalization factors

Ae=1/sqrt(aA/2-sin(2*kine*aA)/(4*kine)+sin(kine*aA)**2/(2*koute))
Be=Ae*sin(kine*aA)*exp(koute*aA) !*(koute/kine)

Ah1=1/sqrt(aA/2-sin(2*kinh1*aA)/(4*kinh1)+sin(kinh1*aA)**2/(2*kouth1))
Bh1=Ah1*sin(kinh1*aA)*exp(kouth1*aA) !*(koute/kine)

Ah2=1/sqrt(aA/2-sin(2*kinh2*aA)/(4*kinh2)+sin(kinh2*aA)**2/(2*kouth2))
Bh2=Ah2*sin(kinh2*aA)*exp(kouth2*aA) !*(koute/kine)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Comlomb correction (analytical) between electron and holes

I1eh1=(Ae**2*Ah1**2*aA/4)*(1-sin(2*alphae)/(2*alphae)-s13adf(2*alphae,ifail)/(2*alphae)+&
         (s13adf(2*alphae-2*alphah1,ifail)+s13adf(2*alphae+2*alphah1,ifail))/(4*alphae)) + &
         (Ae**2*Ah1**2*aA/4)*(1-sin(2*alphah1)/(2*alphah1)-s13adf(2*alphah1,ifail)/(2*alphah1)+&
         (s13adf(2*alphah1-2*alphae,ifail)+s13adf(2*alphah1+2*alphae,ifail))/(4*alphah1)) 

I2eh1=(Ae**2*Bh1**2*aA/2)*(1-sin(2*alphae)/(2*alphae))*eone(2*betah1)+&
         (Ah1**2*Be**2*aA/2)*(1-sin(2*alphah1)/(2*alphah1))*eone(2*betae)

I3eh1=(Be**2*Bh1**2*aA/(2*betae)*(exp(-2*betae)*eone(2*betah1)-eone(2*betae+2*betah1)))+&
         (Be**2*Bh1**2*aA/(2*betah1)*(exp(-2*betah1)*eone(2*betae)-eone(2*betah1+2*betae)))

I1eh2=(Ae**2*Ah2**2*aA/4)*(1-sin(2*alphae)/(2*alphae)-s13adf(2*alphae,ifail)/(2*alphae)+&
         (s13adf(2*alphae-2*alphah2,ifail)+s13adf(2*alphae+2*alphah2,ifail))/(4*alphae)) + &
         (Ae**2*Ah2**2*aA/4)*(1-sin(2*alphah2)/(2*alphah2)-s13adf(2*alphah2,ifail)/(2*alphah1)+&
         (s13adf(2*alphah2-2*alphae,ifail)+s13adf(2*alphah2+2*alphae,ifail))/(4*alphah2)) 

I2eh2=(Ae**2*Bh2**2*aA/2)*(1-sin(2*alphae)/(2*alphae))*eone(2*betah2)+&
         (Ah2**2*Be**2*aA/2)*(1-sin(2*alphah2)/(2*alphah2))*eone(2*betae)

I3eh2=(Be**2*Bh2**2*aA/(2*betae)*(exp(-2*betae)*eone(2*betah2)-eone(2*betae+2*betah2)))+&
         (Be**2*Bh2**2*aA/(2*betah2)*(exp(-2*betah2)*eone(2*betae)-eone(2*betah2+2*betae)))


Eeh1=(elec**2/(4*pi*epsin*eps0))*(I1eh1+I2eh1+I3eh1)/elec

Eeh2=(elec**2/(4*pi*eps*eps0))*(I1eh2+I2eh2+I3eh2)/elec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(20,*) aA, minEe(1)/elec, minEh(1)/elec, minEh(2)/elec
write(11,*) aA, (minEe(1)+minEh(1))/elec+V0eV
write(12,*) aA, (minEe(1)+minEh(2))/elec+V0eV

write(6,*) Eeh1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Comlomb correction (numerical) between electron and holes
!

!oo=15
!
!do r=1,rsteps
!
!aA=stepr0*r
!
!epsinf = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*aA))**1.2)
!!disteh = (36.0/35.0)*aA
!!epsinfrR= 1.0/((1.0/epsinf)-((1.0/epsinf)-(1.0/(epsinf+3.5)))*(1-(exp(-disteh/rhoe)+exp(-disteh/rhoh))/2))
!
!write(23,*) aA , (minEe(1)+minEh(1))/elec+V0eV - &
!                  (elec/(4*pi*epsinf*eps0))*C_num_in_in(oo, Ae, Ah1, kine, kinh1, aA) , &
!                  (minEe(1)+minEh(2))/elec+V0eV - &
!                  (elec/(4*pi*epsinf*eps0))*C_num_in_in(oo, Ae, Ah2, kine, kinh2, aA) !, &
!                   !(elec/(4*pi*eps*eps0))*C_num_in_out(oo, Ae, Bh1, kine, kouth1, r), &
!                   !(elec/(4*pi*eps*eps0))*C_num_out_in(oo, Be, Ah1, koute, kinh1, r), &
!                   !(elec/(4*pi*eps*eps0))*C_num_out_out(oo, Be, Bh1, koute, kouth1, r)
!
!enddo

!Comlomb correction (numerical) between electron and holes with size dependent dielectric constant
!
!do r=1,200
!oo=1000
!
!write(6,*) 1+r*0.1, (elec/(4*pi*eps0))*C_num_in_in_R(oo, Ae(100+r*10), Ah1(100+r*10), kine(100+r*10), kinh1(100+r*10), r)
!                    (elec/(4*pi*eps0))*C_num_in_out_R(oo, Ae(100+r*10), Bh1(100+r*10), kine(100+r*10), kouth1(100+r*10), r), &
!                    (elec/(4*pi*eps0))*C_num_out_in_R(oo, Be(100+r*10), Ah1(100+r*10), koute(100+r*10), kinh1(100+r*10), r), &
!                    (elec/(4*pi*eps0))*C_num_out_out_R(oo, Be(100+r*10), Bh1(100+r*10), koute(100+r*10), kouth1(100+r*10), r)
!
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(13,*) aA, (minEe(1)+minEh(1))/elec+V0eV-Eeh1
write(14,*) aA, (minEe(1)+minEh(2))/elec+V0eV-Eeh2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Wave function QDA and radial probability

do r=1,rsteps

Rine=Ae*sin(kine*r0)/r0
Route=Be*exp(-1*koute*r0)/r0

Rinh1=Ah1*sin(kinh1*r0)/r0
Routh1=Bh1*exp(-1*kouth1*r0)/r0

Rinh2=Ah2*sin(kinh2*r0)/r0
Routh2=Bh2*exp(-1*kouth2*r0)/r0

if (aA .le. aA) then
        RadProbe=Rine
        RadProbh1=Rinh1
        RadProbh2=Rinh2
else if (aA .gt. aA) then
        RadProbe=Route
        RadProbh1=Routh1
        RadProbh2=Routh2
endif

write(17,*) r0, RadProbe, RadProbh1, RadProbh2
write(15,*) r0, aA**2*RadProbe**2/(4*pi) , aA**2*RadProbh1**2/(4*pi), aA**2*RadProbh2**2/(4*pi)

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Normalization integrals

!write(6,*) -Ae(raA)**2*(sin(2*kine(raA)*aA)-2*kine(raA)*aA)/(4*kine(raA)) + Be(raA)**2*exp(-2*koute(raA)*aA)/(2*koute(raA))
!write(6,*) -Ah1(raA)**2*(sin(2*kinh1(raA)*aA)-2*kinh1(raA)*aA)/(4*kinh1(raA)) + Bh1(raA)**2*exp(-2*kouth1(raA)*aA)/(2*kouth1(raA))
!write(6,*) -Ah2(raA)**2*(sin(2*kinh2(raA)*aA)-2*kinh2(raA)*aA)/(4*kinh2(raA)) + Bh2(raA)**2*exp(-2*kouth2(raA)*aA)/(2*kouth2(raA))
!write(6,*) -Ae(raB)**2*(sin(2*kine(raB)*aB)-2*kine(raB)*aB)/(4*kine(raB)) + Be(raB)**2*exp(-2*koute(raB)*aB)/(2*koute(raB))
!write(6,*) -Ah1(raB)**2*(sin(2*kinh1(raB)*aB)-2*kinh1(raB)*aB)/(4*kinh1(raB)) + Bh1(raB)**2*exp(-2*kouth1(raB)*aB)/(2*kouth1(raB))
!write(6,*) -Ah2(raB)**2*(sin(2*kinh2(raB)*aB)-2*kinh2(raB)*aB)/(4*kinh2(raB)) + Bh2(raB)**2*exp(-2*kouth2(raB)*aB)/(2*kouth2(raB))


!InteB=Norm(m,Ae(raB),Be(raB),kine(raB),koute(raB),aB)
!Inth1B=Norm(m,Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),aB)
!Inth2B=Norm(m,Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),aB)

!write(6,*) Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB
!write(6,*) m , Norm_cart(m,Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB)
!write(60,*) 'Rdm', Norm_cart_Rdm(Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB)
!write(60,*) 'MC' , Norm_cart_MC(m,Ae(raA),Be(raA),kine(raA),koute(raA),aA,aB)


!InteA=Norm(m,Ae(raA),Be(raA),kine(raA),koute(raA),aA)
!Inth1A=Norm(m,Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),aA)
!Inth2A=Norm(m,Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),aA)

!write(6,*) InteA, Inth1A, Inth2A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!RDF overlap one dot

!write(50,*) "eA - h1A num", Overlapinin(m,Ah1(raA),Ae(raA),kinh1(raA),kine(raA),aA) +  &
!                           Overlapoutout(m,Bh1(raA),Be(raA),kouth1(raA),koute(raA),aA)
!write(6,*) "eA - h1A ana", abs(Oin(Ae(raA),Ah1(raA),kine(raA),kinh1(raA),aA) , Oout(Be(raA),Bh1(raA),koute(raA),kouth1(raA),aA))
!write(6,*) "eA - h1A ana", abs(Oin(Ae(raA),Ah1(raA),kine(raA),kinh1(raA),aA) - Oout(Be(raA),Bh1(raA),koute(raA),kouth1(raA),aA))
!write(6,*) "eA - h2A ana", abs(Oin(Ae(raA),Ah2(raA),kine(raA),kinh2(raA),aA) , Oout(Be(raA),Bh2(raA),koute(raA),kouth2(raA),aA))
!write(6,*) "eA - h2A ana", abs(Oin(Ae(raA),Ah2(raA),kine(raA),kinh2(raA),aA) - Oout(Be(raA),Bh2(raA),koute(raA),kouth2(raA),aA))
!write(6,*) "eB - h1B ana", abs(Oin(Ae(raB),Ah1(raB),kine(raB),kinh1(raB),aB) , Oout(Be(raB),Bh1(raB),koute(raB),kouth1(raB),aB))
!write(6,*) "eB - h1B ana", abs(Oin(Ae(raB),Ah1(raB),kine(raB),kinh1(raB),aB) - Oout(Be(raB),Bh1(raB),koute(raB),kouth1(raB),aB))
!write(6,*) "eB - h2B ana", abs(Oin(Ae(raB),Ah2(raB),kine(raB),kinh2(raB),aB) , Oout(Be(raB),Bh2(raB),koute(raB),kouth2(raB),aB))
!write(6,*) "eB - h2B ana", abs(Oin(Ae(raB),Ah2(raB),kine(raB),kinh2(raB),aB) - Oout(Be(raB),Bh2(raB),koute(raB),kouth2(raB),aB))
!write(6,*) raA, raB, linker, raA*0.01e-9, raB*0.01e-9

m=10000

!write(6,*) r, aA

!Overlap 

!write(6,*) aA, Overlapoutout(m,Bh1,Be,kouth1,koute,aA)

!Dipole moment h1-e same dot
!write(6,*) aA, TransDipIn(m,Ah1,Ae,kinh1,kine,aA) +&
!                   TransDipOut(m,Bh1,Be,kouth1,koute,aA),&
                   !TransDip_cart(m,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,aA)
!                   TransDip_Ana(Ah1,kinh1,Ae,kine,Bh1,kouth1,Be,koute,aA)

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

!aA = aA
!aB = aA

!write(6,*) aA, linker, TransDip_dimer_MC(Ah1,Bh1,kinh1,kouth1,&
!                                           Ae,Be,kine,koute,aA,aA,linker)


!write(6,*) aA, linker, TransDip_dimer_cart(m,Ah1,Bh1,kinh1,kouth1,&
!                                              Ae,Be,kine,koute,aA,aA,linker)
!                 TransDip_dimer_cart(m,Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),&
!                                              Ae(raB),Be(raB),kine(raB),koute(raB),aA,aA,linker)
!                    TransDip_dimer_cart(m,Ae(raA),Be(raA),kine(raA),koute(raA),&
!                                              Ah1(raB),Bh1(raB),kinh1(raB),kouth1(raB),aA,aA,linker),&
!                    TransDip_dimer_cart(m,Ae(raA),Be(raA),kine(raA),koute(raA),&
!                                              Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),aA,aA,linker)

!enddo
!write(6,*) "         "


!enddo

!Dipole moment h2-e same dot
!write(6,*) aA ,TransDipIn(m,Ah2,Ae,kinh2,kine,aA) +&
!                   TransDipOut(m,Bh2,Be,kouth2,koute,aA),&
!                   TransDip_cart(m,Ah2,Bh2,kinh2,kouth2,Ae,Be,kine,koute,aA)
!                   TransDip_Ana(Ah2,kinh2,Ae,kine,Bh2,kouth2,Be,koute,aA)

!write(23,*) linker,TransDip_dimer_MC(Ah1(raA),Bh1(raA),kinh1(raA),kouth1(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker),&
!                 TransDip_dimer_MC(Ah2(raA),Bh2(raA),kinh2(raA),kouth2(raA),Ae(raB),Be(raB),kine(raB),koute(raB),aA,aB,linker)
!write(6,*) linker , TransDip_dimer_MC(Ah2(raB),Bh2(raB),kinh2(raB),kouth2(raB),Ae(raA),Be(raA),kine(raA),koute(raA),aB,aA,linker)

!enddo

!Dipole moment from EMA equation
!write(6,*) aA, sqrt(((elec**2*hbar**2)/(6*(minEe(1)+minEh(1)+V0-Eeh1*1.60218e-19)**2*m0))*&
!       ((m0/me)-1)*((V0*(V0+(0.42*1.60218e-19)))/(V0+(2/3)*(0.42*1.60218e-19)))) ,&
!        sqrt(((elec**2*hbar**2)/(6*(minEe(1)+minEh(2)+V0-Eeh2*1.60218e-19)**2*m0))*&
!        ((m0/me)-1)*((V0*(V0+(0.42*1.60218e-19)))/(V0+(2/3)*(0.42*1.60218e-19))))

!Oscilator strength
!write(6,*) aA, (minEe(1)+minEh(1)+V0-Eeh1*1.60218e-19), &
!                ((2*me)/(3*elec**2*hbar**2))*(minEe(1)+minEh(1)+V0-Eeh1*1.60218e-19)*&
!           (TransDipIn(m,Ah1,Ae,kinh1,kine,aA) + TransDipOut(m,Bh1,Be,kouth1,koute,aA))**2,&
!              ((2*me)/(3*elec**2*hbar**2))*(minEe(1)+minEh(2)+V0-Eeh1*1.60218e-19)*&
!           (TransDipIn(m,Ah2,Ae,kinh2,kine,aA) + TransDipOut(m,Bh2,Be,kouth2,koute,aA))**2

!molar extinction coefficient
!A=1.25e5*elec**2
!mu=(me*mh)/(me+mh)
!write(6,*) aA,(minEe(1)+minEh(1)+V0-Eeh1*1.60218e-19), &
        !(pi*6.022e23*elec**2*20*1.60218e-19*(2/3)*hbar)/(2*log(10)*me*3.0e8*3.06*eps0*&
          ! (minEe(1)+minEh(1)+V0-Eeh1*1.60218e-19)*sqrt(2*pi)*

!        (A*mu*aA**2)/(sqrt(2*pi)*(minEe(1)+minEh(1)+V0-Eeh1*1.60218e-19)*pi*hbar**2*0.1)
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
write(6,*) Ah2,Bh2,kinh2,kouth2

write(6,*) aA, D12dir(oo,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah2,Bh2,kinh2,kouth2,Ae,Be,kine,koute,aA), &
               D12ex(oo,Ah1,Bh1,kinh1,kouth1,Ae,Be,kine,koute,Ah2,Bh2,kinh2,kouth2,Ae,Be,kine,koute,aA)
!               D12_out_in(oo,Bh1,kouth1,Ae,kine,Bh1,kouth1,Ae,kine,aA) + &
!               D12_in_out(oo,Ah1,kinh1,Be,koute,Ah1,kinh1,Be,koute,aA) + &
!               D12_out_out(oo,Bh1,kouth1,Be,koute,Bh1,kouth1,Be,koute,aA)
!               D12ex_in_in(oo,Ah1,kinh1,Ae,kine,Ah1,kinh1,Ae,kine,aA) !+ & 
               !D12ex_out_in(oo,Bh1,kouth1,Be,koute,Ah1,kinh1,Ae,kine,aA) + &
               !D12ex_in_out(oo,Ah1,kinh1,Ae,kine,Bh1,kouth1,Be,koute,aA) + &
               !D12ex_out_out(oo,Bh1,kouth1,Be,koute,Bh1,kouth1,Be,koute,aA)

!  OA= Ana_in_in(m,Ah1(raA),Ae(raA),kinh1(raA),kine(raA),aA) + & 
!      Ana_out_out(m,Bh1(raA),Be(raA),kouth1(raA),koute(raA),aA) 
!
!  OB= Ana_in_in(m,Ah2(raA),Ae(raA),kinh2(raA),kine(raA),aA) + & 
!      Ana_out_out(m,Bh2(raA),Be(raA),kouth2(raA),koute(raA),aA) 

!write(6,*) Ana_in_in(m,Ah1(raA),Ae(raA),kinh1(raA),kine(raA),aA), & 
!      Ana_out_out(m,Bh1(raA),Be(raA),kouth1(raA),koute(raA),aA) 

!write(6,*) Ana_in_in(m,Ah2(raA),Ae(raA),kinh2(raA),kine(raA),aA), & 
!      Ana_out_out(m,Bh2(raA),Be(raA),kouth2(raA),koute(raA),aA) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(E,diffe,diffh,minEe,minEh)

END
