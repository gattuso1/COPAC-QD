module Variables

use Constants

implicit none

   character*1 :: o_Norm, o_Over, o_Coul, o_DipS, o_Osci, o_Exti, o_DipD, o_DCou
   integer :: ndots, n, nstates, i, j
   real(dp) :: aA, aB, linker, me, mh, eps, epsout, V0, omegaLO, rhoe, rhoh, slope, V0eV, side
   real(dp),allocatable :: aR(:), epsin(:), epsR(:), V0e(:), V0h(:)

contains 

subroutine getVariables

NAMELIST /outputs/ o_Norm,o_Over,o_Coul,o_DipS,o_Osci,o_Exti,o_DipD,o_DCou
NAMELIST /elecSt/ me,mh,eps,epsout,V0eV,omegaLO,slope,side
NAMELIST /syst/ ndots,aA,aB,linker,nstates

open(150,file='WaveFunctiond.def',form='formatted')
read(150,NML=outputs)
read(150,NML=elecSt)
read(150,NML=syst)

allocate(aR(ndots))
allocate(epsin(ndots))
allocate(epsR(ndots))
allocate(V0e(ndots))
allocate(V0h(ndots))

me=me*m0
mh=mh*m0
rhoe  = 1.0/sqrt((2*me*omegaLO)/hbar)
rhoh  = 1.0/sqrt((2*mh*omegaLO)/hbar)
V0 = V0eV*elec

aR(1) = aA
aR(2) = aB

do n=1,ndots
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*aR(n)))**1.2)
epsR(n) = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36/35)*aR(n)/rhoe)+exp(-(36/35)*aR(n)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
enddo

end subroutine getVariables
   
end module Variables
