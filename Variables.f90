module Variables

use Constants

implicit none

   character*1 :: o_Norm, o_Over, o_Coul, o_DipS, o_Osci, o_Exti, o_DipD
   real(dp) :: aA, aB, linker, me, mh, eps, epsout, V0, omegaLO, rhoe, rhoh, epsin, V0e, V0h, V0eV, epsR, slope

contains 

subroutine getVariables

NAMELIST /outputs/ o_Norm,o_Over,o_Coul,o_DipS,o_Osci,o_Exti,o_DipD
NAMELIST /elecSt/ me,mh,eps,epsout,V0eV,omegaLO,slope
NAMELIST /syst/ aA,aB,linker

open(150,file='WaveFunction1.def',form='formatted')
read(150,NML=outputs)
read(150,NML=elecSt)
read(150,NML=syst)

me=me*m0
mh=mh*m0
rhoe  = 1.0/sqrt((2*me*omegaLO)/hbar)
rhoh  = 1.0/sqrt((2*mh*omegaLO)/hbar)
epsin = 1.0 + (eps - 1.0) / (1.0 + (0.75e-9/(2*aA))**1.2)
epsR = 1.0/((1.0/epsin)-((1.0/epsin)-(1.0/(epsin+3.5)))*(1-(exp(-(36/35)*aA/rhoe)+exp(-(36/35)*aA/rhoh))/2))
V0 = V0eV*elec
V0e=-1*(-3.49+2.47*(1d9*2*aA)**(-1.32))*elec
V0h=-1*(-5.23-0.74*(1d9*2*aA)**(-0.95))*elec

end subroutine getVariables
   
end module Variables
