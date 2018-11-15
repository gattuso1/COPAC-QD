module Variables

use Constants

implicit none

   real(dp) :: aA, linker, me, mh, eps, epsout, V0, omegaLO

contains 

subroutine getVariables

NAMELIST /elecSt/ me,mh,eps,epsout,V0,omegaLO
NAMELIST /syst/ aA,linker

open(150,file='WaveFunction1.def',form='formatted')
read(150,NML=elecSt)
read(150,NML=syst)

end subroutine getVariables
   
end module Variables
