module Output

use Constants
use Variables
use Integrals
use Vectors

contains

subroutine makeOutputSingle

open(60,file='Output.txt',form='formatted')

if ( vers .eq. 'singl') then

write(60,*) "Single QD"
write(60,*) 
write(60,*) "The wavefunction of states e, h1 and h2 on QDA has been stored in:"
write(60,*) "The wavefunction of states e, h1 and h2 on QDB has been stored in:"
write(60,*) 
write(60,*) "The initial parameters are:"
write(60,*)
write(60,'("Electron effective mass              :",f8.4," m0")') me
write(60,'("Hole effective mass                  :",f8.4," m0")') mh
write(60,'("Bulk dielectric constant             :",f8.4)')      eps
write(60,'("Ligands dielectric constant          :",f8.4)')      epsout
write(60,'("Bulk band gap                        :",f8.4)')      V0eV
write(60,*)
write(60,*)
write(60,*)
write(60,'(56x,"QDA" , 22x, "QDB")')
write(60,'("Radius                               :",ES25.4E2)') aR(1)
write(60,'("Confining potential electron         :",ES25.4E2)') V0e(1)
write(60,'("Confining potential holes            :",ES25.4E2)') V0h(1)
write(60,'("Dielectric constant                  :",ES25.4E2)') epsin(1)
write(60,'("Transition energy h1e                :",ES25.4E2)') Eeh1(1)/elec
write(60,'("Transition energy h2e                :",ES25.4E2)') Eeh2(1)/elec
write(60,'("Coulomb correction h1e               :",ES25.4E2)') Cb_eh1(1)/elec
write(60,'("Coulomb correction h12               :",ES25.4E2)') Cb_eh2(1)/elec
write(60,'("Norm e                               :",ES25.4E2)') Norm_Ana_e(1)
write(60,'("Norm h1                              :",ES25.4E2)') Norm_Ana_h1(1)
write(60,'("Norm h2                              :",ES25.4E2)') Norm_Ana_h2(1)
write(60,*)
write(60,'("Wavefunction overlap                 :",2f18.6)') OverlapAna_h1e(1), OverlapAna_h2e(1)
write(60,'("Transition dipole moment             :",2f18.6)') TransDip_Ana_h1e(1)/Cm_to_D, TransDip_Ana_h2e(1)/Cm_to_D
write(60,'("Oscillator strength                  :",2f18.6)') Oscillator_Ana_h1e(1), Oscillator_Ana_h2e(1)
write(60,'("Extinction coefficient               :",2f18.6)') ExctCoef_h1e(1), ExctCoef_h2e(1)
write(60,'("Coulomb correction                   :",2f18.6)') Cb_Num_eh1(rmax+1), Cb_Num_eh2(1)
write(60,*) 

endif

end subroutine makeOutputSingle

subroutine makeOutputDimer

open(60,file='Output.txt',form='formatted')

if ( vers .eq. 'dimer') then

if ( aA .eq. aB ) then
write(60,*) "Results for a QD homodimer"
elseif ( aA .ne. aB ) then
write(60,*) "Results for a QD heterodimer"
endif
write(60,*) 
write(60,*) "The wavefunction of states e, h1 and h2 on QDA has been stored in:"
write(60,*) "The wavefunction of states e, h1 and h2 on QDB has been stored in:"
write(60,*) 
write(60,*) "The initial parameters are:"
write(60,*)
write(60,'("Electron effective mass              :",f8.4," m0")') me
write(60,'("Hole effective mass                  :",f8.4," m0")') mh
write(60,'("Bulk dielectric constant             :",f8.4)')      eps
write(60,'("Ligands dielectric constant          :",f8.4)')      epsout
write(60,'("Bulk band gap                        :",f8.4)')      V0eV
write(60,*)
write(60,*)
write(60,*)
write(60,'(56x,"QDA" , 22x, "QDB")')
write(60,'("Radius                               :",ES23.4E2,ES30.4E2)') aR(1) , aR(2)
write(60,'("Confining potential electron         :",ES23.4E2,ES30.4E2)') V0e(1), V0e(2)
write(60,'("Confining potential holes            :",ES23.4E2,ES30.4E2)') V0h(1), V0h(2)
write(60,'("Dielectric constant                  :",ES23.4E2,ES30.4E2)') epsin(1), epsin(2)
write(60,'("Transition energy h1e                :",ES23.4E2,ES30.4E2)') Eeh1(1)/elec, Eeh1(2)/elec
write(60,'("Transition energy h2e                :",ES23.4E2,ES30.4E2)') Eeh2(1)/elec, Eeh2(2)/elec
write(60,'("Coulomb correction h1e               :",ES23.4E2,ES30.4E2)') Cb_eh1(1)/elec, Cb_eh1(2)/elec
write(60,'("Coulomb correction h12               :",ES23.4E2,ES30.4E2)') Cb_eh2(1)/elec, Cb_eh2(2)/elec
write(60,'("Norm e                               :",ES23.4E2,ES30.4E2)') Norm_Ana_e(1), Norm_Ana_e(2)
write(60,'("Norm h1                              :",ES23.4E2,ES30.4E2)') Norm_Ana_h1(1), Norm_Ana_h1(2) 
write(60,'("Norm h2                              :",ES23.4E2,ES30.4E2)') Norm_Ana_h2(1), Norm_Ana_h2(2)
write(60,*)
write(60,'("                                                   /       \                      /      \  ")')
write(60,'("                                                  /         \                    /        \ ")')
write(60,'("                                                 /           \                  /          \")')
write(60,*)
write(60,'(47x,"h1e",12x, "h2e",12x, "h1e",12x, "h2e")')
write(60,'("Wavefunction overlap                 :",4f15.6)') OverlapAna_h1e(1), OverlapAna_h2e(1), &
                                                               OverlapAna_h1e(2), OverlapAna_h2e(2)
write(60,'("Transition dipole moment             :",4ES15.4E2)') 
write(60,'("Transition dipole moment (fit)       :",4ES15.4E2)') TransHam(0,5), TransHam(0,6), TransHam(0,7), TransHam(0,8) 
write(60,'("Oscillator strength                  :",4ES15.4E2)') Oscillator_Ana_h1e(1), Oscillator_Ana_h2e(1), &
                                                                Oscillator_Ana_h1e(2), Oscillator_Ana_h2e(2)
write(60,'("Extinction coefficient               :",4ES15.4E2)') ExctCoef_h1e(1), ExctCoef_h2e(1), ExctCoef_h1e(2), ExctCoef_h2e(2)
write(60,'("Coulomb correction                   :",4ES15.4E2)') Cb_Num_eh1(rmax+1), Cb_Num_eh2(1), Cb_Num_eh1(2), Cb_Num_eh2(2)
write(60,*)
write(60,*)
write(60,*)
write(60,*) "The initial parameters of the dynamic are:"
write(60,*)
write(60,'("Number of states        :",2x,i2)')   nstates  
write(60,'("Number of pulses        :",2x,i2)')   npulses  
if ( npulses .eq. 1 ) then
write(60,'("t0 of pulse             :",ES15.6E2)')   t01      
elseif ( npulses .ge. 2 ) then
write(60,'("t0 of first pulse       :",ES15.6E2)')   t01      
write(60,'("t0 of second pulse      :",ES15.6E2)')   t02      
elseif ( npulses .eq. 3 ) then
write(60,'("t0 of first pulse       :",ES15.6E2)')   t01      
write(60,'("t0 of second pulse      :",ES15.6E2)')   t02      
write(60,'("t0 of third pulse       :",ES15.6E2)')   t03      
endif
write(60,'("Time step               :",ES15.6E2)')   timestep 
write(60,'("Time length of dynamic  :",ES15.6E2)')   totaltime
write(60,'("Omega                   :",ES15.6E2)')   omega    
write(60,'("Phase                   :",ES15.6E2)')   phase   
write(60,'("Width                   :",ES15.6E2)')   width    
write(60,'("Power                   :",ES15.6E2)')   Ed       
write(60,*)
write(60,*)
write(60,*) 
write(60,*) "The 0th order Hamiltonian is:"
write(60,*) 
write(60,'(11x,a2,12x,a2,12x,a2,12x,a2,12x,a2,12x,a2,12x,a2,12x,a2,12x,a2)') "0",  "1", "2", "3", "4", "5", "6", "7", "8"
do i=0,nstates-1
write(60,'(i2,2x,9es14.6e2)') i, (real(Ham(i,j))/elec, j=0,nstates-1)  
enddo

write(60,*) 
write(60,*) "And the transition dipole matrix is:"
write(60,*) 
write(60,'(11x,a2,12x,a2,12x,a2,12x,a2,12x,a2,12x,a2,12x,a2,12x,a2,12x,a2)') "0",  "1", "2", "3", "4", "5", "6", "7", "8"
do i=0,nstates-1
write(60,'(i2,2x,9es14.6e2)') i, (real(TransHam(i,j)/Cm_to_D), j=0,nstates-1)
enddo


endif

end subroutine makeOutputDimer

subroutine makeOutputRange

open(60,file='Output.txt',form='formatted')

write(60,*) "Range of single QD"
write(60,*)
write(60,*) "The wavefunction of states e, h1 and h2 on QDA has been stored in:"
write(60,*) "The wavefunction of states e, h1 and h2 on QDB has been stored in:"
write(60,*)
write(60,*) "The initial parameters are:"
write(60,*)
write(60,'("Electron effective mass              :",f8.4," m0")') me
write(60,'("Hole effective mass                  :",f8.4," m0")') mh
write(60,'("Bulk dielectric constant             :",f8.4)')      eps
write(60,'("Ligands dielectric constant          :",f8.4)')      epsout
write(60,'("Bulk band gap                        :",f8.4)')      V0eV

do n = rmin, rmax
write(11,*) aR(n), linker(n),  Eeh1(n), Eeh2(n)
enddo

end subroutine makeOutputRange

subroutine makeOutputRandm

open(60,file='Output.txt',form='formatted')

open(40,file="Dimer.dat")
open(11,file='Etransitions.dat')

write(40,'("#     Number                  QDA                       QDB                    linker")')
if ( aA .eq. aB) then
write(11,'(a12,18x,a3,23x,a3,20x,a6,22x,a4,22x,a4)') "#     Number", "QDA", "QDB", "linker", "Eh1e", "Eh2e"
else if ( aA .ne. aB) then
write(11,'(a12,18x,a3,23x,a3,20x,a6,22x,a4,21x,a5,21x,a5,21x,a5)') "#     Number", "QDA", "QDB", "linker", "Eh1eA", "Eh2eA", &
                                                                   "Eeh1B", "Eeh2B"
endif

do n = rmin, rmax
if ( aA .eq. aB) then
write(40,*) n, aR(n), aR(n), linker(n)
write(11,*) n, aR(n), aR(n), linker(n), Eeh1(n), Eeh2(n)
else if ( aA .ne. aB) then
write(40,*) n, aR(n), aR(n+nsys), linker(n)
write(11,*) n, aR(n), aR(n+nsys), linker(n), Eeh1(n), Eeh2(n), Eeh1(n+nsys), Eeh2(n+nsys)
endif
enddo

if ( ( vers .eq. 'randm' ) .and. ( aA .eq. aB ) ) then
write(60,*) "Random homodimers"
else if ( ( vers .eq. 'randm' ) .and. ( aA .ne. aB ) ) then
write(60,*) "Random heterodimers"
endif
write(60,*)
write(60,*)
write(60,*) "The process has generated the following files:"
write(60,*)
write(60,'("Dimers.dat            contains radius QDA, radius QDB, linker length")')
write(60,'("Etransitions.dat      contains radius QDA, radius QDB, linker length, Eh1e, Eh2e")')
write(60,'("Ham0.dat              contains Hamiltonian of 0th order")')
write(60,'("Pulse.dat             contains the pulse(s)")')
write(60,'("TransHam.dat          contains the transition dipole moment matrix elements")')
write(60,'("Hamt.dat              contains the time dependent Hamiltonian of dimer n")')
write(60,'("Popc.dat              contains the population evolution of each state |c**2(t)| of dimer n")')
write(60,'("Norm.dat              contains the time dependent norm of dimer n")')
write(60,*)
write(60,*)
write(60,*)
write(60,*) "The initial parameters are:"
write(60,*)
write(60,'("Number of dimers                     :",i5)') nsys
write(60,'("Average radius of QDA                :",f8.4)') aA
write(60,'("Average radius of QDB                :",f8.4)') aB
write(60,'("Average length of linker             :",f8.4)') link
write(60,'("Size dispersion of QDA radius        :",f8.4)') dispQD
write(60,'("Length dispersion of linker          :",f8.4)') displink
write(60,'("Electron effective mass              :",f8.4," m0")') me
write(60,'("Hole effective mass                  :",f8.4," m0")') mh
write(60,'("Bulk dielectric constant             :",f8.4)')      eps
write(60,'("Ligands dielectric constant          :",f8.4)')      epsout
write(60,'("Bulk band gap                        :",f8.4)')      V0eV
write(60,*)
write(60,*)
write(60,*) "The initial parameters of the dynamic are:"
write(60,*)
write(60,'("Number of states        :",2x,i2)')   nstates
write(60,'("Number of pulses        :",2x,i2)')   npulses
if ( npulses .eq. 1 ) then
write(60,'("t0 of pulse             :",ES15.6E2)')   t01
elseif ( npulses .ge. 2 ) then
write(60,'("t0 of first pulse       :",ES15.6E2)')   t01
write(60,'("t0 of second pulse      :",ES15.6E2)')   t02
elseif ( npulses .eq. 3 ) then
write(60,'("t0 of first pulse       :",ES15.6E2)')   t01
write(60,'("t0 of second pulse      :",ES15.6E2)')   t02
write(60,'("t0 of third pulse       :",ES15.6E2)')   t03
endif
write(60,'("Time step               :",ES15.6E2)')   timestep
write(60,'("Time length of dynamic  :",ES15.6E2)')   totaltime
write(60,'("Omega                   :",ES15.6E2)')   omega
write(60,'("Phase                   :",ES15.6E2)')   phase
write(60,'("Width                   :",ES15.6E2)')   width
write(60,'("Power                   :",ES15.6E2)')   Ed
write(60,*)


end subroutine makeOutputRandm

end module Output
