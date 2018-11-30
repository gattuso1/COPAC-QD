module Output

use Constants
use Variables
use Integrals
use Vectors

contains

subroutine makeOutputSingle

if ( vers .eq. 'singl') then
open(60,file='Output_dimer.txt',form='formatted')

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

if ( vers .eq. 'dimer') then
open(60,file='Output_dimer.txt',form='formatted')

write(60,*) "DIMER OF QD"
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
write(60,'("                                                  /          \                  /          \")')
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
write(60,'("Number of pulse         :",2x,i2)')   npulses  
write(60,'("t0 of first pulse       :",ES15.6E2)')   t01      
write(60,'("t1 of first pulse       :",ES15.6E2)')   t02      
write(60,'("t2 of first pulse       :",ES15.6E2)')   t03      
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

subroutine makeOutputRdmho

write(60,*) "Random homodimer"
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

end subroutine makeOutputRdmho

subroutine makeOutputRdmhe

write(60,*) "Random heterodimer"
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

end subroutine makeOutputRdmhe

end module Output
