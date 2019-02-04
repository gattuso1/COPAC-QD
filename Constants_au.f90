module Constants_au

use omp_lib

implicit none 

   integer,  parameter :: dp = SELECTED_REAL_KIND(15,307)
   real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
   real(dp), parameter :: E_au = 5.1422065211d11
   real(dp), parameter :: Dip_au = 8.4783532619d-30
   real(dp), parameter :: a0 = 5.291772109217d-11
   real(dp), parameter :: t_au = 2.41888432650516d-17
   real(dp), parameter :: Energ_au = 4.3597441775d-18
   real(dp), parameter :: h  = 6.62607004d-34
   real(dp), parameter :: hbar_au  = 1
   real(dp), parameter :: hbar = h/(2.0d0*pi)
   real(dp), parameter :: cl  = 299792458
   real(dp), parameter :: eps0  = 8.85418781762d-12
   real(dp), parameter :: elec  = 1.60217662d-19
   real(dp), parameter :: elec_au  = 1
   real(dp), parameter :: m0  = 9.10938356d-31
   real(dp), parameter :: Cm_to_D  = 3.33564d-30
   real(dp), parameter :: D_to_au  = 2.5417462310548435 

end module Constants_au 
