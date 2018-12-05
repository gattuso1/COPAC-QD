module Constants  
implicit none 

   integer,  parameter :: dp = SELECTED_REAL_KIND(16,307)
   real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
   real(dp), parameter :: h  = 6.62607004d-34
   real(dp), parameter :: hbar  = h/(2.0d0*pi)
   real(dp), parameter :: cl  = 299792458
   real(dp), parameter :: eps0  = 8.85418781762d-12
   real(dp), parameter :: elec  = 1.60217662d-19
   real(dp), parameter :: m0  = 9.10938356d-31
   real(dp), parameter :: Cm_to_D  = 3.33564d-30

end module Constants 
