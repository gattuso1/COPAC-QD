module Constants_au
implicit none 

   integer,  parameter :: dp = SELECTED_REAL_KIND(15,307)
   real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
   real(dp), parameter :: E_au = 5.14d11
   real(dp), parameter :: a0 = 5.291772109217d-11
   real(dp), parameter :: t_au = 2.41888432650516d-17
   real(dp), parameter :: Energ_au = 4.3597441775e-18
   real(dp), parameter :: hbar  = 1
   real(dp), parameter :: cl  = 299792458
   real(dp), parameter :: eps0  = 8.85418781762d-12
   real(dp), parameter :: elec  = 1
   real(dp), parameter :: m0  = 1
   real(dp), parameter :: Cm_to_D  = 3.33564d-30

end module Constants_au 
