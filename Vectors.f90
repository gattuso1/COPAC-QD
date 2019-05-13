module Vectors

use Constants_au
use omp_lib

contains

!get vector of module d and random orientation
function vector(d)

      implicit none
      real(kind=8), dimension(3) :: vector
      real(kind=8) :: d, phi, teta, rand1, rand2

      call random_number(rand1)
      call random_number(rand2)
      teta = rand1*2.d0*PI
      phi = acos(2*rand2-1.d0)

      vector(1)=0.d0
      vector(2)=0.d0
      vector(3)=0.d0

      vector(1)= d*sqrt(1-cos(phi)**2)*cos(teta)

      vector(2)= d*sqrt(1-cos(phi)**2)*sin(teta)

      vector(3)= d*cos(phi)

!write(6,*) vector(1), vector(2), vector(3)

end function vector

!get vector of random module and random orientation inside volume ddot
function vectorin(ddot)

      implicit none
      real(kind=8), dimension(3) :: vectorin
      real(kind=8) :: ddot, phi, teta, rand1, rand2, rand3, dval

      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      teta = rand1*2*PI
      phi = acos(2*rand2-1)
      dval = rand3*ddot

      vectorin(1)=0
      vectorin(2)=0
      vectorin(3)=0

      vectorin(1)= dval*sqrt(1-cos(phi)**2)*cos(teta)

      vectorin(2)= dval*sqrt(1-cos(phi)**2)*sin(teta)

      vectorin(3)= dval*cos(phi)

end function vectorin

!get vector of random module and random orientation outside volume ddot
function vectorout(ddot, dmax)

      implicit none
      double precision, external:: s13adf
      real(kind=8), dimension(3) :: vectorout
      real(kind=8) :: ddot,  phi, teta, rand1, rand2, dmax, rand3, PI, dval

      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand3)
      teta = rand1*2*PI
      phi = acos(2*rand2-1)
      dval = ddot + rand3 * (dmax-ddot)

      vectorout(1)=0
      vectorout(2)=0
      vectorout(3)=0

      vectorout(1)= dval*sqrt(1-cos(phi)**2)*cos(teta)

      vectorout(2)= dval*sqrt(1-cos(phi)**2)*sin(teta)

      vectorout(3)= dval*cos(phi)

end function vectorout

end module Vectors
