program UnitPulse
implicit none

real*8 :: w, c, e0, n, e, P, I, Ener, PnJ, l, EVm, EnJ

w = 1.d-4 !m
c = 299792458 !m/s
e0 = 8.8541878d-15
n = 1.33
e = 10.d-9
EVm = 5.0d9
l = 10.d-15
!P = e/1.d-15

!I = 2*P/(3.1416*w*w)

!Ener = sqrt(2*I/(c*e0*n))

I = (EVm**2 * c * e0 * n)/2

P = (3.1416*w*w*I)/2

EnJ = P * l

write(6,*) EVm, EnJ

end 
