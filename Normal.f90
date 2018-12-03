module Normal

use Constants

contains

real(dp) function r8_uniform_01(seed)
  implicit none

  integer(kind=8) k
  integer(kind=4) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01



real(dp) function r8_normal_ab(a,b,seed)
  implicit none

  real(dp) a
  real(dp) b
  real(dp) r1
  real(dp) r2
  integer(kind = 4) seed
  real(dp) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )

  r8_normal_ab = a + b * x

!write(6,*) seed

  return
end

end module Normal
