module Variables

use Constants_au
use Normal

implicit none

   character*5 :: vers
   character*64 :: popc, hmti, norm, tdmM, hmt0, outputdir
   character*1 :: o_Norm, o_Over, o_Coul, o_DipS, o_Osci, o_Exti, o_DipD, dyn, hamilt
   integer :: ndots, n, rmin, rmax, nsys, npulses, nstates, ntime, i, j, t, lwork, info
   integer,allocatable :: seed(:)
   real(dp) :: aA, aB, me, mh, eps, epsout, V0, omegaLO, rhoe, rhoh, slope, V0eV, minr, maxr, rsteps, side
   real(dp) :: dispQD, displink, rdmlinker, rdmQDA, rdmQDB, link, t01, t02, t03, timestep, totaltime, omega, phase, width, Ed
   real(dp) :: pulse1, pulse2, pulse3, test
   real(dp),allocatable :: aR(:), aRA(:), aRB(:), epsin(:), epsR(:), V0e(:), V0h(:), linker(:)
   real(dp),allocatable :: epsinA(:), epsinB(:), epsRA(:), epsRB(:), V0eA(:), V0eB(:), V0hA(:), V0hB(:)
   real(dp),allocatable :: Eeh1(:), Eeh2(:), Cb_eh1(:), Cb_eh2(:), Norm_Ana_e(:), Norm_Ana_h1(:), Norm_Ana_h2(:)
   real(dp),allocatable :: OverlapAna_h1e(:), OverlapAna_h2e(:), Cb_Num_eh1(:), Cb_Num_eh1_eh2(:), Cb_Num_eh2(:)
   real(dp),allocatable :: minEe(:,:),minEh(:,:), TransDip_Num_h1e(:), TransDip_Num_h2e(:), work(:), lambda(:)
   real(dp),allocatable :: TransDip_Ana_h1e(:), TransDip_Ana_h2e(:), Oscillator_Ana_h1e(:), Oscillator_Ana_h2e(:)
   real(dp),allocatable :: ExctCoef_h1e(:), ExctCoef_h2e(:), Ham(:,:), E0(:), c0(:), TransHam(:,:), Hamt(:,:,:), c(:,:)
   complex(dp) :: ct1, ct2, ct3, ct4, xt01, xt02, xt03, xhbar, im, xwidth, xomega , xEd, xh, xphase, xtime
   complex(dp),allocatable :: xHam(:,:) , xHamt(:,:,:), xTransHam(:,:), xE0(:), xHamtk2(:,:,:), xHamtk3(:,:,:), xHamtk4(:,:,:)
   complex(dp),allocatable :: xc0(:), xc(:,:), xc_ei(:,:), xcnew(:,:), k1(:), k2(:), k3(:) , k4(:), xHam_ei(:,:)

contains 

subroutine getVariables

NAMELIST /version/ vers
NAMELIST /outputs/ o_Norm,o_Over,o_Coul,o_DipS,o_Osci,o_Exti,o_DipD,dyn,hamilt
NAMELIST /elecSt/ me,mh,eps,epsout,V0eV,omegaLO,slope,side
NAMELIST /pulses/ nstates,npulses,t01,t02,t03,timestep,totaltime,omega,phase,width,Ed
NAMELIST /syst_single/ aA
NAMELIST /syst_dimer/ aA,aB,link
NAMELIST /syst_range/ rsteps,minr,maxr,link
NAMELIST /syst_random/ nsys,aA,aB,link,displink,dispQD

open(150,file='QD_quest.def',form='formatted')
read(150,NML=version)
read(150,NML=outputs)
read(150,NML=elecSt)

im         = dcmplx(0.0d0,1.0d0)
me         = me*m0
mh         = mh*m0
rhoe       = 1.0d0/sqrt((2.d0*me*omegaLO)/hbar)
rhoh       = 1.0d0/sqrt((2.d0*mh*omegaLO)/hbar)
V0         = V0eV*elec

if ( dyn .eq. 'y' ) then

rewind 150
read(150,NML=pulses)

xhbar      = dcmplx(hbar,0.0d0)
xh         = dcmplx(timestep,0.0d0)
xt01       = dcmplx(t01,0.0d0)
xt02       = dcmplx(t02,0.0d0)
xt03       = dcmplx(t03,0.0d0)
xphase     = dcmplx(phase,0.0d0)
xomega     = dcmplx(omega,0.0d0)
xwidth     = dcmplx(width,0.0d0)
xEd        = dcmplx(Ed,0.0d0)
ntime      = nint(totaltime/timestep)
phase      = pi

allocate(Hamt(0:nstates-1,0:nstates-1,0:ntime))
allocate(xHamt(0:nstates-1,0:nstates-1,0:ntime))
allocate(xHamtk2(0:nstates-1,0:nstates-1,0:ntime))
allocate(xHamtk3(0:nstates-1,0:nstates-1,0:ntime))
allocate(xHamtk4(0:nstates-1,0:nstates-1,0:ntime))
allocate(c(0:nstates-1,0:ntime))
allocate(xc(0:nstates-1,0:ntime))
allocate(xc_ei(0:nstates-1,0:ntime))
allocate(xcnew(0:nstates-1,ntime+1))
allocate(TransHam(0:nstates-1,0:nstates-1))
allocate(xTransHam(0:nstates-1,0:nstates-1))
allocate(Ham(0:nstates-1,0:nstates-1))
allocate(xHam(0:nstates-1,0:nstates-1))
allocate(xHam_ei(0:nstates-1,0:nstates-1))
allocate(c0(0:nstates-1))
allocate(k1(0:nstates-1))
allocate(k2(0:nstates-1))
allocate(k3(0:nstates-1))
allocate(k4(0:nstates-1))

if ( npulses .eq. 3) then
pulse1 = 1.d0 
pulse2 = 1.d0
pulse3 = 1.d0
elseif ( npulses .eq. 2) then
pulse1 = 1.d0
pulse2 = 1.d0
pulse3 = 0.d0
elseif ( npulses .eq. 1) then
pulse1 = 1.d0
pulse2 = 0.d0
pulse3 = 0.d0
elseif ( npulses .eq. 0) then
pulse1 = 0.d0
pulse2 = 0.d0
pulse3 = 0.d0
endif

endif

if ( vers .eq. 'randm' ) then
rewind 150
read(150,NML=syst_random)
endif


if ( vers .eq. 'singl' ) then

write(6,*) "You are requesting me to tackle a single QD"

allocate(aR(1))
allocate(linker(1))
allocate(epsin(1))
allocate(epsR(1))
allocate(V0e(1))
allocate(V0h(1))

rewind 150
read(150,NML=syst_single)

rmin= 1
rmax = 1

aR(1) = aA
epsin(1) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(1)))**1.2)
epsR(1) = 1.0/((1.0/epsin(1))-((1.0/epsin(1))-(1.0/(epsin(1)+3.5)))*(1-(exp(-(36/35)*aR(1)/rhoe)+exp(-(36/35)*aR(1)/rhoh))/2))
V0e(1)=-1*(-3.49+2.47*(1d9*2*aR(1))**(-1.32))*elec
V0h(1)=-1*(-5.23-0.74*(1d9*2*aR(1))**(-0.95))*elec

else if ( vers .eq. 'dimer' ) then

allocate(aR(2))
allocate(linker(2))
allocate(epsin(2))
allocate(epsR(2))
allocate(V0e(2))
allocate(V0h(2))

rewind 150
read(150,NML=syst_dimer)

if ( aA .eq. aB ) then
write(6,*) "You are requesting me to tackle a homodimer QD"
else if ( aA .ne. aB ) then
write(6,*) "You are requesting me to tackle a heterodimer QD"
endif

aR(1) = aA
aR(2) = aB
linker(:) = link

rmin = 1
rmax = 2

do n = 1,2
epsin(n) = 1.d0 + (eps - 1.d0) / (1.d0 + (0.75d-9/(2*aR(n)))**1.2d0)
epsR(n) = 1.d0/((1.d0/epsin(n))-((1.d0/epsin(n))-(1.d0/(epsin(n)+3.5d0)))*&
(1.d0-(exp(-(36.d0/35.d0)*aR(n)/rhoe)+exp(-(36.d0/35.d0)*aR(n)/rhoh))/2.d0))
V0e(n)=-1.d0*(-3.49d0+2.47d0*(1d9*2.d0*aR(n))**(-1.32d0))*elec
V0h(n)=-1.d0*(-5.23d0-0.74d0*(1d9*2.d0*aR(n))**(-0.95d0))*elec
enddo

else if ( vers .eq. 'range' ) then

write(6,*) "You are requesting me to tackle a range of QD"

rewind 150
read(150,NML=syst_range)

rmin = minr/rsteps
rmax = maxr/rsteps
nsys= maxr/rsteps - minr/rsteps

allocate(aR(rmin:rmax))
allocate(linker(rmin:rmax))
allocate(epsin(rmin:rmax))
allocate(epsR(rmin:rmax))
allocate(V0e(rmin:rmax))
allocate(V0h(rmin:rmax))

linker(:) = link

do n = rmin,rmax
aR(n)= n*rsteps
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n)))**1.2)
epsR(n) = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36/35)*aR(n)/rhoe)+exp(-(36/35)*aR(n)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
enddo

else if ( ( vers .eq. 'randm' ) .and. ( aA .eq. aB ) ) then

write(6,*) "You are requesting me to tackle a random set of homodimer QD"

allocate(aR(nsys))
allocate(linker(nsys))
allocate(epsin(nsys))
allocate(epsR(nsys))
allocate(V0e(nsys))
allocate(V0h(nsys))

  call random_seed(size = n)
  allocate(seed(n))
  call random_seed(get=seed)

rmin = 1
rmax = nsys

do n = 1,nsys
aR(n)= r8_NORMAL_AB(aA, dispQD*1d-9, seed(1))
linker(n) = r8_NORMAL_AB(link, displink*1d-9, seed(2))
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n)))**1.2)
epsR(n) = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36/35)*aR(n)/rhoe)+exp(-(36/35)*aR(n)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
!write(6,*) aR(n), linker(n)
enddo

else if ( ( vers .eq. 'randm' ) .and. ( aA .ne. aB ) ) then

write(6,*) "You are requesting me to tackle a random set of heterodimer QD"

allocate(aR(2*nsys))
allocate(linker(2*nsys))
allocate(epsin(2*nsys))
allocate(epsR(2*nsys))
allocate(V0e(2*nsys))
allocate(V0h(2*nsys))

linker = link

  call random_seed(size = n)
  allocate(seed(n))
  call random_seed(get=seed)

rmin = 1
rmax = 2*nsys

do n = 1,nsys
aR(n) = r8_NORMAL_AB(aA,dispQD*1d-9,seed(1))
aR(n+nsys) = r8_NORMAL_AB(aB,dispQD*1d-9,seed(2))
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n)))**1.2)
epsin(n+nsys) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n+nsys)))**1.2)
epsR(n)= 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36/35)*aR(n)/rhoe)+exp(-(36/35)*aR(n)/rhoh))/2))
epsR(n+nsys)= 1.0/((1.0/epsin(n+nsys))-((1.0/epsin(n+nsys))-(1.0/(epsin(n+nsys)+3.5)))*&
                  (1-(exp(-(36/35)*aR(n+nsys)/rhoe)+exp(-(36/35)*aR(n+nsys)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0e(n+nsys)=-1*(-3.49+2.47*(1d9*2*aR(n+nsys))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
V0h(n+nsys)=-1*(-5.23-0.74*(1d9*2*aR(n+nsys))**(-0.95))*elec
write(6,*) aR(n), aR(n+nsys)
enddo

endif 

end subroutine getVariables
   
end module Variables
