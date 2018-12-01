module Variables

use Constants

implicit none

   character*5 :: vers
   character*1 :: o_Norm, o_Over, o_Coul, o_DipS, o_Osci, o_Exti, o_DipD, dyn
   integer :: ndots, n, rmin, rmax, nsys, npulses, nstates, ntime, i, j, t
   real(dp) :: aA, aB, me, mh, eps, epsout, V0, omegaLO, rhoe, rhoh, slope, V0eV, minr, maxr, rsteps, side
   real(dp) :: dispQD, displink, rdmlinker, rdmQDA, rdmQDB, link, t01, t02, t03, timestep, totaltime, omega, phase, width, Ed
   real(dp) :: pulse1, pulse2, pulse3
   real(dp),allocatable :: aR(:), aRA(:), aRB(:), epsin(:), epsR(:), V0e(:), V0h(:), linker(:)
   real(dp),allocatable :: epsinA(:), epsinB(:), epsRA(:), epsRB(:), V0eA(:), V0eB(:), V0hA(:), V0hB(:)
   real(dp),allocatable :: Eeh1(:), Eeh2(:), Cb_eh1(:), Cb_eh2(:), Norm_Ana_e(:), Norm_Ana_h1(:), Norm_Ana_h2(:)
   real(dp),allocatable :: OverlapAna_h1e(:), OverlapAna_h2e(:), Cb_Num_eh1(:), Cb_Num_eh1_eh2(:), Cb_Num_eh2(:)
   real(dp),allocatable :: minEe(:,:),minEh(:,:)
   real(dp),allocatable :: TransDip_Ana_h1e(:), TransDip_Ana_h2e(:), Oscillator_Ana_h1e(:), Oscillator_Ana_h2e(:)
   real(dp),allocatable :: ExctCoef_h1e(:), ExctCoef_h2e(:), Ham(:,:), E0(:), c0(:), TransHam(:,:), Hamt(:,:,:), c(:,:)
   complex(KIND=8) :: ct1, ct2, ct3, ct4, xt01, xt02, xt03, xhbar, im, xwidth, xomega , xEd, xh, xphase, xtime
   complex(KIND=8),allocatable :: xHam(:,:) , xHamt(:,:,:), xTransHam(:,:), xE0(:), xHamtk2(:,:,:), xHamtk3(:,:,:), xHamtk4(:,:,:)
   complex(KIND=8),allocatable :: xc0(:), xc(:,:), xcnew(:,:), k1(:), k2(:), k3(:) , k4(:)

contains 

subroutine getVariables

NAMELIST /version/ vers
NAMELIST /outputs/ o_Norm,o_Over,o_Coul,o_DipS,o_Osci,o_Exti,o_DipD,dyn
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

im         = dcmplx(0.0,1.0)
me         = me*m0
mh         = mh*m0
rhoe       = 1.0/sqrt((2*me*omegaLO)/hbar)
rhoh       = 1.0/sqrt((2*mh*omegaLO)/hbar)
V0         = V0eV*elec

if ( dyn .eq. 'y' ) then

rewind 150
read(150,NML=pulses)



xhbar      = dcmplx(hbar,0.0)
xh         = dcmplx(timestep,0.0)
xt01       = dcmplx(t01,0.0)
xt02       = dcmplx(t02,0.0)
xt03       = dcmplx(t03,0.0)
xphase     = dcmplx(phase,0.0)
xomega     = dcmplx(omega,0.0)
xwidth     = dcmplx(width,0.0)
xEd        = dcmplx(Ed,0.0)
ntime      = nint(totaltime/timestep)
phase      = pi

allocate(Hamt(0:nstates-1,0:nstates-1,0:ntime))
allocate(xHamt(0:nstates-1,0:nstates-1,0:ntime))
allocate(xHamtk2(0:nstates-1,0:nstates-1,0:ntime))
allocate(xHamtk3(0:nstates-1,0:nstates-1,0:ntime))
allocate(xHamtk4(0:nstates-1,0:nstates-1,0:ntime))
allocate(c(0:nstates-1,0:ntime))
allocate(xc(0:nstates-1,0:ntime))
allocate(xcnew(0:nstates-1,ntime+1))
allocate(TransHam(0:nstates-1,0:nstates-1))
allocate(xTransHam(0:nstates-1,0:nstates-1))
allocate(Ham(0:nstates-1,0:nstates-1))
allocate(xHam(0:nstates-1,0:nstates-1))
allocate(c0(0:nstates-1))
allocate(k1(0:nstates-1))
allocate(k2(0:nstates-1))
allocate(k3(0:nstates-1))
allocate(k4(0:nstates-1))

if ( npulses .eq. 3) then
pulse1 = 1.0
pulse2 = 1.0
pulse3 = 1.0
elseif ( npulses .eq. 2) then
pulse1 = 1.0
pulse2 = 1.0
pulse3 = 0.0
elseif ( npulses .eq. 1) then
pulse1 = 1.0
pulse2 = 0.0
pulse3 = 0.0
endif

endif


if ( vers .eq. 'singl' ) then

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

aR(1) = aA
aR(2) = aB
linker(:) = link

rmin = 1
rmax = 2

do n = 1,2
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n)))**1.2)
epsR(n) = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36/35)*aR(n)/rhoe)+exp(-(36/35)*aR(n)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
enddo

else if ( vers .eq. 'range' ) then

rewind 150
read(150,NML=syst_range)

rmin = minr/rsteps
rmax = maxr/rsteps
ndots= maxr/rsteps - minr/rsteps

allocate(aR(rmax-rmin))
allocate(linker(rmax-rmin))
allocate(epsin(rmax-rmin))
allocate(epsR(rmax-rmin))
allocate(V0e(rmax-rmin))
allocate(V0h(rmax-rmin))

linker(:) = link

do n = rmin,rmax
aR(n)= n*rsteps
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n)))**1.2)
epsR(n) = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36/35)*aR(n)/rhoe)+exp(-(36/35)*aR(n)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
enddo

else if ( ( vers .eq. 'randm' ) .and. ( aA .eq. aB ) ) then

rewind 150
read(150,NML=syst_random)

allocate(aR(nsys))
allocate(linker(nsys))
allocate(epsin(nsys))
allocate(epsR(nsys))
allocate(V0e(nsys))
allocate(V0h(nsys))

call random_seed() 

rmin = 1
rmax = nsys

do n = 1,nsys
call random_number(rdmQDA)
call random_number(rdmlinker)
rdmQDA = 2*(rdmQDA-0.5)
rdmlinker = 2*(rdmlinker-0.5)
aR(n)= rdmQDA*dispQD*aA+aA
linker(n) = rdmlinker*dispQD*link+link
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n)))**1.2)
epsR(n) = 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36/35)*aR(n)/rhoe)+exp(-(36/35)*aR(n)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
enddo

else if ( ( vers .eq. 'randm' ) .and. ( aA .eq. aB ) ) then

rewind 150
read(150,NML=syst_random)

allocate(aR(2*nsys))
allocate(linker(2*nsys))
allocate(epsin(2*nsys))
allocate(epsR(2*nsys))
allocate(V0e(2*nsys))
allocate(V0h(2*nsys))

call random_seed()
linker(:) = link

rmin = 1
rmax = nsys

do n = 1,2*nsys
   if ( n .le. nsys ) then 
      call random_number(rdmQDA)
      rdmQDA    = 2*(rdmQDA-0.5)
      aR(n)     = (rdmQDA)*dispQD*aA+aA
   elseif ( n .gt. nsys ) then
      call random_number(rdmQDB)
      rdmQDB    = 2*(rdmQDB-0.5)
      aR(n+nsys)= (rdmQDB)*dispQD*aB+aB
   endif
epsin(n) = 1.0 + (eps - 1.0) / (1.0 + (0.75d-9/(2*aR(n)))**1.2)
epsR(n)= 1.0/((1.0/epsin(n))-((1.0/epsin(n))-(1.0/(epsin(n)+3.5)))*(1-(exp(-(36/35)*aR(n)/rhoe)+exp(-(36/35)*aR(n)/rhoh))/2))
V0e(n)=-1*(-3.49+2.47*(1d9*2*aR(n))**(-1.32))*elec
V0h(n)=-1*(-5.23-0.74*(1d9*2*aR(n))**(-0.95))*elec
enddo

endif 

end subroutine getVariables
   
end module Variables
