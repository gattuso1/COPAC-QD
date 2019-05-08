program Absorption

implicit none

   integer,  parameter :: dp = SELECTED_REAL_KIND(15,307)
   real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
   character*128 :: line, linex, liney, linez, Re_c_ei, Im_c_ei
   integer :: t,t2, tmax,wmin,wmax,error,nn,i,j,t0ei,nsamples,nstates,n, en_id, en, k, pol, npol
   real(dp) :: normt, normw, time0ei, timestep, w1,w2, wstep, w, elec, h, pow_ave, dummy, en_i, en_f, en_step, tottime, dummy2
   real(dp) :: integPol
   real(dp),allocatable :: time(:), cohe(:), pow_gaus(:), E(:,:), raA(:), raB(:), Ediff(:,:,:), c_ei(:,:)
   real(dp),allocatable :: maxid(:), TransDipx(:,:,:), TransDipy(:,:,:), TransDipz(:,:,:), k_1(:), k_2(:), k_3(:)
   real(dp),allocatable :: Dcenter(:,:), l1(:), l2(:), l3(:), Energy(:,:)
   real(dp),allocatable :: TransDip_eix(:,:,:), TransDip_eiy(:,:,:), TransDip_eiz(:,:,:), pulse(:)
   real(dp),allocatable :: Re_c(:,:,:), Im_c(:,:,:), sum_pow(:,:), Ham_ei(:,:,:), pow_gaus2(:), p1(:), p2(:), p3(:)

nstates=9
h  = 6.62607004d-34
elec  = 1.60217662d-19

nn=5000

open(20,file='Etransitions-he_ei.dat')
open(21,file='TransMat_ei.dat')
open(24,file='Absorption.dat')

allocate(TransDip_eix(nn,0:nstates-1,0:nstates-1))
allocate(Energy(nn,0:nstates-1))

read(21,'(A128)')

do n=1,nn

read(20,*) dummy , dummy2 , Energy(n,0), Energy(n,1),&
             Energy(n,2),Energy(n,3),Energy(n,4),Energy(n,5),Energy(n,6),Energy(n,7),Energy(n,8)

read(21,'(A128)')
read(21,*) (TransDip_eix(n,0,j),j=0,nstates-1)
read(21,'(A128)')
read(21,'(A128)')
read(21,'(A128)')
read(21,'(A128)')
read(21,'(A128)')
read(21,'(A128)')
read(21,'(A128)')
read(21,'(A128)')

do i=1,nstates-1
write(6,*) Energy(n,i), (TransDip_eix(n,0,i))**2 ,i
enddo

enddo

end
