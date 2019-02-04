module Make_Ham

use omp_lib
use Constants_au
use Variables_au
use Integrals

real(dp), parameter :: a11_1d_ho = 0.0191697d0  
real(dp), parameter :: a11_2d_ho = 0.242775d0  
real(dp), parameter :: a11_3d_ho = 1.28322d0 
real(dp), parameter :: a11_1e_ho = 0.0131599d0       
real(dp), parameter :: a11_2e_ho = 0.199444d0        
real(dp), parameter :: a11_3e_ho = 1.10133d0
real(dp), parameter :: a22_1d_ho = 0.0178298d0  
real(dp), parameter :: a22_2d_ho = 0.239631d0   
real(dp), parameter :: a22_3d_ho = 1.25828d0
real(dp), parameter :: a22_1e_ho = 0.00477726d0      
real(dp), parameter :: a22_2e_ho = 0.0379395d0       
real(dp), parameter :: a22_3e_ho = 1.66235d0
real(dp), parameter :: a12_1d_ho = -0.00288509
real(dp), parameter :: a12_2d_ho = 0.0260066
real(dp), parameter :: a12_3d_ho = 0.57151
real(dp), parameter :: a12_1e_ho = -0.00651569   
real(dp), parameter :: a12_2e_ho = 0.0530662     
real(dp), parameter :: a12_3e_ho = 1.71364
real(dp), parameter :: a55_1d_ho = 0.019678d0
real(dp), parameter :: a55_2d_ho = 0.26841d0
real(dp), parameter :: a55_3d_ho = 1.22231d+09
real(dp), parameter :: a55_4d_ho = 4.6499d+08
real(dp), parameter :: a66_1d_ho = 0.0196255d0
real(dp), parameter :: a66_2d_ho = 0.270825d0
real(dp), parameter :: a66_3d_ho = 1.21959d+09
real(dp), parameter :: a66_4d_ho = 4.64529d+08
real(dp), parameter :: a13_1d_ho = -6.07575d-08
real(dp), parameter :: a13_2d_ho = 1.94165d-10
real(dp), parameter :: a13_3d_ho = 6.65937d0
real(dp), parameter :: a13_4d_ho = 1.34041d+10
real(dp), parameter :: a15_1d_ho = -3.03287d-06
real(dp), parameter :: a15_2d_ho = 3.28361d-10
real(dp), parameter :: a15_3d_ho = 4.53023d0
real(dp), parameter :: a15_4d_ho = 5.84667d+09
real(dp), parameter :: a17_1d_ho = -1.18454d-08
real(dp), parameter :: a17_2d_ho = 1.93051d-10
real(dp), parameter :: a17_3d_ho = 4.90209d0
real(dp), parameter :: a17_4d_ho = 1.05339d+10
real(dp), parameter :: a24_1d_ho = 2.01287d-08
real(dp), parameter :: a24_2d_ho = 3.13757d-10
real(dp), parameter :: a24_3d_ho = 7.90071d0
real(dp), parameter :: a24_4d_ho = 1.42154d+10
real(dp), parameter :: a26_1d_ho = -4.09828d-06
real(dp), parameter :: a26_2d_ho = 3.25598d-10
real(dp), parameter :: a26_3d_ho = 4.51295d0
real(dp), parameter :: a26_4d_ho = 5.72301d+09
real(dp), parameter :: a28_1d_ho = -3.33925d-07
real(dp), parameter :: a28_2d_ho = 2.72707d-10
real(dp), parameter :: a28_3d_ho = 5.15739d0
real(dp), parameter :: a28_4d_ho = 8.8487d+09
real(dp), parameter :: a68_1d_ho = 4.1131d-08
real(dp), parameter :: a68_2d_ho = 3.2414d-10
real(dp), parameter :: a68_3d_ho = 8.09373d0
real(dp), parameter :: a68_4d_ho = 1.42896d+10
real(dp), parameter :: a56_1d_ho = -0.000176992d0
real(dp), parameter :: a56_2d_ho = 3.09602d-14
real(dp), parameter :: a56_3d_ho = 0.663239d0
real(dp), parameter :: a56_4d_ho = 1.32932d+08
real(dp), parameter :: a57_1d_ho = -9.32088d-09
real(dp), parameter :: a57_2d_ho = 2.35068d-10
real(dp), parameter :: a57_3d_ho = 7.27929d0
real(dp), parameter :: a57_4d_ho = 1.57266d+10
real(dp), parameter :: a25_1d_ho = 1.8281d-06
real(dp), parameter :: a25_2d_ho = 3.02171d-10
real(dp), parameter :: a25_3d_ho = 6.27232d0
real(dp), parameter :: a25_4d_ho = 1.89898d+10
real(dp), parameter :: a18_1d_ho = -2.4777d-07
real(dp), parameter :: a18_2d_ho = 2.29293d-10
real(dp), parameter :: a18_3d_ho = 5.01759d0
real(dp), parameter :: a18_4d_ho = 9.33159d+09
real(dp), parameter :: a27_1d_ho = -1.20211d-07
real(dp), parameter :: a27_2d_ho = 2.29099d-10
real(dp), parameter :: a27_3d_ho = 5.02518d0
real(dp), parameter :: a27_4d_ho = 9.60451d+09
real(dp), parameter :: a16_1d_ho = -7.52249d-07
real(dp), parameter :: a16_2d_ho = 4.62924d-11
real(dp), parameter :: a16_3d_ho = 3.28286d0
real(dp), parameter :: a16_4d_ho = 3.92155d+09
real(dp), parameter :: a14_1d_ho = 2.21733d-08
real(dp), parameter :: a14_2d_ho = 2.79711d-10
real(dp), parameter :: a14_3d_ho = 7.71956d0
real(dp), parameter :: a14_4d_ho = 1.51704d+10
real(dp), parameter :: a58_1d_ho = -4.36445d-08
real(dp), parameter :: a58_2d_ho = 2.60451d-10
real(dp), parameter :: a58_3d_ho = 7.39341d0
real(dp), parameter :: a58_4d_ho = 1.41517d+10
real(dp), parameter :: a13_1e_ho = 0.0172999d0
real(dp), parameter :: a13_2e_ho = 0.167727d0
real(dp), parameter :: a13_3e_ho = 1.06735d+09
real(dp), parameter :: a13_4e_ho = 4.93919d+08
real(dp), parameter :: a14_1e_ho = -0.00397311d0
real(dp), parameter :: a14_2e_ho = -0.105774d0
real(dp), parameter :: a14_3e_ho = 1.38848d+09
real(dp), parameter :: a14_4e_ho = 3.93184d+08
real(dp), parameter :: a24_1e_ho = 0.000964533d0
real(dp), parameter :: a24_2e_ho = 0.0633182d0
real(dp), parameter :: a24_3e_ho = 1.71083d+09
real(dp), parameter :: a24_4e_ho = 3.59354d+08
real(dp), parameter :: a15_1e_ho = -3.07142d-07
real(dp), parameter :: a15_2e_ho = 2.92033d-10
real(dp), parameter :: a15_3e_ho = 5.14692d0
real(dp), parameter :: a15_4e_ho = 6.99873d+09
real(dp), parameter :: a17_1e_ho = -4.63489d-07
real(dp), parameter :: a17_2e_ho = 2.48465d-10
real(dp), parameter :: a17_3e_ho = 4.72745d0
real(dp), parameter :: a17_4e_ho = 7.36577d+09
real(dp), parameter :: a55_1e_ho = -3.37156d-08
real(dp), parameter :: a55_2e_ho = 2.66642d-10
real(dp), parameter :: a55_3e_ho = 7.74613d0
real(dp), parameter :: a55_4e_ho = 1.31425d+10
real(dp), parameter :: a66_1e_ho = -1.60013d-08
real(dp), parameter :: a66_2e_ho = 3.05203d-10
real(dp), parameter :: a66_3e_ho = 7.65275d0
real(dp), parameter :: a66_4e_ho = 1.30185d+10
real(dp), parameter :: a26_1e_ho = -3.0698d-07
real(dp), parameter :: a26_2e_ho = 2.81903d-10
real(dp), parameter :: a26_3e_ho = 5.35289d0
real(dp), parameter :: a26_4e_ho = 6.74812d+09
real(dp), parameter :: a28_1e_ho = -2.97071d-07
real(dp), parameter :: a28_2e_ho = 2.61368d-10
real(dp), parameter :: a28_3e_ho = 5.14267d0
real(dp), parameter :: a28_4e_ho = 7.10254d+09
real(dp), parameter :: a68_1e_ho = -2.58585d-08
real(dp), parameter :: a68_2e_ho = 3.06916d-10
real(dp), parameter :: a68_3e_ho = 7.70315d0
real(dp), parameter :: a68_4e_ho = 1.31844d+10
real(dp), parameter :: a57_1e_ho = -2.78302d-08
real(dp), parameter :: a57_2e_ho = 2.13605d-10
real(dp), parameter :: a57_3e_ho = 6.92294d0
real(dp), parameter :: a57_4e_ho = 1.48671d+10
real(dp), parameter :: a25_1e_ho = -1.21351d-07
real(dp), parameter :: a25_2e_ho = 2.66787d-10
real(dp), parameter :: a25_3e_ho = 5.52544d0
real(dp), parameter :: a25_4e_ho = 7.04239d+09
real(dp), parameter :: a27_1e_ho = -4.37725d-07
real(dp), parameter :: a27_2e_ho = 2.31345d-10
real(dp), parameter :: a27_3e_ho = 5.15123d0
real(dp), parameter :: a27_4e_ho = 7.02434d+09
real(dp), parameter :: a16_1e_ho = -1.48139d-06
real(dp), parameter :: a16_2e_ho = 2.86902d-10
real(dp), parameter :: a16_3e_ho = 4.72314d0
real(dp), parameter :: a16_4e_ho = 6.71625d+09
real(dp), parameter :: a18_1e_ho = -1.33651d-06
real(dp), parameter :: a18_2e_ho = 2.70946d-10
real(dp), parameter :: a18_3e_ho = 4.60265d0
real(dp), parameter :: a18_4e_ho = 7.01259d+09
real(dp), parameter :: a58_1e_ho = 7.84145d-08
real(dp), parameter :: a58_2e_ho = 3.26014d-10
real(dp), parameter :: a58_3e_ho = 8.50289d0
real(dp), parameter :: a58_4e_ho = 1.4545d+10
real(dp), parameter :: a56_1e_ho = -2.51904d-08
real(dp), parameter :: a56_2e_ho = 2.89506d-10
real(dp), parameter :: a56_3e_ho = 7.77886d0
real(dp), parameter :: a56_4e_ho = 1.31881d+10

real(dp) :: a11_1d_he
real(dp) :: a11_2d_he
real(dp) :: a11_3d_he
real(dp) :: a11_1e_he
real(dp) :: a11_2e_he
real(dp) :: a11_3e_he
real(dp) :: a22_1d_he
real(dp) :: a22_2d_he
real(dp) :: a22_3d_he
real(dp) :: a22_1e_he
real(dp) :: a22_2e_he
real(dp) :: a22_3e_he
real(dp) :: a12_1d_he 
real(dp) :: a12_2d_he 
real(dp) :: a12_3d_he 
real(dp) :: a12_1e_he 
real(dp) :: a12_2e_he 
real(dp) :: a12_3e_he 
real(dp) :: a13_1d_he 
real(dp) :: a13_2d_he 
real(dp) :: a13_3d_he 
real(dp) :: a13_4d_he 
real(dp) :: a15_1d_he 
real(dp) :: a15_2d_he 
real(dp) :: a15_3d_he 
real(dp) :: a15_4d_he 
real(dp) :: a17_1d_he 
real(dp) :: a17_2d_he 
real(dp) :: a17_3d_he 
real(dp) :: a17_4d_he 
real(dp) :: a24_1d_he 
real(dp) :: a24_2d_he 
real(dp) :: a24_3d_he 
real(dp) :: a24_4d_he 
real(dp) :: a26_1d_he 
real(dp) :: a26_2d_he 
real(dp) :: a26_3d_he 
real(dp) :: a26_4d_he 
real(dp) :: a28_1d_he 
real(dp) :: a28_2d_he 
real(dp) :: a28_3d_he 
real(dp) :: a28_4d_he 
real(dp) :: a35_1d_he 
real(dp) :: a35_2d_he 
real(dp) :: a35_3d_he 
real(dp) :: a35_4d_he 
real(dp) :: a37_1d_he 
real(dp) :: a37_2d_he 
real(dp) :: a37_3d_he 
real(dp) :: a37_4d_he 
real(dp) :: a46_1d_he 
real(dp) :: a46_2d_he 
real(dp) :: a46_3d_he 
real(dp) :: a46_4d_he 
real(dp) :: a48_1d_he 
real(dp) :: a48_2d_he 
real(dp) :: a48_3d_he 
real(dp) :: a48_4d_he 
real(dp) :: a55_1d_he 
real(dp) :: a55_2d_he 
real(dp) :: a55_3d_he 
real(dp) :: a55_4d_he 
real(dp) :: a57_1d_he 
real(dp) :: a57_2d_he 
real(dp) :: a57_3d_he 
real(dp) :: a57_4d_he 
real(dp) :: a66_1d_he 
real(dp) :: a66_2d_he 
real(dp) :: a66_3d_he 
real(dp) :: a66_4d_he 
real(dp) :: a68_1d_he 
real(dp) :: a68_2d_he 
real(dp) :: a68_3d_he 
real(dp) :: a68_4d_he 
real(dp) :: a77_1d_he 
real(dp) :: a77_2d_he 
real(dp) :: a77_3d_he 
real(dp) :: a77_4d_he 
real(dp) :: a88_1d_he 
real(dp) :: a88_2d_he 
real(dp) :: a88_3d_he 
real(dp) :: a88_4d_he 
real(dp) :: a14_1d_he 
real(dp) :: a14_2d_he 
real(dp) :: a14_3d_he 
real(dp) :: a14_4d_he 
real(dp) :: a16_1d_he 
real(dp) :: a16_2d_he 
real(dp) :: a16_3d_he 
real(dp) :: a16_4d_he 
real(dp) :: a18_1d_he 
real(dp) :: a18_2d_he 
real(dp) :: a18_3d_he 
real(dp) :: a18_4d_he 
real(dp) :: a23_1d_he 
real(dp) :: a23_2d_he 
real(dp) :: a23_3d_he 
real(dp) :: a23_4d_he 
real(dp) :: a25_1d_he 
real(dp) :: a25_2d_he 
real(dp) :: a25_3d_he 
real(dp) :: a25_4d_he 
real(dp) :: a27_1d_he 
real(dp) :: a27_2d_he 
real(dp) :: a27_3d_he 
real(dp) :: a27_4d_he 
real(dp) :: a36_1d_he 
real(dp) :: a36_2d_he 
real(dp) :: a36_3d_he 
real(dp) :: a36_4d_he 
real(dp) :: a38_1d_he 
real(dp) :: a38_2d_he 
real(dp) :: a38_3d_he 
real(dp) :: a38_4d_he 
real(dp) :: a45_1d_he 
real(dp) :: a45_2d_he 
real(dp) :: a45_3d_he 
real(dp) :: a45_4d_he 
real(dp) :: a47_1d_he 
real(dp) :: a47_2d_he 
real(dp) :: a47_3d_he 
real(dp) :: a47_4d_he 
real(dp) :: a56_1d_he 
real(dp) :: a56_2d_he 
real(dp) :: a56_3d_he 
real(dp) :: a56_4d_he 
real(dp) :: a58_1d_he 
real(dp) :: a58_2d_he 
real(dp) :: a58_3d_he 
real(dp) :: a58_4d_he 
real(dp) :: a67_1d_he 
real(dp) :: a67_2d_he 
real(dp) :: a67_3d_he 
real(dp) :: a67_4d_he 
real(dp) :: a78_1d_he 
real(dp) :: a78_2d_he 
real(dp) :: a78_3d_he 
real(dp) :: a78_4d_he 
real(dp) :: a13_1e_he 
real(dp) :: a13_2e_he 
real(dp) :: a13_3e_he 
real(dp) :: a13_4e_he 
real(dp) :: a15_1e_he 
real(dp) :: a15_2e_he 
real(dp) :: a15_3e_he 
real(dp) :: a15_4e_he 
real(dp) :: a17_1e_he 
real(dp) :: a17_2e_he 
real(dp) :: a17_3e_he 
real(dp) :: a17_4e_he 
real(dp) :: a24_1e_he 
real(dp) :: a24_2e_he 
real(dp) :: a24_3e_he 
real(dp) :: a24_4e_he 
real(dp) :: a26_1e_he 
real(dp) :: a26_2e_he 
real(dp) :: a26_3e_he 
real(dp) :: a26_4e_he 
real(dp) :: a28_1e_he 
real(dp) :: a28_2e_he 
real(dp) :: a28_3e_he 
real(dp) :: a28_4e_he 
real(dp) :: a35_1e_he 
real(dp) :: a35_2e_he 
real(dp) :: a35_3e_he 
real(dp) :: a35_4e_he 
real(dp) :: a37_1e_he 
real(dp) :: a37_2e_he 
real(dp) :: a37_3e_he 
real(dp) :: a37_4e_he 
real(dp) :: a46_1e_he 
real(dp) :: a46_2e_he 
real(dp) :: a46_3e_he 
real(dp) :: a46_4e_he 
real(dp) :: a48_1e_he 
real(dp) :: a48_2e_he 
real(dp) :: a48_3e_he 
real(dp) :: a48_4e_he 
real(dp) :: a55_1e_he 
real(dp) :: a55_2e_he 
real(dp) :: a55_3e_he 
real(dp) :: a55_4e_he 
real(dp) :: a57_1e_he 
real(dp) :: a57_2e_he 
real(dp) :: a57_3e_he 
real(dp) :: a57_4e_he 
real(dp) :: a66_1e_he 
real(dp) :: a66_2e_he 
real(dp) :: a66_3e_he 
real(dp) :: a66_4e_he 
real(dp) :: a68_1e_he 
real(dp) :: a68_2e_he 
real(dp) :: a68_3e_he 
real(dp) :: a68_4e_he 
real(dp) :: a77_1e_he 
real(dp) :: a77_2e_he 
real(dp) :: a77_3e_he 
real(dp) :: a77_4e_he 
real(dp) :: a88_1e_he 
real(dp) :: a88_2e_he 
real(dp) :: a88_3e_he 
real(dp) :: a88_4e_he 
real(dp) :: a14_1e_he 
real(dp) :: a14_2e_he 
real(dp) :: a14_3e_he 
real(dp) :: a14_4e_he 
real(dp) :: a16_1e_he 
real(dp) :: a16_2e_he 
real(dp) :: a16_3e_he 
real(dp) :: a16_4e_he 
real(dp) :: a18_1e_he 
real(dp) :: a18_2e_he 
real(dp) :: a18_3e_he 
real(dp) :: a18_4e_he 
real(dp) :: a23_1e_he 
real(dp) :: a23_2e_he 
real(dp) :: a23_3e_he 
real(dp) :: a23_4e_he 
real(dp) :: a25_1e_he 
real(dp) :: a25_2e_he 
real(dp) :: a25_3e_he 
real(dp) :: a25_4e_he 
real(dp) :: a27_1e_he 
real(dp) :: a27_2e_he 
real(dp) :: a27_3e_he 
real(dp) :: a27_4e_he 
real(dp) :: a36_1e_he 
real(dp) :: a36_2e_he 
real(dp) :: a36_3e_he 
real(dp) :: a36_4e_he 
real(dp) :: a38_1e_he 
real(dp) :: a38_2e_he 
real(dp) :: a38_3e_he 
real(dp) :: a38_4e_he 
real(dp) :: a45_1e_he 
real(dp) :: a45_2e_he 
real(dp) :: a45_3e_he 
real(dp) :: a45_4e_he 
real(dp) :: a47_1e_he 
real(dp) :: a47_2e_he 
real(dp) :: a47_3e_he 
real(dp) :: a47_4e_he 
real(dp) :: a56_1e_he 
real(dp) :: a56_2e_he 
real(dp) :: a56_3e_he 
real(dp) :: a56_4e_he 
real(dp) :: a58_1e_he 
real(dp) :: a58_2e_he 
real(dp) :: a58_3e_he 
real(dp) :: a58_4e_he 
real(dp) :: a67_1e_he 
real(dp) :: a67_2e_he 
real(dp) :: a67_3e_he 
real(dp) :: a67_4e_he 
real(dp) :: a78_1e_he 
real(dp) :: a78_2e_he 
real(dp) :: a78_3e_he 
real(dp) :: a78_4e_he 

contains

subroutine make_Ham_ho

Ham(0,0) = 0.0d0

Ham(1,1) = minEe(1,n) + minEh(1,n) + V0 - 1.0 * elec*(a11_1d_ho + a11_2d_ho / ((aR(n)*1d9)**a11_3d_ho)) !&
!              + 2 * elec*(a11_1e_he + a11_2e_he / ((aR(n)*1d9)**a11_3e_he))

Ham(2,2) = minEe(1,n) + minEh(2,n) + V0 - 1.0 * elec*(a22_1d_ho + a22_2d_ho / ((aR(n)*1d9)**a22_3d_ho)) !&
!              + 2 * elec*(a22_1e_he + a22_2e_he / ((aR(n)*1d9)**a22_3e_he))

Ham(3,3) = minEe(1,n) + minEh(1,n) + V0 - 1.0 * elec*(a11_1d_ho + a11_2d_ho / ((aR(n)*1d9)**a11_3d_ho)) !&
!              + 2 * elec*(a11_1e_he + a11_2e_he / ((aR(n+nsys)*1d9)**a11_3e_he))

Ham(4,4) = minEe(1,n) + minEh(2,n) + V0 - 1.0 * elec*(a22_1d_ho + a22_2d_ho / ((aR(n)*1d9)**a22_3d_ho)) !&
!              + 2 * elec*(a22_1e_he + a22_2e_he / ((aR(n+nsys)*1d9)**a22_3e_he))

Ham(5,5) = minEe(1,n) + minEh(1,n) + V0 &
                               - elec*(a55_1d_ho + a55_2d_ho * exp(-1.0d0*(a55_3d_ho*aR(n) + a55_4d_ho*linker(n)))) !&
!                             + 2*elec*(a55_1e_ho + a55_2e_ho * exp(-1.0d0*(a55_3e_ho*aR(n) + a55_4e_ho*linker(n))))

Ham(6,6) = minEe(1,n) + minEh(2,n) + V0 - elec*(a66_1d_ho + a66_2d_ho * exp(-1.0d0*(a66_3d_ho*aR(n) + a66_4d_ho*linker(n)))) !&
!                                  + 2*elec*(a66_1e_ho + a66_2e_ho * exp(-1.0d0*(a66_3e_ho*aR(n) + a66_4e_ho*linker(n))))
Ham(7,7) = Ham(5,5)
Ham(8,8) = Ham(6,6)
Ham(0,1) = 0.0d0
Ham(0,2) = 0.0d0
Ham(0,3) = 0.0d0
Ham(0,4) = 0.0d0
Ham(0,5) = 0.0d0
Ham(0,6) = 0.0d0
Ham(0,7) = 0.0d0
Ham(0,8) = 0.0d0
Ham(1,0) = Ham(0,1)
Ham(2,0) = Ham(0,2)
Ham(3,0) = Ham(0,3)
Ham(4,0) = Ham(0,4)
Ham(5,0) = Ham(0,5)
Ham(6,0) = Ham(0,6)
Ham(7,0) = Ham(0,7)
Ham(8,0) = Ham(0,8)

Ham(1,2) = - 1.0d0 * elec*(a12_1d_ho + a12_2d_ho / ((aR(n)*1d9)**a12_3d_ho)) &
             + elec*(a12_1e_ho + a12_2e_ho / ((aR(n)*1d9)**a12_3e_ho))

!write(20,*) aR(n), (a12_1d_ho + a12_2d_ho / exp(aR(n)*a12_3d_ho)), (a12_1e_ho + a12_2e_ho / exp(aR(n)*a12_3e_ho))

Ham(1,3) = - 1.0d0 * elec*(a13_1d_ho + a13_2d_ho * exp(-1.0d0*(a13_3d_ho*aR(n) + a13_4d_ho*linker(n)))) &
             + elec*(a13_1e_ho + a13_2e_ho * exp(-1.0d0*(a13_3e_ho*aR(n) + a13_4e_ho*linker(n))))

Ham(1,4) = - 1.0d0 * elec*(a14_1d_ho + a14_2d_ho * exp(-1.0d0*(a14_3d_ho*aR(n) + a14_4d_ho*linker(n)))) &
             + elec*(a14_1e_ho + a14_2e_ho * exp(-1.0d0*(a14_3e_ho*aR(n) + a14_4e_ho*linker(n))))

Ham(1,5) =  - 1.0d0 * elec*(a15_1d_ho + (a15_2d_ho/aR(n))**a15_3d_ho * exp(-1.0d0*a15_4d_ho*linker(n))) &
             + elec*(a15_1e_ho + (a15_2e_ho/aR(n))**a15_3e_ho * exp(-1.0d0*a15_4e_ho*linker(n)))

Ham(1,6) = - 1.0d0 *( -1.0d0 * elec * (a16_1d_ho + (a16_2d_ho/aR(n))**a16_3d_ho * exp(-1.0d0*a16_4d_ho*linker(n))) &
                     + elec * (a16_1e_ho + (a16_2e_ho/aR(n))**a16_3e_ho * exp(-1.0d0*a16_4e_ho*linker(n))))

Ham(1,7) = - 1.0d0 * elec*(a17_1d_ho + (a17_2d_ho/aR(n))**a17_3d_ho * exp(-1.0d0*a17_4d_ho*linker(n))) &
             + elec*(a17_1e_ho + (a17_2e_ho/aR(n))**a17_3e_ho * exp(-1.0d0*a17_4e_ho*linker(n)))

Ham(1,8) =  - 1.0d0 *( -1.0d0 * elec * (a18_1d_ho + (a18_2d_ho/aR(n))**a18_3d_ho * exp(-1.0d0*a18_4d_ho*linker(n))) &
                      +  elec * (a18_1e_ho + (a18_2e_ho/aR(n))**a18_3e_ho * exp(-1.0d0*a18_4e_ho*linker(n))))

Ham(2,3) = Ham(1,4)

Ham(2,4) = - 1.0d0 * elec*(a24_1d_ho + a24_2d_ho * exp(-1.0d0*(a24_3d_ho*aR(n) + a24_4d_ho*linker(n)))) &
             + elec*(a24_1e_ho + a24_2e_ho * exp(-1.0d0*(a24_3e_ho*aR(n) + a24_4e_ho*linker(n))))

Ham(2,5) =  - 1.0d0 *( -1.0d0 * elec*(a25_1d_ho + (a25_2d_ho/aR(n))**a25_3d_ho * exp(-1.0d0*a25_4d_ho*linker(n))) &
                      + elec*(a25_1e_ho + (a25_2e_ho/aR(n))**a25_3e_ho * exp(-1.0d0*a25_4e_ho*linker(n))))

Ham(2,6) =  - 1.0d0 * elec*(a26_1d_ho + (a26_2d_ho/aR(n))**a26_3d_ho * exp(-1.0d0*a26_4d_ho*linker(n))) &
              + elec*(a26_1e_ho + (a26_2e_ho/aR(n))**a26_3e_ho * exp(-1.0d0*a26_4e_ho*linker(n)))

Ham(2,7) =  - 1.0d0 *( -1.0d0 * elec * (a27_1d_ho + (a27_2d_ho/aR(n))**a27_3d_ho * exp(-1.0d0*a27_4d_ho*linker(n))) &
                      + elec * (a27_1e_ho + (a27_2e_ho/aR(n))**a27_3e_ho * exp(-1.0d0*a27_4e_ho*linker(n))))

Ham(2,8) =  - 1.0d0 * elec*(a28_1d_ho + (a28_2d_ho/aR(n))**a28_3d_ho * exp(-1.0d0*a28_4d_ho*linker(n))) &
              + elec*(a28_1e_ho + (a28_2e_ho/aR(n))**a28_3e_ho * exp(-1.0d0*a28_4e_ho*linker(n)))

Ham(3,4) = Ham(1,2)
Ham(3,5) = Ham(1,7)
Ham(3,6) = Ham(1,8)
Ham(3,7) = Ham(1,5)
Ham(3,8) = Ham(1,6)
Ham(4,5) = Ham(2,7)
Ham(4,6) = Ham(2,8)
Ham(4,7) = Ham(2,5)
Ham(4,8) = Ham(2,6)

Ham(5,6) = - 1.0d0 * elec*(a56_1d_ho + (a56_2d_ho/aR(n))**a56_3d_ho * exp(-1.0d0*a56_4d_ho*linker(n))) &
             + elec*(a56_1e_ho + (a56_2e_ho/aR(n))**a56_3e_ho * exp(-1.0d0*a56_4e_ho*linker(n)))

Ham(5,7) = - 1.0d0 * elec*(a57_1d_ho + (a57_2d_ho/aR(n))**a57_3d_ho * exp(-1.0d0*a57_4d_ho*linker(n))) &
             + elec*(a57_1e_ho + (a57_2e_ho/aR(n))**a57_3e_ho * exp(-1.0d0*a57_4e_ho*linker(n)))

Ham(5,8) = - 1.0d0 * elec*(a58_1d_ho + (a58_2d_ho/aR(n))**a58_3d_ho * exp(-1.0d0*a58_4d_ho*linker(n))) &
             + elec*(a58_1e_ho + (a58_2e_ho/aR(n))**a58_3e_ho * exp(-1.0d0*a58_4e_ho*linker(n)))

Ham(6,7) = Ham(5,8)

Ham(6,8) = - 1.0d0 * elec*(a68_1d_ho + (a68_2d_ho/aR(n))**a68_3d_ho * exp(-1.0d0*a68_4d_ho*linker(n))) &
             + elec*(a68_1e_ho + (a68_2e_ho/aR(n))**a68_3e_ho * exp(-1.0d0*a68_4e_ho*linker(n)))

Ham(7,8) = Ham(5,6)

Ham(2,1) = Ham(1,2)
Ham(3,1) = Ham(1,3)
Ham(4,1) = Ham(1,4)
Ham(5,1) = Ham(1,5)
Ham(6,1) = Ham(1,6)
Ham(7,1) = Ham(1,7)
Ham(8,1) = Ham(1,8)
Ham(3,2) = Ham(2,3)
Ham(4,2) = Ham(2,4)
Ham(5,2) = Ham(2,5)
Ham(6,2) = Ham(2,6)
Ham(7,2) = Ham(2,7)
Ham(8,2) = Ham(2,8)
Ham(4,3) = Ham(3,4)
Ham(5,3) = Ham(3,5)
Ham(6,3) = Ham(3,6)
Ham(7,3) = Ham(3,7)
Ham(8,3) = Ham(3,8)
Ham(5,4) = Ham(4,5)
Ham(6,4) = Ham(4,6)
Ham(7,4) = Ham(4,7)
Ham(8,4) = Ham(4,8)
Ham(6,5) = Ham(5,6)
Ham(7,5) = Ham(5,7)
Ham(8,5) = Ham(5,8)
Ham(7,6) = Ham(6,7)
Ham(8,6) = Ham(6,8)
Ham(8,7) = Ham(7,8)

!Ham = 0.1d-19
Ham = Ham/Energ_au
TransHam(0,1) = TransDip_Ana_h1e(n)
TransHam(0,2) = TransDip_Ana_h2e(n)
TransHam(0,3) = TransHam(0,1)
TransHam(0,4) = TransHam(0,2)
TransHam(0,5) = TransDip_Fit_h1e_ho(aR(n),linker(n))
TransHam(0,6) = TransDip_Fit_h2e_ho(aR(n),linker(n))
TransHam(0,7) = TransHam(0,5)
TransHam(0,8) = TransHam(0,6)
!TransHam = 0.2d-33
TransHam = TransHam/Dip_au
do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo
do i=0,nstates-1
do j=0,nstates-1
if ( ( Ham(i,j) .ne. Ham(j,i) ) .or. ( TransHam(i,j) .ne. TransHam(j,i) ) ) then
write(6,*) i, j, "HAMILTONIAN NON HERMITIAN"
endif
enddo
enddo

end subroutine make_Ham_ho

subroutine make_Ham_he

if ( idlink .eq. 20 ) then
include 'Parameters-dir-02.f90'
include 'Parameters-ex-02.f90'
elseif ( idlink .eq. 55 ) then
include 'Parameters-dir-55.f90'
include 'Parameters-ex-55.f90'
endif

if (vers .eq. 'dimer') then

Ham     = 0.0
Ham_0   = 0.0
Ham_dir = 0.0
Ham_ex  = 0.0

Ham_0(1)     = minEe(1,1) + minEh(1,1) + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aR(1)*1d9)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(1)*1d9)**a11_3e_ho))

Ham_0(2)     = minEe(1,1) + minEh(2,1) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aR(1)*1d9)**a22_3d_ho)) 
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(1)*1d9)**a22_3e_ho))

Ham_0(3)     = minEe(1,2) + minEh(1,2) + V0 
Ham_dir(3,3) = elec*(a11_1d_ho + a11_2d_ho / ((aR(2)*1d9)**a11_3d_ho)) 
Ham_ex(3,3)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(2)*1d9)**a11_3e_ho))

Ham_0(4)     = minEe(1,2) + minEh(2,2) + V0 
Ham_dir(4,4) = elec*(a22_1d_ho + a22_2d_ho / ((aR(2)*1d9)**a22_3d_ho)) 
Ham_ex(4,4)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(2)*1d9)**a22_3e_ho))

Ham_0(5)     = minEe(1,2) + minEh(1,1) + V0 
Ham_dir(5,5) = elec*(a55_1d_he + a55_2d_he / ((aR(1)*1d9)**a55_3d_he * (aR(2)*1d9)**a55_4d_he)) 
Ham_ex(5,5)  = elec*(a55_1e_he + a55_2e_he / ((aR(1)*1d9)**a55_3e_he * (aR(2)*1d9)**a55_4e_he))
                                  
Ham_0(6)     = minEe(1,2) + minEh(2,1) + V0 
Ham_dir(6,6) = elec*(a66_1d_he + a66_2d_he / ((aR(1)*1d9)**a66_3d_he * (aR(2)*1d9)**a66_4d_he)) 
Ham_ex(6,6)  = elec*(a66_1e_he + a66_2e_he / ((aR(1)*1d9)**a66_3e_he * (aR(2)*1d9)**a66_4e_he)) 
                                  
!Ham(7,7) = minEe(1,1) + minEh(1,2) + V0 &
! - 1.0 * elec*(a77_1d_he + a77_2d_he / ((aR(1)*1d9)**a77_3d_he * (aR(2)*1d9)**a77_4d_he)) &
!  +  elec*(a77_1e_he + a77_2e_he / ((aR(1)*1d9)**a77_3e_he * (aR(2)*1d9)**a77_4e_he))
Ham_0(7)     = minEe(1,1) + minEh(1,2) + V0 
Ham_dir(7,7) = elec*(a55_1d_he + a55_2d_he / ((aR(2)*1d9)**a55_3d_he * (aR(1)*1d9)**a55_4d_he)) 
Ham_ex(7,7)  = elec*(a55_1e_he + a55_2e_he / ((aR(2)*1d9)**a55_3e_he * (aR(1)*1d9)**a55_4e_he))

!Ham(8,8) = minEe(1,1) + minEh(2,2) + V0 &
! - 1.0 * elec*(a88_1d_he + a88_2d_he / ((aR(1)*1d9)**a88_3d_he * (aR(2)*1d9)**a88_4d_he)) &
!  +  elec*(a88_1e_he + a88_2e_he / ((aR(1)*1d9)**a88_3e_he * (aR(2)*1d9)**a88_4e_he))
Ham_0(8)     = minEe(1,1) + minEh(2,2) + V0 
Ham_dir(8,8) = elec*(a66_1d_he + a66_2d_he / ((aR(2)*1d9)**a66_3d_he * (aR(1)*1d9)**a66_4d_he)) 
Ham_ex(8,8)  = elec*(a66_1e_he + a66_2e_he / ((aR(2)*1d9)**a66_3e_he * (aR(1)*1d9)**a66_4e_he))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aR(1)*1d9)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(1)*1d9)**a12_3e_ho))

Ham_dir(1,3) = elec*(a13_1d_he + a13_2d_he / ((aR(1)*1d9)**a13_3d_he * (aR(2)*1d9)**a13_4d_he)) 
Ham_ex(1,3)  = elec*(a13_1e_he + a13_2e_he / ((aR(1)*1d9)**a13_3e_he * (aR(2)*1d9)**a13_4e_he))

Ham_dir(1,4) = elec*(a14_1d_he + a14_2d_he / ((aR(1)*1d9)**a14_3d_he * (aR(2)*1d9)**a14_4d_he)) 
Ham_ex(1,4)  = elec*(a14_1e_he + a14_2e_he / ((aR(1)*1d9)**a14_3e_he * (aR(2)*1d9)**a14_4e_he))

Ham_dir(1,5) = elec*(a15_1d_he + a15_2d_he / ((aR(1)*1d9)**a15_3d_he * (aR(2)*1d9)**a15_4d_he)) 
Ham_ex(1,5)  = elec*(a15_1e_he + a15_2e_he / ((aR(1)*1d9)**a15_3e_he * (aR(2)*1d9)**a15_4e_he))

Ham_dir(1,6) = elec*(a16_1d_he + a16_2d_he / ((aR(1)*1d9)**a16_3d_he * (aR(2)*1d9)**a16_4d_he)) 
Ham_ex(1,6)  = elec*(a16_1e_he + a16_2e_he / ((aR(1)*1d9)**a16_3e_he * (aR(2)*1d9)**a16_4e_he))

Ham_dir(1,7) = elec*(a17_1d_he + a17_2d_he / ((aR(1)*1d9)**a17_3d_he * (aR(2)*1d9)**a17_4d_he)) 
Ham_ex(1,7)  = elec*(a17_1e_he + a17_2e_he / ((aR(1)*1d9)**a17_3e_he * (aR(2)*1d9)**a17_4e_he))

Ham_dir(1,8) = elec*(a18_1d_he + a18_2d_he / ((aR(1)*1d9)**a18_3d_he * (aR(2)*1d9)**a18_4d_he)) 
Ham_ex(1,8)  = elec*(a18_1e_he + a18_2e_he / ((aR(1)*1d9)**a18_3e_he * (aR(2)*1d9)**a18_4e_he))

!Ham(2,3) = - 1.0 * elec*(a23_1d_he + a23_2d_he / ((aR(1)*1d9)**a23_3d_he * (aR(2)*1d9)**a23_4d_he)) &
!             +  elec*(a23_1e_he + a23_2e_he / ((aR(1)*1d9)**a23_3e_he * (aR(2)*1d9)**a23_4e_he))
Ham_dir(2,3) = elec*(a14_1d_he + a14_2d_he / ((aR(2)*1d9)**a14_3d_he * (aR(1)*1d9)**a14_4d_he)) 
Ham_ex(2,3)  = elec*(a14_1e_he + a14_2e_he / ((aR(2)*1d9)**a14_3e_he * (aR(1)*1d9)**a14_4e_he))

Ham_dir(2,4) = elec*(a24_1d_he + a24_2d_he / ((aR(1)*1d9)**a24_3d_he * (aR(2)*1d9)**a24_4d_he)) 
Ham_ex(2,4)  = elec*(a24_1e_he + a24_2e_he / ((aR(1)*1d9)**a24_3e_he * (aR(2)*1d9)**a24_4e_he))

Ham_dir(2,5) = elec*(a25_1d_he + a25_2d_he / ((aR(1)*1d9)**a25_3d_he * (aR(2)*1d9)**a25_4d_he)) 
Ham_ex(2,5)  = elec*(a25_1e_he + a25_2e_he / ((aR(1)*1d9)**a25_3e_he * (aR(2)*1d9)**a25_4e_he))

Ham_dir(2,6) = elec*(a26_1d_he + a26_2d_he / ((aR(1)*1d9)**a26_3d_he * (aR(2)*1d9)**a26_4d_he)) 
Ham_ex(2,6)  = elec*(a26_1e_he + a26_2e_he / ((aR(1)*1d9)**a26_3e_he * (aR(2)*1d9)**a26_4e_he))

Ham_dir(2,7) = elec*(a27_1d_he + a27_2d_he / ((aR(1)*1d9)**a27_3d_he * (aR(2)*1d9)**a27_4d_he)) 
Ham_ex(2,7)  = elec*(a27_1e_he + a27_2e_he / ((aR(1)*1d9)**a27_3e_he * (aR(2)*1d9)**a27_4e_he))

Ham_dir(2,8) = elec*(a28_1d_he + a28_2d_he / ((aR(1)*1d9)**a28_3d_he * (aR(2)*1d9)**a28_4d_he)) 
Ham_ex(2,8)  = elec*(a28_1e_he + a28_2e_he / ((aR(1)*1d9)**a28_3e_he * (aR(2)*1d9)**a28_4e_he))

Ham_dir(3,4) = elec*(a12_1d_ho + a12_2d_ho / ((aR(2)*1d9)**a12_3d_ho)) 
Ham_ex(3,4)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(2)*1d9)**a12_3e_ho))

!Ham(3,5) =  - 1.0 * elec*(a35_1d_he + a35_2d_he / ((aR(1)*1d9)**a35_3d_he * (aR(2)*1d9)**a35_4d_he)) &
!              +  elec*(a35_1e_he + a35_2e_he / ((aR(1)*1d9)**a35_3e_he * (aR(2)*1d9)**a35_4e_he))
Ham_dir(3,5) = elec*(a17_1d_he + a17_2d_he / ((aR(2)*1d9)**a17_3d_he * (aR(1)*1d9)**a17_4d_he)) 
Ham_ex(3,5)  = elec*(a17_1e_he + a17_2e_he / ((aR(2)*1d9)**a17_3e_he * (aR(1)*1d9)**a17_4e_he))

!Ham(3,6) =  - 1.0 * elec*(a36_1d_he + a36_2d_he / ((aR(1)*1d9)**a36_3d_he * (aR(2)*1d9)**a36_4d_he)) &
!              +  elec*(a36_1e_he + a36_2e_he / ((aR(1)*1d9)**a36_3e_he * (aR(2)*1d9)**a36_4e_he))
Ham_dir(3,6) = elec*(a18_1d_he + a18_2d_he / ((aR(2)*1d9)**a18_3d_he * (aR(1)*1d9)**a18_4d_he)) 
Ham_ex(3,6)  = elec*(a18_1e_he + a18_2e_he / ((aR(2)*1d9)**a18_3e_he * (aR(1)*1d9)**a18_4e_he))

!Ham(3,7) =  - 1.0 * elec*(a37_1d_he + a37_2d_he / ((aR(1)*1d9)**a37_3d_he * (aR(2)*1d9)**a37_4d_he)) &
!              +  elec*(a37_1e_he + a37_2e_he / ((aR(1)*1d9)**a37_3e_he * (aR(2)*1d9)**a37_4e_he))
Ham_dir(3,7) = elec*(a15_1d_he + a15_2d_he / ((aR(2)*1d9)**a15_3d_he * (aR(1)*1d9)**a15_4d_he)) 
Ham_ex(3,7)  = elec*(a15_1e_he + a15_2e_he / ((aR(2)*1d9)**a15_3e_he * (aR(1)*1d9)**a15_4e_he))

!Ham(3,8) =  - 1.0 * elec*(a38_1d_he + a38_2d_he / ((aR(1)*1d9)**a38_3d_he * (aR(2)*1d9)**a38_4d_he)) &
!              +  elec*(a38_1e_he + a38_2e_he / ((aR(1)*1d9)**a38_3e_he * (aR(2)*1d9)**a38_4e_he))
Ham_dir(3,8) = elec*(a16_1d_he + a16_2d_he / ((aR(2)*1d9)**a16_3d_he * (aR(1)*1d9)**a16_4d_he)) 
Ham_ex(3,8)  = elec*(a16_1e_he + a16_2e_he / ((aR(2)*1d9)**a16_3e_he * (aR(1)*1d9)**a16_4e_he))

!Ham(4,5) =  - 1.0 * elec*(a45_1d_he + a45_2d_he / ((aR(1)*1d9)**a45_3d_he * (aR(2)*1d9)**a45_4d_he)) &
!              +  elec*(a45_1e_he + a45_2e_he / ((aR(1)*1d9)**a45_3e_he * (aR(2)*1d9)**a45_4e_he))
Ham_dir(4,5) = elec*(a27_1d_he + a27_2d_he / ((aR(2)*1d9)**a27_3d_he * (aR(1)*1d9)**a27_4d_he)) 
Ham_ex(4,5)  = elec*(a27_1e_he + a27_2e_he / ((aR(2)*1d9)**a27_3e_he * (aR(1)*1d9)**a27_4e_he))

!Ham(4,6) =  - 1.0 * elec*(a46_1d_he + a46_2d_he / ((aR(1)*1d9)**a46_3d_he * (aR(2)*1d9)**a46_4d_he)) &
!              +  elec*(a46_1e_he + a46_2e_he / ((aR(1)*1d9)**a46_3e_he * (aR(2)*1d9)**a46_4e_he))
Ham_dir(4,6) = elec*(a28_1d_he + a28_2d_he / ((aR(2)*1d9)**a28_3d_he * (aR(1)*1d9)**a28_4d_he)) 
Ham_ex(4,6)  = elec*(a28_1e_he + a28_2e_he / ((aR(2)*1d9)**a28_3e_he * (aR(1)*1d9)**a28_4e_he))

!Ham(4,7) =  - 1.0 * elec*(a47_1d_he + a47_2d_he / ((aR(1)*1d9)**a47_3d_he * (aR(2)*1d9)**a47_4d_he)) &
!              +  elec*(a47_1e_he + a47_2e_he / ((aR(1)*1d9)**a47_3e_he * (aR(2)*1d9)**a47_4e_he))
Ham_dir(4,7) = elec*(a25_1d_he + a25_2d_he / ((aR(2)*1d9)**a25_3d_he * (aR(1)*1d9)**a25_4d_he)) 
Ham_ex(4,7)  = elec*(a25_1e_he + a25_2e_he / ((aR(2)*1d9)**a25_3e_he * (aR(1)*1d9)**a25_4e_he))

!Ham(4,8) =  - 1.0 * elec*(a48_1d_he + a48_2d_he / ((aR(1)*1d9)**a48_3d_he * (aR(2)*1d9)**a48_4d_he)) &
!              +  elec*(a48_1e_he + a48_2e_he / ((aR(1)*1d9)**a48_3e_he * (aR(2)*1d9)**a48_4e_he))
Ham_dir(4,8) = elec*(a26_1d_he + a26_2d_he / ((aR(2)*1d9)**a26_3d_he * (aR(1)*1d9)**a26_4d_he)) 
Ham_ex(4,8)  = elec*(a26_1e_he + a26_2e_he / ((aR(2)*1d9)**a26_3e_he * (aR(1)*1d9)**a26_4e_he))

Ham_dir(5,6) = elec*(a56_1d_he + a56_2d_he / ((aR(1)*1d9)**a56_3d_he * (aR(2)*1d9)**a56_4d_he)) 
Ham_ex(5,6)  = elec*(a56_1e_he + a56_2e_he / ((aR(1)*1d9)**a56_3e_he * (aR(2)*1d9)**a56_4e_he))

Ham_dir(5,7) = elec*(a57_1d_he + a57_2d_he / ((aR(1)*1d9)**a57_3d_he * (aR(2)*1d9)**a57_4d_he))  
Ham_ex(5,7)  = elec*(a57_1e_he + a57_2e_he / ((aR(1)*1d9)**a57_3e_he * (aR(2)*1d9)**a57_4e_he))   

Ham_dir(5,8) = elec*(a58_1d_he + a58_2d_he / ((aR(1)*1d9)**a58_3d_he * (aR(2)*1d9)**a58_4d_he))  
Ham_ex(5,8)  = elec*(a58_1e_he + a58_2e_he / ((aR(1)*1d9)**a58_3e_he * (aR(2)*1d9)**a58_4e_he))  

!Ham(6,7) = - 1.0 * elec*(a67_1d_he + a67_2d_he / ((aR(1)*1d9)**a67_3d_he * (aR(2)*1d9)**a67_4d_he)) & 
!             +  elec*(a67_1e_he + a67_2e_he / ((aR(1)*1d9)**a67_3e_he * (aR(2)*1d9)**a67_4e_he))   
Ham_dir(6,7) = elec*(a58_1d_he + a58_2d_he / ((aR(2)*1d9)**a58_3d_he * (aR(1)*1d9)**a58_4d_he))  
Ham_ex(6,7)  = elec*(a58_1e_he + a58_2e_he / ((aR(2)*1d9)**a58_3e_he * (aR(1)*1d9)**a58_4e_he))   

Ham_dir(6,8) = elec*(a68_1d_he + a68_2d_he / ((aR(1)*1d9)**a68_3d_he * (aR(2)*1d9)**a68_4d_he))  
Ham_ex(6,8)  = elec*(a68_1e_he + a68_2e_he / ((aR(1)*1d9)**a68_3e_he * (aR(2)*1d9)**a68_4e_he))   

!Ham(7,8) = - 1.0 * elec*(a78_1d_he + a78_2d_he / ((aR(1)*1d9)**a78_3d_he * (aR(2)*1d9)**a78_4d_he)) &
!             +  elec*(a78_1e_he + a78_2e_he / ((aR(1)*1d9)**a78_3e_he * (aR(2)*1d9)**a78_4e_he))
Ham_dir(7,8) = elec*(a56_1d_he + a56_2d_he / ((aR(2)*1d9)**a56_3d_he * (aR(1)*1d9)**a56_4d_he)) 
Ham_ex(7,8)  = elec*(a56_1e_he + a56_2e_he / ((aR(2)*1d9)**a56_3e_he * (aR(1)*1d9)**a56_4e_he))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1.d0 * Ham_dir(i,j) + Ham_ex(i,j)
    endif
  enddo
enddo

do i=1,nstates-1
  do j=i+1,nstates-1
    Ham(j,i) = Ham(i,j)
    Ham_dir(j,i) = Ham_dir(i,j)
    Ham_ex(j,i) = Ham_ex(i,j)
  enddo
enddo

else if ( vers .eq. 'randm') then

Ham     = 0.0
Ham_0   = 0.0
Ham_dir = 0.0
Ham_ex  = 0.0

Ham_0(1)   = minEe(1,n) + minEh(1,n)  + V0 
Ham_dir(1,1) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n)*1d9)**a11_3d_ho)) 
Ham_ex(1,1)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n)*1d9)**a11_3e_ho))

Ham_0(2)   = minEe(1,n) + minEh(2,n) + V0 
Ham_dir(2,2) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n)*1d9)**a22_3d_ho))
Ham_ex(2,2)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n)*1d9)**a22_3e_ho))

Ham_0(3)   = minEe(1,n+nsys) + minEh(1,n+nsys) + V0 
Ham_dir(3,3) = elec*(a11_1d_ho + a11_2d_ho / ((aR(n+nsys)*1d9)**a11_3d_ho))
Ham_ex(3,3)  = elec*(a11_1e_ho + a11_2e_ho / ((aR(n+nsys)*1d9)**a11_3e_ho))

Ham_0(4)   = minEe(1,n+nsys) + minEh(2,n+nsys) + V0 
Ham_dir(4,4) = elec*(a22_1d_ho + a22_2d_ho / ((aR(n+nsys)*1d9)**a22_3d_ho))
Ham_ex(4,4)  = elec*(a22_1e_ho + a22_2e_ho / ((aR(n+nsys)*1d9)**a22_3e_ho))

Ham_0(5)   = minEe(1,n+nsys) + minEh(1,n) + V0 
Ham_dir(5,5) = elec*(a55_1d_he + a55_2d_he / ((aR(n)*1d9)**a55_3d_he * (aR(n+nsys)*1d9)**a55_4d_he)) 
Ham_ex(5,5)  = elec*(a55_1e_he + a55_2e_he / ((aR(n)*1d9)**a55_3e_he * (aR(n+nsys)*1d9)**a55_4e_he))
                                  
Ham_0(6)   = minEe(1,n+nsys) + minEh(2,n) + V0 
Ham_dir(6,6) = elec*(a66_1d_he + a66_2d_he / ((aR(n)*1d9)**a66_3d_he * (aR(n+nsys)*1d9)**a66_4d_he)) 
Ham_ex(6,6)  = elec*(a66_1e_he + a66_2e_he / ((aR(n)*1d9)**a66_3e_he * (aR(n+nsys)*1d9)**a66_4e_he)) 
                                  
!Ham(7,7) = minEe(1,n) + minEh(1,n+nsys) + V0 &
!  - 1.0 * elec*(a77_1d_he + a77_2d_he / ((aR(n)*1d9)**a77_3d_he * (aR(n+nsys)*1d9)**a77_4d_he)) &
!  +  elec*(a77_1e_he + a77_2e_he / ((aR(n)*1d9)**a77_3e_he * (aR(n+nsys)*1d9)**a77_4e_he))
Ham_0(7)   = minEe(1,n) + minEh(1,n+nsys) + V0 
Ham_dir(7,7) = elec*(a55_1d_he + a55_2d_he / ((aR(n+nsys)*1d9)**a55_3d_he * (aR(n)*1d9)**a55_4d_he)) 
Ham_ex(7,7)  = elec*(a55_1e_he + a55_2e_he / ((aR(n+nsys)*1d9)**a55_3e_he * (aR(n)*1d9)**a55_4e_he))

!Ham(8,8) = minEe(1,n) + minEh(2,n+nsys) + V0 &
! - 1.0 * elec*(a88_1d_he + a88_2d_he / ((aR(n)*1d9)**a88_3d_he * (aR(n+nsys)*1d9)**a88_4d_he)) &
!  +  elec*(a88_1e_he + a88_2e_he / ((aR(n)*1d9)**a88_3e_he * (aR(n+nsys)*1d9)**a88_4e_he))
Ham_0(8)   = minEe(1,n) + minEh(2,n+nsys) + V0 
Ham_dir(8,8) = elec*(a66_1d_he + a66_2d_he / ((aR(n+nsys)*1d9)**a66_3d_he * (aR(n)*1d9)**a66_4d_he)) 
Ham_ex(8,8)  = elec*(a66_1e_he + a66_2e_he / ((aR(n+nsys)*1d9)**a66_3e_he * (aR(n)*1d9)**a66_4e_he))

Ham_dir(1,2) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n)*1d9)**a12_3d_ho)) 
Ham_ex(1,2)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n)*1d9)**a12_3e_ho))

Ham_dir(1,3) = elec*(a13_1d_he + a13_2d_he / ((aR(n)*1d9)**a13_3d_he * (aR(n+nsys)*1d9)**a13_4d_he))
Ham_ex(1,3)  = elec*(a13_1e_he + a13_2e_he / ((aR(n)*1d9)**a13_3e_he * (aR(n+nsys)*1d9)**a13_4e_he))

Ham_dir(1,4) = elec*(a14_1d_he + a14_2d_he / ((aR(n)*1d9)**a14_3d_he * (aR(n+nsys)*1d9)**a14_4d_he))
Ham_ex(1,4)  = elec*(a14_1e_he + a14_2e_he / ((aR(n)*1d9)**a14_3e_he * (aR(n+nsys)*1d9)**a14_4e_he))

Ham_dir(1,5) = elec*(a15_1d_he + a15_2d_he / ((aR(n)*1d9)**a15_3d_he * (aR(n+nsys)*1d9)**a15_4d_he))
Ham_ex(1,5)  = elec*(a15_1e_he + a15_2e_he / ((aR(n)*1d9)**a15_3e_he * (aR(n+nsys)*1d9)**a15_4e_he))

Ham_dir(1,6) = elec*(a16_1d_he + a16_2d_he / ((aR(n)*1d9)**a16_3d_he * (aR(n+nsys)*1d9)**a16_4d_he))
Ham_ex(1,6)  = elec*(a16_1e_he + a16_2e_he / ((aR(n)*1d9)**a16_3e_he * (aR(n+nsys)*1d9)**a16_4e_he))

Ham_dir(1,7) = elec*(a17_1d_he + a17_2d_he / ((aR(n)*1d9)**a17_3d_he * (aR(n+nsys)*1d9)**a17_4d_he))
Ham_ex(1,7)  = elec*(a17_1e_he + a17_2e_he / ((aR(n)*1d9)**a17_3e_he * (aR(n+nsys)*1d9)**a17_4e_he))

Ham_dir(1,8) = elec*(a18_1d_he + a18_2d_he / ((aR(n)*1d9)**a18_3d_he * (aR(n+nsys)*1d9)**a18_4d_he))
Ham_ex(1,8)  = elec*(a18_1e_he + a18_2e_he / ((aR(n)*1d9)**a18_3e_he * (aR(n+nsys)*1d9)**a18_4e_he))

!Ham_dir(2,3) = - 1.0 * elec*(a23_1d_he + a23_2d_he / ((aR(n)*1d9)**a23_3d_he * (aR(n+nsys)*1d9)**a23_4d_he)) &
!             +  elec*(a23_1e_he + a23_2e_he / ((aR(n)*1d9)**a23_3e_he * (aR(n+nsys)*1d9)**a23_4e_he))
Ham_dir(2,3) = elec*(a14_1d_he + a14_2d_he / ((aR(n+nsys)*1d9)**a14_3d_he * (aR(n)*1d9)**a14_4d_he))
Ham_ex(2,3)  = elec*(a14_1e_he + a14_2e_he / ((aR(n+nsys)*1d9)**a14_3e_he * (aR(n)*1d9)**a14_4e_he))

Ham_dir(2,4) = elec*(a24_1d_he + a24_2d_he / ((aR(n)*1d9)**a24_3d_he * (aR(n+nsys)*1d9)**a24_4d_he))
Ham_ex(2,4)  = elec*(a24_1e_he + a24_2e_he / ((aR(n)*1d9)**a24_3e_he * (aR(n+nsys)*1d9)**a24_4e_he))

Ham_dir(2,5) = elec*(a25_1d_he + a25_2d_he / ((aR(n)*1d9)**a25_3d_he * (aR(n+nsys)*1d9)**a25_4d_he))
Ham_ex(2,5)  = elec*(a25_1e_he + a25_2e_he / ((aR(n)*1d9)**a25_3e_he * (aR(n+nsys)*1d9)**a25_4e_he))

Ham_dir(2,6) = elec*(a26_1d_he + a26_2d_he / ((aR(n)*1d9)**a26_3d_he * (aR(n+nsys)*1d9)**a26_4d_he))
Ham_ex(2,6)  = elec*(a26_1e_he + a26_2e_he / ((aR(n)*1d9)**a26_3e_he * (aR(n+nsys)*1d9)**a26_4e_he))

Ham_dir(2,7) = elec*(a27_1d_he + a27_2d_he / ((aR(n)*1d9)**a27_3d_he * (aR(n+nsys)*1d9)**a27_4d_he))
Ham_ex(2,7)  = elec*(a27_1e_he + a27_2e_he / ((aR(n)*1d9)**a27_3e_he * (aR(n+nsys)*1d9)**a27_4e_he))

Ham_dir(2,8) = elec*(a28_1d_he + a28_2d_he / ((aR(n)*1d9)**a28_3d_he * (aR(n+nsys)*1d9)**a28_4d_he))
Ham_ex(2,8)  = elec*(a28_1e_he + a28_2e_he / ((aR(n)*1d9)**a28_3e_he * (aR(n+nsys)*1d9)**a28_4e_he))

Ham_dir(3,4) = elec*(a12_1d_ho + a12_2d_ho / ((aR(n+nsys)*1d9)**a12_3d_ho)) 
Ham_ex(3,4)  = elec*(a12_1e_ho + a12_2e_ho / ((aR(n+nsys)*1d9)**a12_3e_ho))

!Ham_dir(3,5) =  - 1.0 * elec*(a35_1d_he + a35_2d_he / ((aR(n)*1d9)**a35_3d_he * (aR(n+nsys)*1d9)**a35_4d_he)) &
!              +  elec*(a35_1e_he + a35_2e_he / ((aR(n)*1d9)**a35_3e_he * (aR(n+nsys)*1d9)**a35_4e_he))
Ham_dir(3,5) = elec*(a17_1d_he + a17_2d_he / ((aR(n+nsys)*1d9)**a17_3d_he * (aR(n)*1d9)**a17_4d_he)) 
Ham_ex(3,5)  = elec*(a17_1e_he + a17_2e_he / ((aR(n+nsys)*1d9)**a17_3e_he * (aR(n)*1d9)**a17_4e_he))

!Ham_dir(3,6) =  - 1.0 * elec*(a36_1d_he + a36_2d_he / ((aR(n)*1d9)**a36_3d_he * (aR(n+nsys)*1d9)**a36_4d_he)) &
!              +  elec*(a36_1e_he + a36_2e_he / ((aR(n)*1d9)**a36_3e_he * (aR(n+nsys)*1d9)**a36_4e_he))
Ham_dir(3,6) = elec*(a18_1d_he + a18_2d_he / ((aR(n+nsys)*1d9)**a18_3d_he * (aR(n)*1d9)**a18_4d_he)) 
Ham_ex(3,6)  = elec*(a18_1e_he + a18_2e_he / ((aR(n+nsys)*1d9)**a18_3e_he * (aR(n)*1d9)**a18_4e_he))

!Ham_dir(3,7) =  - 1.0 * elec*(a37_1d_he + a37_2d_he / ((aR(n)*1d9)**a37_3d_he * (aR(n+nsys)*1d9)**a37_4d_he)) &
!              +  elec*(a37_1e_he + a37_2e_he / ((aR(n)*1d9)**a37_3e_he * (aR(n+nsys)*1d9)**a37_4e_he))
Ham_dir(3,7) = elec*(a15_1d_he + a15_2d_he / ((aR(n+nsys)*1d9)**a15_3d_he * (aR(n)*1d9)**a15_4d_he)) 
Ham_ex(3,7)  = elec*(a15_1e_he + a15_2e_he / ((aR(n+nsys)*1d9)**a15_3e_he * (aR(n)*1d9)**a15_4e_he))

!Ham_dir(3,8) =  - 1.0 * elec*(a38_1d_he + a38_2d_he / ((aR(n)*1d9)**a38_3d_he * (aR(n+nsys)*1d9)**a38_4d_he)) &
!              +  elec*(a38_1e_he + a38_2e_he / ((aR(n)*1d9)**a38_3e_he * (aR(n+nsys)*1d9)**a38_4e_he))
Ham_dir(3,8) = elec*(a16_1d_he + a16_2d_he / ((aR(n+nsys)*1d9)**a16_3d_he * (aR(n)*1d9)**a16_4d_he)) 
Ham_ex(3,8)  = elec*(a16_1e_he + a16_2e_he / ((aR(n+nsys)*1d9)**a16_3e_he * (aR(n)*1d9)**a16_4e_he))

!Ham_dir(4,5) =  - 1.0 * elec*(a45_1d_he + a45_2d_he / ((aR(n)*1d9)**a45_3d_he * (aR(n+nsys)*1d9)**a45_4d_he)) &
!              +  elec*(a45_1e_he + a45_2e_he / ((aR(n)*1d9)**a45_3e_he * (aR(n+nsys)*1d9)**a45_4e_he))
Ham_dir(4,5) = elec*(a27_1d_he + a27_2d_he / ((aR(n+nsys)*1d9)**a27_3d_he * (aR(n)*1d9)**a27_4d_he)) 
Ham_ex(4,5)  = elec*(a27_1e_he + a27_2e_he / ((aR(n+nsys)*1d9)**a27_3e_he * (aR(n)*1d9)**a27_4e_he))

!Ham_dir(4,6) =  - 1.0 * elec*(a46_1d_he + a46_2d_he / ((aR(n)*1d9)**a46_3d_he * (aR(n+nsys)*1d9)**a46_4d_he)) &
!              +  elec*(a46_1e_he + a46_2e_he / ((aR(n)*1d9)**a46_3e_he * (aR(n+nsys)*1d9)**a46_4e_he))
Ham_dir(4,6) = elec*(a28_1d_he + a28_2d_he / ((aR(n+nsys)*1d9)**a28_3d_he * (aR(n)*1d9)**a28_4d_he)) 
Ham_ex(4,6)  = elec*(a28_1e_he + a28_2e_he / ((aR(n+nsys)*1d9)**a28_3e_he * (aR(n)*1d9)**a28_4e_he))

!Ham_dir(4,7) =  - 1.0 * elec*(a47_1d_he + a47_2d_he / ((aR(n)*1d9)**a47_3d_he * (aR(n+nsys)*1d9)**a47_4d_he)) &
!              +  elec*(a47_1e_he + a47_2e_he / ((aR(n)*1d9)**a47_3e_he * (aR(n+nsys)*1d9)**a47_4e_he))
Ham_dir(4,7) = elec*(a25_1d_he + a25_2d_he / ((aR(n+nsys)*1d9)**a25_3d_he * (aR(n)*1d9)**a25_4d_he)) 
Ham_ex(4,7)  = elec*(a25_1e_he + a25_2e_he / ((aR(n+nsys)*1d9)**a25_3e_he * (aR(n)*1d9)**a25_4e_he))

!Ham_dir(4,8) =  - 1.0 * elec*(a48_1d_he + a48_2d_he / ((aR(n)*1d9)**a48_3d_he * (aR(n+nsys)*1d9)**a48_4d_he)) &
!              +  elec*(a48_1e_he + a48_2e_he / ((aR(n)*1d9)**a48_3e_he * (aR(n+nsys)*1d9)**a48_4e_he))
Ham_dir(4,8) = elec*(a26_1d_he + a26_2d_he / ((aR(n+nsys)*1d9)**a26_3d_he * (aR(n)*1d9)**a26_4d_he)) 
Ham_ex(4,8)  = elec*(a26_1e_he + a26_2e_he / ((aR(n+nsys)*1d9)**a26_3e_he * (aR(n)*1d9)**a26_4e_he))

Ham_dir(5,6) = elec*(a56_1d_he + a56_2d_he / ((aR(n)*1d9)**a56_3d_he * (aR(n+nsys)*1d9)**a56_4d_he)) 
Ham_ex(5,6)  = elec*(a56_1e_he + a56_2e_he / ((aR(n)*1d9)**a56_3e_he * (aR(n+nsys)*1d9)**a56_4e_he))

Ham_dir(5,7) = elec*(a57_1d_he + a57_2d_he / ((aR(n)*1d9)**a57_3d_he * (aR(n+nsys)*1d9)**a57_4d_he))  
Ham_ex(5,7)  = elec*(a57_1e_he + a57_2e_he / ((aR(n)*1d9)**a57_3e_he * (aR(n+nsys)*1d9)**a57_4e_he))   

Ham_dir(5,8) = elec*(a58_1d_he + a58_2d_he / ((aR(n)*1d9)**a58_3d_he * (aR(n+nsys)*1d9)**a58_4d_he))  
Ham_ex(5,8)  = elec*(a58_1e_he + a58_2e_he / ((aR(n)*1d9)**a58_3e_he * (aR(n+nsys)*1d9)**a58_4e_he))   

!Ham_dir(6,7) = - 1.0 * elec*(a67_1d_he + a67_2d_he / ((aR(n)*1d9)**a67_3d_he * (aR(n+nsys)*1d9)**a67_4d_he)) & 
!             +  elec*(a67_1e_he + a67_2e_he / ((aR(n)*1d9)**a67_3e_he * (aR(n+nsys)*1d9)**a67_4e_he))   
Ham_dir(6,7) = elec*(a58_1d_he + a58_2d_he / ((aR(n+nsys)*1d9)**a58_3d_he * (aR(n)*1d9)**a58_4d_he))  
Ham_ex(6,7)  = elec*(a58_1e_he + a58_2e_he / ((aR(n+nsys)*1d9)**a58_3e_he * (aR(n)*1d9)**a58_4e_he))   

Ham_dir(6,8) = elec*(a68_1d_he + a68_2d_he / ((aR(n)*1d9)**a68_3d_he * (aR(n+nsys)*1d9)**a68_4d_he)) 
Ham_ex(6,8)  = elec*(a68_1e_he + a68_2e_he / ((aR(n)*1d9)**a68_3e_he * (aR(n+nsys)*1d9)**a68_4e_he))   

!Ham_dir(7,8) = - 1.0 * elec*(a78_1d_he + a78_2d_he / ((aR(n)*1d9)**a78_3d_he * (aR(n+nsys)*1d9)**a78_4d_he)) &
!             +  elec*(a78_1e_he + a78_2e_he / ((aR(n)*1d9)**a78_3e_he * (aR(n+nsys)*1d9)**a78_4e_he))
Ham_dir(7,8) = elec*(a56_1d_he + a56_2d_he / ((aR(n+nsys)*1d9)**a56_3d_he * (aR(n)*1d9)**a56_4d_he)) 
Ham_ex(7,8)  = elec*(a78_1e_he + a78_2e_he / ((aR(n+nsys)*1d9)**a78_3e_he * (aR(n)*1d9)**a78_4e_he))

do i=1,nstates-1
  do j=1,nstates-1
    if ( i .eq. j ) then
    Ham(i,j) = Ham_0(i) - Ham_dir(i,j) + Ham_ex(i,j)
    elseif ( i .ne. j ) then
    Ham(i,j) = -1.d0 * Ham_dir(i,j) + Ham_ex(i,j)
    endif
  enddo
enddo

do i=1,nstates-1
  do j=i+1,nstates-1
    Ham(j,i) = Ham(i,j)
    Ham_dir(j,i) = Ham_dir(i,j)
    Ham_ex(j,i) = Ham_ex(i,j)
  enddo
enddo

endif

Ham = Ham/Energ_au
Ham_dir = Ham_dir/Energ_au
Ham_ex = Ham_ex/Energ_au
TransHam(0,1) = TransDip_Ana_h1e(n)
TransHam(0,2) = TransDip_Ana_h2e(n)
TransHam(0,3) = TransDip_Ana_h1e(n+nsys)
TransHam(0,4) = TransDip_Ana_h2e(n+nsys)
TransHam(0,5) = TransDip_Fit_h1e_he(aR(n+nsys),aR(n))
TransHam(0,6) = TransDip_Fit_h2e_he(aR(n+nsys),aR(n))
TransHam(0,7) = TransDip_Fit_h1e_he(aR(n),aR(n+nsys))
TransHam(0,8) = TransDip_Fit_h2e_he(aR(n),aR(n+nsys))
TransHam = TransHam*D_to_au

do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

end subroutine make_Ham_he

subroutine make_Ham_fineSt

Ham = 0.d0

Ham(1,1)   = Eeh1(n)  - Dso1/3.d0 - 3.d0*Kas
Ham(2,2)   = Eeh1(n)  - Dso1/3.d0 - 3.d0*Kcs
Ham(3,3)   = Eeh1(n)  - Kas
Ham(4,4)   = Eeh1(n)  - 3.d0*Kas
Ham(5,5)   = Eeh1(n)  - 3.d0*Kbs + Dxf
Ham(6,6)   = Eeh1(n)  - 3.d0*Kbs + Dxf
Ham(7,7)   = Eeh1(n)  - Kcs
Ham(8,8)   = Eeh1(n)  - 3.d0*Kcs
Ham(9,9)   = Eeh1(n)  + Dso1/3.d0 - 3.d0*Kas
Ham(10,10) = Eeh1(n)  - Kbs + Dxf
Ham(11,11) = Eeh1(n)  - 3.d0*Kbs + Dxf
Ham(12,12) = Eeh1(n)  + Dso1/3.d0 - 3.d0*Kcs

Ham(3,4) = -1.d0*Dso1/3.d0
Ham(4,3) = Ham(3,4)
Ham(3,5) = Ham(3,4)
Ham(5,3) = Ham(3,4)
Ham(6,7) = Ham(3,4)
Ham(7,6) = Ham(3,4)
Ham(6,8) = Ham(3,4)
Ham(8,6) = Ham(3,4)
Ham(9,10) = Ham(3,4)
Ham(10,9) = Ham(3,4)
Ham(9,11) = Ham(3,4)
Ham(11,9) = Ham(3,4)
Ham(10,12) = Ham(3,4)
Ham(12,10) = Ham(3,4)

Ham(4,5) = Dso1/3.d0
Ham(5,4) = Ham(4,5)
Ham(7,8) = Ham(4,5)
Ham(8,7) = Ham(4,5)
Ham(11,12) = Ham(4,5)
Ham(12,11) = Ham(4,5)

Ham(13,13)   = Eeh2(n)  - Dso2/3.d0 - 3.d0*Kas
Ham(14,14)   = Eeh2(n)  - Dso2/3.d0 - 3.d0*Kcs
Ham(15,15)   = Eeh2(n)  - Kas
Ham(16,16)   = Eeh2(n)  - 3.d0*Kas
Ham(17,17)   = Eeh2(n)  - 3*Kbs + Dxf
Ham(18,18)   = Eeh2(n)  - 3*Kbs + Dxf
Ham(19,19)   = Eeh2(n)  - Kcs
Ham(20,20)   = Eeh2(n)  - 3.d0*Kcs
Ham(21,21)   = Eeh2(n)  + Dso2/3.d0 - 3.d0*Kas
Ham(22,22)   = Eeh2(n)  - Kbs + Dxf
Ham(23,23)   = Eeh2(n)  - 3.d0*Kbs + Dxf
Ham(24,24)   = Eeh2(n)  + Dso2/3.d0 - 3.d0*Kcs

Ham(15,16) = -1.d0*Dso2/3.d0
Ham(16,15) = Ham(15,16)
Ham(15,17) = Ham(15,16)
Ham(17,15) = Ham(15,16)
Ham(18,19) = Ham(15,16)
Ham(19,18) = Ham(15,16)
Ham(18,20) = Ham(15,16)
Ham(20,18) = Ham(15,16)
Ham(21,22) = Ham(15,16)
Ham(22,21) = Ham(15,16)
Ham(21,23) = Ham(15,16)
Ham(23,21) = Ham(15,16)
Ham(22,24) = Ham(15,16)
Ham(24,22) = Ham(15,16)

Ham(16,17) = Dso2/3
Ham(17,16) = Ham(16,17)
Ham(19,20) = Ham(16,17)
Ham(20,19) = Ham(16,17)
Ham(23,24) = Ham(16,17)
Ham(24,23) = Ham(16,17)

Ham = Ham/Energ_au

end subroutine make_Ham_fineSt

subroutine make_TransHam_0_fineSt

TransHam0 = 0.d0

TransHam0(0,3)  = abs(TransDip_Ana_h1e(n))
TransHam0(0,7)  = abs(TransDip_Ana_h1e(n)) 
TransHam0(0,10) = abs(TransDip_Ana_h1e(n))
TransHam0(0,15) = abs(TransDip_Ana_h2e(n))
TransHam0(0,19) = abs(TransDip_Ana_h2e(n))
TransHam0(0,22) = abs(TransDip_Ana_h2e(n))

TransHam0 = TransHam0/Dip_au

do i=0,nstates-1
TransHam0(i,0) = TransHam0(0,i)
enddo

end subroutine make_TransHam_0_fineSt

subroutine make_TransHam_ei_fineSt

TransHam = 0.d0

do i=0,nstates-1
do j=0,nstates-1
TransHam(0,i) = TransHam(0,i) +  TransHam0(0,j) * Ham_ei(j,i)
enddo
enddo

TransHam = abs(TransHam)

write(6,*) (abs(TransHam(0,i)*Dip_au/Cm_to_D), i=0,nstates-1)

do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

Ham = 0.d0

do i=0,nstates-1
Ham(i,i) = lambda(i)
enddo

end subroutine make_TransHam_ei_fineSt

end module Make_Ham
