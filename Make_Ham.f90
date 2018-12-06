module Make_Ham

use omp_lib
use Constants_au
use Variables_au
use Integrals

real(dp), parameter :: a12_1d_ho = 0.00876695d0
real(dp), parameter :: a12_2d_ho = 0.0346284d0
real(dp), parameter :: a12_3d_ho = 8.94293d+08
real(dp), parameter :: a12_1e_ho = 0.000588561d0
real(dp), parameter :: a12_2e_ho = 0.231987d0
real(dp), parameter :: a12_3e_ho = 1.64474d+09
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
Ham(1,1) = Eeh1(1) - Cb_eh1(1)
Ham(2,2) = Eeh2(1) - Cb_eh2(1)                                                                         
Ham(3,3) = Ham(1,1)
Ham(4,4) = Ham(2,2)
Ham(5,5) = Eeh1(1) - elec*(a55_1d_ho + a55_2d_ho * exp(-1.0d0*(a55_3d_ho*aR(1) + a55_4d_ho*link))) &
                                  + 2*elec*(a55_1e_ho + a55_2e_ho * exp(-1.0d0*(a55_3e_ho*aR(1) + a55_4e_ho*link)))
Ham(6,6) = Eeh2(2) - elec*(a66_1d_ho + a66_2d_ho * exp(-1.0d0*(a66_3d_ho*aR(1) + a66_4d_ho*link))) &
                                  + 2*elec*(a66_1e_ho + a66_2e_ho * exp(-1.0d0*(a66_3e_ho*aR(1) + a66_4e_ho*link)))
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

Ham(1,2) = - 1.0d0 * elec*(a12_1d_ho + a12_2d_ho / exp(aR(1)*a12_3d_ho)) &
             - 2.0d0 * elec*(a12_1e_ho + a12_2e_ho / exp(aR(1)*a12_3e_ho))

Ham(1,3) = - 1.0d0 * elec*(a13_1d_ho + a13_2d_ho * exp(-1.0d0*(a13_3d_ho*aR(1) + a13_4d_ho*link))) &
             + 2.0d0 * elec*(a13_1e_ho + a13_2e_ho * exp(-1.0d0*(a13_3e_ho*aR(1) + a13_4e_ho*link)))

Ham(1,4) = - 1.0d0 * elec*(a14_1d_ho + a14_2d_ho * exp(-1.0d0*(a14_3d_ho*aR(1) + a14_4d_ho*link))) &
             + 2.0d0 * elec*(a14_1e_ho + a14_2e_ho * exp(-1.0d0*(a14_3e_ho*aR(1) + a14_4e_ho*link)))

Ham(1,5) = - 1.0d0 * elec*(a15_1d_ho + (a15_2d_ho/aR(1))**a15_3d_ho * exp(-1.0d0*a15_4d_ho*link)) &
             + 2.0d0 * elec*(a15_1e_ho + (a15_2e_ho/aR(1))**a15_3e_ho * exp(-1.0d0*a15_4e_ho*link))

Ham(1,6) = - 1.0d0 *( -1.0d0 * elec * (a16_1d_ho + (a16_2d_ho/aR(1))**a16_3d_ho * exp(-1.0d0*a16_4d_ho*link)) &
                     + 2.0d0 * elec * (a16_1e_ho + (a16_2e_ho/aR(1))**a16_3e_ho * exp(-1.0d0*a16_4e_ho*link)))

Ham(1,7) = - 1.0d0 * elec*(a17_1d_ho + (a17_2d_ho/aR(1))**a17_3d_ho * exp(-1.0d0*a17_4d_ho*link)) &
             + 2.0d0 * elec*(a17_1e_ho + (a17_2e_ho/aR(1))**a17_3e_ho * exp(-1.0d0*a17_4e_ho*link))

Ham(1,8) =  - 1.0d0 *( -1.0d0 * elec * (a18_1d_ho + (a18_2d_ho/aR(1))**a18_3d_ho * exp(-1.0d0*a18_4d_ho*link)) &
                      + 2.0d0 * elec * (a18_1e_ho + (a18_2e_ho/aR(1))**a18_3e_ho * exp(-1.0d0*a18_4e_ho*link)))

Ham(2,3) = Ham(1,4)

Ham(2,4) = - 1.0d0 * elec*(a24_1d_ho + a24_2d_ho * exp(-1.0d0*(a24_3d_ho*aR(1) + a24_4d_ho*link))) &
             + 2.0d0 * elec*(a24_1e_ho + a24_2e_ho * exp(-1.0d0*(a24_3e_ho*aR(1) + a24_4e_ho*link)))

Ham(2,5) =  - 1.0d0 *( -1.0d0 * elec*(a25_1d_ho + (a25_2d_ho/aR(1))**a25_3d_ho * exp(-1.0d0*a25_4d_ho*link)) &
                      + 2.0d0 * elec*(a25_1e_ho + (a25_2e_ho/aR(1))**a25_3e_ho * exp(-1.0d0*a25_4e_ho*link)))

Ham(2,6) =  - 1.0d0 * elec*(a26_1d_ho + (a26_2d_ho/aR(1))**a26_3d_ho * exp(-1.0d0*a26_4d_ho*link)) &
              + 2.0d0 * elec*(a26_1e_ho + (a26_2e_ho/aR(1))**a26_3e_ho * exp(-1.0d0*a26_4e_ho*link))

Ham(2,7) =  - 1.0d0 *( -1.0d0 * elec * (a27_1d_ho + (a27_2d_ho/aR(1))**a27_3d_ho * exp(-1.0d0*a27_4d_ho*link)) &
                      + 2.0d0 * elec * (a27_1e_ho + (a27_2e_ho/aR(1))**a27_3e_ho * exp(-1.0d0*a27_4e_ho*link)))

Ham(2,8) =  - 1.0d0 * elec*(a28_1d_ho + (a28_2d_ho/aR(1))**a28_3d_ho * exp(-1.0d0*a28_4d_ho*link)) &
              + 2.0d0 * elec*(a28_1e_ho + (a28_2e_ho/aR(1))**a28_3e_ho * exp(-1.0d0*a28_4e_ho*link))

Ham(3,4) = Ham(1,2)
Ham(3,5) = Ham(1,7)
Ham(3,6) = Ham(1,8)
Ham(3,7) = Ham(1,5)
Ham(3,8) = Ham(1,6)
Ham(4,5) = Ham(2,7)
Ham(4,6) = Ham(2,8)
Ham(4,7) = Ham(2,5)
Ham(4,8) = Ham(2,6)

Ham(5,6) = - 1.0d0 * elec*(a56_1d_ho + (a56_2d_ho/aR(1))**a56_3d_ho * exp(-1.0d0*a56_4d_ho*link)) &
             + 2.0d0 * elec*(a56_1e_ho + (a56_2e_ho/aR(1))**a56_3e_ho * exp(-1.0d0*a56_4e_ho*link))

Ham(5,7) = - 1.0d0 * elec*(a57_1d_ho + (a57_2d_ho/aR(1))**a57_3d_ho * exp(-1.0d0*a57_4d_ho*link)) &
             + 2.0d0 * elec*(a57_1e_ho + (a57_2e_ho/aR(1))**a57_3e_ho * exp(-1.0d0*a57_4e_ho*link))

Ham(5,8) = - 1.0d0 * elec*(a58_1d_ho + (a58_2d_ho/aR(1))**a58_3d_ho * exp(-1.0d0*a58_4d_ho*link)) &
             + 2.0d0 * elec*(a58_1e_ho + (a58_2e_ho/aR(1))**a58_3e_ho * exp(-1.0d0*a58_4e_ho*link))

Ham(6,7) = Ham(5,8)

Ham(6,8) = - 1.0d0 * elec*(a68_1d_ho + (a68_2d_ho/aR(1))**a68_3d_ho * exp(-1.0d0*a68_4d_ho*link)) &
             + 2.0d0 * elec*(a68_1e_ho + (a68_2e_ho/aR(1))**a68_3e_ho * exp(-1.0d0*a68_4e_ho*link))

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
TransHam(0,5) = TransDip_Fit_h1e_ho(aR(n),link)
TransHam(0,6) = TransDip_Fit_h2e_ho(aR(n),link)
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

if ( link .eq. 0.2e-9 ) then
a12_1d_he = 0.00876695
a12_2d_he = 0.0346284
a12_3d_he = 8.94293d+08
a12_1e_he = 0.000588561
a12_2e_he = 0.231987
a12_3e_he = 1.64474d+09
a13_1d_he = -5.7308d-10
a13_2d_he = 1.16996d-06
a13_3d_he = 3.64279
a13_4d_he = 3.80378
a15_1d_he = 3.32236d-06
a15_2d_he = 0.00204328
a15_3d_he = 2.77545
a15_4d_he = 1.85776
a17_1d_he = -3.56487d-08
a17_2d_he = 4.00073d-05
a17_3d_he = 2.68943
a17_4d_he = 2.0951
a24_1d_he = 2.14157d-09
a24_2d_he = 6.46666d-06
a24_3d_he = 3.95207
a24_4d_he = 4.00739
a26_1d_he = 4.13952d-06
a26_2d_he = 0.00206339
a26_3d_he = 2.80285
a26_4d_he = 1.85298
a28_1d_he = 5.16485d-08
a28_2d_he = 0.000216119
a28_3d_he = 2.85638
a28_4d_he = 2.27227
a35_1d_he = -4.63653d-08
a35_2d_he = 3.91962d-05
a35_3d_he = 2.68723
a35_4d_he = 2.0612
a37_1d_he = 3.69519d-06
a37_2d_he = 0.00204119
a37_3d_he = 2.77624
a37_4d_he = 1.86059
a46_1d_he = -3.1446d-07
a46_2d_he = 0.000210886
a46_3d_he = 2.80464
a46_4d_he = 2.20824
a48_1d_he = 6.30777d-06
a48_2d_he = 0.00211127
a48_3d_he = 2.87692
a48_4d_he = 1.85623
a55_1d_he = 0.0123067
a55_2d_he = 0.0834401
a55_3d_he = 0.751985
a55_4d_he = 0.822935
a57_1d_he = -2.18594d-10
a57_2d_he = 1.20859d-06
a57_3d_he = 3.72691
a57_4d_he = 3.81947
a66_1d_he = 0.0120069
a66_2d_he = 0.084655
a66_3d_he = 0.75655
a66_4d_he = 0.813347
a68_1d_he = 2.26161d-09
a68_2d_he = 6.64211d-06
a68_3d_he = 3.9102
a68_4d_he = 4.10851
a77_1d_he = 0.0131285
a77_2d_he = 0.0829806
a77_3d_he = 0.772324
a77_4d_he = 0.839392
a88_1d_he = 0.0110918
a88_2d_he = 0.0853286
a88_3d_he = 0.75124
a88_4d_he = 0.781124
a14_1d_he = -4.51135d-10
a14_2d_he = 2.84653d-06
a14_3d_he = 3.65651
a14_4d_he = 4.12674
a16_1d_he = 4.50443d-07
a16_2d_he = 2.55061d-05
a16_3d_he = 4.46738
a16_4d_he = 0.606741
a18_1d_he = 3.16132d-08
a18_2d_he = 9.86635d-05
a18_3d_he = 2.73099
a18_4d_he = 2.33682
a23_1d_he = 1.43653d-09
a23_2d_he = 2.99651d-06
a23_3d_he = 4.03847
a23_4d_he = 3.89306
a25_1d_he = 3.80318d-07
a25_2d_he = 3.90141d-05
a25_3d_he = 4.43569
a25_4d_he = 1.3408
a27_1d_he = -1.38977d-07
a27_2d_he = 9.05307d-05
a27_3d_he = 2.81297
a27_4d_he = 2.05185
a36_1d_he = -7.52562d-08
a36_2d_he = 9.46348d-05
a36_3d_he = 2.64458
a36_4d_he = 2.2994
a38_1d_he = -1.87329d-07
a38_2d_he = 3.50221d-05
a38_3d_he = 3.41031
a38_4d_he = 1.25104
a45_1d_he = -1.62686d-07
a45_2d_he = 9.0782d-05
a45_3d_he = 2.84149
a45_4d_he = 2.02033
a47_1d_he = 1.54893d-07
a47_2d_he = 2.28436d-05
a47_3d_he = 4.12838
a47_4d_he = 0.493992
a56_1d_he = 0.000250184
a56_2d_he = 0.00138359
a56_3d_he = 2.70362
a56_4d_he = -0.0663319
a58_1d_he = 1.50601d-09
a58_2d_he = 2.85302d-06
a58_3d_he = 3.72181
a58_4d_he = 4.19241
a67_1d_he = 5.44859d-10
a67_2d_he = 2.90479d-06
a67_3d_he = 3.98366
a67_4d_he = 3.83348
a78_1d_he = 0.000183386
a78_2d_he = 0.00107079
a78_3d_he = 2.28933
a78_4d_he = -0.466982
a13_1e_he = 0.0129891
a13_2e_he = 0.059974
a13_3e_he = 0.745223
a13_4e_he = 0.746591
a15_1e_he = 7.81765d-07
a15_2e_he = 0.000436249
a15_3e_he = 2.94985
a15_4e_he = 2.23921
a17_1e_he = 8.32306d-07
a17_2e_he = 0.000325249
a17_3e_he = 3.0237
a17_4e_he = 1.92494
a24_1e_he = 0.000165069
a24_2e_he = 0.0120975
a24_3e_he = 1.09414
a24_4e_he = 1.0831
a26_1e_he = -3.83814d-08
a26_2e_he = 0.000298365
a26_3e_he = 3.27584
a26_4e_he = 2.06783
a28_1e_he = 8.50008d-08
a28_2e_he = 0.000253351
a28_3e_he = 3.3815
a28_4e_he = 1.89094
a35_1e_he = 8.11209d-07
a35_2e_he = 0.000323139
a35_3e_he = 3.01545
a35_4e_he = 1.91859
a37_1e_he = 7.11155d-07
a37_2e_he = 0.000434065
a37_3e_he = 2.93827
a37_4e_he = 2.23353
a46_1e_he = 3.42974d-07
a46_2e_he = 0.00025947
a46_3e_he = 3.47192
a46_4e_he = 1.9116
a48_1e_he = -4.01333d-07
a48_2e_he = 0.000290705
a48_3e_he = 3.16886
a48_4e_he = 2.05449
a55_1e_he = 2.18641d-09
a55_2e_he = 2.66753d-06
a55_3e_he = 3.7282
a55_4e_he = 4.73166
a57_1e_he = 5.7703d-10
a57_2e_he = 2.01368d-06
a57_3e_he = 3.99743
a57_4e_he = 4.00244
a66_1e_he = -1.90295d-09
a66_2e_he = 8.64309d-06
a66_3e_he = 3.58118
a66_4e_he = 4.22789
a68_1e_he = 2.24073d-09
a68_2e_he = 8.0202d-06
a68_3e_he = 3.92005
a68_4e_he = 4.06342
a77_1e_he = 1.44698d-09
a77_2e_he = 2.5936d-06
a77_3e_he = 3.61919
a77_4e_he = 4.71832
a88_1e_he = -3.75173d-10
a88_2e_he = 8.84303d-06
a88_3e_he = 3.5768
a88_4e_he = 4.2987
a14_1e_he = 0.00100523
a14_2e_he = 0.0286904
a14_3e_he = 0.571985
a14_4e_he = 1.11753
a16_1e_he = 5.20996d-08
a16_2e_he = 0.000733493
a16_3e_he = 2.81241
a16_4e_he = 2.05309
a18_1e_he = 9.14038d-07
a18_2e_he = 0.000618908
a18_3e_he = 2.90339
a18_4e_he = 1.91276
a23_1e_he = 0.00072624
a23_2e_he = 0.0288513
a23_3e_he = 1.08402
a23_4e_he = 0.561365
a25_1e_he = -5.80522d-08
a25_2e_he = 0.000160715
a25_3e_he = 3.26161
a25_4e_he = 2.1577
a27_1e_he = 1.47805d-07
a27_2e_he = 0.000132624
a27_3e_he = 3.4939
a27_4e_he = 1.90809
a36_1e_he = 1.11843d-06
a36_2e_he = 0.00062124
a36_3e_he = 2.90829
a36_4e_he = 1.92675
a38_1e_he = 1.85001d-08
a38_2e_he = 0.000731357
a38_3e_he = 2.80983
a38_4e_he = 2.04924
a45_1e_he = 1.43075d-07
a45_2e_he = 0.000133518
a45_3e_he = 3.48729
a45_4e_he = 1.92688
a47_1e_he = 1.78003d-07
a47_2e_he = 0.000164668
a47_3e_he = 3.28665
a47_4e_he = 2.24905
a56_1e_he = 1.19765d-09
a56_2e_he = 4.7171d-06
a56_3e_he = 3.62877
a56_4e_he = 4.45655
a58_1e_he = 8.52181d-10
a58_2e_he = 4.03096d-06
a58_3e_he = 3.89968
a58_4e_he = 4.1013
a67_1e_he = -4.70846d-10
a67_2e_he = 4.045d-06
a67_3e_he = 4.05712
a67_4e_he = 3.87227
a78_1e_he = 3.51426d-10
a78_2e_he = 4.75278d-06
a78_3e_he = 3.67051
a78_4e_he = 4.39769
else if ( link .eq. 0.55e-9 ) then
a13_1d_he = 1.57399e-12
a13_2d_he = 3.50893e-09
a13_3d_he = 4.03868
a13_4d_he = 4.14056
a15_1d_he = 6.85547e-07
a15_2d_he = 0.000244362
a15_3d_he = 3.04719
a15_4d_he = 2.21621
a17_1d_he = -3.59973e-09
a17_2d_he = 8.81978e-07
a17_3d_he = 2.41012
a17_4d_he = 2.05063
a24_1d_he = 3.77683e-11
a24_2d_he = 3.12744e-08
a24_3d_he = 4.65005
a24_4d_he = 4.72889
a26_1d_he = 6.00368e-07
a26_2d_he = 0.000245588
a26_3d_he = 3.02745
a26_4d_he = 2.21944
a28_1d_he = 7.68638e-09
a28_2d_he = 7.54373e-06
a28_3d_he = 2.8655
a28_4d_he = 2.54412
a35_1d_he = -2.43665e-09
a35_2d_he = 8.98824e-07
a35_3d_he = 2.44482
a35_4d_he = 2.08731
a37_1d_he = 5.90618e-07
a37_2d_he = 0.000241336
a37_3d_he = 3.0211
a37_4d_he = 2.2027
a46_1d_he = 9.16578e-09
a46_2d_he = 7.65946e-06
a46_3d_he = 2.89934
a46_4d_he = 2.55858
a48_1d_he = 6.41298e-07
a48_2d_he = 0.000244868
a48_3d_he = 3.02887
a48_4d_he = 2.22031
a55_1d_he = 0.010562
a55_2d_he = 0.0720194
a55_3d_he = 0.678201
a55_4d_he = 0.751055
a57_1d_he = 2.4545e-12
a57_2d_he = 3.61305e-09
a57_3d_he = 4.07425
a57_4d_he = 4.20705
a66_1d_he = 0.0113674
a66_2d_he = 0.0728677
a66_3d_he = 0.711906
a66_4d_he = 0.770273
a68_1d_he = 2.21986e-11
a68_2d_he = 2.97398e-08
a68_3d_he = 4.50655
a68_4d_he = 4.60787
a77_1d_he = 0.0101749
a77_2d_he = 0.0723521
a77_3d_he = 0.674328
a77_4d_he = 0.738953
a88_1d_he = 0.0117021
a88_2d_he = 0.0727629
a88_3d_he = 0.719738
a88_4d_he = 0.785395
a14_1d_he = 7.21212e-12
a14_2d_he = 1.05695e-08
a14_3d_he = 3.96831
a14_4d_he = 4.82051
a16_1d_he = 3.32166e-08
a16_2d_he = 4.42872e-06
a16_3d_he = 4.53249
a16_4d_he = 1.68844
a18_1d_he = 2.55479e-09
a18_2d_he = 2.91375e-06
a18_3d_he = 2.50374
a18_4d_he = 2.69826
a23_1d_he = 1.33306e-11
a23_2d_he = 1.15501e-08
a23_3d_he = 4.8955
a23_4d_he = 4.14756
a25_1d_he = 2.94051e-08
a25_2d_he = 3.82774e-06
a25_3d_he = 4.39842
a25_4d_he = 1.41442
a27_1d_he = 6.30465e-10
a27_2d_he = 2.63317e-06
a27_3d_he = 2.96484
a27_4d_he = 2.08035
a36_1d_he = 6.34424e-09
a36_2d_he = 2.93056e-06
a36_3d_he = 2.55331
a36_4d_he = 2.72335
a38_1d_he = -3.21317e-08
a38_2d_he = 2.62097e-06
a38_3d_he = 3.66472
a38_4d_he = 0.74023
a45_1d_he = 3.60765e-09
a45_2d_he = 2.60853e-06
a45_3d_he = 3.00693
a45_4d_he = 2.08212
a47_1d_he = 2.96099e-08
a47_2d_he = 3.36313e-06
a47_3d_he = 4.41206
a47_4d_he = 1.31914
a56_1d_he = 0.000131706
a56_2d_he = 0.000768637
a56_3d_he = 1.75453
a56_4d_he = -0.622038
a58_1d_he = 9.97983e-12
a58_2d_he = 1.07316e-08
a58_3d_he = 3.99161
a58_4d_he = 4.874
a67_1d_he = 8.23298e-12
a67_2d_he = 1.11974e-08
a67_3d_he = 4.75997
a67_4d_he = 4.08435
a78_1d_he = -0.000136197
a78_2d_he = 0.00133936
a78_3d_he = 1.11817
a78_4d_he = -0.104085

a13_1e_he = 0.0109327
a13_2e_he = 0.0519701
a13_3e_he = 0.659029
a13_4e_he = 0.658513
a15_1e_he = 1.47982e-07
a15_2e_he = 3.67145e-05
a15_3e_he = 2.94317
a15_4e_he = 3.07214
a17_1e_he = 1.06945e-07
a17_2e_he = 2.28155e-05
a17_3e_he = 3.70037
a17_4e_he = 1.9026
a24_1e_he = 0.000136337
a24_2e_he = 0.0106178
a24_3e_he = 1.03583
a24_4e_he = 1.05693
a26_1e_he = 2.96311e-08
a26_2e_he = 2.57749e-05
a26_3e_he = 3.26181
a26_4e_he = 2.74179
a28_1e_he = 2.30387e-08
a28_2e_he = 1.89098e-05
a28_3e_he = 3.93186
a28_4e_he = 1.92073
a35_1e_he = 1.09057e-07
a35_2e_he = 2.28292e-05
a35_3e_he = 3.71505
a35_4e_he = 1.89522
a37_1e_he = 1.42111e-07
a37_2e_he = 3.65006e-05
a37_3e_he = 2.9313
a37_4e_he = 3.05999
a46_1e_he = 5.48839e-08
a46_2e_he = 1.98823e-05
a46_3e_he = 4.07912
a46_4e_he = 1.98258
a48_1e_he = 3.63558e-08
a48_2e_he = 2.62514e-05
a48_3e_he = 3.32534
a48_4e_he = 2.73036
a55_1e_he = 2.42958e-11
a55_2e_he = 1.83338e-08
a55_3e_he = 3.66039
a55_4e_he = 6.74665
a57_1e_he = 1.38706e-11
a57_2e_he = 1.12304e-08
a57_3e_he = 4.76099
a57_4e_he = 4.84575
a66_1e_he = 5.84948e-11
a66_2e_he = 6.85974e-08
a66_3e_he = 3.70117
a66_4e_he = 6.00316
a68_1e_he = 3.15193e-11
a68_2e_he = 5.10034e-08
a68_3e_he = 4.67405
a68_4e_he = 4.6257
a77_1e_he = 1.93027e-11
a77_2e_he = 1.85891e-08
a77_3e_he = 3.59972
a77_4e_he = 6.772
a88_1e_he = 3.94163e-11
a88_2e_he = 6.77891e-08
a88_3e_he = 3.67444
a88_4e_he = 5.88743
a14_1e_he = 0.000550669
a14_2e_he = 0.0250181
a14_3e_he = 0.502764
a14_4e_he = 1.03466
a16_1e_he = 1.61622e-07
a16_2e_he = 6.43726e-05
a16_3e_he = 2.86348
a16_4e_he = 2.76026
a18_1e_he = 1.52405e-07
a18_2e_he = 4.70051e-05
a18_3e_he = 3.47505
a18_4e_he = 1.98549
a23_1e_he = 4.2061e-05
a23_2e_he = 0.0252244
a23_3e_he = 0.998988
a23_4e_he = 0.467619
a25_1e_he = 3.32503e-08
a25_2e_he = 1.38297e-05
a25_3e_he = 3.29087
a25_4e_he = 2.99276
a27_1e_he = 3.15513e-08
a27_2e_he = 9.54817e-06
a27_3e_he = 4.28083
a27_4e_he = 1.8885
a36_1e_he = 1.53191e-07
a36_2e_he = 4.69363e-05
a36_3e_he = 3.48391
a36_4e_he = 1.97577
a38_1e_he = 1.59786e-07
a38_2e_he = 6.40722e-05
a38_3e_he = 2.85755
a38_4e_he = 2.75703
a45_1e_he = 2.62878e-08
a45_2e_he = 9.42823e-06
a45_3e_he = 4.28118
a45_4e_he = 1.82196
a47_1e_he = 2.87981e-08
a47_2e_he = 1.37012e-05
a47_3e_he = 3.24459
a47_4e_he = 2.9942
a56_1e_he = 2.69684e-11
a56_2e_he = 3.44717e-08
a56_3e_he = 3.60976
a56_4e_he = 6.25379
a58_1e_he = 2.1819e-11
a58_2e_he = 2.37258e-08
a58_3e_he = 4.36133
a58_4e_he = 5.09821
a67_1e_he = 2.9366e-11
a67_2e_he = 2.44913e-08
a67_3e_he = 5.00022
a67_4e_he = 4.60152
a78_1e_he = 2.93389e-11
a78_2e_he = 3.42756e-08
a78_3e_he = 3.58614
a78_4e_he = 6.26156
endif

Ham(0,0) = 0.0
Ham(1,1) = Eeh1(1) - Cb_eh1(1)
Ham(2,2) = Eeh1(2) - Cb_eh1(2)
Ham(3,3) = Ham(1,1)
Ham(4,4) = Ham(2,2)
Ham(5,5) = minEe(1,2) + minEh(1,1) + V0 - 1.0 * elec*(a55_1d_he + a55_2d_he / ((aR(1)*1d9)**a55_3d_he * (aR(2)*1d9)**a55_4d_he)) &
                                      + 2 * elec*(a55_1e_he + a55_2e_he / ((aR(1)*1d9)**a55_3e_he * (aR(2)*1d9)**a55_4e_he))
                                  
Ham(6,6) = minEe(1,2) + minEh(1,2) + V0 - 1.0 * elec*(a66_1d_he + a66_2d_he / ((aR(1)*1d9)**a66_3d_he * (aR(2)*1d9)**a66_4d_he)) &
                                      + 2 * elec*(a66_1e_he + a66_2e_he / ((aR(1)*1d9)**a66_3e_he * (aR(2)*1d9)**a66_4e_he)) 
                                  
Ham(7,7) = minEe(1,1) + minEh(1,1) + V0 - 1.0 * elec*(a77_1d_he + a77_2d_he / ((aR(1)*1d9)**a77_3d_he * (aR(2)*1d9)**a77_4d_he)) &
                                      + 2 * elec*(a77_1e_he + a77_2e_he / ((aR(1)*1d9)**a77_3e_he * (aR(2)*1d9)**a77_4e_he))

Ham(8,8) = minEe(1,1) + minEh(2,2) + V0 - 1.0 * elec*(a88_1d_he + a88_2d_he / ((aR(1)*1d9)**a88_3d_he * (aR(2)*1d9)**a88_4d_he)) &
                                      + 2 * elec*(a88_1e_he + a88_2e_he / ((aR(1)*1d9)**a88_3e_he * (aR(2)*1d9)**a88_4e_he))

Ham(0,1) = 0.0 
Ham(0,2) = 0.0 
Ham(0,3) = 0.0 
Ham(0,4) = 0.0 
Ham(0,5) = 0.0 
Ham(0,6) = 0.0 
Ham(0,7) = 0.0 
Ham(0,8) = 0.0 
Ham(1,0) = Ham(0,1)
Ham(2,0) = Ham(0,2)
Ham(3,0) = Ham(0,3)
Ham(4,0) = Ham(0,4)
Ham(5,0) = Ham(0,5)
Ham(6,0) = Ham(0,6)
Ham(7,0) = Ham(0,7)
Ham(8,0) = Ham(0,8)

Ham(1,2) = - 1.0 * elec*(a12_1d_he + a12_2d_he / exp((aR(1)*1d9)*a12_3d_he)) &
             + 2 * elec*(a12_1e_he + a12_2e_he / exp((aR(1)*1d9)*a12_3e_he))

Ham(1,3) = - 1.0 * elec*(a13_1d_he + a13_2d_he / ((aR(1)*1d9)**a13_3d_he * (aR(2)*1d9)**a13_4d_he)) &
             + 2 * elec*(a13_1e_he + a13_2e_he / ((aR(1)*1d9)**a13_3e_he * (aR(2)*1d9)**a13_4e_he))

Ham(1,4) = - 1.0 * elec*(a14_1d_he + a14_2d_he / ((aR(1)*1d9)**a14_3d_he * (aR(2)*1d9)**a14_4d_he)) &
             + 2 * elec*(a14_1e_he + a14_2e_he / ((aR(1)*1d9)**a14_3e_he * (aR(2)*1d9)**a14_4e_he))

Ham(1,5) = - 1.0 * elec*(a15_1d_he + a15_2d_he / ((aR(1)*1d9)**a15_3d_he * (aR(2)*1d9)**a15_4d_he)) &
             + 2 * elec*(a15_1e_he + a15_2e_he / ((aR(1)*1d9)**a15_3e_he * (aR(2)*1d9)**a15_4e_he))

Ham(1,6) = - 1.0 *( -1.0 * elec*(a16_1d_he + a16_2d_he / ((aR(1)*1d9)**a16_3d_he * (aR(2)*1d9)**a16_4d_he)) &
                     + 2 * elec*(a16_1e_he + a16_2e_he / ((aR(1)*1d9)**a16_3e_he * (aR(2)*1d9)**a16_4e_he)))

Ham(1,7) = - 1.0 * elec*(a17_1d_he + a17_2d_he / ((aR(1)*1d9)**a17_3d_he * (aR(2)*1d9)**a17_4d_he)) &
             + 2 * elec*(a17_1e_he + a17_2e_he / ((aR(1)*1d9)**a17_3e_he * (aR(2)*1d9)**a17_4e_he))

Ham(1,8) =  - 1.0 *( -1.0 * elec*(a18_1d_he + a18_2d_he / ((aR(1)*1d9)**a18_3d_he * (aR(2)*1d9)**a18_4d_he)) &
                      + 2 * elec*(a18_1e_he + a18_2e_he / ((aR(1)*1d9)**a18_3e_he * (aR(2)*1d9)**a18_4e_he)))

Ham(2,3) = - 1.0 * elec*(a23_1d_he + a23_2d_he / ((aR(1)*1d9)**a23_3d_he * (aR(2)*1d9)**a23_4d_he)) &
             + 2 * elec*(a23_1e_he + a23_2e_he / ((aR(1)*1d9)**a23_3e_he * (aR(2)*1d9)**a23_4e_he))

Ham(2,4) = - 1.0 * elec*(a24_1d_he + a24_2d_he / ((aR(1)*1d9)**a24_3d_he * (aR(2)*1d9)**a24_4d_he)) &
             + 2 * elec*(a24_1e_he + a24_2e_he / ((aR(1)*1d9)**a24_3e_he * (aR(2)*1d9)**a24_4e_he))

Ham(2,5) =  - 1.0 *( -1.0 * elec*(a25_1d_he + a25_2d_he / ((aR(1)*1d9)**a25_3d_he * (aR(2)*1d9)**a25_4d_he)) &
                      + 2 * elec*(a25_1e_he + a25_2e_he / ((aR(1)*1d9)**a25_3e_he * (aR(2)*1d9)**a25_4e_he)))

Ham(2,6) =  - 1.0 * elec*(a26_1d_he + a26_2d_he / ((aR(1)*1d9)**a26_3d_he * (aR(2)*1d9)**a26_4d_he)) &
              + 2 * elec*(a26_1e_he + a26_2e_he / ((aR(1)*1d9)**a26_3e_he * (aR(2)*1d9)**a26_4e_he))

Ham(2,7) =  - 1.0 *( -1.0 * elec*(a27_1d_he + a27_2d_he / ((aR(1)*1d9)**a27_3d_he * (aR(2)*1d9)**a27_4d_he)) &
                      + 2 * elec*(a27_1e_he + a27_2e_he / ((aR(1)*1d9)**a27_3e_he * (aR(2)*1d9)**a27_4e_he)))

Ham(2,8) =  - 1.0 * elec*(a28_1d_he + a28_2d_he / ((aR(1)*1d9)**a28_3d_he * (aR(2)*1d9)**a28_4d_he)) &
              + 2 * elec*(a28_1e_he + a28_2e_he / ((aR(1)*1d9)**a28_3e_he * (aR(2)*1d9)**a28_4e_he))

Ham(3,4) = - 1.0 * elec*(a12_1d_he + a12_2d_he / exp((aR(2)*1d9)*a12_3d_he)) &
             + 2 * elec*(a12_1e_he + a12_2e_he / exp((aR(2)*1d9)*a12_3e_he))

Ham(3,5) =  - 1.0 * elec*(a35_1d_he + a35_2d_he / ((aR(1)*1d9)**a35_3d_he * (aR(2)*1d9)**a35_4d_he)) &
              + 2 * elec*(a35_1e_he + a35_2e_he / ((aR(1)*1d9)**a35_3e_he * (aR(2)*1d9)**a35_4e_he))

Ham(3,6) =  - 1.0 * elec*(a36_1d_he + a36_2d_he / ((aR(1)*1d9)**a36_3d_he * (aR(2)*1d9)**a36_4d_he)) &
              + 2 * elec*(a36_1e_he + a36_2e_he / ((aR(1)*1d9)**a36_3e_he * (aR(2)*1d9)**a36_4e_he))

Ham(3,7) =  - 1.0 * elec*(a37_1d_he + a37_2d_he / ((aR(1)*1d9)**a37_3d_he * (aR(2)*1d9)**a37_4d_he)) &
              + 2 * elec*(a37_1e_he + a37_2e_he / ((aR(1)*1d9)**a37_3e_he * (aR(2)*1d9)**a37_4e_he))

Ham(3,8) =  - 1.0 * elec*(a38_1d_he + a38_2d_he / ((aR(1)*1d9)**a38_3d_he * (aR(2)*1d9)**a38_4d_he)) &
              + 2 * elec*(a38_1e_he + a38_2e_he / ((aR(1)*1d9)**a38_3e_he * (aR(2)*1d9)**a38_4e_he))

Ham(4,5) =  - 1.0 * elec*(a45_1d_he + a45_2d_he / ((aR(1)*1d9)**a45_3d_he * (aR(2)*1d9)**a45_4d_he)) &
              + 2 * elec*(a45_1e_he + a45_2e_he / ((aR(1)*1d9)**a45_3e_he * (aR(2)*1d9)**a45_4e_he))

Ham(4,6) =  - 1.0 * elec*(a46_1d_he + a46_2d_he / ((aR(1)*1d9)**a46_3d_he * (aR(2)*1d9)**a46_4d_he)) &
              + 2 * elec*(a46_1e_he + a46_2e_he / ((aR(1)*1d9)**a46_3e_he * (aR(2)*1d9)**a46_4e_he))

Ham(4,7) =  - 1.0 * elec*(a47_1d_he + a47_2d_he / ((aR(1)*1d9)**a47_3d_he * (aR(2)*1d9)**a47_4d_he)) &
              + 2 * elec*(a47_1e_he + a47_2e_he / ((aR(1)*1d9)**a47_3e_he * (aR(2)*1d9)**a47_4e_he))

Ham(4,8) =  - 1.0 * elec*(a48_1d_he + a48_2d_he / ((aR(1)*1d9)**a48_3d_he * (aR(2)*1d9)**a48_4d_he)) &
              + 2 * elec*(a48_1e_he + a48_2e_he / ((aR(1)*1d9)**a48_3e_he * (aR(2)*1d9)**a48_4e_he))

Ham(5,6) = - 1.0 * elec*(a56_1d_he + a56_2d_he / ((aR(1)*1d9)**a56_3d_he * (aR(2)*1d9)**a56_4d_he)) &
             + 2 * elec*(a56_1e_he + a56_2e_he / ((aR(1)*1d9)**a56_3e_he * (aR(2)*1d9)**a56_4e_he))

Ham(5,7) = - 1.0 * elec*(a57_1d_he + a57_2d_he / ((aR(1)*1d9)**a57_3d_he * (aR(2)*1d9)**a57_4d_he)) & 
             + 2 * elec*(a57_1e_he + a57_2e_he / ((aR(1)*1d9)**a57_3e_he * (aR(2)*1d9)**a57_4e_he))   

Ham(5,8) = - 1.0 * elec*(a58_1d_he + a58_2d_he / ((aR(1)*1d9)**a58_3d_he * (aR(2)*1d9)**a58_4d_he)) & 
             + 2 * elec*(a58_1e_he + a58_2e_he / ((aR(1)*1d9)**a58_3e_he * (aR(2)*1d9)**a58_4e_he))   

Ham(6,7) = - 1.0 * elec*(a67_1d_he + a67_2d_he / ((aR(1)*1d9)**a67_3d_he * (aR(2)*1d9)**a67_4d_he)) & 
             + 2 * elec*(a67_1e_he + a67_2e_he / ((aR(1)*1d9)**a67_3e_he * (aR(2)*1d9)**a67_4e_he))   

Ham(6,8) = - 1.0 * elec*(a68_1d_he + a68_2d_he / ((aR(1)*1d9)**a68_3d_he * (aR(2)*1d9)**a68_4d_he)) & 
             + 2 * elec*(a68_1e_he + a68_2e_he / ((aR(1)*1d9)**a68_3e_he * (aR(2)*1d9)**a68_4e_he))   

Ham(7,8) = - 1.0 * elec*(a78_1d_he + a78_2d_he / ((aR(1)*1d9)**a78_3d_he * (aR(2)*1d9)**a78_4d_he)) &
             + 2 * elec*(a78_1e_he + a78_2e_he / ((aR(1)*1d9)**a78_3e_he * (aR(2)*1d9)**a78_4e_he))

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

Ham = Ham/Energ_au
TransHam(0,1) = TransDip_Ana_h1e(n)
TransHam(0,2) = TransDip_Ana_h2e(n)
TransHam(0,3) = TransDip_Ana_h1e(n+nsys)
TransHam(0,4) = TransDip_Ana_h2e(n+nsys)
TransHam(0,5) = TransDip_Fit_h1e_he(aR(n),aR(n+nsys))
TransHam(0,6) = TransDip_Fit_h2e_he(aR(n),aR(n+nsys))
TransHam(0,7) = TransDip_Fit_h1e_he(aR(n),aR(n+nsys))
TransHam(0,8) = TransDip_Fit_h2e_he(aR(n),aR(n+nsys))
TransHam = TransHam/Dip_au
do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

end subroutine make_Ham_he

subroutine make_Ham_fineSt

Ham(1,1)   = Eeh1(1) - Cb_eh1(1) - Dso/3 - 3*Kas
Ham(2,2)   = Eeh1(1) - Cb_eh1(1) - Dso/3 - 3*Kcs
Ham(3,3)   = Eeh1(1) - Cb_eh1(1) - Kas
Ham(4,4)   = Eeh1(1) - Cb_eh1(1) - 3*Kas
Ham(5,5)   = Eeh1(1) - Cb_eh1(1) - 3*Kbs + Dxf
Ham(6,6)   = Eeh1(1) - Cb_eh1(1) - 3*Kbs + Dxf
Ham(7,7)   = Eeh1(1) - Cb_eh1(1) - Kcs
Ham(8,8)   = Eeh1(1) - Cb_eh1(1) - 3*Kcs
Ham(9,9)   = Eeh1(1) - Cb_eh1(1) + Dso/3 - 3*Kas
Ham(10,10) = Eeh1(1) - Cb_eh1(1) - Kbs + Dxf
Ham(11,11) = Eeh1(1) - Cb_eh1(1) - 3*Kbs + Dxf
Ham(12,12) = Eeh1(1) - Cb_eh1(1) + Dso/3 - 3*Kcs

Ham(3,4) = -1.0*Dso/3
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

Ham(4,5) = Dso/3
Ham(5,4) = Ham(4,5)
Ham(7,8) = Ham(4,5)
Ham(8,7) = Ham(4,5)
Ham(11,12) = Ham(4,5)
Ham(12,11) = Ham(4,5)

Ham(13,13)   = Eeh2(1) - Cb_eh2(1) - Dso/3 - 3*Kas
Ham(14,14)   = Eeh2(1) - Cb_eh2(1) - Dso/3 - 3*Kcs
Ham(15,15)   = Eeh2(1) - Cb_eh2(1) - Kas
Ham(16,16)   = Eeh2(1) - Cb_eh2(1) - 3*Kas
Ham(17,17)   = Eeh2(1) - Cb_eh2(1) - 3*Kbs + Dxf
Ham(18,18)   = Eeh2(1) - Cb_eh2(1) - 3*Kbs + Dxf
Ham(19,19)   = Eeh2(1) - Cb_eh2(1) - Kcs
Ham(20,20)   = Eeh2(1) - Cb_eh2(1) - 3*Kcs
Ham(21,21)   = Eeh2(1) - Cb_eh2(1) + Dso/3 - 3*Kas
Ham(22,22)   = Eeh2(1) - Cb_eh2(1) - Kbs + Dxf
Ham(23,23)   = Eeh2(1) - Cb_eh2(1) - 3*Kbs + Dxf
Ham(24,24)   = Eeh2(1) - Cb_eh2(1) + Dso/3 - 3*Kcs

Ham(15,16) = -1.0*Dso/3
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

Ham(16,17) = Dso/3
Ham(17,16) = Ham(16,17)
Ham(19,20) = Ham(16,17)
Ham(20,19) = Ham(16,17)
Ham(23,24) = Ham(16,17)
Ham(24,23) = Ham(16,17)

TransHam(0,3)  = 10.07* 3.33564e-30   
TransHam(0,4)  = 10.07* 3.33564e-30 
TransHam(0,6)  = 17.48* 3.33564e-30 
TransHam(0,7)  = 17.48* 3.33564e-30 
TransHam(0,8)  = 19.50* 3.33564e-30 
TransHam(0,10) = 14.63* 3.33564e-30 
TransHam(0,11) = 14.63* 3.33564e-30 
TransHam(0,12) = 15.53* 3.33564e-30 
TransHam(0,15) = 6.76 * 3.33564e-30 
TransHam(0,16) = 6.76 * 3.33564e-30 
TransHam(0,18) = 11.73* 3.33564e-30 
TransHam(0,19) = 11.73* 3.33564e-30 
TransHam(0,20) = 13.08* 3.33564e-30 
TransHam(0,22) = 9.82 * 3.33564e-30 
TransHam(0,23) = 9.82 * 3.33564e-30 
TransHam(0,24) = 10.42* 3.33564e-30 

do i=0,nstates-1
TransHam(i,0) = TransHam(0,i)
enddo

end subroutine make_Ham_fineSt

end module Make_Ham
