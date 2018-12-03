module Make_Ham

use Constants
use Variables

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

real(dp), parameter :: a12_1d_he = 0.00876695
real(dp), parameter :: a12_2d_he = 0.0346284
real(dp), parameter :: a12_3d_he = 8.94293d+08
real(dp), parameter :: a12_1e_he = 0.000588561
real(dp), parameter :: a12_2e_he = 0.231987
real(dp), parameter :: a12_3e_he = 1.64474d+09
real(dp), parameter :: a13_1d_he = -5.7308d-10
real(dp), parameter :: a13_2d_he = 1.16996d-06
real(dp), parameter :: a13_3d_he = 3.64279
real(dp), parameter :: a13_4d_he = 3.80378
real(dp), parameter :: a15_1d_he = 3.32236d-06
real(dp), parameter :: a15_2d_he = 0.00204328
real(dp), parameter :: a15_3d_he = 2.77545
real(dp), parameter :: a15_4d_he = 1.85776
real(dp), parameter :: a17_1d_he = -3.56487d-08
real(dp), parameter :: a17_2d_he = 4.00073d-05
real(dp), parameter :: a17_3d_he = 2.68943
real(dp), parameter :: a17_4d_he = 2.0951
real(dp), parameter :: a24_1d_he = 2.14157d-09
real(dp), parameter :: a24_2d_he = 6.46666d-06
real(dp), parameter :: a24_3d_he = 3.95207
real(dp), parameter :: a24_4d_he = 4.00739
real(dp), parameter :: a26_1d_he = 4.13952d-06
real(dp), parameter :: a26_2d_he = 0.00206339
real(dp), parameter :: a26_3d_he = 2.80285
real(dp), parameter :: a26_4d_he = 1.85298
real(dp), parameter :: a28_1d_he = 5.16485d-08
real(dp), parameter :: a28_2d_he = 0.000216119
real(dp), parameter :: a28_3d_he = 2.85638
real(dp), parameter :: a28_4d_he = 2.27227
real(dp), parameter :: a35_1d_he = -4.63653d-08
real(dp), parameter :: a35_2d_he = 3.91962d-05
real(dp), parameter :: a35_3d_he = 2.68723
real(dp), parameter :: a35_4d_he = 2.0612
real(dp), parameter :: a37_1d_he = 3.69519d-06
real(dp), parameter :: a37_2d_he = 0.00204119
real(dp), parameter :: a37_3d_he = 2.77624
real(dp), parameter :: a37_4d_he = 1.86059
real(dp), parameter :: a46_1d_he = -3.1446d-07
real(dp), parameter :: a46_2d_he = 0.000210886
real(dp), parameter :: a46_3d_he = 2.80464
real(dp), parameter :: a46_4d_he = 2.20824
real(dp), parameter :: a48_1d_he = 6.30777d-06
real(dp), parameter :: a48_2d_he = 0.00211127
real(dp), parameter :: a48_3d_he = 2.87692
real(dp), parameter :: a48_4d_he = 1.85623
real(dp), parameter :: a55_1d_he = 0.0123067
real(dp), parameter :: a55_2d_he = 0.0834401
real(dp), parameter :: a55_3d_he = 0.751985
real(dp), parameter :: a55_4d_he = 0.822935
real(dp), parameter :: a57_1d_he = -2.18594d-10
real(dp), parameter :: a57_2d_he = 1.20859d-06
real(dp), parameter :: a57_3d_he = 3.72691
real(dp), parameter :: a57_4d_he = 3.81947
real(dp), parameter :: a66_1d_he = 0.0120069
real(dp), parameter :: a66_2d_he = 0.084655
real(dp), parameter :: a66_3d_he = 0.75655
real(dp), parameter :: a66_4d_he = 0.813347
real(dp), parameter :: a68_1d_he = 2.26161d-09
real(dp), parameter :: a68_2d_he = 6.64211d-06
real(dp), parameter :: a68_3d_he = 3.9102
real(dp), parameter :: a68_4d_he = 4.10851
real(dp), parameter :: a77_1d_he = 0.0131285
real(dp), parameter :: a77_2d_he = 0.0829806
real(dp), parameter :: a77_3d_he = 0.772324
real(dp), parameter :: a77_4d_he = 0.839392
real(dp), parameter :: a88_1d_he = 0.0110918
real(dp), parameter :: a88_2d_he = 0.0853286
real(dp), parameter :: a88_3d_he = 0.75124
real(dp), parameter :: a88_4d_he = 0.781124
real(dp), parameter :: a14_1d_he = -4.51135d-10
real(dp), parameter :: a14_2d_he = 2.84653d-06
real(dp), parameter :: a14_3d_he = 3.65651
real(dp), parameter :: a14_4d_he = 4.12674
real(dp), parameter :: a16_1d_he = 4.50443d-07
real(dp), parameter :: a16_2d_he = 2.55061d-05
real(dp), parameter :: a16_3d_he = 4.46738
real(dp), parameter :: a16_4d_he = 0.606741
real(dp), parameter :: a18_1d_he = 3.16132d-08
real(dp), parameter :: a18_2d_he = 9.86635d-05
real(dp), parameter :: a18_3d_he = 2.73099
real(dp), parameter :: a18_4d_he = 2.33682
real(dp), parameter :: a23_1d_he = 1.43653d-09
real(dp), parameter :: a23_2d_he = 2.99651d-06
real(dp), parameter :: a23_3d_he = 4.03847
real(dp), parameter :: a23_4d_he = 3.89306
real(dp), parameter :: a25_1d_he = 3.80318d-07
real(dp), parameter :: a25_2d_he = 3.90141d-05
real(dp), parameter :: a25_3d_he = 4.43569
real(dp), parameter :: a25_4d_he = 1.3408
real(dp), parameter :: a27_1d_he = -1.38977d-07
real(dp), parameter :: a27_2d_he = 9.05307d-05
real(dp), parameter :: a27_3d_he = 2.81297
real(dp), parameter :: a27_4d_he = 2.05185
real(dp), parameter :: a36_1d_he = -7.52562d-08
real(dp), parameter :: a36_2d_he = 9.46348d-05
real(dp), parameter :: a36_3d_he = 2.64458
real(dp), parameter :: a36_4d_he = 2.2994
real(dp), parameter :: a38_1d_he = -1.87329d-07
real(dp), parameter :: a38_2d_he = 3.50221d-05
real(dp), parameter :: a38_3d_he = 3.41031
real(dp), parameter :: a38_4d_he = 1.25104
real(dp), parameter :: a45_1d_he = -1.62686d-07
real(dp), parameter :: a45_2d_he = 9.0782d-05
real(dp), parameter :: a45_3d_he = 2.84149
real(dp), parameter :: a45_4d_he = 2.02033
real(dp), parameter :: a47_1d_he = 1.54893d-07
real(dp), parameter :: a47_2d_he = 2.28436d-05
real(dp), parameter :: a47_3d_he = 4.12838
real(dp), parameter :: a47_4d_he = 0.493992
real(dp), parameter :: a56_1d_he = 0.000250184
real(dp), parameter :: a56_2d_he = 0.00138359
real(dp), parameter :: a56_3d_he = 2.70362
real(dp), parameter :: a56_4d_he = -0.0663319
real(dp), parameter :: a58_1d_he = 1.50601d-09
real(dp), parameter :: a58_2d_he = 2.85302d-06
real(dp), parameter :: a58_3d_he = 3.72181
real(dp), parameter :: a58_4d_he = 4.19241
real(dp), parameter :: a67_1d_he = 5.44859d-10
real(dp), parameter :: a67_2d_he = 2.90479d-06
real(dp), parameter :: a67_3d_he = 3.98366
real(dp), parameter :: a67_4d_he = 3.83348
real(dp), parameter :: a78_1d_he = 0.000183386
real(dp), parameter :: a78_2d_he = 0.00107079
real(dp), parameter :: a78_3d_he = 2.28933
real(dp), parameter :: a78_4d_he = -0.466982
real(dp), parameter :: a13_1e_he = 0.0129891
real(dp), parameter :: a13_2e_he = 0.059974
real(dp), parameter :: a13_3e_he = 0.745223
real(dp), parameter :: a13_4e_he = 0.746591
real(dp), parameter :: a15_1e_he = 7.81765d-07
real(dp), parameter :: a15_2e_he = 0.000436249
real(dp), parameter :: a15_3e_he = 2.94985
real(dp), parameter :: a15_4e_he = 2.23921
real(dp), parameter :: a17_1e_he = 8.32306d-07
real(dp), parameter :: a17_2e_he = 0.000325249
real(dp), parameter :: a17_3e_he = 3.0237
real(dp), parameter :: a17_4e_he = 1.92494
real(dp), parameter :: a24_1e_he = 0.000165069
real(dp), parameter :: a24_2e_he = 0.0120975
real(dp), parameter :: a24_3e_he = 1.09414
real(dp), parameter :: a24_4e_he = 1.0831
real(dp), parameter :: a26_1e_he = -3.83814d-08
real(dp), parameter :: a26_2e_he = 0.000298365
real(dp), parameter :: a26_3e_he = 3.27584
real(dp), parameter :: a26_4e_he = 2.06783
real(dp), parameter :: a28_1e_he = 8.50008d-08
real(dp), parameter :: a28_2e_he = 0.000253351
real(dp), parameter :: a28_3e_he = 3.3815
real(dp), parameter :: a28_4e_he = 1.89094
real(dp), parameter :: a35_1e_he = 8.11209d-07
real(dp), parameter :: a35_2e_he = 0.000323139
real(dp), parameter :: a35_3e_he = 3.01545
real(dp), parameter :: a35_4e_he = 1.91859
real(dp), parameter :: a37_1e_he = 7.11155d-07
real(dp), parameter :: a37_2e_he = 0.000434065
real(dp), parameter :: a37_3e_he = 2.93827
real(dp), parameter :: a37_4e_he = 2.23353
real(dp), parameter :: a46_1e_he = 3.42974d-07
real(dp), parameter :: a46_2e_he = 0.00025947
real(dp), parameter :: a46_3e_he = 3.47192
real(dp), parameter :: a46_4e_he = 1.9116
real(dp), parameter :: a48_1e_he = -4.01333d-07
real(dp), parameter :: a48_2e_he = 0.000290705
real(dp), parameter :: a48_3e_he = 3.16886
real(dp), parameter :: a48_4e_he = 2.05449
real(dp), parameter :: a55_1e_he = 2.18641d-09
real(dp), parameter :: a55_2e_he = 2.66753d-06
real(dp), parameter :: a55_3e_he = 3.7282
real(dp), parameter :: a55_4e_he = 4.73166
real(dp), parameter :: a57_1e_he = 5.7703d-10
real(dp), parameter :: a57_2e_he = 2.01368d-06
real(dp), parameter :: a57_3e_he = 3.99743
real(dp), parameter :: a57_4e_he = 4.00244
real(dp), parameter :: a66_1e_he = -1.90295d-09
real(dp), parameter :: a66_2e_he = 8.64309d-06
real(dp), parameter :: a66_3e_he = 3.58118
real(dp), parameter :: a66_4e_he = 4.22789
real(dp), parameter :: a68_1e_he = 2.24073d-09
real(dp), parameter :: a68_2e_he = 8.0202d-06
real(dp), parameter :: a68_3e_he = 3.92005
real(dp), parameter :: a68_4e_he = 4.06342
real(dp), parameter :: a77_1e_he = 1.44698d-09
real(dp), parameter :: a77_2e_he = 2.5936d-06
real(dp), parameter :: a77_3e_he = 3.61919
real(dp), parameter :: a77_4e_he = 4.71832
real(dp), parameter :: a88_1e_he = -3.75173d-10
real(dp), parameter :: a88_2e_he = 8.84303d-06
real(dp), parameter :: a88_3e_he = 3.5768
real(dp), parameter :: a88_4e_he = 4.2987
real(dp), parameter :: a14_1e_he = 0.00100523
real(dp), parameter :: a14_2e_he = 0.0286904
real(dp), parameter :: a14_3e_he = 0.571985
real(dp), parameter :: a14_4e_he = 1.11753
real(dp), parameter :: a16_1e_he = 5.20996d-08
real(dp), parameter :: a16_2e_he = 0.000733493
real(dp), parameter :: a16_3e_he = 2.81241
real(dp), parameter :: a16_4e_he = 2.05309
real(dp), parameter :: a18_1e_he = 9.14038d-07
real(dp), parameter :: a18_2e_he = 0.000618908
real(dp), parameter :: a18_3e_he = 2.90339
real(dp), parameter :: a18_4e_he = 1.91276
real(dp), parameter :: a23_1e_he = 0.00072624
real(dp), parameter :: a23_2e_he = 0.0288513
real(dp), parameter :: a23_3e_he = 1.08402
real(dp), parameter :: a23_4e_he = 0.561365
real(dp), parameter :: a25_1e_he = -5.80522d-08
real(dp), parameter :: a25_2e_he = 0.000160715
real(dp), parameter :: a25_3e_he = 3.26161
real(dp), parameter :: a25_4e_he = 2.1577
real(dp), parameter :: a27_1e_he = 1.47805d-07
real(dp), parameter :: a27_2e_he = 0.000132624
real(dp), parameter :: a27_3e_he = 3.4939
real(dp), parameter :: a27_4e_he = 1.90809
real(dp), parameter :: a36_1e_he = 1.11843d-06
real(dp), parameter :: a36_2e_he = 0.00062124
real(dp), parameter :: a36_3e_he = 2.90829
real(dp), parameter :: a36_4e_he = 1.92675
real(dp), parameter :: a38_1e_he = 1.85001d-08
real(dp), parameter :: a38_2e_he = 0.000731357
real(dp), parameter :: a38_3e_he = 2.80983
real(dp), parameter :: a38_4e_he = 2.04924
real(dp), parameter :: a45_1e_he = 1.43075d-07
real(dp), parameter :: a45_2e_he = 0.000133518
real(dp), parameter :: a45_3e_he = 3.48729
real(dp), parameter :: a45_4e_he = 1.92688
real(dp), parameter :: a47_1e_he = 1.78003d-07
real(dp), parameter :: a47_2e_he = 0.000164668
real(dp), parameter :: a47_3e_he = 3.28665
real(dp), parameter :: a47_4e_he = 2.24905
real(dp), parameter :: a56_1e_he = 1.19765d-09
real(dp), parameter :: a56_2e_he = 4.7171d-06
real(dp), parameter :: a56_3e_he = 3.62877
real(dp), parameter :: a56_4e_he = 4.45655
real(dp), parameter :: a58_1e_he = 8.52181d-10
real(dp), parameter :: a58_2e_he = 4.03096d-06
real(dp), parameter :: a58_3e_he = 3.89968
real(dp), parameter :: a58_4e_he = 4.1013
real(dp), parameter :: a67_1e_he = -4.70846d-10
real(dp), parameter :: a67_2e_he = 4.045d-06
real(dp), parameter :: a67_3e_he = 4.05712
real(dp), parameter :: a67_4e_he = 3.87227
real(dp), parameter :: a78_1e_he = 3.51426d-10
real(dp), parameter :: a78_2e_he = 4.75278d-06
real(dp), parameter :: a78_3e_he = 3.67051
real(dp), parameter :: a78_4e_he = 4.39769

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

end subroutine make_Ham_ho

subroutine make_Ham_he

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

end subroutine make_Ham_he

end module Make_Ham
