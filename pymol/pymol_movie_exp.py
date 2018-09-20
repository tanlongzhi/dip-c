set matrix_mode, 1
set movie_panel, 0
set scene_buttons, 1
set cache_frames, 1
viewport 1080, 1080

set ray_shadows,0
as sticks, all
set_bond stick_radius, 0.5, all
spectrum b, rainbow, all, 1, 23

extract 1_mat, chain "1(mat)"
extract 1_pat, chain "1(pat)"
extract 10_mat, chain "10(mat)"
extract 10_pat, chain "10(pat)"
extract 11_mat, chain "11(mat)"
extract 11_pat, chain "11(pat)"
extract 12_mat, chain "12(mat)"
extract 12_pat, chain "12(pat)"
extract 13_mat, chain "13(mat)"
extract 13_pat, chain "13(pat)"
extract 14_mat, chain "14(mat)"
extract 14_pat, chain "14(pat)"
extract 15_mat, chain "15(mat)"
extract 15_pat, chain "15(pat)"
extract 16_mat, chain "16(mat)"
extract 16_pat, chain "16(pat)"
extract 17_mat, chain "17(mat)"
extract 17_pat, chain "17(pat)"
extract 18_mat, chain "18(mat)"
extract 18_pat, chain "18(pat)"
extract 19_mat, chain "19(mat)"
extract 19_pat, chain "19(pat)"
extract 2_mat, chain "2(mat)"
extract 2_pat, chain "2(pat)"
extract 20_mat, chain "20(mat)"
extract 20_pat, chain "20(pat)"
extract 21_mat, chain "21(mat)"
extract 21_pat, chain "21(pat)"
extract 22_mat, chain "22(mat)"
extract 22_pat, chain "22(pat)"
extract 3_mat, chain "3(mat)"
extract 3_pat, chain "3(pat)"
extract 4_mat, chain "4(mat)"
extract 4_pat, chain "4(pat)"
extract 5_mat, chain "5(mat)"
extract 5_pat, chain "5(pat)"
extract 6_mat, chain "6(mat)"
extract 6_pat, chain "6(pat)"
extract 7_mat, chain "7(mat)"
extract 7_pat, chain "7(pat)"
extract 8_mat, chain "8(mat)"
extract 8_pat, chain "8(pat)"
extract 9_mat, chain "9(mat)"
extract 9_pat, chain "9(pat)"
extract X_mat, chain "X(mat)"
extract X_pat, chain "X(pat)"

set_view (\
     0.750269890,    0.324408263,    0.576068759,\
     0.147000015,    0.767668426,   -0.623760700,\
    -0.644581974,    0.552670181,    0.528270006,\
     0.000002384,    0.000001788, -1049.974853516,\
    -6.021775723,   -6.870912552,   -2.041615725,\
   812.712585449, 1287.237182617,  -20.000000000 )

mset 1 x240

   
frame 1
mview store, object=1_mat
mview store, object=1_pat
mview store, object=10_mat
mview store, object=10_pat
mview store, object=11_mat
mview store, object=11_pat
mview store, object=12_mat
mview store, object=12_pat
mview store, object=13_mat
mview store, object=13_pat
mview store, object=14_mat
mview store, object=14_pat
mview store, object=15_mat
mview store, object=15_pat
mview store, object=16_mat
mview store, object=16_pat
mview store, object=17_mat
mview store, object=17_pat
mview store, object=18_mat
mview store, object=18_pat
mview store, object=19_mat
mview store, object=19_pat
mview store, object=2_mat
mview store, object=2_pat
mview store, object=20_mat
mview store, object=20_pat
mview store, object=21_mat
mview store, object=21_pat
mview store, object=22_mat
mview store, object=22_pat
mview store, object=3_mat
mview store, object=3_pat
mview store, object=4_mat
mview store, object=4_pat
mview store, object=5_mat
mview store, object=5_pat
mview store, object=6_mat
mview store, object=6_pat
mview store, object=7_mat
mview store, object=7_pat
mview store, object=8_mat
mview store, object=8_pat
mview store, object=9_mat
mview store, object=9_pat
mview store, object=X_mat
mview store, object=X_pat
mview store

frame 61
mview store, object=1_mat
mview store, object=1_pat
mview store, object=10_mat
mview store, object=10_pat
mview store, object=11_mat
mview store, object=11_pat
mview store, object=12_mat
mview store, object=12_pat
mview store, object=13_mat
mview store, object=13_pat
mview store, object=14_mat
mview store, object=14_pat
mview store, object=15_mat
mview store, object=15_pat
mview store, object=16_mat
mview store, object=16_pat
mview store, object=17_mat
mview store, object=17_pat
mview store, object=18_mat
mview store, object=18_pat
mview store, object=19_mat
mview store, object=19_pat
mview store, object=2_mat
mview store, object=2_pat
mview store, object=20_mat
mview store, object=20_pat
mview store, object=21_mat
mview store, object=21_pat
mview store, object=22_mat
mview store, object=22_pat
mview store, object=3_mat
mview store, object=3_pat
mview store, object=4_mat
mview store, object=4_pat
mview store, object=5_mat
mview store, object=5_pat
mview store, object=6_mat
mview store, object=6_pat
mview store, object=7_mat
mview store, object=7_pat
mview store, object=8_mat
mview store, object=8_pat
mview store, object=9_mat
mview store, object=9_pat
mview store, object=X_mat
mview store, object=X_pat
turn y,90
mview store

frame 121
translate [-50.573105313,47.7330591577,74.1278178996], object=1_mat, camera=0
translate [-13.1357309253,-94.8084330769,8.4931842669], object=1_pat, camera=0
translate [-93.9998927298,-52.2533840041,32.845810202], object=10_mat, camera=0
translate [-113.242247665,-3.56242013626,72.2642629726], object=10_pat, camera=0
translate [44.5965755417,-7.73747765276,-118.886541216], object=11_mat, camera=0
translate [-96.6557544499,-4.58341825817,-52.2020253907], object=11_pat, camera=0
translate [-102.955348824,35.9164344531,84.6534045198], object=12_mat, camera=0
translate [-86.0586318589,-55.3600178366,-66.4917392689], object=12_pat, camera=0
translate [55.6850210667,7.17220601691,90.6350121511], object=13_mat, camera=0
translate [117.926172832,-13.9399517924,11.4593568631], object=13_pat, camera=0
translate [-35.0864626096,-36.7379814144,129.836291657], object=14_mat, camera=0
translate [-90.1410724534,40.8081236724,-29.9233384234], object=14_pat, camera=0
translate [48.308362928,-90.2308074972,-33.5707222744], object=15_mat, camera=0
translate [1.07732898254,-77.9521886966,-31.7786429694], object=15_pat, camera=0
translate [25.6957738104,-8.05464153745,80.9118773749], object=16_mat, camera=0
translate [-40.2681584777,-79.7462161508,22.7855555632], object=16_pat, camera=0
translate [-46.0352393275,-15.7579122789,8.21833279353], object=17_mat, camera=0
translate [35.0864257972,6.30814060123,-101.527383148], object=17_pat, camera=0
translate [-55.4934022548,118.607642906,59.8029037006], object=18_mat, camera=0
translate [13.0141117432,-39.2105476411,-53.3637078487], object=18_pat, camera=0
translate [-57.3322422559,24.2791605626,-77.9856617662], object=19_mat, camera=0
translate [-17.7641256774,61.3981098935,-9.02152763102], object=19_pat, camera=0
translate [26.0589048771,110.427908577,-7.35546026234], object=2_mat, camera=0
translate [63.9435239283,39.5529058449,0.910299436047], object=2_pat, camera=0
translate [37.1388001155,-56.702553347,100.405922904], object=20_mat, camera=0
translate [-14.905450053,-20.2673895134,79.1241648251], object=20_pat, camera=0
translate [9.88389916627,-19.5165763475,138.911799454], object=21_mat, camera=0
translate [-3.80690283991,-43.8533420936,-33.9092683513], object=21_pat, camera=0
translate [-99.5684227941,-4.63323814867,24.7105534268], object=22_mat, camera=0
translate [1.44705107983,-12.2571639626,-27.571309466], object=22_pat, camera=0
translate [107.907500507,-41.0521568758,-26.6700557391], object=3_mat, camera=0
translate [73.7661493548,-82.1811618648,-73.9075983871], object=3_pat, camera=0
translate [-8.58885996012,54.8838225305,116.89380896], object=4_mat, camera=0
translate [-28.4138353718,-115.462028775,-58.972958804], object=4_pat, camera=0
translate [-38.4464491394,-26.1697174238,79.4459578982], object=5_mat, camera=0
translate [-19.4534943505,-63.0789501673,-125.641992792], object=5_pat, camera=0
translate [-25.9582754118,75.6152241819,-85.8921827626], object=6_mat, camera=0
translate [58.7253214897,87.0750566429,-75.6465115294], object=6_pat, camera=0
translate [-96.4747825099,85.1103285974,13.5695788907], object=7_mat, camera=0
translate [-17.7388207673,119.268544049,22.0593360286], object=7_pat, camera=0
translate [48.9872383503,67.0899777755,84.6164676867], object=8_mat, camera=0
translate [-9.57549499521,27.208814551,-125.310392466], object=8_pat, camera=0
translate [81.6684592598,27.6542751367,-81.8300056854], object=9_mat, camera=0
translate [80.0937574521,-29.8435114742,-15.3388418139], object=9_pat, camera=0
translate [98.9599356232,-25.9348903078,71.8748956543], object=X_mat, camera=0
translate [53.653125317,-91.2443213672,65.9293798679], object=X_pat, camera=0
mview store, object=1_mat
mview store, object=1_pat
mview store, object=10_mat
mview store, object=10_pat
mview store, object=11_mat
mview store, object=11_pat
mview store, object=12_mat
mview store, object=12_pat
mview store, object=13_mat
mview store, object=13_pat
mview store, object=14_mat
mview store, object=14_pat
mview store, object=15_mat
mview store, object=15_pat
mview store, object=16_mat
mview store, object=16_pat
mview store, object=17_mat
mview store, object=17_pat
mview store, object=18_mat
mview store, object=18_pat
mview store, object=19_mat
mview store, object=19_pat
mview store, object=2_mat
mview store, object=2_pat
mview store, object=20_mat
mview store, object=20_pat
mview store, object=21_mat
mview store, object=21_pat
mview store, object=22_mat
mview store, object=22_pat
mview store, object=3_mat
mview store, object=3_pat
mview store, object=4_mat
mview store, object=4_pat
mview store, object=5_mat
mview store, object=5_pat
mview store, object=6_mat
mview store, object=6_pat
mview store, object=7_mat
mview store, object=7_pat
mview store, object=8_mat
mview store, object=8_pat
mview store, object=9_mat
mview store, object=9_pat
mview store, object=X_mat
mview store, object=X_pat
turn y,180
mview store

frame 240
mview store, object=1_mat
mview store, object=1_pat
mview store, object=10_mat
mview store, object=10_pat
mview store, object=11_mat
mview store, object=11_pat
mview store, object=12_mat
mview store, object=12_pat
mview store, object=13_mat
mview store, object=13_pat
mview store, object=14_mat
mview store, object=14_pat
mview store, object=15_mat
mview store, object=15_pat
mview store, object=16_mat
mview store, object=16_pat
mview store, object=17_mat
mview store, object=17_pat
mview store, object=18_mat
mview store, object=18_pat
mview store, object=19_mat
mview store, object=19_pat
mview store, object=2_mat
mview store, object=2_pat
mview store, object=20_mat
mview store, object=20_pat
mview store, object=21_mat
mview store, object=21_pat
mview store, object=22_mat
mview store, object=22_pat
mview store, object=3_mat
mview store, object=3_pat
mview store, object=4_mat
mview store, object=4_pat
mview store, object=5_mat
mview store, object=5_pat
mview store, object=6_mat
mview store, object=6_pat
mview store, object=7_mat
mview store, object=7_pat
mview store, object=8_mat
mview store, object=8_pat
mview store, object=9_mat
mview store, object=9_pat
mview store, object=X_mat
mview store, object=X_pat
turn y,360
mview store

mview reinterpolate, power=1