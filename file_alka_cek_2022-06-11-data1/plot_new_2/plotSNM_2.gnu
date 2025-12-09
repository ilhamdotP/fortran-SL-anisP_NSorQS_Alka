set term postscript portrait enh color  "Times-Roman" 16
set output "plotSNM_2.eps"

set size 1.0,1.0

set multiplot

set lmargin 10
set rmargin 8

set origin 0.0,0.5
set size 1,0.5
set xlabel '{/Symbol r}_N/{/Symbol r}_0' font ",20"
set ylabel 'P(MeV fm^{-3})' font ",20"
set key t l vertical Left reverse enhanced autotitle nobox font ",13"
set xrange [1:5]
set xtics 1,1,5 offset 0,0.2
set yrange [0:260]
set ytics 0,50,200 

plot "daniel.dat" u 2:3 with filledcurve closed fc "grey" fs solid 0.5 border lc "black" dt 3 lw 2 notitle,\
	 "new_result/basic_BSR12_SNM.dat" u 1:5 w l lw 4 dt 1 lc "red" t "SNM basic",\
     "new_result/basic_BSR12_SNM_GUP10.dat" u 1:5 w l lw 4 dt 2 lc "grey10" t "{/Symbol b}=10^{-7}",\
     "new_result/basic_BSR12_SNM_GUP20.dat" u 1:5 w l lw 4 dt 3 lc "grey20" t "{/Symbol b}=2x10^{-7}",\
	 "new_result/basic_BSR12_SNM_GUP30.dat" u 1:5 w l lw 4 dt 4 lc "grey30" t "{/Symbol b}=3x10^{-7}",\
	 "new_result/basic_BSR12_SNM_GUP40.dat" u 1:5 w l lw 4 dt 5 lc "grey40" t "{/Symbol b}=4x10^{-7}",\
#	 "new_result/basic_BSR12_SNM_GUP50.dat" u 1:5 w l lw 4 dt 6 lc "grey50" t "{/Symbol b}=5x10^{-7}",\
     
     
     
set origin 0.0,0.0
set size 1,0.5
set xlabel '{/Symbol r}_N (fm^{-3})' font ",20"
set ylabel 'E/N(MeV)' font ",20"
set tics font ",15"
set xrange [0.1:0.4]
set xtics 0.1,0.1,0.4 offset 0,0.2
set yrange [-20:0]
set ytics -20,5,0 offset 0,0.2
unset key

plot "EOS1.pngX.dat" with filledcurve closed fc "grey" fs pattern 2 border lc "black" dt 3 lw 3,\
     "EOS2.pngX.dat" w l lc "black" dt 3 lw 2,\
	 "new_result/basic_BSR12_SNM.dat" u 2:3 w l lw 4 dt 1 lc "red" t "SNM basic",\
     "new_result/basic_BSR12_SNM_GUP10.dat" u 2:3 w l lw 4 dt 2 lc "grey10" t "SNM GUP {/Symbol b}=10^{-7}",\
     "new_result/basic_BSR12_SNM_GUP20.dat" u 2:3 w l lw 4 dt 3 lc "grey20" t "SNM GUP {/Symbol b}=2x10^{-7}",\
	 "new_result/basic_BSR12_SNM_GUP30.dat" u 2:3 w l lw 4 dt 4 lc "grey30" t "SNM GUP {/Symbol b}=3x10^{-7}",\
	 "new_result/basic_BSR12_SNM_GUP40.dat" u 2:3 w l lw 4 dt 5 lc "grey40" t "SNM GUP {/Symbol b}=4x10^{-7}",\
#	 "new_result/basic_BSR12_SNM_GUP50.dat" u 2:3 w l lw 4 dt 6 lc "grey50" t "SNM GUP {/Symbol b}=5x10^{-7}",\
     "Brown2014.dat" u 1:2:($2-$3):($2+$4) w yerrorlines pt 7 lc "red"
     
      
      
unset multiplot

