set term postscript portrait enh color  "Times-Roman" 16
set output "plotSNM.eps"

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
     "old_result/basic_BSR12_SNM_GUP15.dat" u 1:5 w l lw 4 t "{/Symbol b}=1.5x10^{-7}",\
     "old_result/basic_BSR12_SNM_GUP10.dat" u 1:5 w l lw 4 t "{/Symbol b}=10^{-7}",\
     "old_result/basic_BSR12_SNM_GUP5.dat" u 1:5 w l lw 4 t "{/Symbol b}=0.5x10^{-7}",\
     "old_result/basic_BSR12_SNM_GUP1.dat" u 1:5 w l lw 4 t "{/Symbol b}=10^{-8}",\
     "old_result/basic_BSR12_SNM_GUP0.dat" u 1:5 w l lw 4 t "{/Symbol b}=0",\
     "old_result/basic_BSR12_SNM.dat" u 1:5 w l lw 4 dt 2 t "SNM basic",\
     
     
     
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
     "old_result/basic_BSR12_SNM_GUP15.dat" u 2:3 w l lw 4 t "SNM GUP {/Symbol b}=1.5x10^{-7}",\
     "old_result/basic_BSR12_SNM_GUP10.dat" u 2:3 w l lw 4 t "SNM GUP {/Symbol b}=10^{-7}",\
     "old_result/basic_BSR12_SNM_GUP5.dat" u 2:3 w l lw 4 t "SNM GUP {/Symbol b}=0.5x10^{-7}",\
     "old_result/basic_BSR12_SNM_GUP1.dat" u 2:3 w l lw 4 t "SNM GUP {/Symbol b}=10^{-8}",\
     "old_result/basic_BSR12_SNM_GUP0.dat" u 2:3 w l lw 4 t "SNM GUP {/Symbol b}=0",\
     "old_result/basic_BSR12_SNM.dat" u 2:3 w l lw 4 dt 2 t "SNM basic",\
     "EOS2.pngX.dat" w l lc "black" dt 3 lw 2,\
     "Brown2014.dat" u 1:2:($2-$3):($2+$4) w yerrorlines pt 7 lc "red"
      
      
unset multiplot

