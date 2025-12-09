set term postscript portrait enh color  "Times-Roman" 16
set output "plotPNM.eps"

set size 1.0,1.0

set multiplot

set lmargin 10
set rmargin 8

set origin 0.0,0.5
set size 1,0.5
set xlabel '{/Symbol r}_N/{/Symbol r}_0' font ",20"
set ylabel 'P(MeV fm^{-3})' font ",20"
set key t l vertical Left reverse enhanced autotitle nobox font ",11"
set xrange [1:5]
set xtics 1,1,5 offset 0,0.2
set yrange [0:260]
set ytics 0,50,200 

plot "daniel.dat" u 2:4 with filledcurve closed fc "grey" fs solid 0.5 border lc "black" dt 3 lw 2 notitle,\
     "daniel.dat" u 2:5 with filledcurve closed fc "grey" fs pattern 2 border lc "black" dt 3 lw 2 notitle,\
     "basic_BSR12_PNM_GUP10.dat" u 1:5 w l lw 4 t "{/Symbol b}=10^{-7}",\
     "basic_BSR12_PNM_GUP1.dat" u 1:5 w l lw 4 t "{/Symbol b}=10^{-8}",\
     "basic_BSR12_PNM_GUP0.dat" u 1:5 w l lw 4 t "{/Symbol b}=0",\
     "basic_BSR12_PNM.dat" u 1:5 w l lw 4 dt 2 t "PNM basic",\
     "daniel.dat" u 2:4 w l lc "blue" dt 3 lw 2 notitle,\
     
     
     
set origin 0.0,0.0
set size 1,0.5
set xlabel '{/Symbol r}_N (fm^{-3})' font ",20"
set ylabel 'E/A(MeV)' font ",20"
set tics font ",15"
set xrange [0:0.25]
set xtics 0,0.05,0.25 offset 0,0.2
set yrange [0:30]
set ytics 0,6,30 offset 0,0.2
unset key

plot "1fUMhz.png.dat" with filledcurve closed fc "yellow" fs pattern 2,\
     "1fUMhz.png.dat" with line lc "black" dt 3 lw 3,\
     "basic_BSR12_PNM_GUP10.dat" u 2:3 w l lw 4 t "PNM GUP {/Symbol b}=10^{-7}",\
     "basic_BSR12_PNM_GUP1.dat" u 2:3 w l lw 4 t "PNM GUP {/Symbol b}=10^{-7}",\
     "basic_BSR12_PNM_GUP0.dat" u 2:3 w l lw 4 t "PNM GUP {/Symbol b}=0",\
     "basic_BSR12_PNM.dat" u 2:3 w l lw 4 dt 2 t "PNM basic"
      
      
unset multiplot

