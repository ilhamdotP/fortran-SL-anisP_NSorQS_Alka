set multiplot
set origin 0, 0
set size 1, 1
set xlabel 'P(MeV/fm^3)'
set ylabel '{/Symbol e}(MeV/fm^3)'
unset grid

set key t l vertical Left reverse enhanced autotitle nobox font ",9"

#
plot [0:800][0:1800] 'EOS_NSM_FSU2HZ2.dat' u 3:4 w l lw 4 t 'FSU 2Hz', 'EOS_NSM_FSU2HZ2_GUP.dat' u 3:4 w l lw 4 t '+GUP {/Symbol b} = 0', 'EOS_NSM_FSU2HZ2_GUP1D-7.dat' u 3:4 w l lw 4 t '+GUP {/Symbol b} = 10^{-7}', 'EOS_NSM_FSU2HZ2_GUP1D-8.dat' u 3:4 w l lw 4 t '+GUP {/Symbol b} = 10^{-8}', 'EOS_NSM_FSU2HZ2_GUP1D-7a.dat' u 3:4 w l lw 4 dt 2 lc 'blue' t '+GUP {/Symbol b} = 10^{-7}(analytic P)', 'EOS_NSM_FSU2HZ2_GUP1D-8a.dat' u 3:4 w l lw 4 dt 1 lc 'black' t '+GUP {/Symbol b} = 10^{-8}(analytic P)'


set origin 0.45, 0.2
set size 0.5, 0.6
unset key
unset xlabel
unset ylabel
set grid

plot[-1:1][0:150] 'EOS_NSM_FSU2HZ2.dat' u 3:4 w l lw 1.5, 'EOS_NSM_FSU2HZ2_GUP.dat' u 3:4 w l lw 1.5, 'EOS_NSM_FSU2HZ2_GUP1D-7.dat' u 3:4 w l lw 1.5 , 'EOS_NSM_FSU2HZ2_GUP1D-8.dat' u 3:4 w l lw 1.5 , 'EOS_NSM_FSU2HZ2_GUP1D-7a.dat' u 3:4 w l lw 1.5 lc 'blue' dt 1, 'EOS_NSM_FSU2HZ2_GUP1D-8a.dat' u 3:4 w l lw 1.5 lc 'black' dt 1

unset multiplot
