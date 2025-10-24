set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'system_plot.png'
set multiplot layout 2,2 title 'Система нелинейных уравнений' font 'Arial,14'
set xlabel 'x1'
set ylabel 'x2'

set title 'f1(x1,x2) = x1^2 + x2^2 - 5sin(x1) - 9 = 0'
set contour base
set cntrparam levels discrete 0
unset surface
set view map
splot 'system_plot.dat' using 1:2:3 with lines title 'f1=0' lw 2

set title 'f2(x1,x2) = 2x1^2 + 2x1x2 - 3x2^2 - 4x1 - x2cos(x1) + 3 = 0'
splot 'system_plot.dat' using 1:2:4 with lines title 'f2=0' lw 2

set title 'Пересечения: f1=0 и f2=0 (корни системы)'
splot 'system_plot.dat' using 1:2:3 with lines title 'f1=0' lw 2, \
      '' using 1:2:4 with lines title 'f2=0' lw 2

set title '3D поверхности f1 и f2'
set surface
unset contour
set view 60,30
splot 'system_plot.dat' using 1:2:3 with lines title 'f1(x1,x2)', \
      '' using 1:2:4 with lines title 'f2(x1,x2)'

unset multiplot
