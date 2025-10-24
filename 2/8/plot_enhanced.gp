set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'
set output 'enhanced_system_plot.png'
set multiplot layout 2,2 title 'Анализ системы нелинейных уравнений' font 'Arial,16'
set xlabel 'x₁'
set ylabel 'x₂'
set xrange [-2.5:3.5]
set yrange [-3:3]
set grid

set title 'f₁(x₁,x₂) = x₁² + x₂² - 5sin(x₁) - 9 = 0'
set contour base
set cntrparam levels discrete 0
unset surface
set view map
splot 'system_data.txt' using 1:2:3 with lines lw 2 lc rgb 'blue' notitle

set title 'f₂(x₁,x₂) = 2x₁² + 2x₁x₂ - 3x₂² - 4x₁ - x₂cos(x₁) + 3 = 0'
splot 'system_data.txt' using 1:2:4 with lines lw 2 lc rgb 'red' notitle

set title 'Корни системы (пересечения кривых)'
splot 'system_data.txt' using 1:2:3 with lines lw 2 lc rgb 'blue' title 'f₁=0', \
      '' using 1:2:4 with lines lw 2 lc rgb 'red' title 'f₂=0', \
      'system_roots.txt' using 1:2:(0) with points pt 7 ps 2 lc rgb 'green' title 'Корни'

set title '3D поверхности (нулевая плоскость)'
set contour base
set cntrparam levels discrete 0
set surface
set view 60,30,1,1
splot 'system_data.txt' using 1:2:3 with lines lw 1 lc rgb 'blue' title 'f₁', \
      '' using 1:2:4 with lines lw 1 lc rgb 'red' title 'f₂', \
      0 with lines lw 1 lc rgb 'black' title 'z=0'

unset multiplot
