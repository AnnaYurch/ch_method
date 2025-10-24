set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'system_plot.png'
set title 'Система нелинейных уравнений'
set xlabel 'x₁'
set ylabel 'x₂'
set xrange [-2.5:3.5]
set yrange [-3:3]
set grid
set key outside top center horizontal
set zeroaxis lt -1 lc 'black'

set contour base
set cntrparam levels discrete 0
unset surface
set view map
splot 'system_data.txt' using 1:2:3 with lines lw 3 lc rgb 'blue' title 'f₁(x₁,x₂) = x₁² + x₂² - 5sin(x₁) - 9 = 0', \
      '' using 1:2:4 with lines lw 3 lc rgb 'red' title 'f₂(x₁,x₂) = 2x₁² + 2x₁x₂ - 3x₂² - 4x₁ - x₂cos(x₁) + 3 = 0', \
      'system_roots.txt' using 1:2:(0) with points pt 7 ps 2 lc rgb 'green' title 'Найденные корни (3 шт.)'
