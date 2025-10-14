set terminal pngcairo size 1000,600 enhanced font 'Arial,12'
set output 'equation_plot.png'
set title 'График функции f(x) = 2x^{3} - 5x^{2} - sin(3x+2) + 4'
set xlabel 'x'
set ylabel 'f(x)'
set grid
set key top left
set zeroaxis lt -1 lc 'black'
plot 'function_data.txt' with lines lw 2 lc rgb 'blue' title 'f(x)', \
     'roots.txt' with points pt 7 ps 1.5 lc rgb 'red' title 'Найденные корни (3 шт.)'
