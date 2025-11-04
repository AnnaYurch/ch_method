set terminal pngcairo size 1200,900 enhanced font 'Arial,12'
set output 'separate_derivatives.png'
set multiplot layout 3,1 title 'Функция и ее производные' font 'Arial,16'
set title 'Исходная функция y = f(x)'
set xlabel 'x'
set ylabel 'y'
set grid
set xrange [0.1:3.0]
plot 'separate_data.txt' index 0 with lines lw 3 lc rgb 'red' title 'f(x)'
set title 'Первая производная y'' = f''(x)'
set xlabel 'x'
set ylabel 'y'''
set grid
set xrange [0.1:3.0]
plot 'separate_data.txt' index 1 with lines lw 3 lc rgb 'blue' title 'f''(x)'
set title 'Вторая производная y'''' = f''''(x)'
set xlabel 'x'
set ylabel 'y''''
set grid
set xrange [0.1:3.0]
plot 'separate_data.txt' index 2 with lines lw 3 lc rgb 'green' title 'f''''(x)'
unset multiplot
