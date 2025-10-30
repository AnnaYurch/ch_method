set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'classic_derivatives.png'
set title 'Функция и ее производные' font 'Arial,16'
set xlabel 'x' font 'Arial,14'
set ylabel 'y' font 'Arial,14'
set grid
set key top left box
set xrange [0.1:3.0]
set yrange [-4:8]
set style line 1 lc rgb 'red' lw 3
set style line 2 lc rgb 'blue' lw 3
set style line 3 lc rgb 'green' lw 3
plot 'classic_data.txt' index 0 with lines ls 1 title 'y = f(x)', \
     'classic_data.txt' index 1 with lines ls 2 title 'y'' = f''(x)', \
     'classic_data.txt' index 2 with lines ls 3 title 'y'''' = f''''(x)'
