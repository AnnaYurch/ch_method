set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'second_derivative.png'
set title 'Вторая производная кубического сплайна' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'S''''(x)' font 'Arial,12'
set grid
set key top left box
set xrange [-0.0:3.9]
set yrange [-82.9:63.9]
set style line 1 lc rgb 'purple' lw 3
set style line 2 lc rgb 'orange' pt 9 ps 2
set style line 3 lc rgb 'black' pt 7 ps 1.5
plot 'second_derivative_data.txt' index 1 with lines ls 1 title 'Вторая производная', \
     'second_derivative_data.txt' index 2 with points ls 2 title 'x*=1.285 (d²S/dx²=26.025)', \
     'second_derivative_data.txt' index 0 with points ls 3 title 'Исходные точки'
