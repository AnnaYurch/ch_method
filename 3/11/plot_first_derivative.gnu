set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'first_derivative.png'
set title 'Первая производная кубического сплайна' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'S''(x)' font 'Arial,12'
set grid
set key top left box
set xrange [-0.0:3.9]
set yrange [-11.3:8.7]
set style line 1 lc rgb 'blue' lw 3
set style line 2 lc rgb 'red' pt 9 ps 2
set style line 3 lc rgb 'black' pt 7 ps 1.5
plot 'derivative_data.txt' index 1 with lines ls 1 title 'Первая производная S'(x)', \
     'derivative_data.txt' index 2 with points ls 2 title 'x*=1.285 (S'=-0.442)', \
     'derivative_data.txt' index 0 with points ls 3 title 'Исходные точки'
