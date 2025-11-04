set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'spline_derivative2.png'
set title 'Вторая производная кубического сплайна' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'S''(x)' font 'Arial,12'
set grid
set key top left box
set xrange [-0.0:3.9]
set yrange [-71.1:52.1]
set style line 1 lc rgb 'black' pt 7 ps 1.5
set style line 2 lc rgb 'purple' lw 2
set style line 3 lc rgb 'orange' pt 9 ps 2
plot 'spline_data.txt' index 0 with points ls 1 title 'Исходные точки', \
     'spline_data.txt' index 3 with lines ls 2 title 'Вторая производная', \
     '+' using (1.285):(0) with points ls 3 title 'x*=1.285'
