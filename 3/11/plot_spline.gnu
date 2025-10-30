set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'spline_graph.png'
set title 'Естественный кубический сплайн дефекта 1' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'y' font 'Arial,12'
set grid
set key top left box
set xrange [-0.0:3.9]
set yrange [0.2:4.2]
set style line 1 lc rgb 'black' pt 7 ps 1.5
set style line 2 lc rgb 'red' lw 2
set style line 3 lc rgb 'blue' pt 9 ps 2
plot 'spline_data.txt' index 0 with points ls 1 title 'Исходные точки', \
     'spline_data.txt' index 1 with lines ls 2 title 'Кубический сплайн', \
     'spline_data.txt' index 2 with points ls 3 title 'x*=1.285 (S=0.543)'
