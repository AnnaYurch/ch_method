set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'newton_graph.png'
set title 'Интерполяция Ньютона (вторая формула) - Вариант 44' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'y' font 'Arial,12'
set grid
set key top left box
set xrange [-1.5:2.5]
set yrange [-3:2]
set style line 1 lc rgb 'black' pt 7 ps 1.5
set style line 2 lc rgb 'red' lw 2
set style line 3 lc rgb 'blue' lw 2
set style line 4 lc rgb 'dark-red' pt 2 ps 2
set style line 5 lc rgb 'dark-blue' pt 4 ps 2
set style line 6 lc rgb 'green' pt 9 ps 2
plot 'newton_data.txt' index 0 with points ls 1 title 'Все исходные точки', \
     'newton_data.txt' index 1 with lines ls 2 title 'P2(x) (2-я степень)', \
     'newton_data.txt' index 2 with lines ls 3 title 'P3(x) (3-я степень)', \
     'newton_data.txt' index 3 with points ls 4 title 'Узлы P2', \
     'newton_data.txt' index 4 with points ls 5 title 'Узлы P3', \
     'newton_data.txt' index 5 with points ls 6 title 'x*=1.708'
