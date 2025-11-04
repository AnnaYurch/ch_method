set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'newton_graph.png'
set title 'Интерполяция Ньютона с выбором узлов - Вариант 44' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'y' font 'Arial,12'
set grid
set key top left box
set xrange [-1.5:2.5]
set yrange [-3:2]
plot 'newton_data.txt' index 0 with points pt 7 ps 1.5 lc rgb 'black' title 'Все исходные точки', \
     'newton_data.txt' index 1 with lines lw 2 lc rgb 'red' title 'P2(x) (2-я степень)', \
     'newton_data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'P3(x) (3-я степень)', \
     'newton_data.txt' index 3 with points pt 2 ps 2 lc rgb 'dark-red' title 'Узлы P2', \
     'newton_data.txt' index 4 with points pt 4 ps 2 lc rgb 'dark-blue' title 'Узлы P3', \
     'newton_data.txt' index 5 with points pt 9 ps 3 lc rgb 'green' title 'x*=1.708'
