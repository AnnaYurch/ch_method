set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'lagrange_graph.png'
set title 'Интерполяция Лагранжа - Вариант 44' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'y' font 'Arial,12'
set grid
set key top left box
set xrange [-1.5:2.5]
set yrange [-3:2]
plot 'lagrange_data.txt' index 0 with points pt 7 ps 1.5 lc rgb 'black' title 'Все исходные точки', \
     'lagrange_data.txt' index 1 with lines lw 2 lc rgb 'red' title 'L2(x) (2-я степень)', \
     'lagrange_data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'L3(x) (3-я степень)', \
     'lagrange_data.txt' index 3 with points pt 2 ps 2 lc rgb 'dark-red' title 'Узлы L2', \
     'lagrange_data.txt' index 4 with points pt 4 ps 2 lc rgb 'dark-blue' title 'Узлы L3', \
     'lagrange_data.txt' index 5 with points pt 9 ps 3 lc rgb 'green' title 'x* = 1.708'
