set terminal pngcairo size 1400,1000 enhanced font 'Arial,12'
set output 'first_derivative_methods_h_0.150.png'
set title 'ПЕРВАЯ ПРОИЗВОДНАЯ - Сравнение методов (h=0.150)' font 'Arial,16'
set xlabel 'x' font 'Arial,14'
set ylabel 'f''(x)' font 'Arial,14'
set grid linecolor rgb '#dddddd' linewidth 1
set key outside center top horizontal box font 'Arial,11'
set xrange [0.1:3.0]
set style line 1 lc rgb '#000000' lw 4 title 'Аналитическая'
set style line 2 lc rgb '#FF0000' lw 2 dt 1 title 'Правая разность'
set style line 3 lc rgb '#0000FF' lw 2 dt 1 title 'Левая разность'
set style line 4 lc rgb '#00AA00' lw 3 dt 1 title 'Центральная разность'
set style line 5 lc rgb '#FF00FF' lw 2 dt 2 title '3-точ. вперед'
set style line 6 lc rgb '#FFA500' lw 2 dt 3 title '3-точ. назад'
set style line 7 lc rgb '#8B4513' lw 2 dt 4 title '4-точ.'
plot 'visible_methods_data.txt' index 0 with lines ls 1, \
     'visible_methods_data.txt' index 1 with lines ls 2, \
     'visible_methods_data.txt' index 2 with lines ls 3, \
     'visible_methods_data.txt' index 3 with lines ls 4, \
     'visible_methods_data.txt' index 4 with lines ls 5, \
     'visible_methods_data.txt' index 5 with lines ls 6, \
     'visible_methods_data.txt' index 6 with lines ls 7
