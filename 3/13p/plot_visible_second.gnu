set terminal pngcairo size 1400,1000 enhanced font 'Arial,12'
set output 'second_derivative_methods_h_0.150.png'
set title 'ВТОРАЯ ПРОИЗВОДНАЯ - Сравнение методов (h=0.150)' font 'Arial,16'
set xlabel 'x' font 'Arial,14'
set ylabel 'f''''(x)' font 'Arial,14'
set grid linecolor rgb '#dddddd' linewidth 1
set key outside center top horizontal box font 'Arial,11'
set xrange [0.1:3.0]
set style line 1 lc rgb '#000000' lw 4 title 'Аналитическая'
set style line 2 lc rgb '#FF0000' lw 2 dt 1 title 'Правая разность'
set style line 3 lc rgb '#0000FF' lw 2 dt 1 title 'Левая разность'
set style line 4 lc rgb '#00AA00' lw 3 dt 1 title 'Центральная разность'
set style line 5 lc rgb '#FF00FF' lw 2 dt 2 title '5-точ.'
plot 'visible_methods_data.txt' index 7 with lines ls 1, \
     'visible_methods_data.txt' index 8 with lines ls 2, \
     'visible_methods_data.txt' index 9 with lines ls 3, \
     'visible_methods_data.txt' index 10 with lines ls 4, \
     'visible_methods_data.txt' index 11 with lines ls 5
