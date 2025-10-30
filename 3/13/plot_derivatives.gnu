set terminal pngcairo size 1600,1200 enhanced font 'Arial,10'
set output 'derivative_comparison.png'
set multiplot layout 2,2 title 'Сравнение численных и аналитических производных (h=0.150)' font 'Arial,14'
set title 'Исходная функция'
set xlabel 'x'
set ylabel 'f(x)'
set grid
plot 'derivative_data.txt' index 0 with lines lw 2 lc rgb 'black' title 'f(x)'
set title 'Первая производная'
set xlabel 'x'
set ylabel 'f''(x)'
set grid
plot 'derivative_data.txt' index 1 with lines lw 2 lc rgb 'blue' title 'Аналитическая', \
     'derivative_data.txt' index 2 with lines lw 2 lc rgb 'red' dt 2 title 'Численная (центральная разность)'
set title 'Вторая производная'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
plot 'derivative_data.txt' index 3 with lines lw 2 lc rgb 'green' title 'Аналитическая', \
     'derivative_data.txt' index 4 with lines lw 2 lc rgb 'orange' dt 2 title 'Численная (центральная разность)'
set title 'Абсолютные ошибки'
set xlabel 'x'
set ylabel 'Ошибка'
set grid
set logscale y
plot 'derivative_data.txt' index 5 with lines lw 2 lc rgb 'red' title 'Ошибка первой производной', \
     'derivative_data.txt' index 6 with lines lw 2 lc rgb 'blue' title 'Ошибка второй производной'
unset multiplot
