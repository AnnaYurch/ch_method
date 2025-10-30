set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'numerical_comparison.png'
set multiplot layout 2,1 title 'Сравнение аналитических и численных производных (h=0.150)' font 'Arial,14'
set title 'Первая производная'
set xlabel 'x'
set ylabel 'f''(x)'
set grid
set key top left
plot 'numerical_data.txt' index 0 with lines lw 2 lc rgb 'blue' title 'Аналитическая', \
     'numerical_data.txt' index 1 with lines lw 2 lc rgb 'red' dt 2 title 'Численная (центральная разность)'
set title 'Вторая производная'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key top left
plot 'numerical_data.txt' index 2 with lines lw 2 lc rgb 'green' title 'Аналитическая', \
     'numerical_data.txt' index 3 with lines lw 2 lc rgb 'orange' dt 2 title 'Численная (центральная разность)'
unset multiplot
