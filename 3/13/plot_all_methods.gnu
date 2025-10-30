set terminal pngcairo size 1600,1200 enhanced font 'Arial,10'
set output 'all_methods_comparison.png'
set multiplot layout 2,1 title 'Сравнение всех методов численного дифференцирования (h=0.075)' font 'Arial,14'
set title 'ПЕРВАЯ ПРОИЗВОДНАЯ - Все методы'
set xlabel 'x'
set ylabel 'f''(x)'
set grid
set key outside center top horizontal
plot 'all_methods_data.txt' index 0 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'all_methods_data.txt' index 1 with lines lw 2 lc rgb 'red' title 'Правая разность', \
     'all_methods_data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'Левая разность', \
     'all_methods_data.txt' index 3 with lines lw 2 lc rgb 'green' title 'Центральная разность', \
     'all_methods_data.txt' index 4 with lines lw 2 lc rgb 'orange' title '3-точ. вперед', \
     'all_methods_data.txt' index 5 with lines lw 2 lc rgb 'purple' title '3-точ. назад', \
     'all_methods_data.txt' index 6 with lines lw 2 lc rgb 'brown' title '4-точ.'
set title 'ВТОРАЯ ПРОИЗВОДНАЯ - Все методы'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key outside center top horizontal
plot 'all_methods_data.txt' index 7 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'all_methods_data.txt' index 8 with lines lw 2 lc rgb 'red' title 'Правая разность', \
     'all_methods_data.txt' index 9 with lines lw 2 lc rgb 'blue' title 'Левая разность', \
     'all_methods_data.txt' index 10 with lines lw 2 lc rgb 'green' title 'Центральная разность', \
     'all_methods_data.txt' index 11 with lines lw 2 lc rgb 'orange' title '5-точ.'
unset multiplot
