set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'
set output 'second_derivative_comparison.png'
set multiplot layout 2,2 title 'Вторая производная: сравнение методов (h=0.075)' font 'Arial,14'
set title 'Сравнение методов второй производной'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key top left
plot 'second_derivative_data.txt' index 0 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'second_derivative_data.txt' index 1 with lines lw 2 lc rgb 'red' title 'Правая разность', \
     'second_derivative_data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'Левая разность', \
     'second_derivative_data.txt' index 3 with lines lw 2 lc rgb 'green' title 'Центральная разность', \
     'second_derivative_data.txt' index 4 with lines lw 2 lc rgb 'orange' title 'Пятиточечная'
set title 'Абсолютные ошибки методов (линейная шкала)'
set xlabel 'x'
set ylabel 'Ошибка'
set grid
set key top left
plot 'second_derivative_data.txt' index 5 with lines lw 2 lc rgb 'red' title 'Ошибка правой разности', \
     'second_derivative_data.txt' index 6 with lines lw 2 lc rgb 'blue' title 'Ошибка левой разности', \
     'second_derivative_data.txt' index 7 with lines lw 2 lc rgb 'green' title 'Ошибка центральной разности', \
     'second_derivative_data.txt' index 8 with lines lw 2 lc rgb 'orange' title 'Ошибка пятиточечной'
set title 'Абсолютные ошибки методов (логарифмическая шкала)'
set xlabel 'x'
set ylabel 'Ошибка'
set grid
set logscale y
set key top left
plot 'second_derivative_data.txt' index 5 with lines lw 2 lc rgb 'red' title 'Ошибка правой разности', \
     'second_derivative_data.txt' index 6 with lines lw 2 lc rgb 'blue' title 'Ошибка левой разности', \
     'second_derivative_data.txt' index 7 with lines lw 2 lc rgb 'green' title 'Ошибка центральной разности', \
     'second_derivative_data.txt' index 8 with lines lw 2 lc rgb 'orange' title 'Ошибка пятиточечной'
set title 'Лучшие методы: центральная и пятиточечная'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key top left
unset logscale y
plot 'second_derivative_data.txt' index 0 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'second_derivative_data.txt' index 3 with lines lw 2 lc rgb 'green' title 'Центральная разность', \
     'second_derivative_data.txt' index 4 with lines lw 2 lc rgb 'orange' title 'Пятиточечная'
unset multiplot
