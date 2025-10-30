set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'individual_second_derivatives.png'
set title 'Вторая производная: отдельные методы (h=0.075)' font 'Arial,14'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key top left
plot 'individual_second_derivatives.txt' using 1:2 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     '' using 1:3 with lines lw 2 lc rgb 'red' title 'Правая разность', \
     '' using 1:4 with lines lw 2 lc rgb 'blue' title 'Левая разность', \
     '' using 1:5 with lines lw 2 lc rgb 'green' title 'Центральная разность', \
     '' using 1:6 with lines lw 2 lc rgb 'orange' title 'Пятиточечная'
