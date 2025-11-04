set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'simple_second_derivatives.png'
set title 'Вторая производная: сравнение методов (h=0.075)' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'f''''(x)' font 'Arial,12'
set grid
set key top left box
set xrange [0.1:3.0]
plot 'simple_second_derivatives.txt' using 1:2 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     '' using 1:3 with lines lw 2 lc rgb 'red' title 'Правая разность', \
     '' using 1:4 with lines lw 2 lc rgb 'blue' title 'Левая разность', \
     '' using 1:5 with lines lw 2 lc rgb 'green' title 'Центральная разность', \
     '' using 1:6 with lines lw 2 lc rgb 'orange' title 'Пятиточечная'
