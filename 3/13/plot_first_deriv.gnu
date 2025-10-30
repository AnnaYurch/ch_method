set terminal png size 1200,800 enhanced font 'Arial,10'
set output 'first_derivative_plot.png'
set title 'Первая производная: сравнение методов (h=0.075)' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'f'(x)' font 'Arial,12'
set grid
set key top left box
set xrange [0.25:2.85]
plot 'first_deriv_data.txt' using 1:2 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     '' using 1:3 with lines lw 2 lc rgb 'red' title 'Правая разность', \
     '' using 1:4 with lines lw 2 lc rgb 'blue' title 'Левая разность', \
     '' using 1:5 with lines lw 2 lc rgb 'green' title 'Центральная разность', \
     '' using 1:6 with lines lw 2 lc rgb 'orange' title '3-точ. вперед', \
     '' using 1:7 with lines lw 2 lc rgb 'purple' title '3-точ. назад', \
     '' using 1:8 with lines lw 2 lc rgb 'brown' title '4-точ.'
