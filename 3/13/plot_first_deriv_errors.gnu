set terminal png size 1200,800 enhanced font 'Arial,10'
set output 'first_derivative_errors.png'
set title 'Ошибки методов первой производной (h=0.075)' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'Абсолютная ошибка' font 'Arial,12'
set grid
set key top left box
set logscale y
set xrange [0.25:2.85]
plot 'first_deriv_errors_data.txt' using 1:2 with lines lw 2 lc rgb 'red' title 'Ошибка правой разности', \
     '' using 1:3 with lines lw 2 lc rgb 'blue' title 'Ошибка левой разности', \
     '' using 1:4 with lines lw 2 lc rgb 'green' title 'Ошибка центральной разности', \
     '' using 1:5 with lines lw 2 lc rgb 'orange' title 'Ошибка 3-точ. вперед', \
     '' using 1:6 with lines lw 2 lc rgb 'purple' title 'Ошибка 3-точ. назад', \
     '' using 1:7 with lines lw 2 lc rgb 'brown' title 'Ошибка 4-точ.'
