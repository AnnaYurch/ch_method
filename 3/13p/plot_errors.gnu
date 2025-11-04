set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'errors_plot.png'
set title 'Ошибки методов второй производной (h=0.075)' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'Абсолютная ошибка' font 'Arial,12'
set grid
set key top left box
set logscale y
set xrange [0.1:3.0]
plot 'errors_data.txt' using 1:2 with lines lw 2 lc rgb 'red' title 'Ошибка правой разности', \
     '' using 1:3 with lines lw 2 lc rgb 'blue' title 'Ошибка левой разности', \
     '' using 1:4 with lines lw 2 lc rgb 'green' title 'Ошибка центральной разности', \
     '' using 1:5 with lines lw 2 lc rgb 'orange' title 'Ошибка пятиточечной'
