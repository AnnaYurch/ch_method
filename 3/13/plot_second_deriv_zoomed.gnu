set terminal pngcairo size 1800,600 enhanced font 'Arial,12'
set output 'second_derivative_zoomed.png'
set multiplot layout 1,3 title 'Увеличенные фрагменты второй производной (h=0.150)' font 'Arial,16'
set title 'Диапазон [0.5:1.5]'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key top left box
set xrange [0.5:1.5]
plot 'second_deriv_zoomed.txt' index 0 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'second_deriv_zoomed.txt' index 1 with lines lw 2 lc rgb 'green' dt 2 title 'Центральная разность', \
     'second_deriv_zoomed.txt' index 2 with lines lw 2 lc rgb 'purple' dt 3 title 'Пятиточечная'
set title 'Диапазон [1.0:2.0]'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key top left box
set xrange [1.0:2.0]
plot 'second_deriv_zoomed.txt' index 3 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'second_deriv_zoomed.txt' index 4 with lines lw 2 lc rgb 'green' dt 2 title 'Центральная разность', \
     'second_deriv_zoomed.txt' index 5 with lines lw 2 lc rgb 'purple' dt 3 title 'Пятиточечная'
set title 'Диапазон [1.5:2.5]'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key top left box
set xrange [1.5:2.5]
plot 'second_deriv_zoomed.txt' index 6 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'second_deriv_zoomed.txt' index 7 with lines lw 2 lc rgb 'green' dt 2 title 'Центральная разность', \
     'second_deriv_zoomed.txt' index 8 with lines lw 2 lc rgb 'purple' dt 3 title 'Пятиточечная'
unset multiplot
