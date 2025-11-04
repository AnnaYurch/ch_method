set terminal pngcairo size 2000,1200 enhanced font 'Arial,12'
set output 'detailed_derivatives.png'
set multiplot layout 2,1 title 'Сравнение методов численного дифференцирования (h=0.150)' font 'Arial,16'
set title 'ПЕРВЫЕ ПРОИЗВОДНЫЕ (область без краевых эффектов)' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'f''(x)' font 'Arial,12'
set grid
set key outside right top vertical box
set xrange [0.47:2.9]
set yrange [-3:3]
plot 'detailed_derivatives.txt' index 0 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'detailed_derivatives.txt' index 1 with lines lw 2 lc rgb 'red' title 'Правая разность', \
     'detailed_derivatives.txt' index 2 with lines lw 2 lc rgb 'blue' title 'Левая разность', \
     'detailed_derivatives.txt' index 3 with lines lw 2 lc rgb 'green' title 'Центральная', \
     'detailed_derivatives.txt' index 4 with lines lw 2 lc rgb 'purple' title '3-точ. вперед', \
     'detailed_derivatives.txt' index 5 with lines lw 2 lc rgb 'orange' title '3-точ. назад', \
     'detailed_derivatives.txt' index 6 with lines lw 2 lc rgb 'brown' title '4-точ.'
set title 'ВТОРЫЕ ПРОИЗВОДНЫЕ (область без краевых эффектов)' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'f''''(x)' font 'Arial,12'
set grid
set key outside right top vertical box
set xrange [0.55:2.7]
set yrange [-4:7]
plot 'detailed_derivatives.txt' index 7 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     'detailed_derivatives.txt' index 8 with lines lw 2 lc rgb 'red' title 'Правая разность', \
     'detailed_derivatives.txt' index 9 with lines lw 2 lc rgb 'blue' title 'Левая разность', \
     'detailed_derivatives.txt' index 10 with lines lw 2 lc rgb 'green' title 'Центральная', \
     'detailed_derivatives.txt' index 11 with lines lw 2 lc rgb 'purple' title '5-точ.'
unset multiplot
