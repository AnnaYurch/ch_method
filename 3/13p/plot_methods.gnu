set terminal pngcairo size 1600,800 enhanced font 'Arial,10'
set output 'methods_comparison.png'
set multiplot layout 1,2 title 'Сравнение методов численного дифференцирования (h=0.150)' font 'Arial,14'
set title 'Методы первой производной'
set xlabel 'x'
set ylabel 'f''(x)'
set grid
set key outside
plot 'methods_data.txt' index 0 using 1:2 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     '' index 0 using 1:3 with points pt 1 title 'Правая разность', \
     '' index 0 using 1:4 with points pt 2 title 'Левая разность', \
     '' index 0 using 1:5 with points pt 3 title 'Центральная', \
     '' index 0 using 1:6 with points pt 4 title '3-точ. вперед', \
     '' index 0 using 1:7 with points pt 5 title '3-точ. назад', \
     '' index 0 using 1:8 with points pt 6 title '4-точ.'
set title 'Методы второй производной'
set xlabel 'x'
set ylabel 'f''''(x)'
set grid
set key outside
plot 'methods_data.txt' index 1 using 1:2 with lines lw 3 lc rgb 'black' title 'Аналитическая', \
     '' index 1 using 1:3 with points pt 1 title 'Правая разность', \
     '' index 1 using 1:4 with points pt 2 title 'Левая разность', \
     '' index 1 using 1:5 with points pt 3 title 'Центральная', \
     '' index 1 using 1:6 with points pt 4 title '5-точ.'
unset multiplot
