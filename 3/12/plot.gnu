set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'mnk_graph.png'
set title 'Метод Наименьших Квадратов - Вариант 44' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'y' font 'Arial,12'
set grid
set key top left box
plot 'data.txt' using 1:2 with points pt 7 ps 1.5 lc rgb 'black' title 'Исходные данные', \
     'data.txt' using 1:3 with lines lw 2 lc rgb 'red' title 'F1(x) = -1.148 + 0.422 x', \
     'data.txt' using 1:4 with lines lw 2 lc rgb 'blue' title 'F2(x) = -1.288 + 0.835 x + -0.078 x^2', \
     'data.txt' using 1:5 with lines lw 2 lc rgb 'green' title 'F3(x) = -1.727 + -0.049 x + 0.519 x^2 + -0.076 x^3'
