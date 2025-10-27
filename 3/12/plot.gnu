set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'mnk_graph.png'
set title 'Метод Наименьших Квадратов - Вариант 44' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'y' font 'Arial,12'
set grid
set key top left box
plot 'data.txt' index 0 with points pt 7 ps 1.5 lc rgb 'black' title 'Исходные данные', \
     'data.txt' index 1 with lines lw 2 lc rgb 'red' title 'F1(x)', \
     'data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'F2(x)', \
     'data.txt' index 3 with lines lw 2 lc rgb 'green' title 'F3(x)'
