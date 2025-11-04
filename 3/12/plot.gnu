set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'mnk_graph.png'
set title 'Метод Наименьших Квадратов - Аппроксимация многочленами 1-6 степени' font 'Arial,14'
set xlabel 'x' font 'Arial,12'
set ylabel 'y' font 'Arial,12'
set grid
set key top left box
plot 'data.txt' index 0 with points pt 7 ps 1.5 lc rgb 'black' title 'Исходные данные', \
     'data.txt' index 1 with lines lw 2 lc rgb 'red' title 'F1(x) (1 степень)', \
     'data.txt' index 2 with lines lw 2 lc rgb 'blue' title 'F2(x) (2 степень)', \
     'data.txt' index 3 with lines lw 2 lc rgb 'green' title 'F3(x) (3 степень)', \
     'data.txt' index 4 with lines lw 2 lc rgb 'purple' title 'F4(x) (4 степень)', \
     'data.txt' index 5 with lines lw 2 lc rgb 'orange' title 'F5(x) (5 степень)', \
     'data.txt' index 6 with lines lw 2 lc rgb 'brown' title 'F6(x) (6 степень)'
