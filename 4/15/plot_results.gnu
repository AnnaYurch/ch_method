set terminal pngcairo size 1200,800 enhanced font 'Verdana,12'
set output 'ode_solution.png'
set title 'Решение ОДУ: y'' = y^2 - y*sin(x) + cos(x)' font 'Verdana,14'
set xlabel 'x' font 'Verdana,12'
set ylabel 'y(x)' font 'Verdana,12'
set xrange [0:4]
set yrange [-1.2:1.2]
set grid
set key outside right top vertical box
set key spacing 1.5

plot 'solutions.dat' using 1:2 with lines linewidth 3 linecolor rgb 'black' title 'Точное решение (sin(x))', \
     '' using 1:3 with points pointtype 7 pointsize 1.2 linecolor rgb 'red' title 'Явный Эйлер', \
     '' using 1:4 with points pointtype 9 pointsize 1.2 linecolor rgb 'blue' title 'Неявный Эйлер', \
     '' using 1:5 with points pointtype 11 pointsize 1.2 linecolor rgb 'green' title 'Хойн', \
     '' using 1:6 with points pointtype 13 pointsize 1.2 linecolor rgb 'purple' title 'Эйлера-Коши', \
     '' using 1:7 with points pointtype 15 pointsize 1.2 linecolor rgb 'orange' title 'Рунге-Кутта 4', \
     '' using 1:8 with points pointtype 17 pointsize 1.2 linecolor rgb 'brown' title 'Трапеций'
