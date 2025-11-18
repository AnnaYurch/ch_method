set terminal pngcairo size 1200,800 enhanced font 'Arial,10'
set multiplot layout 2,1 title 'Решение ОДУ: y'' = y^2 - y*sin(x) + cos(x)' font 'Arial,12'

set title 'Численные решения и точное решение'
set xlabel 'x'
set ylabel 'y(x)'
set grid
set key outside right top
plot 'solutions.dat' using 1:2 with lines linewidth 2 title 'Точное решение (sin(x))', '' using 1:3 with points pointtype 1 linecolor rgb 'red' title 'Явный Эйлер', '' using 1:4 with points pointtype 2 linecolor rgb 'blue' title 'Неявный Эйлер', '' using 1:5 with points pointtype 3 linecolor rgb 'green' title 'Хойн', '' using 1:6 with points pointtype 4 linecolor rgb 'purple' title 'Эйлера-Коши', '' using 1:7 with points pointtype 5 linecolor rgb 'orange' title 'Рунге-Кутта 4', '' using 1:8 with points pointtype 6 linecolor rgb 'brown' title 'Трапеций'

set title 'Погрешности численных методов'
set xlabel 'x'
set ylabel 'Погрешность |y_{числ} - y_{точн}|'
set grid
set logscale y
set key outside right top
plot 'errors.dat' using 1:2 with lines linewidth 2 title 'Ошибка Явный Эйлер', '' using 1:3 with lines linewidth 2 linecolor rgb 'blue' title 'Ошибка Неявный Эйлер', '' using 1:4 with lines linewidth 2 linecolor rgb 'green' title 'Ошибка Хойн', '' using 1:5 with lines linewidth 2 linecolor rgb 'purple' title 'Ошибка Эйлера-Коши', '' using 1:6 with lines linewidth 2 linecolor rgb 'orange' title 'Ошибка Рунге-Кутта 4', '' using 1:7 with lines linewidth 2 linecolor rgb 'brown' title 'Ошибка Трапеций'

unset multiplot
set output 'results_comparison.png'
replot
