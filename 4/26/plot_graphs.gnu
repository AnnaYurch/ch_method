# График решения
set terminal pngcairo size 800,600
set output 'solution.png'
set title 'Решение уравнения'
set xlabel 'x'
set ylabel 'y(x)'
plot 'solution_plot.dat' using 1:2 with linespoints title 'Численное', \
     'solution_plot.dat' using 1:3 with lines title 'Аналитическое'

# График погрешности  
set output 'error.png'
set title 'Погрешность решения'
set ylabel 'Погрешность'
plot 'error_plot.dat' using 1:2 with linespoints title 'Ошибка'

# График метода стрельбы
set output 'shooting.png'
set title 'Метод стрельбы'
set xlabel 'α'
set ylabel 'y(4.5)'
plot 'shooting_plot.dat' using 1:2 with linespoints title 'y(4.5;α)', \
     'shooting_plot.dat' using 1:3 with lines title 'Целевое значение'

# Все графики вместе
set terminal pngcairo size 1200,900
set output 'all_graphs.png'
set multiplot layout 2,2
set title 'Решение краевой задачи методом стрельбы'

set title 'Решение'
plot 'solution_plot.dat' using 1:2 with linespoints title 'Численное', \
     'solution_plot.dat' using 1:3 with lines title 'Аналитическое'

set title 'Погрешность'
plot 'error_plot.dat' using 1:2 with linespoints title 'Ошибка'

set title 'Метод стрельбы'
plot 'shooting_plot.dat' using 1:2 with linespoints title 'y(4.5;α)', \
     'shooting_plot.dat' using 1:3 with lines title 'Цель'

set title 'Производная'
plot 'solution_plot.dat' using 1:4 with lines title "y'(x)"

unset multiplot