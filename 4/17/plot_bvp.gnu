set terminal pngcairo size 1200,800 enhanced font 'Arial,12'
set output 'bvp_solution.png'
set title 'Решение краевой задачи №44'
set xlabel 'x'
set ylabel 'y(x)'
set grid
set key top left box
set xrange [-1:1]

plot 'bvp_solutions.dat' using 1:2 with lines linewidth 3 title 'Аналитическое решение', \
     '' using 1:3 with points pointtype 7 pointsize 1.5 title 'Численное (1 порядок)', \
     '' using 1:4 with points pointtype 9 pointsize 1.5 title 'Численное (2 порядок)'
