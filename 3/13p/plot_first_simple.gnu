set terminal pngcairo size 800,600
set output 'first_derivative_simple.png'
set title 'First Derivative Methods (h=0.150)'
set xlabel 'x'
set ylabel 'f'(x)'
set grid
set key top left
plot 'first_deriv_simple.txt' using 1:2 with lines lw 3 title 'Analytic', \
     '' using 1:3 with lines lw 2 title 'Right Difference', \
     '' using 1:4 with lines lw 2 title 'Left Difference', \
     '' using 1:5 with lines lw 2 title 'Central Difference'
