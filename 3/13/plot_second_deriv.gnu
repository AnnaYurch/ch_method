set terminal png size 1000,600
set output 'second_derivative_plot.png'
set title 'Second Derivative Comparison (h=0.075)'
set xlabel 'x'
set ylabel 'f''(x)'
set grid
set key top left
set xrange [0.25:2.85]
plot 'second_deriv_data.txt' using 1:2 with lines lw 3 title 'Analytic', \
     '' using 1:3 with lines lw 2 title 'Right Difference', \
     '' using 1:4 with lines lw 2 title 'Left Difference', \
     '' using 1:5 with lines lw 2 title 'Central Difference', \
     '' using 1:6 with lines lw 2 title 'Five Point'
