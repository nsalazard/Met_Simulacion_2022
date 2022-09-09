set term pdf
set out 'plot.pdf'
set xlabel 'Xrot'
set ylabel 'Yrot'
set title 'Tiempo de CPU vs Blocksize'
plot "c.txt" u 1:2 w lp t "Sol", "c.txt" u 3:4 w lp t "Jupiter", "c.txt" u 5:6 w lp t "Troyanos"