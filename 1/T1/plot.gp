set term pdf
set out 'c.pdf'
set xlabel 'Xrot'
set ylabel 'Yrot'
set xrange[-1100:1100]
set yrange[-1100:1100]
set size ratio -1
set title 'Trayectoria de planetas'
plot "c.txt" u 1:2 w l t "Sol", "c.txt" u 3:4 w l t "Jupiter", "c.txt" u 5:6 w l t "Troyanos"
