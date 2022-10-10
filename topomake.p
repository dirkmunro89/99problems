set term eps crop
set output 'topology.eps'
set view map
set size ratio -1

unset key
set palette rgbformulae -3, -3, -3

x1=-0.49
x2=xi - 0.55
y1=-0.49
y2=yi - 0.50

set xrange [x1:x2] noreverse nowriteback
set yrange [y1:y2] noreverse nowriteback
unset xtics
unset ytics
unset colorbox

set offset 0,0

splot 'topology.dat' matrix with image
unset multiplot
