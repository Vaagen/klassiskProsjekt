gnuplot << EOF

a = 1.1 # range of axis

set grid
show grid
set size square
set xrange [-a:a]
set yrange [-a:a]
set xlabel 'Distance relative to starting radius of particle 1'
#set ylabel 'Distance relative to starting radius of particle 1'
plot '+' using (0):(0) title 'Central mass', 'tertiarypos1.txt' title 'Trajectory of particle 1' with line , 'tertiarypos2.txt' title 'Trajectory of particle 2' with line
set term aqua 1
plot '+' using (0):(0) title 'Central mass' , 'tertiarypos1.txt' title 'Trajectory of particle 1' with line , 'tertiarypos1EXACT.txt' title 'Trajectory of particle 1 calculated with exact solution'
set term aqua 2
set xlabel 'Time'
set ylabel 'Energy'
set autoscale
set xrange [0:12]
set yrange [-0.0227792:-0.0227778]
plot 'tertiaryenergy.txt' title 'Total energy of system' with line
set term aqua 3
set ylabel 'Angular momentum'
plot 'tertiaryangularmomentum.txt' title 'Total angular momentum of system' with line
EOF