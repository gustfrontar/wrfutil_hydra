set terminal png size 1280,640
set output "Energy.png"

set xlabel "Step"
set ylabel "Energy[J]"

set multiplot layout 3,2 title "Energy conservation check"

plot [] [-5e14:5e13] "energy.dat" using 1:2       ti "Total Energy(ENGT)"     with linespoints pointsize 0.2
plot [] [-6e14:5e13] "energy.dat" using 1:($3+$5) ti "ENGP+ENGI"              with linespoints pointsize 0.2
plot [] [-1e13:2e14] "energy.dat" using 1:3       ti "Potential Energy(ENGP)" with linespoints pointsize 0.2
plot [] [0e0:8e12] "energy.dat" using 1:4       ti "Kinetic Energy(ENGK)"   with linespoints pointsize 0.2
plot [] [-8e14:1.5e14] "energy.dat" using 1:5       ti "Internal Energy(ENGI)"  with linespoints pointsize 0.2

unset multiplot