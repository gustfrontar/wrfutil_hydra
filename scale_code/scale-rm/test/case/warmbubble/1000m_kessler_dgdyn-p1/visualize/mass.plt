set terminal png size 1280,640
set output "Mass.png"

set xlabel "Step"
set ylabel "Mass[kg]"

set multiplot layout 2,2 title "Mass conservation check"

plot [] [-4e0:0] "mass.dat" using 1:2 ti "Dry   Mass(QDRY)" with linespoints pointsize 0.2
plot [] [-4e8:5e7] "mass.dat" using 1:3 ti "Water Mass(QTOT)" with linespoints pointsize 0.2
plot [] [-1:1] "mass.dat" using 1:4 ti "Evaporation"      with linespoints pointsize 0.2
plot [] [0:4e8] "mass.dat" using 1:5 ti "Precipitation"    with linespoints pointsize 0.2

unset multiplot
