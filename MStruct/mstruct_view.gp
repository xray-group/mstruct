# very simple script to plot data
# 1 ... 2Theta (deg)
# 2 ... Iobs
# 3 ... Ical
# 3 ... Iobs - Icalc

print "plotting file: $0.dat"

set autoscale xy

# Sqrt scale ?
psqrt = ("$1" eq "sqrt" || "$1" eq "Sqrt" || "$1" eq "SQRT")

#set style line 1 lt rgb "cyan" lw 3 pt 6

# linear plot
if (!psqrt) plot "$0.dat" u 1:2 title 'measured' w p lt rgb "red", "$0.dat" u 1:3 title 'calculated' w l lt rgb "blue", "$0.dat" u 1:4 title 'difference' w l lt rgb "green';

# sqrt plot
if (psqrt) plot "$0.dat" u ($$1):(sqrt($$2)) title 'measured' w p lt rgb "red", "$0.dat" u ($$1):(sqrt($$3)) title 'calculated' w l lt rgb "blue", "$0.dat" u ($$1):(sgn($$4)*sqrt(abs($$4))) title 'difference' w l lt rgb "green";

set xlabel "2Theta (deg)"

# linear plot
if (!psqrt) set ylabel "Intensity (arbitrary units)";

# sqrt plot
if (psqrt) set ylabel "Square Root of Intensity (arbitrary units)";

replot

# to zoom use mouse right button
# to un-zoom type: set autoscale xy; replot

