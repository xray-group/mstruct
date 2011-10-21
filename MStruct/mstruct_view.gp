# very simple script to plot data
# 1 ... 2Theta (deg)
# 2 ... Iobs
# 3 ... Ical
# 3 ... Iobs - Icalc

print "plotting file: $0.dat"

set autoscale xy
plot "$0.dat" u 1:2 w p, "$0.dat" u 1:3 w l, "$0.dat" u 1:4 w l

set xlabel "2Theta (deg)"
set ylabel "Intensity (arbitrary units)"

# to zoom use mouse right button
# to un-zoom type: set autoscale xy; replot

