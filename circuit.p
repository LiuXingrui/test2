# gnuplot column and row index starts with 1, but the element numbering in an array starts with 0
set format x "10^{%L}"
set format y "10^{%L}"
set term pngcairo dashed

set grid
 
 #    set   autoscale                        # scale axes automatically
           set log x
  set log y                         # remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      
      set xlabel "Gate Depolarizing Rate"
      set ylabel "Logical Error Rate"

# Enable second y-axis (right-hand side)

     set key outside bottom
#  set label "Yield Point" at 0.003,260
#     set arrow from 0.0028,250 to 0.003,280
#    set xr [0.0:0.022]
#     set yr [0:325]
#    set origin 0.0,0.5


set title "GB code circuit-level decoding BP+OSD  "
set output "circuit_GB.png"
    
	plot   "CBP_GB_w4" using 2:5 title "GB n=22 wt=4 d=4 BP" with linespoints lw 3 lc rgb "purple" ,\
	"CBP_GB_w4" using 2:4 title "GB n=22 wt=4 d=4 BP+OSD" with linespoints lw 3 dashtype 2 lc rgb "purple" ,\
	 "CBP_GB_w6" using 2:5 title "GB n=22 wt=6 d=5 BP" with linespoints lw 3 lc rgb "red" ,\
	"CBP_GB_w6" using 2:4 title "GB n=22 wt=6 d=5 BP+OSD" with linespoints lw 3 dashtype 2 lc rgb "red" ,\
	 "CBP_GB_w8" using 2:5 title "GB n=22 wt=8 d=6 BP" with linespoints lw 3 lc rgb "blue" ,\
	"CBP_GB_w8" using 2:4 title "GB n=22 wt=8 d=6 BP+OSD" with linespoints lw 3 dashtype 2 lc rgb "blue" ,\
	 x  title "y=x"  lc rgb "yellow",\
	 x*x title "y=x^2" lc rgb "black",\
	 x*x*x title "y=x^3" lc rgb "orange"
	 


	  



	  

