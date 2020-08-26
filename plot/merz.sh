

    #setting
    # to run this file: type <gnuplot "filename" >

    set terminal png size 1200,800
    set output  "mass.png"
    set title "HMC"

    set xlabel "m_q"
    set ylabel "m_20,m_31"
    set grid

  
    plot "fort.20" using 1:2:3 with errorbars,\
         "fort.20" using 1:2 with linespoints,\
         "fort.31" using 1:2:3 with errorbars,\
         "fort.31" using 1:2 with linespoints,

 




