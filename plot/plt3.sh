


    #!/bin/bash
    #setting
    # to run this file: type <gnuplot "filename" >

    set terminal png size 1200,800
    set output  "LFthree.png"
    set title "HMC: Leap Frog,L=10,LF=50"
    set autoscale
    set xlabel "Frequency"
    set ylabel "Action"
    set grid
    set nokey

    plot "fort.200" using 1:2 with linespoints,\
         "fort.200" using 1:3 with linespoints,\
         "fort.200" using 1:4 with linespoints,\
         "fort.200" using 1:5 with linespoints



   #pause -1 "Hit any key to continue"


