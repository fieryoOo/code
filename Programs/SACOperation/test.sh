#!/bin/bash

for sacf in `ls sacs/ft*SAC`; do
   ./chhdr $sacf KEVNM WAHAHA
#${SACHOME}/bin/sac /home/tianye/usr/macros/sacinit.m << END
#RH $sacf
#CH KEVNM WAHAHA
#WH OVER
#quit
#END
done
