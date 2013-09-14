#foreach month ( 2000.JAN 2000.FEB 2000.MAR 2000.APR 2000.MAY 2000.JUN )
foreach month ( 2000.JAN 2000.FEB 2000.MAR 2000.APR )
   echo $month
   cd $month
   foreach day ( `ls -d ${month}.*` )
      echo $day
      cd $day
      foreach sacf ( `ls ${day}.*.SAC` )
         set outf = 'ft_'$sacf
         set stlo = `echo $sacf | cut -d. -f4 | sed s/'STA'/''/ | awk '{print $1*2}'`
sac << END
r $sacf
bp c 0.01 0.3 n 4 p 2
w temp1.sac
cut 1000 85000
r temp1.sac
ch STLA 0
ch STLO $stlo
w $outf
quit
END
      end
      rm -f temp1.sac
      cd ..
   end
   cd ..
end
