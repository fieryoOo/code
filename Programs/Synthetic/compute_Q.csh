set model = Q100
set dep = 5
set dir = test_alpha
mkdir -p $dir
cp $model $dir
cd $dir
foreach azi (35.8 90. 123.9 257.8)
   ### produce synthetic sac files at different distances
   set sigdir = SynSig_${model}_${dep}
   mkdir -p $sigdir
   foreach dis ( 100 200 300 400 500 600 800 1000 1200 1400 1600 1800 2000 )
      if( ! -e ${model}_${dep}/${dis}.grn.0 ) /home/tianye/Software/fk/fk.pl -M${model}/${dep}/k -N512/1 $dis
      set oname = ${sigdir}/dis${dis}_azi${azi}.z
      if( ! -e $oname )/home/tianye/Software/fk/syn -M1.0/30/60 -D1 -A${azi} -O${oname} -G${model}_${dep}/${dis}.grn.0
   end #dis

### read amplitudes from sac files at each period and fit for alpha
   set outf = per_alpha_azi${azi}.txt
   set tmpf = amp.tmp
   rm -f $outf
   foreach per ( 4. 6. 9. 12. 15. 20. 25. 30. 38. 46. 54. 65. 75. 85. 95.)
      set fl = `echo $per | awk '{print 0.8/$1}'`
      set fh = `echo $per | awk '{print 1.2/$1}'`
      rm -f $tmpf
      foreach dis ( 100 200 300 400 500 600 800 1000 1200 1400 1600 1800 2000 )
         if( `echo $per $dis | awk '{if($1*10>$2){print 1} else {print 0}}'` ) continue
         set sacf = ${sigdir}/dis${dis}_azi${azi}.z
sac << END
r $sacf
mul 1.0e25
bp c $fl $fh n 4 p 2
envelope
w temp1.sac
quit
END
         set amp = `saclst DEPMAX f temp1.sac | awk '{print $2}'`
         if( `echo $dis | awk '{if($1<1000.){print 1} else {print 0}}'` ) continue
         echo $dis $amp | awk '{print $1,log($2)+0.5*log($1)}' >> $tmpf
      end # dis
      ~/code/Programs/FIT/least_squares_line $tmpf 0 | awk -v per=$per 'BEGIN{pi=3.14159265359}{print per,-$1,-pi/per/3./$1}' >> $outf
   end #per
end #azi
rm -f $tmpf temp1.sac
