#foreach per ( 2. 4. 6. 9. 12. 15. 20. 25. 30. 38. 46. 54. 65. 75. 85. 95.)
foreach per ( 15 25 50 )
   mkdir -p $per
   set fl = `echo $per | awk '{print 0.8/$1}'`
   set fh = `echo $per | awk '{print 1.2/$1}'`
   foreach sta ( STA0 STA1 STA2 STA3 STA4 STA5 )
      set outf = ${per}/${sta}_amp.txt
      rm -f $outf
      foreach file ( `ls ${sta}/COR_${sta}_STA?.SAC` )
sac << END
r $file
bp c $fl $fh n 4 p 2
envelope
w temp1.sac
cut 0 e
r temp1.sac
w over
quit
END
         set dis = `saclst dist f temp1.sac | awk '{print $2}'`
         set amp = `saclst depmax f temp1.sac | awk '{print $2}'`
         echo $dis $amp | awk '{print $1,$2,log($2)+0.5*log($1)}' >> $outf
      end
   end
end
