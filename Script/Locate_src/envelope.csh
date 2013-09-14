#set dir='STACK_OCT2DEC'
#cp station.lst $dir
#cp event.dat $dir
#cd $dir
#foreach file (`awk '$2<30 {print $1}' file_daynum.lst`)
#rm -f $file
#end
#awk '$2>=30' file_daynum.lst > file_daynum_30.lst

foreach file ( `awk 'NR>22100 {print $1"@"$2}' file_daynum.lst` )
set f1=`echo $file | cut -d@ -f1`
if( ! -e $f1)continue
set f2='/home/tianye/data_locate_src/'$f1'.env'
set num=`echo $file | cut -d@ -f2`
echo 'Working on file '$f1'...'
sac << END
r $f1
bp c 0.07 0.09 n 4 p 2
envelope
mul 1e6
div $num
w $f2
quit
END
end
