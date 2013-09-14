#!/bin/csh

if ($#argv != 6) then
  echo "USAGE: "$0" [dir_to_be_moved] [dir_target] [merge orig_SACs?] [merge ft_SACs?] [merge am&phs?] [merge COR directories?]"
  exit 1
endif

set dirf = $argv[1]
set dirt = $argv[2]
set do_orSAC = $argv[3]
set do_ftSAC = $argv[4]
set do_amphs = $argv[5]
set do_CORRs = $argv[6]

# check/compare month names
if( ! -e $dirf || ! -e $dirt ) then
echo "origin or/and target directory not accessable!"
exit 0
endif
set month = `echo $dirf | awk -F/ '{print $NF}'`
set month2 = `echo $dirt | awk -F/ '{print $NF'}`
if( $month != $month2 ) then
echo "No, you shouldn't  merge 2 different months: "$month" -> "$month2"!!!"
exit 0
endif

# prompt to continue?
echo "Warning: files with same file_name will be overwritten\!\! Ready to continue? "
set flag = $<
set flag = `echo $flag | head -c1`
if( $flag != 'y' && $flag != 'Y' ) exit 0
echo "Started merging for month "$month"..."

# do one day at a time
echo "Moving sigle-station files for each day..."
cd $dirf
foreach day (`ls -d $month'.'*`)
if( `ls $day | wc -l` == 0 ) then
rm -rf $day
continue
else if( ! -e $dirt'/'$day ) then
mv $day $dirt'/'
continue
endif
cd $day
if( $do_orSAC ) then
( mv $day*'SAC' $dirt'/'$day ) >& /dev/null
endif
if( $do_ftSAC ) then
( mv 'ft_'$day*'SAC' $dirt'/'$day ) >& /dev/null
( mv 'ft_'$day*'SAC_rec'* $dirt'/'$day ) >& /dev/null
endif
if( $do_amphs ) then
( mv 'ft_'$day*'SAC.am' $dirt'/'$day ) >& /dev/null
( mv 'ft_'$day*'SAC.ph' $dirt'/'$day ) >& /dev/null
endif
cd ..
if( `ls $day | wc -l` == 0 ) rm -rf $day
end #day

# work on the COR directory
   # Cor_dayflag.lst
if( $do_CORRs ) then
if( ! -e COR ) then
echo "No COR directory! Skipped!"
else
echo "Moving CCs and merging Cor_dayflag.lst ..."
cd COR
foreach CC ( `awk '{print $1"@"$2}' $dirt'/COR/Cor_dayflag.lst'` )
set file = `echo $CC | cut -d@ -f1`
set rec = `echo $CC | cut -d@ -f2`
if( `awk -v file=$file '$1==file' Cor_dayflag.lst | wc -l` > 0 ) continue
echo $file"   "$rec >> Cor_dayflag.lst
end
mv Cor_dayflag.lst $dirt'/COR/Cor_dayflag.lst'
   # Cor_dayflag.lst_old
cat Cor_dayflag.lst_old >> $dirt'/COR/Cor_dayflag.lst_old'
rm -f Cor_dayflag.lst_old
   # each station directory
foreach dir (`ls -d */`)
if( ! -e $dirt'/COR/'$dir ) then
mv $dir $dirt'/COR/'
else
( mv $dir'/'* $dirt'/COR/'$dir'/' ) >& /dev/null
if( `ls $dir | wc -l` == 0 ) rm -rf $dir
endif
end #dir
   # check empty of COR directory
cd ..
if( `ls COR | wc -l` == 0 ) rm -rf COR
endif # do_CORRs

# all done! check empty of the top directory
cd ..
if( `ls $month | wc -l` == 0 ) then
rm -rf $month
echo "Merges all done!"
else
echo "Warning: originated directory isn't empty after merging!"
endif
