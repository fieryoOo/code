#set srcdir = Programs_furtado/ # trailing / is important!
#set desdir = Programs/ # trailing / is important!
set srcdir = Script_furtado/
set desdir = Script/
#rsync -abvizP --suffix=_tobechecked $srcdir $desdir
# -u to ignore older version of files from srcdir

### handle duplicated files ###
set fcheck = merge_check.txt
rm -f $fcheck
foreach ofile ( `find $desdir -type f -name "*_tobechecked"` )
set nfile = `echo $ofile | sed s/'_tobechecked'/''/`
set Ndiff = `sdiff -s -w 148 $nfile $ofile | awk 'substr($1,0,1)!="#"' | wc -l`
if( $Ndiff == 0 )then
   rm -f $ofile
   continue
endif
echo '########################## '$nfile' # '$ofile' ##########################' >> $fcheck
sdiff -s -w 148 $nfile $ofile | awk 'substr($1,0,1)!="#"' >> $fcheck
end
