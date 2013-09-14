set here = `pwd`
set bp = 5to150_new
set mainpass=/utera/tianye/data_check_amp
#set mainpass=/home/jiayi/Tianshan/data/Data_phamp/
foreach year ( 2007 2008 )
	
	if ( $year == 2007 ) then
	cd ${mainpass}
	\rm *.lst
	foreach month(AUG SEP OCT NOV DEC)
	cd ${mainpass}
	cd ${year}${month}/
	cd ${bp}/COR/
#	pwd > ${year}${month}.lst
	set dir = `pwd`
	echo ${dir}"/" > ${year}${month}.lst
	ls COR*SAC >> ${year}${month}.lst
	mv ${year}${month}.lst ${mainpass} 
	end #foreach month 2007
	
	else if ( $year == 2008 ) then
	foreach month(JAN FEB MAR APR MAY JUN JUL AUG SEP) 
	cd ${mainpass}
	cd ${year}${month}/
	cd ${bp}/COR/
#        pwd > ${year}${month}.lst
	set dir = `pwd`
	echo ${dir}"/" > ${year}${month}.lst
        ls COR*SAC >> ${year}${month}.lst
        mv ${year}${month}.lst ${mainpass}
	
	end #foreach month 2008
	endif
end#foreach yea
cd ${mainpass}
ls -1 20*.lst > all.lst
cd $here
