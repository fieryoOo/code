#!/bin/csh
if ($#argv != 2)then
echo "Usage: inter_sta_dis.csh [sta_lst_one] [sta_lst_two]"
exit
endif

rm -f inter_sta_dis.txt
foreach sta1 (`more $argv[1] | sed s/' '/'@'/g`)
set sta_name1=`echo $sta1 | cut -d@ -f1`
set long1=`echo $sta1 | cut -d@ -f2`
set lati1=`echo $sta1 | cut -d@ -f3`

foreach sta2 (`more $argv[2] | sed s/' '/'@'/g`)
set sta_name2=`echo $sta2 | cut -d@ -f1`
set long2=`echo $sta2 | cut -d@ -f2`
set lati2=`echo $sta2 | cut -d@ -f3`
set dis=`echo $long1 $lati1 $long2 $lati2 | awk '{printf "%.10f\n",sqrt(($1-$3)^2+($2-$4)^2)}'`
echo $sta_name1" ("$long1 $lati1") " $sta_name2" ("$long2 $lati2") dis: "$dis >> inter_sta_dis.txt
end
end
