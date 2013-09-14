#!/bin/csh
#echo on
if ($#argv != 2) then
  echo "USAGE: correct_2pi_curvature [region_infile] [sac_path]"
  exit 1
endif

cd $argv[2]

set REG = `more $argv[1]`
foreach per ( 10 15 20 25 )
#@ dis = 6 * $per
cd Ph_Amp_Map_$per"sec"

foreach sta (`ls *_center_ph_amp_map | cut -d_ -f1`)
#foreach sta (M12A)
/home/tianye/code/Programs/ASN/correct_2pi_v1 $sta"_center_ph_amp_map" $per 50
/home/tianye/code/Script/GMT/C_plot_travel_positive $sta"_center_ph_amp_map_v1" $argv[1] 0.1
mv $sta"_center_ph_amp_map_v1.HD" $sta"_center_ph_map_v1.HD"
/home/tianye/code/Script/GMT/C_plot_ASN_am  $sta"_center_ph_amp_map_v1" $argv[1]
/home/tianye/code/Programs/ASN/correct_curvature_v1 $sta $per 50 0.02 0.05
#end

rm -f `wc $sta"_center_ph_amp_map_v2" | grep "0       0       0" | awk '{print $4}'`
if( ! -e $sta"_center_ph_amp_map_v2" ) continue
/home/tianye/code/Script/GMT/C_plot_travel_positive $sta"_center_ph_amp_map_v2" $argv[1] 0.1
mv $sta"_center_ph_amp_map_v2.HD" $sta"_center_ph_map_v2.HD"
/home/tianye/code/Script/GMT/C_plot_travel_T0.2 $sta"_center_ph_amp_map_v2" $argv[1]
/home/tianye/code/Script/GMT/C_plot_ASN_am $sta"_center_ph_amp_map_v2" $argv[1]

end
cd ..
end
