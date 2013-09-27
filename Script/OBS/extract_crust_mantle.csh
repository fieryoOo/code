foreach file ( MC.1.0.5.mod MC.1.1.0.mod MC.1.1.5.mod MC.1.2.0.mod MC.1.2.5.mod MC.1.3.0.mod )
   set age = `echo $file | cut -d. -f3,4`
   awk 'BEGIN{depcur=-1; depbase=0; nl=0;}{if($1-depcur==0){nl++; depbase=$1}; depcur=$1;if(nl==2){print $1-depbase,$2,600}}' $file > $file'_crust_Qmiu'
   awk 'BEGIN{depcur=-1; depbase=0; nl=0;}{if($1-depcur==0){nl++; depbase=$1}; depcur=$1;if(nl==3){print $1-depbase,$2}}' $file > $file'_mantle'
end #file
