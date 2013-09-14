
foreach per (8 10 12 14)
rm -f 'grv_map_US_global'$per'sec.txt.HD'
set lati=-89
while($lati<90)
echo $lati
set long=0

while($long<360)
if( $long >= 200 && $long <= 300 && $lati <= 80 && $lati >= 0 )then
awk -v long=$long -v lati=$lati '$1==long && $2==lati' 'grv_map_US_'$per'sec.txt.HD' >> 'grv_map_US_global'$per'sec.txt.HD'
else
echo $long'\t'$lati'\t0' >> 'grv_map_US_global'$per'sec.txt.HD'
endif
@ long += 1
end
 @ lati += 1
end

end
