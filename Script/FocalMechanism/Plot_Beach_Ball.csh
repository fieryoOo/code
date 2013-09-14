psmeca -H1 -P -R-120/-100/30/50 -JM15c -Sm5.0c -U -K << eof > plotfile.ps
longitude latitude depth f1 f2 f3 f4 f5 f6 ex newlon newlat title
-115 40 16.2 -8.19 1.46 6.73 0.610 1.31 3.14 24 -175 -15 Wells_true
eof

psmeca -H1 -P -R-120/-100/30/50 -JM15c -Sm5.0c -U -O << eof >> plotfile.ps
longitude latitude depth f1 f2 f3 f4 f5 f6 ex newlon newlat title
-105 40 16.2 -6.86 1.46 5.07 0.0 0.937 1.99 24 -175 -15 Wells_recovered
eof

