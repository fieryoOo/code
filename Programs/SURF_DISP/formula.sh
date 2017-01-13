#!/bin/bash

fin=
fou=
awk '{if($2==0){print $1,1.45,0,0.77,$3}else{if($3==80){ratio=2.}else{ratio=1.73} print $1,$2*ratio,$2,$2*ratio*0.32 + 0.77,$3}}' $fin > $fout
