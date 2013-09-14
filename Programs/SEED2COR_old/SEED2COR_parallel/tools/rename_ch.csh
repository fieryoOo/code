#!/bin/csh

if ($#argv != 3) then
  echo "USAGE: "$0" [dir] [old ch name] [new ch name]"
  exit 1
endif

set chold = $argv[2]
set chnew = $argv[3]
set month = `echo $argv[1] | awk -F/ '{print $NF}'`
cd $argv[1]
echo "Started renaming for month "$month"..."
foreach day (`ls -d $month'.'*`)
   cd $day
   echo $day
   foreach file (`ls *${day}*'.'${chold}'.SAC'*`)
      set nfile = `echo $file | sed s/$chold/$chnew/`
      mv $file $nfile
   end #file
   cd ..
end #day
