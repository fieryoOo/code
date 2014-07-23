#!/bin/bash

IFS=$'\n'
for file in `ls Listing\ *.cpp`; do
	newname=`echo $file | awk '{print $1"_"$2}'`
	echo $file "->" $newname
	mv $file $newname
done
