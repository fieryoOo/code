#!/bin/bash

if [ $# != 1 ] && [ $# != 2 ]; then
   echo "Usage: "$0" [lib name] [-llib]"
   exit
fi

libname=`echo $1 | awk -F/ '{print $1}'`
if [ ! -e $libname ]; then
   echo "Directory "$libname" not found"
   exit
fi


cd $libname
fsrc=${libname}.cpp
if [ ! -e $fsrc ]; then
   echo "source file "$fsrc" not found"
   exit
fi

flibsrc=lib${libname}.so.1.0
if [ $# == 2 ]; then
	extralib=$2
else
	extralib=""
fi
rm -f $flibsrc ${libname}.o
gcc -std=c++11 -O3 -fPIC -fopenmp -c $fsrc
gcc -shared -Wl,-soname,lib${libname}.so.1 $extralib -o $flibsrc ${libname}.o

cd ../../lib
flib=lib${libname}.so
flib1=lib${libname}.so.1
flibsrc=../sources/${libname}/lib${libname}.so.1.0
rm -f $flib $flib1
ln -sf ${flibsrc} $flib
ln -sf ${flibsrc} $flib1
