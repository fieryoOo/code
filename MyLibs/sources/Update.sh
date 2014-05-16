#!/bin/bash

if [ $# != 1 ]; then
   echo "Usage: "$0" [lib name]"
   exit
fi

libname=$1
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
rm -f $flibsrc ${libname}.o
gcc -std=c++0x -O3 -fPIC -c $fsrc
gcc -shared -Wl,-soname,lib${libname}.so.1 -o $flibsrc ${libname}.o

cd ../../lib
flib=lib${libname}.so
flib1=lib${libname}.so.1
flibsrc=../sources/${libname}/lib${libname}.so.1.0
rm -f $flib $flib1
ln -sf ${flibsrc} $flib
ln -sf ${flibsrc} $flib1
