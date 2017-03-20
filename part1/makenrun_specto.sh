#!/bin/bash

make clean
make

# SPECTOGRAM - For all pictures
for n in pictures/*
do
    if [[ ! -d "$n" ]]
	then continue
    fi
    s=$( echo $n | cut -d '/' -f 2)
    ns=$n/$s.png   
    ./watermark 1.1 $ns spectogram_$s.png
done
