#!/bin/bash

make clean
make

# WATERMARK CHECK - For all pictures
for n in pictures/*
do
    if [[ ! -d "$n" ]]
	then continue
    fi
    s=$( echo $n | cut -d '/' -f 2)
    for ns in $n/watermark_*
    do
	t=$( echo $ns | cut -d '_' -f 3)
	val=$( echo $t | cut -d '.' -f 1)
	
	RANDOM=$$
	c=0
	echo "Checking "$( echo $ns | cut -d '/' -f 3)
	for i in {0..100}
    	do
	    r=$RANDOM
	    if [ "$r" -eq "$i" ]
	        then c=1
	    fi
	    ./watermark 1.3 $ns "checked_"$s"_"$r".png" check $r
        done
	if [ "$c" -eq 0 ]
	    then ./watermark 1.3 $ns "checked_"$s"_"$val".png" check $val
	fi
	echo
	echo
    done
done
