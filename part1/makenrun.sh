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

# CLEANING - Just for noise1.png
./watermark 1.2 pictures/noise1/noise1.png cleaned_noise1.png

# WATERMARK ADD - For all pictures
RANDOM=$$
for n in pictures/*
do
    if [[ ! -d "$n" ]]
	then continue
    fi
    s=$( echo $n | cut -d '/' -f 2)
    for ns in $n/watermark_*
    do
	rm $ns
    done
    ns=$n/$s.png
    r=$RANDOM
    ./watermark 1.3 $ns "watermark_"$s"_"$r".png" add $r
done

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
	for i in {0..10}
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
