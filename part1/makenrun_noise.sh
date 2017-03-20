#!/bin/bash

make clean
make

# CLEANING - Just for noise1.png
./watermark 1.2 pictures/noise1/noise1.png cleaned_noise1.png
