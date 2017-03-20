cmake_minimum_required(VERSION 3.6)
project(miajiang_jphender_snehvish_a1)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
    part1/DTwoDimArray.h
    part1/fft.h
    part1/Makefile
    part1/noise1.png
    part1/SImage.h
    part1/SImageIO.h
    part1/watermark.cpp
    part2/detect.cpp
    part2/DTwoDimArray.h
    part2/Informatics-2.png
    part2/Informatics.png
    part2/Makefile
    part2/Plaza.png
    part2/SImage.h
    part2/SImageIO.h
    part2/SRSC-2.png
    part2/SRSC.png
    README.md)

add_executable(miajiang_jphender_snehvish_a1 ${SOURCE_FILES})