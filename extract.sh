#!/bin/bash

filename=${1:-"tests/4_226.fts"}
filename=$(readlink -f "$filename")

if [ $(xpaaccess ds9) == 'no' ]; then
    ds9 &
    sleep 2
fi

# Display FITS
xpaset -p ds9 file fits '"'$filename'"'

if [ ! -f $filename.psf ]; then
    # Create the PSF and store it to filename.psf
    ./extract -psfex $filename
fi

# Find inital peaks
./extract -simple $filename $filename.txt

cat $filename.txt | awk '{print "point("($1+1)", "($2+1)") # point=x"}' |xpaset ds9 regions
