# extract

**extract** is source extraction and photometry code for Mini-MegaTORTORA camera

    ./extract -help

Currently it may perform one of the following tasks.

PSF extraction
==============

Spatially variable PSF extraction using **psfex** code

    ./extract -psfex $filename

By default output is stored to *$filename.psf*. Output file may also be specified on command line.

Simple peak extraction
======================

Simple extraction and characterization of peaks in the image.

    ./extract -simple $filename $outname

Ideally, it should extract all the peaks for the following PSF-fitting and properly handle spatially variable background etc.

First commit contains the code to extract all local peaks (after simple smoothing with 1-pixel Gaussian).

Second commit is an attempt to implement more clever uplink-based algorithm which may not only locate local peaks but also to merge and/or deblend spatially close objects. It is unfinished but somehow works.

PSF photometry
==============

PSF photometry

    ./extract $filename $outname

It performs peaks extraction from previous step and then tries to fit the PSF model to every peak.

PSF model will be created from the image, or loaded from *$filename.psf* if exists, or loaded from the PSF file spectfied on command line.

The code fits PSF models for spatially close sources simultaneosly, and for spatially distant - sequentially. It performs fit quality checks and rejects bad fits. The background is fitted independently for every object and therefore may not be flat.

Misc
====

    ./extract.sh $filename

This simple script loads the file into DS9, extracts the initial peaks and then marks their locations on the image. It may be easily modified to perform the PSF photometry and then display its results.
