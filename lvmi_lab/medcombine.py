#!/usr/bin/env python


import numpy as np
from astropy.io import fits
import argparse
import ccdproc
from astropy.nddata import CCDData


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Median combine images')
    parser.add_argument('files', type=str, nargs='*', help="File name")

    args = parser.parse_args()
    dats = []

    for file in args.files:
        ff = fits.open(file)
        dats.append(CCDData(ff[0].data, unit='adu'))
    
    c = ccdproc.Combiner(dats)

    ff[0].data = c.median_combine()

    ff.writeto("med.fits")
    

