#!/usr/bin/env python


import numpy as np
from astropy.io import fits
import argparse
from astropy.nddata import CCDData


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Median combine images')
    parser.add_argument('files', type=str, nargs='*', help="File name")

    args = parser.parse_args()
    dats = []

    for file in args.files:
        ff = fits.open(file)
        dats.append(ff[0].data)
    

    ff[0].data = np.median(dats, 0)

    ff.writeto("med.fits")
    

