#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import argparse



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Diff')
    parser.add_argument('a', type=str, help="File name")
    parser.add_argument('b', type=str, help="File name")


    args = parser.parse_args()

    hdua = fits.open(args.a)
    hdub = fits.open(args.b)
    dat = hdua[0].data.astype(np.float) - hdub[0].data.astype(np.float)

    hdus = fits.PrimaryHDU(dat)

    hdus.writeto("diff.fits")
