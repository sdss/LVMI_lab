#!/usr/bin/env python 

import argparse
import numpy as np

from lvmi_lab import xcor_frames
from astropy.io import fits





if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Measure XCor of SH Pair')
    parser.add_argument('file1', type=str, nargs=1, help="File name")
    parser.add_argument('file2', type=str, nargs=1, help="File name")


    args = parser.parse_args()


    print(args.file1, args.file2)
    hdu1 = fits.open(args.file1[0])
    hdu2 = fits.open(args.file2[0])

    if hdu1[0].header["HARTMANN"] != '1 0':
        print("LEFT IS NOT LEFT!!!")
    if hdu2[0].header["HARTMANN"] != '0 1':
        print("RIGHT IS NOT RIGHT")

    res = xcor_frames(hdu1[0].data, hdu2[0].data)
    #return x_shifts, y_shifts, best, Warnings
    x,y,best,warn = res

    conv = 12/.2 # micron / pix
    x *= conv
    y *= conv
    print("values")
    best /= np.nanmax(best)
    print(np.array_str(best,precision=2))
    bad = best < 1e-4
    x[bad] = np.nan
    y[bad] = np.nan
    print("X Array [µm]")
    print(np.array_str(x, precision=2))
    print("Y Array [µm]")
    print(np.array_str(y, precision=2))

    print(" X: %3.2f µm Y: %3.2f µm" % (np.nanmedian(x), np.nanmedian(y)))

