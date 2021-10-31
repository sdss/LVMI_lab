#!/usr/bin/env python 

from lvmi_lab import xcor_frames
import numpy as np
from astropy.io import fits
import argparse





if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Measure XCor of SH Pair')
    parser.add_argument('file1', type=str, nargs=1, help="File name")
    parser.add_argument('file2', type=str, nargs=1, help="File name")


    args = parser.parse_args()


    print(args.file1, args.file2)
    hdu1 = fits.open(args.file1[0])
    hdu2 = fits.open(args.file2[0])

    res = xcor_frames(hdu1[0].data, hdu2[0].data)
    x,y,err = res

    print(err)
    conv = 12/.2 # micron / pix
    x *= conv
    y *= conv
    print("X Array [µm]")
    print(np.array_str(x, precision=2))
    print("Y Array [µm]")
    print(np.array_str(y, precision=2))

    print(" X: %3.2f µm Y: %3.2f µm" % (np.nanmedian(x), np.nanmedian(y)))

