#!/usr/bin/env python 

import argparse
import numpy as np

from lvmi_lab import xcor_frames
from astropy.io import fits


r1 = "r1"
z1 = "z1"

bpm_paths = {r1: "/home/npk/code/LVMI_lab/lvmi_lab/bpm-r1.fits", 
             z1: "/home/npk/code/LVMI_lab/lvmi_lab/bpm-z1.fits"}



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Measure XCor of SH Pair')
    parser.add_argument('file1', type=str, nargs=1, help="File name")
    parser.add_argument('file2', type=str, nargs=1, help="File name")



    args = parser.parse_args()

    print(args)

    print(args.file1, args.file2)
    hdu1 = fits.open(args.file1[0])
    hdu2 = fits.open(args.file2[0])

    if hdu1[0].header["HARTMANN"] != '1 0':
        print("LEFT IS NOT LEFT!!! -- LIKELY SIGN IS WRONG!!!")
    if hdu2[0].header["HARTMANN"] != '0 1':
        print("RIGHT IS NOT RIGHT -- LIKELY SIGN IS WRONG!!!")
        
    
    d1 = hdu1[0].data
    d2 = hdu2[0].data
    
    bpm = []
    if True:
        if r1 in args.file1[0]:
            print("Using BPM: %s" % bpm_paths[r1])
            bpm = fits.open(bpm_paths[r1])[0].data
            bpm = np.where(bpm==1)

        if z1 in args.file1[0]:
            print("Using BPM: %s" % bpm_paths[z1])
            bpm = fits.open(bpm_paths[z1])[0].data
            bpm = np.where(bpm==1)
    
    d1[bpm] = 0
    d2[bpm] = 0
        
    res = xcor_frames(d1, d2)
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
    print("X Array [micron]")
    print(np.array_str(x, precision=2))
    print("Y Array [micron]")
    print(np.array_str(y, precision=2))

    print(" X: %3.2f micron  Y: %3.2f micron" % (np.nanmedian(x), np.nanmedian(y)))

