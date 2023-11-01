#!/usr/bin/env python 

import argparse
import numpy as np
import json


from lvmi_lab import xcor_frames, get_positions_on_ccd, hartman_focus_by_peak_finding
from astropy.io import fits


from sklearn.linear_model import HuberRegressor



def regress(x, y):
    """ Returns a linear fitting function """
    X = np.zeros((x.shape[0],1))
    X[:,0] = x
    huber = HuberRegressor(epsilon=1)
    huber.fit(X,y)

    return np.poly1d([huber.coef_[0], huber.intercept_])



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

    f1 = args.file1[0]
    f2 = args.file2[0]
    print(f1,f2)

    hdu1 = fits.open(f1)
    hdu2 = fits.open(f2)

    if hdu1[0].header["HARTMANN"] != '1 0':
        print("LEFT IS NOT LEFT!!! -- LIKELY SIGN IS WRONG!!!")
    if hdu2[0].header["HARTMANN"] != '0 1':
        print("RIGHT IS NOT RIGHT -- LIKELY SIGN IS WRONG!!!")
        
    OBSTIME = hdu1[0].header["OBSTIME"]
    d1 = hdu1[0].data
    d2 = hdu2[0].data

    tl, tr, ox, oy = hartman_focus_by_peak_finding(d1,d2)
    x,y = tl.data.T[0], tl.data.T[1]
    ox *= -12/.2

    from pylab import *
    l,r = map(lambda x: x.rstrip(".fits").split("-")[-1], [f1,f2])
    flav = f1.split("/")[-1].split("-")[2]
    framenum = l.split(".")[0]

    fig = figure(figsize=(12,6))
    subplot(2,2,1)
    scatter(x,y,c=ox) ; colorbar()
    clim(-60,60)
    title("%s: %s/%s" % (flav, l, r))
    

    lims = [-80,80]
    #### TOP RIGHT
    subplot(2,2,2)
    XR = np.arange(0,4096,300)
    plot(y, ox, '.') ; grid(True)
    reg = regress(y, ox)
    plot(XR, reg(XR), lw=3)
    axhline(np.median(ox))
    xlabel("Defocus [micron]")
    title("Y Slope is %3.1f [µm/4096 pix]" % (reg.coef[0]*4096))
    yslope = reg.coef[0] * 4096
    ylim(*lims)


    ##### BOTTOM LEFT
    subplot(2,2,3)
    XR = np.arange(0,4096,300)
    plot(x, ox,'.') ; grid(True)
    reg = regress(x, ox)
    plot(XR, reg(XR), lw=3)
    axhline(np.median(ox))
    title("X Slope is %3.1f [µm/4096 pix]" % (reg.coef[0]*4096))
    xslope = reg.coef[0]*4096

    ylabel("Defocus [micron]")
    ylim(*lims)


    subplot(2,2,4)
    hist(ox, range=[-30,30], bins=50)
    xlabel("Defocus [µm]")
    q10,q50,q90 = np.quantile(ox, [0.1, 0.5, 0.9])
    title("90percent within %3.1f to %3.1f" % (q10,q90))


    tight_layout()
    outname = "/Users/npk/Dropbox/REPOS/LVM-focus-measurements/results/%s-%s-%s-fig.pdf" % (flav,l,r)
    print( "*I***")
    print(outname)
    savefig(outname)

    import os
    #os.system("open %s" % outname)
    #os.system("ds9 %s -zscale %s -zscale -single &" % (f1, f2))

    fname = "/Users/npk/Dropbox/REPOS/LVM-focus-measurements/results/results.json"
    try:
        with open(fname, 'r') as f:
            all_results = json.load(f)
    except FileNotFoundError:
        all_results = {}
    
    key = "%s-%s" % (framenum, flav)
    all_results[key] = [OBSTIME, xslope, yslope, q10, q50, q90]
    print(all_results)

    with open(fname, "w") as f:
        json.dump(all_results, f)

    