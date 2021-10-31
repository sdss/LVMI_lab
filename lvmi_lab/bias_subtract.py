import numpy as np
from astropy.io import fits
import argparse

def floatcompress(data, ndig=12):
    '''Adapted from Finkbeiner IDL routine floatcompress'''

    t = data.dtype
    if not ((t == 'float32') or (t == 'float64')):
         error("Only works on floating point numbers")
         raise Exception("Only works on floating point numbers")

    wzer = np.where(data == 0)
    data[wzer] = 1.0

    log2 = np.ceil(np.log(np.abs(data)) / np.log(2.0))
    mant = np.round(data/2.0**(log2 - ndig))/2.0**ndig
    out = mant*2.0**log2

    out[wzer] = 0.0
    return out


def subtract_overscan(dat):
    """ Subtract LVM overscan """



    os1 = np.mean(dat[:,2040:2060], axis=1)
    os2 = np.mean(dat[:,2061:2080], axis=1)

    A = np.tile(os1, (2040,1))
    B = np.tile(os2, (2040,1))

    D  = np.concatenate((dat[:,0:2040], dat[:,2080:]), axis=1)
    OS = np.concatenate((A,B),axis=0).T

    return np.float32(D-OS)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Bias subtract LVM data')
    parser.add_argument('files', type=str, nargs="*", help="File name")
    parser.add_argument('--no_overwrite', action='store_false', help="Do not overwrite files")
    parser.add_argument('--prefix', default="bs_", help="Prefix on output")


    args = parser.parse_args()

    if args.no_overwrite: overwrite = True
    else: overwrite = False


    for file in args.files:
        hdus = fits.open(file)
        dat = hdus[0].data
        res = subtract_overscan(dat)
        hdus[0].data = res
        hdus[0].header["DEBIAS"] = ("Yes", "Frame was debiased")

        hdus.writeto(args.prefix + file.rstrip(".gz"), overwrite=overwrite)
