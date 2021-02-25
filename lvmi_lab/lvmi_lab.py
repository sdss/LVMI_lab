import numpy as np
import numpy.random
import astropy.io.fits as fits
from scipy import ndimage as ndi
from scipy.signal import convolve, convolve2d

from multiprocessing import Pool


def generate_fake_frame(RN=3, gain=2, n_pix=4096, adc_bits=16, add=0,
    header_add={}):
    """ Return simulated HDU Frame
    
    Args:
        RN: Readnoise in e-
        gain: Gain in e-/ADU 
        n_pix: Number of pixels on a side.
        adc_bits: Number of bits in the ADC
        add: Add this value (e-)
        header_add: Add these header values
    
    Returns:
        PrimaryHDU with frame in ADU

     """

    IMG = (add + np.random.randn(n_pix, n_pix) * RN) / gain
    if adc_bits == 16:
        IMG = IMG.astype(np.int16)
    else:
        print("Raise")
        raise NotImplementedError()

    HDU = fits.PrimaryHDU(IMG)

    HDU.header["GAIN"] = gain
    HDU.header["RN"] = RN
    for k,v in header_add.items():
        HDU.header[k] = v

    return HDU


#PSF comes from zemax simulations
PSF = """  0.0000E+00	  0.0000E+00	  0.0000E+00	  0.0000E+00	  0.0000E+00	  0.0000E+00
  0.0000E+00	  0.0000E+00	  2.0672E+00	  2.7464E+01	  4.7250E+00	  0.0000E+00
  0.0000E+00	  0.0000E+00	  4.3529E+02	  6.3403E+02	  1.9491E+01	  0.0000E+00
  0.0000E+00	  0.0000E+00	  6.0332E+02	  6.8630E+02	  9.4500E+00	  0.0000E+00
  0.0000E+00	  0.0000E+00	  3.1598E+01	  3.2780E+01	  0.0000E+00	  0.0000E+00
  0.0000E+00	  0.0000E+00	  0.0000E+00	  2.9531E-01	  0.0000E+00	  0.0000E+00
"""


def really_fake_line_frame(n_fibers=10, n_pix=4096, nlines=100, counts=30000):
    """ Makes a really fake "emission line" spectrum """

    delta = n_pix/(n_fibers+1)
    img = np.zeros((n_pix,n_pix))

    w=3
    psf=np.fromstring(PSF, sep="\t").reshape(2*w,2*w)
    psf /= np.max(psf)
    
    psf *= counts

    dlines = n_pix/nlines
    poss = np.arange(delta, n_pix-delta, delta).astype(int)
    line_locs = np.arange(dlines,n_pix-dlines,delta).astype(int)
    for pos in poss:
        for line_loc in line_locs:
            s1 = slice(pos-w,pos+w)
            s2 = slice(line_loc-w,line_loc+w)
            img[s1,s2] = np.random.poisson(psf)
    
    return img


def shift_frame(fm, dx, dy):

    return ndi.shift(fm, (dx, dy))
    return sincshift2d(fm, dx, dy)


def piston_focus(fm, microns):
    """ Shift the PSF by microns of defocus.

    Note I use measurements from zemax that 0.2 pix = 12 micron

    Args:
        fm: Frame
        microns: Defocus amount in microns"""

    pixel_motion = microns * 0.016667

    return shift_frame(fm, pixel_motion, 0) 


def generate_shifted_frame(microns, **kwargs):
    lines = really_fake_line_frame(**kwargs)
    lines = piston_focus(lines, microns, **kwargs)

    focus_actuators = {"focus1": microns, "focus2": microns, "focus3": microns}

    hdu = generate_fake_frame(add=lines, header_add=focus_actuators, **kwargs)
    return hdu

#- Utility functions for sinc shifting pixelated PSFs. From specter repo
def _sincfunc(x, dx, dampfac=3.25):
    """sinc helper function for sincshift()"""
    if dx != 0.0:
        xx = (x+dx)*np.pi  #- cache shifted array for 30% faster evals
        return np.exp( -(xx/(dampfac*np.pi))**2 ) * np.sin(xx) / xx
    else:
        xx = np.zeros(len(x))
        xx[len(x)//2] = 1.0
        return xx
    

# From specter repo
def sincshift2d(image, dx, dy, sincrad=10, dampfac=3.25):
    """
    Return image shifted by dx, dy using full 2D sinc interpolation
    """
    s = np.arange(-sincrad, sincrad+1.0)
    sincx = _sincfunc(s, -dx, dampfac=dampfac)
    sincy = _sincfunc(s, -dy, dampfac=dampfac)
    kernel = np.outer(sincy, sincx)
    newimage = convolve2d(image, kernel, mode='same')
    return newimage


def xcor_frames_ascent_helper(ijd, threshold=0, iter_max=500):
    """ Cross correlation - ascent algorithm - helper function 
    
    This bit uses an "ascent" algorithm to find the maximum.
    
    args:
        ijd: i, j, data
            data is a dictionary that contains all the key local 
            variables from the xcor_frames function
        threshold: Threshold value for minimum improvement 
        iter_max: Maximum number of iterations to allow in the
            ascent algorithm
    
    Returns a tuple with:
        i,j: Sames are input
        (best i, best j): Best location in the ascent space 
            for cross correlation
        res[subsample, subsample]: Resulting sparsley populated 2d array
        shiftsize: Same as input
        Warnings: String of any warnings. Contains "Success"
            on success.
    """

    i,j,data = ijd
    for k,v in data.items():
        globals()[k] = v

    sl = (slice(i*cut_into, i*cut_into+cut_into), 
            slice(j*cut_into, j*cut_into+cut_into))
    stepsize = pm_pixels*2/(subsample-1)
    shiftsize = np.arange(-pm_pixels, pm_pixels+stepsize, stepsize)

    res = np.zeros((subsample, subsample))
    res[:] = np.nan

    start_i,start_j= subsample//2, subsample//2

    prev = 0
    improvement = np.inf
    niter = 0
    while (improvement > threshold) and (niter < iter_max):
        for di in [-1,0,1]:
            for dj in [-1,0,1]:
                ix = di + 1 + start_i
                jx = dj + 1 + start_j
                try: 
                    if not np.isfinite(res[ix,jx]):
                        a = shift_frame(A[sl], shiftsize[ix], shiftsize[jx])
                        xcor = np.sum(a*B[sl])
                        res[ix, jx] = xcor
                except IndexError:
                    return (i,j,np.array((start_i,start_j)),res,shiftsize,"Hit Edge")
        
        r = np.nanargmax(res)
        start_i, start_j = np.unravel_index(r, res.shape)

        improvement = res[start_i,start_j] - prev
        prev = res[start_i,start_j]
        niter += 1
    
    if niter >= iter_max:
        return (i,j,np.array((start_i,start_j)),res,shiftsize,"Exceeded iteration limit")

    sl = slice(start_i-1,start_i+2)
    xcom = np.nansum(res[sl,j] * shiftsize[sl])/np.sum(res[sl,j])
    ycom = np.nansum(res[sl,j] * shiftsize[sl])/np.sum(res[sl,j])
    sl = slice(start_j-1,start_j+2)
    xcom = np.nansum(res[start_i,sl] * shiftsize[sl])/np.sum(res[start_i,sl])
    ycom = np.nansum(res[start_i,sl] * shiftsize[sl])/np.sum(res[start_i,sl])
    
    return (i,j,np.array((xcom,ycom)),res,shiftsize,"Success")

def xcor_frames_helper(ijd):

    i,j,data = ijd
    for k,v in data.items():
        globals()[k] = v

    sl = (slice(i*cut_into, i*cut_into+cut_into), 
            slice(j*cut_into, j*cut_into+cut_into))
    stepsize = pm_pixels*2/(subsample-1)
    todo = np.arange(-pm_pixels, pm_pixels+stepsize, stepsize)

    res = np.zeros((subsample, subsample))
    for ix, sx in enumerate(todo):
        for jx,sy  in enumerate(todo):
            a = shift_frame(A[sl], sx, sy)
            res[ix,jx] = np.sum(a*B[sl])
    
    return (i,j,ndi.center_of_mass(res),res,todo)


def xcor_frames(A, B, pm_pixels=1.0, subsample=200, cut_into=1024):
    """ Return the spatial cross correlation of two frames over +-pm_pixes

    Based on simulations, we want a precision of 0.02 pixels (50 subsamples)
    per pixel of cross correlation. So if pm_pixels is 1.0 then 100 subsamples
    are needed. This is really slow, so a compromise may be needed.

    Args:
        A: First frame
        B: Second frame
        pm_pixels: XCor of this many pixels (e.g. 1 means +- 1 pixel)
        subsample: How many cross correlation steps to do
        cut_into: Slice the frame into sections that are this big

    Returns:
        x_shifts[# cuts, # cuts]: Measured pixel offsets between A/B in X
        y_shifts[# cuts, # cuts]: Measured pixel offsets between A/B in X
        warnings: String of warnings

    Note that one will tend to fit a plane to x_shifts and y_shifts in order
    to measure the piston and tilt of the CCD.

"""

    run = A.shape[0]//cut_into
    res = np.zeros((run, run, subsample, subsample))

    
    dat = {"A": A, "B": B, "pm_pixels": pm_pixels, "subsample": subsample, \
        "cut_into": cut_into}

    todo = [(i,j,dat) for i in range(run) for j in range(run)]
    p = Pool()
    #res = list(map(xcor_frames_ascent_helper, todo))
    res = p.map(xcor_frames_ascent_helper, todo)
    p.close()

    # "res" contains the results of all the threads, and needs to be repackaged
    # res is (i,j,((xcom,ycom)),res,shiftsize,"Success")

    x_shifts = np.zeros((run,run))
    y_shifts = np.zeros((run,run))
    x_shifts[:] = np.nan
    y_shifts[:] = np.nan


    Warnings = ""
    for r in res:
        if r[5] != "Success": 
            Warnings = "%s\n%s" % (Warnings, r[5])
        x_shifts[r[0],r[1]] = r[2][0]
        y_shifts[r[0],r[1]] = r[2][1]


    #from IPython import embed; embed()

    return x_shifts, y_shifts, Warnings








