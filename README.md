# LVMI_lab
 LVM-I Laboratory Work 

## LOGIN
ssh as yourself into the server:
ssh XXXXXXXXXXXXXX (ask Nick or Pavan for IP)


## DATA
Go to the data directory
cd /data1/obsld01/npk/lvm

## Sync the data
./SYNC

## cd to the directory
cd 59521

## CODE
Setup aliases etc..
`source ~npk/virtualenv/astropy/bin/activate.csh`

## bias subtract
`~npk/code/LVMI_lab/lvmi_lab/bias_subtract.py *gz` 

### Note, you can bias subtract on individual files, or individual channels. Say you want to do infrared only:
`~npk/code/LVMI_lab/lvmi_lab/bias_subtract.py *z1*gz`


## To assess focus between two files
`~npk/code/LVMI_lab/lvmi_lab/sh_pair.py bs_sdR-s-z1-00001042.fits bs_sdR-s-z1-00001041.fits`

This will spit out three matrices. The first tells you how good the xcor is. The second is the X shift, the third is the Y shift (ignore this one, it will always be 0.34 micron). Units are in µm of focus shift that would lead to this amount of offset in pixel space (if you're off by 0.2 pixels that's 12 µm of defocus, according to my calculations).

```
X Array [micron]
[[-60.77 -60.77 -56.65 -66.95]
 [ -3.09   2.4    4.46  -9.27]
 [  4.46   3.78   1.72  -3.78]
 [  9.27   3.09  -0.34  -1.03]]
Y Array [micron]
[[0.34 0.34 0.34 0.34]
 [0.34 0.34 0.34 0.34]
 [0.34 0.34 0.34 0.34]
 [0.34 0.34 0.34 0.34]]
 X: -0.69 micron  Y: 0.34 micron
X Tilt [mrad]
[[1.3 1.1 1.0 1.2]
 [0.3 0.0 -0.1 0.2]
 [0.2 -0.0 -0.1 0.1]
 [0.0 0.0 0.0 0.0]]
Y Tilt [mrad]
[[nan -0.0 0.2 -0.1]
 [nan 0.6 0.3 -0.1]
 [nan -0.1 -0.1 -0.2]
 [nan -0.7 -0.4 -0.2]]
```

Look at the "X Array" This is the total amount of defocus across the device. You want to get rows 2,3, and 4 to be below 10 microns. Now look at the X tilt and Y tilt. These are the tilts in accordance with the drawing that you made. Call me if you need any further help.


## lvmi_lab

This module (probably should be renamed) currently deals with Hartmann door focus analysis. The purpose of this code is to take in fits files and convert them to focus piston offsets.


![image](https://user-images.githubusercontent.com/3804541/144501541-35c59628-3fa9-484c-be22-687172906f17.png)
