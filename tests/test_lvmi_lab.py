from lvmi_lab import __version__
import lvmi_lab
import pytest
import numpy as np


def test_version():
    assert __version__ == '0.1.0'

def test_fake_frame():

    test = {"hello": "world"}
    ff = lvmi_lab.generate_fake_frame(header_add=test)

    assert ff.shape == (4096, 4096)
    assert ff.header["hello"] == "world"

    rn = np.std(ff.data)

    assert abs(rn - 1.15) < .01


def test_shifted_frame():

    fm10 = lvmi_lab.generate_shifted_frame(10)
    fm0 = lvmi_lab.generate_shifted_frame(0)

    assert np.abs(np.max(fm0.data) - 15000) < 500
    fm10.writeto("shift_10.fits", overwrite=True)
    fm0.writeto("noshift.fits", overwrite=True)


def test_fake_frame_adc():
    with pytest.raises(NotImplementedError):
        lvmi_lab.generate_fake_frame(adc_bits=18)


def test_xcor_frames():
    A = lvmi_lab.really_fake_line_frame()
    B = lvmi_lab.shift_frame(A, 0.5, 0.5)

    x_shifts, y_shifts, warn = lvmi_lab.xcor_frames(A,B)

    assert np.mean(abs(x_shifts - 0.5))/0.5 < 0.01
    assert np.mean(abs(y_shifts - 0.5))/0.5 < 0.01
    assert warn == ""


