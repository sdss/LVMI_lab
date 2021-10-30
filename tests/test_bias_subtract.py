from lvmi_lab import __version__
import lvmi_lab
import pytest
import numpy as np


def test_subtract_overscan():

    dat = np.zeros((4120, 4080))

    res = lvmi_lab.subtract_overscan(dat)

    mx = np.max(np.abs(res))

    assert(mx < 0.0001)



