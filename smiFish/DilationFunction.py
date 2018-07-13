"""This function perform dilation on nuclei. The overlapping is handled."""

import numpy as np
# import pyqtgraph as pg
from scipy.ndimage.morphology import binary_dilation


class DilationFunction:
    def __init__(self, in_vabls):

        inds      =  in_vabls[0]
        nucs_dil  =  in_vabls[1]

        msk  =  np.zeros(nucs_dil.shape, dtype=np.int)
        for i in inds:
            smpl  =   (nucs_dil == i) * 1
            msk   +=  i * binary_dilation(smpl, iterations=4) * (1 - np.sign(msk))

        self.msk  =  msk
        # pg.image(msk)
