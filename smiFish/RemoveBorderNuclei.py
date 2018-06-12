"""This function removes the nuclei regions on the border."""


import numpy as np


class RemoveBorderNuclei:
    def __init__(self, nucs_dil):

        mask              =  np.zeros(nucs_dil.shape)
        mask[1:-1, 1:-1]  =  1

        idxs2rem  =  np.unique((1 - mask) * nucs_dil)[1:]
        for k in idxs2rem:
            nucs_dil  *=  (1 - (nucs_dil == k))

        self.nucs_dil  =  nucs_dil
