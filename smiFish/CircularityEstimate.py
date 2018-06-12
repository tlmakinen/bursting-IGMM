"""This function calculates the circularity, defined as the ratio between
4*np.pi*A and the p**2 (A is the surface and p is the perimeter), of all
the connected objects present in the image-matrix labels. This function
a 2D matrix with the value of the circularity (ies) and the label of the
connected object
"""


import numpy as np
from skimage import measure


class CircularityEstimate:
    def __init__(self, labels):
        self.labels  =  labels

        rgp_left   =  measure.regionprops(self.labels)
        circ       =  np.zeros((2, len(rgp_left)))

        for k in range(len(rgp_left)):
            # print k
            circ[:, k]  =  4 * np.pi * rgp_left[k]['Area'] / rgp_left[k]['Perimeter'] ** 2, rgp_left[k]['Label']

        self.circ  =  circ