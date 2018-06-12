"""This function detects spots in the 4D (z-x-y) stack."""


import numpy as np
from scipy.ndimage import filters
from scipy.stats import norm
from skimage.morphology import label


class SpotsDetection3D:
    def __init__(self, green4d, ker_gauss, thr_val):

        zlen, xlen, ylen  =  green4d.shape

        spts_g        =  filters.gaussian_filter(green4d, ker_gauss)
        spts_gl       =  filters.laplace(spts_g.astype(np.float))
        (mu, sigma)   =  norm.fit(np.abs(spts_gl))                              # histogram is fitted with a Gaussian function
        spts_gl_thr   =  np.abs(spts_gl) > mu + thr_val * sigma                 # thresholding on the histogram
        spts_gl_lbls  =  label(spts_gl_thr)                                     # labelling

        self.spts_gl_lbls  =  spts_gl_lbls
