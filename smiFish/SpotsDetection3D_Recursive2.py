"""This function detects spots in the 3D (z-x-y) stack."""


import numpy as np
from scipy.ndimage import filters
from scipy.stats import norm
from skimage.morphology import label


class SpotsDetection3D_Recursive2:
    def __init__(self, green4d, kern_gauss, num_thrs):

        spts_g        =  filters.gaussian_filter(green4d.astype(np.float), kern_gauss)
        print("Hello")
        spts_gl       =  filters.laplace(spts_g.astype(np.float))
        (mu, sigma)   =  norm.fit(np.abs(spts_gl))                              # histogram is fitted with a Gaussian function

        # thr_vals  =  np.linspace(np.abs(spts_gl.min()), np.abs(spts_gl.max()), num_thrs)
        thr_vals  =  np.arange(int(mu + 5 * sigma), int(mu + 5 * sigma) + num_thrs)
        spts_num  =  np.zeros(thr_vals.shape)
        for k in range(thr_vals.size):
            print(k, ' ', thr_vals.size)
            # spts_gl_thr   =  np.abs(spts_gl) > mu + thr_vals[k] * sigma            # thresholding on the histogram
            spts_gl_thr   =  np.abs(spts_gl) > thr_vals[k]                      # thresholding on the histogram
            spts_gl_lbls  =  label(spts_gl_thr, connectivity=1)                   # labelling
            spts_num[k]   =  spts_gl_lbls.max()

        self.spts_num  =  spts_num
