"""This function detects spots in the 3D (z-x-y) stack."""


import numpy as np
from scipy.ndimage import filters
from scipy.stats import norm
from skimage.morphology import label
import multiprocessing


class SpotsDetection3D_Recursive:
    def __init__(self, green4d, kern_gauss, thr_vals):

        spts_g        =  filters.gaussian_filter(green4d, kern_gauss)
        spts_gl       =  filters.laplace(spts_g.astype(np.float))
        (mu, sigma)   =  norm.fit(np.abs(spts_gl))                              # histogram is fitted with a Gaussian function
        print("One")
        cpu_ow   =  multiprocessing.cpu_count()
        t_chops  =  np.array_split(thr_vals, cpu_ow)
        print("Two")
        job_args  =  []
        for k in range(cpu_ow):
            job_args.append([spts_gl, mu, sigma, t_chops[k]])
        print("Three")
        pool     =  multiprocessing.Pool()
        results  =  pool.map(SpotsDetection3D_ForRecursive, job_args)
        pool.close()
        print("Four")
        spts_num  =  np.array([])
        for j in range(len(results)):
            spts_num = np.append(spts_num, results[j].spts_num)

        self.spts_num  =  spts_num



class SpotsDetection3D_ForRecursive:
    def __init__(self, input):

        spts_gl   =  input[0]
        mu        =  input[1]
        sigma     =  input[2]
        thr_vals  =  input[3]

        spts_num  =  np.zeros(thr_vals.shape)
        for k in range(thr_vals.size):
            spts_gl_thr   =  np.abs(spts_gl) > mu + thr_vals[k] * sigma            # thresholding on the histogram
            spts_gl_lbls  =  label(spts_gl_thr)                                 # labelling
            spts_num[k]   =  spts_gl_lbls.max()

        self.spts_num  =  spts_num










# class SpotsDetection3D:
#     def __init__(self, green4d, ker_gauss, thr_val):
#
#         zlen, xlen, ylen  =  green4d.shape
#
#         spts_g        =  filters.gaussian_filter(green4d, ker_gauss)
#         spts_gl       =  filters.laplace(spts_g.astype(np.float))
#         (mu, sigma)   =  norm.fit(np.abs(spts_gl))                              # histogram is fitted with a Gaussian function
#         spts_gl_thr   =  np.abs(spts_gl) > mu + thr_val * sigma                 # thresholding on the histogram
#         spts_gl_lbls  =  label(spts_gl_thr)                                     # labelling
#
#         self.spts_gl_lbls  =  spts_gl_lbls
