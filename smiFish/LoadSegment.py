"""This function loads raw data and the .xls file output of FISH-Quant software.

Raw data are 3D (z-stack) and we perform segmentation slice. The .xls file contains
the coordinate of the center of mas of all the spots detected by FQ, so we know if
spots are inside or outside nucleus.
"""

import numpy as np
import tifffile
from skimage.filters import threshold_local
from skimage.morphology import remove_small_objects, label
from skimage.measure import regionprops
from scipy.ndimage.morphology import binary_fill_holes, binary_dilation, binary_erosion
from skimage.color import label2rgb


import Exc2Spots


class LoadSegment3D:
    def __init__(self, raw_data_fname, xls_data_fname, num_itrs):

        # raw_data               =  tifffile.imread(raw_data_fname)[0, :, :, :, :]
        # zlen, chs, xlen, ylen  =  raw_data.shape
        #
        # for c in range(chs):
        #     for z in range(raw_data.shape[0]):
        #         raw_data[z, c, :, :]  =  np.rot90(raw_data[z, c, :, ::-1])
        raw_data          =  tifffile.imread(raw_data_fname)
        zlen, xlen, ylen  =  raw_data.shape

        for z in range(raw_data.shape[0]):
            raw_data[z, :, :]  =  np.rot90(raw_data[z, :, ::-1])

        nucs        =  raw_data[:, :, :]
        block_size  =  xlen / 4
        if block_size % 2 == 0:
            block_size  +=  1

        j          =  np.argmax(nucs.sum(2).sum(1))
        adapt_thr  =  threshold_local(nucs[j, :, :], block_size)                # adptive thresholding, chenges in each block of the image (adapt_thr is an image)
        nucs_thr   =  nucs > adapt_thr

        nucs_fin      =  np.zeros((zlen, xlen + 2 * num_itrs, ylen + 2 * num_itrs), dtype=np.int)
        xlen2, ylen2  =  nucs_fin.shape[1:]

        nucs_fin[:, num_itrs:xlen2 - num_itrs, num_itrs:xlen2 - num_itrs]  =  nucs_thr
        msk_er                                                             =  np.ones((xlen2, ylen2), dtype=np.int)
        msk_er[num_itrs:xlen2 - num_itrs, num_itrs:ylen2 - num_itrs]       =  0

        for z in range(zlen):
            nucs_fin[z, :, :]  =   binary_dilation(nucs_fin[z, :, :], iterations=num_itrs)
            nucs_fin[z, :, :]  =   binary_fill_holes(nucs_fin[z, :, :])
            nucs_fin[z, :, :]  +=  msk_er
            nucs_fin[z, :, :]  =   binary_erosion(nucs_fin[z, :, :], iterations=num_itrs)

        nucs_fin  =  nucs_fin[:, num_itrs:xlen2 - num_itrs, num_itrs:xlen2 - num_itrs]

        nucs_lbls  =  np.zeros(nucs.shape, dtype=np.int)
        for z in range(zlen):
            nucs_lbls[z, :, :]  =  label(nucs_fin[z, :, :], connectivity=1)
            nucs_lbls[z, :, :]  =  remove_small_objects(nucs_lbls[z, :, :], 150)


        spts       =  np.zeros(nucs.shape)
        spts_ctrs  =  Exc2Spots.Exc2Spots3D(xls_data_fname).spts_ctrs


        for k in range(spts_ctrs.shape[1]):
            if spts_ctrs[2, k] < zlen:
                spts[spts_ctrs[2, k], spts_ctrs[0, k], spts_ctrs[1, k]]  =  1

        # spts_int  =  spts - spts_ext
        spts_ext       =  (spts * (1 - np.sign(nucs_lbls))).astype(np.int)
        spts_ext_ctrs  =  np.zeros((3, spts_ext.sum()), dtype=np.int)
        counter        =  0
        for z in range(zlen):
            bffr_spts  =  label(spts_ext[z, :, :], connectivity=1).astype(np.int)
            bffr_rgp   =  regionprops(bffr_spts)
            if len(bffr_rgp) > 0:
                for i_bffr in range(len(bffr_rgp)):
                    spts_ext_ctrs[:, counter]  =   bffr_rgp[i_bffr]['centroid'][0], bffr_rgp[i_bffr]['centroid'][1], z
                    counter                    +=  1

        mycmap        =  np.fromfile('mycmap.bin', 'uint16').reshape((10000, 3)) / 255.0
        nucs_lbls_3c  =  label2rgb(nucs_lbls, bg_label=0, bg_color=[0, 0, 0], colors=mycmap)

        nucs_lbls_3c[:, :, :, 0]  *=  (1 - spts)
        nucs_lbls_3c[:, :, :, 1]  *=  (1 - spts)
        nucs_lbls_3c[:, :, :, 2]  *=  (1 - spts)
        nucs_lbls_3c[:, :, :, 1]  +=  spts


        self.nucs           =  nucs
        self.nucs_lbls      =  nucs_lbls
        self.spts_ext       =  spts_ext
        self.spts_ext_ctrs  =  spts_ext_ctrs
        self.nucs_lbls_3c   =  nucs_lbls_3c
