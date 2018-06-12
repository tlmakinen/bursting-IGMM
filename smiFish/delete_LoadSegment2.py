"""This function."""


import numpy as np
import tifffile
from skimage.filters import gaussian, threshold_otsu, median
from skimage.morphology import label, remove_small_objects, remove_small_holes
from skimage.measure import regionprops
from skimage.color import label2rgb
from scipy.optimize import curve_fit
from scipy.ndimage.morphology import binary_fill_holes, binary_dilation, binary_erosion
from scipy.ndimage import distance_transform_edt
from skimage.morphology import watershed
from skimage.feature import peak_local_max

import Exc2Spots


class LoadSegment2D:
    def __init__(self, raw_data_fname, xls_data_fname):

        raw_data               =  tifffile.imread(raw_data_fname)[0, :, :, :, :]
        zlen, chs, xlen, ylen  =  raw_data.shape

        for c in range(chs):
            for t in range(raw_data.shape[0]):
                raw_data[t, c, :, :]  =  np.rot90(raw_data[t, c, :, ::-1])

        nucs  =  raw_data[:, 1, :, :]

        raw_data_mip  =  np.zeros((xlen, ylen, 3))

        for c in range(chs):
            for x in range(xlen):
                raw_data_mip[x, :, c]  =  raw_data[:, c, x, :].max(0)

        nucs_f    =  gaussian(nucs, 5)
        nucs_mip  =  np.zeros((xlen, ylen))
        for x in range(xlen):
            nucs_mip[x, :]  =  nucs_f[:, x, :].max(0)

        val       =  threshold_otsu(nucs_mip)
        nucs_det  =  nucs_mip > val

        mycmap  =  np.fromfile('mycmap.bin', 'uint16').reshape((10000, 3)) / 255.0

        nucs_labls     =  label(nucs_det, connectivity=1)
        nucs_labls     =  remove_small_objects(nucs_labls, 50)
        nucs_labls_3c  =  label2rgb(nucs_labls, bg_color=[0, 0, 0], bg_label=0, colors=mycmap)

        spts       =  np.zeros(nucs_det.shape)
        spts_ctrs  =  Exc2Spots.Exc2Spots2D(xls_data_fname).spts_ctrs
        for k in range(spts_ctrs.shape[1]):
            spts[spts_ctrs[0, k], spts_ctrs[1, k]] = 1


        nucs_labls_3c[:, :, 0]  *=  (1 - spts)
        nucs_labls_3c[:, :, 1]  *=  (1 - spts)
        nucs_labls_3c[:, :, 2]  *=  (1 - spts)
        nucs_labls_3c[:, :, 1]  +=  spts

        self.raw_data_mip   =  raw_data_mip
        self.nucs_det       =  nucs_det
        self.nucs_labls     =  nucs_labls
        self.nucs_labls_3c  =  nucs_labls_3c



class LoadSegment3D_bis:
    def __init__(self, raw_data_fname, xls_data_fname, thr_rate):

        raw_data               =  tifffile.imread(raw_data_fname)[0, :, :, :, :]
        zlen, chs, xlen, ylen  =  raw_data.shape

        for c in range(chs):
            for z in range(raw_data.shape[0]):
                raw_data[z, c, :, :]  =  np.rot90(raw_data[z, c, :, ::-1])

        nucs         =  raw_data[:, 0, :, :]
        bbinns       =  np.arange(nucs.min(), nucs.max(), 150)
        hh           =  np.histogram(nucs, bins=bbinns)

        def gaus(x, a, x0, sigma):
            return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

        mean   =  hh[0].mean()
        sigma  =  hh[0].std()

        popt, pcov  =  curve_fit(gaus, hh[1][1:], hh[0], p0=[1, mean, sigma])
        # yy          =  gaus(hh[1][1:], *popt)
        soglia      =  np.abs(popt[1]) + np.abs(popt[2])
        thr_change  =  nucs > soglia
        nucs2       =  nucs * (1 - thr_change) + thr_change * soglia

        nucs2_gf  =  np.zeros(nucs2.shape)
        for z in range(zlen):
            nucs2_gf[z, :, :]  =  gaussian(nucs2[z, :, :], 4)


        # thr_rate   =  1
        val        =  threshold_otsu(nucs2_gf)
        nucs2_thr  =  nucs2_gf > val
        nucs_fin   =  np.zeros(nucs.shape)
        er_iter    =  5
        for z in range(zlen):
            nucs_fin[z, :, :]  =  binary_dilation(nucs2_thr[z, :, :], iterations=5)
            nucs_fin[z, :, :]  =  binary_fill_holes(nucs_fin[z, :, :])
            nucs_fin[z, :, :]  =  binary_erosion(nucs_fin[z, :, :], iterations=er_iter)

        msk_brd                                                     =   np.ones((nucs.shape))
        msk_brd[:, er_iter:xlen - er_iter, er_iter:ylen - er_iter]  =   0
        nucs_fin                                                   +=  msk_brd * nucs2_thr

        raw_det              =  np.zeros(np.append(nucs.shape, 3))                                                      # this is just for visual purposes
        raw_det[:, :, :, 0]  =  nucs.astype(np.float) / nucs.max()
        raw_det[:, :, :, 1]  =  0.1 * nucs_fin

        nucs_lbls  =  np.zeros(nucs.shape, dtype=np.int)
        for z in range(zlen):
            nucs_lbls[z, :, :]  =  label(nucs_fin[z, :, :], connectivity=1)
            nucs_lbls[z, :, :]  =  remove_small_objects(nucs_lbls[z, :, :], 150)


        areas     =  np.zeros((zlen, nucs_lbls.max()))
        nuc_tags  =  np.zeros((zlen, nucs_lbls.max()))
        # rgp    =  regionprops(nucs_lbls.astype(np.int))
        # areas  =  np.zeros(len(rgp))
        for z in range(zlen):
            rgp  =  regionprops(nucs_lbls[z, :, :])
            for i in range(len(rgp)):
                areas[z, i]     =  rgp[i]['area']
                nuc_tags[z, i]  =  rgp[i]['label']


        ar_vct  =  areas.reshape(areas.size)
        ar_vct  =  np.delete(ar_vct, np.where(ar_vct == 0))

        areas_av        =  ar_vct.mean()
        [z_zrs, j_zrs]  =  np.where(areas > 2*areas_av)
        rad_av          =  np.round(np.sqrt(areas_av)).astype(int)

        patch  =  np.zeros(nucs.shape, dtype=np.int)
        if len(j_zrs) > 0:
            for jj in range(z_zrs.size):
                smpl        =  (nucs_lbls[z_zrs[jj], :, :] == nuc_tags[z_zrs[jj], j_zrs[jj]])
                print jj, smpl.sum()/areas_av
                # pg.image(smpl, title=str(jj))
                distance    =  distance_transform_edt(smpl)
                local_maxi  =  peak_local_max(distance, indices=False, footprint=np.ones((rad_av, rad_av)))
                markers     =  label(local_maxi)
                lab_water   =  watershed(-distance, markers, mask=smpl)
                pg.image(lab_water, title=str(jj))
                patch[z_zrs[jj], :, :]  +=  patch[z_zrs[jj], :, :] * (1 - np.sign(smpl*1)) + lab_water + np.sign(lab_water)
                patch[z_zrs[jj], :, :]  =   label(patch[z_zrs[jj], :, :])

        nucs_lbls  *=  1 - np.sign(patch)
        nucs_lbls  +=  patch

        for z in range(zlen):
            nucs_lbls[z, :, :]  =  label(nucs_lbls[z, :, :], connectivity=1)


        spts       =  np.zeros(nucs.shape)
        spts_ctrs  =  Exc2Spots.Exc2Spots3D(xls_data_fname, zlen).spts_ctrs


        for k in range(spts_ctrs.shape[1]):
            if spts_ctrs[2, k] < zlen + 1:
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

        self.nucs           =  nucs
        self.nucs_lbls      =  nucs_lbls
        self.spts_ext       =  spts_ext
        self.spts_ext_ctrs  =  spts_ext_ctrs




class LoadSegment3D:
    def __init__(self, raw_data_fname, xls_data_fname, thr_rate):

        raw_data               =  tifffile.imread(raw_data_fname)[0, :, :, :, :]
        zlen, chs, xlen, ylen  =  raw_data.shape

        for c in range(chs):
            for z in range(raw_data.shape[0]):
                raw_data[z, c, :, :]  =  np.rot90(raw_data[z, c, :, ::-1])

        nucs      =  raw_data[:, 1, :, :]
        nucs_f    =  gaussian(nucs, 1)
        j         =  np.argmax(nucs_f.sum(2).sum(1))
        val       =  threshold_otsu(nucs_f[j, :, :])
        bkg       =  (nucs_f[j, :, :] < val) * nucs_f[j, :, :]
        bkg       =  bkg.reshape(bkg.size)
        bkg       =  np.delete(bkg, np.where(bkg == 0))
        nucs_det  =  nucs_f > bkg.mean() + thr_rate * bkg.std()

        for z in range(zlen):
            nucs_det[z, :, :]  =  median(nucs_det[z, :, :], selem=np.ones((5, 5)))                                    # median filter on thresholded image to get rid of defects
            nucs_det[z, :, :]  =  remove_small_holes(nucs_det[z, :, :], 400)                                            # remove small holes and small objects
            nucs_det[z, :, :]  =  remove_small_objects(nucs_det[z, :, :], 200)

        raw_det              =  np.zeros(np.append(nucs.shape, 3))                                                      # this is just for visual purposes
        raw_det[:, :, :, 0]  =  nucs.astype(np.float) / nucs.max()
        raw_det[:, :, :, 1]  =  0.1 * nucs_det

        spts       =  np.zeros(nucs.shape)
        spts_ctrs  =  Exc2Spots.Exc2Spots3D(xls_data_fname, zlen).spts_ctrs


        for k in range(spts_ctrs.shape[1]):
            if spts_ctrs[2, k] < zlen + 1:
                spts[spts_ctrs[2, k], spts_ctrs[0, k], spts_ctrs[1, k]]  =  1

        # spts_int  =  spts - spts_ext
        spts_ext       =  (spts * (1 - nucs_det)).astype(np.int)
        spts_ext_ctrs  =  np.zeros((3, spts_ext.sum()), dtype=np.int)
        counter        =  0
        for z in range(zlen):
            bffr_spts  =  label(spts_ext[z, :, :], connectivity=1).astype(np.int)
            bffr_rgp   =  regionprops(bffr_spts)
            if len(bffr_rgp) > 0:
                for i_bffr in range(len(bffr_rgp)):
                    spts_ext_ctrs[:, counter]  =   bffr_rgp[i_bffr]['centroid'][0], bffr_rgp[i_bffr]['centroid'][1], z
                    counter                    +=  1


        nucs_det_lbls  =  np.zeros(nucs.shape, dtype=np.int)
        for z in range(zlen):
            nucs_det_lbls[z, :, :]  =  label(nucs_det[z, :, :], connectivity=1)

        mycmap            =  np.fromfile('mycmap.bin', 'uint16').reshape((10000, 3)) / 255.0
        nucs_det_lbls_3c  =  label2rgb(nucs_det_lbls, bg_label=0, bg_color=[0, 0, 0], colors=mycmap)

        nucs_det_lbls_3c[:, :, :, 0]  *=  (1 - spts)
        nucs_det_lbls_3c[:, :, :, 1]  *=  (1 - spts)
        nucs_det_lbls_3c[:, :, :, 2]  *=  (1 - spts)
        nucs_det_lbls_3c[:, :, :, 1]  +=  spts


        self.nucs              =  nucs
        self.nucs_det_lbls     =  nucs_det_lbls
        self.spts_ext          =  spts_ext
        self.spts_ext_ctrs     =  spts_ext_ctrs
        self.nucs_det_lbls_3c  =  nucs_det_lbls_3c

        # self.spts           =  spts
        # self.spts_ctrs      =  spts_ctrs
