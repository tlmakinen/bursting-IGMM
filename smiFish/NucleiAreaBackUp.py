import numpy as np
import multiprocessing
# import time
from scipy.ndimage.morphology import binary_dilation, binary_erosion
from skimage.morphology import label, remove_small_objects
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import regionprops
from scipy.ndimage import distance_transform_edt
from skimage.morphology import watershed
from skimage.feature import peak_local_max



class NucleiArea:
    def __init__(self, nucs_lbls, spts_ctrs, spts):

        zlen, xlen, ylen  =  nucs_lbls.shape

        nucs_er  =  np.zeros(nucs_lbls.shape, dtype=int)
        for z in range(zlen):
            nucs_er[z, :, :]  =  label(binary_erosion(nucs_lbls[z, :, :], iterations=6), connectivity=1)
            nucs_er[z, :, :]  =  remove_small_objects(nucs_er[z, :, :], 120)


        j_max     =  np.argmax(np.sign(nucs_er).sum(2).sum(1))
        nucs_prj  =  np.copy(nucs_er[j_max, :, :])
        new_idx   =  1000
        # t1 = time.time()
        for z in range(1, zlen):
            print z
            idxs  =  np.unique(nucs_er[z, :, :])[1:]
            for k in idxs:
                smpl  =  ((nucs_er[z, :, :] == k) * nucs_prj).reshape((xlen * ylen))

                if smpl.sum() == 0:
                    nucs_prj  +=  new_idx * (nucs_er[z, :, :] == k)
                    new_idx   +=  1

                else:
                    smpl      =   np.delete(smpl, np.where(smpl == 0)[0])
                    smpl      =   np.median(smpl).astype(int)
                    nucs_prj  +=  smpl * (nucs_er[z, :, :] == k) * (1 - np.sign(nucs_prj))
        # deltaT1 = time.time() - t1

        nucs_prj  =  np.copy(nucs_er[j_max, :, :])
        # # def f_prj(tags, nucs2add, nucs_prj):
        # def f_prj(args):
        #     new_idx  =  args[2].max() + 1
        #     # bffr     =  np.zeros(args[1].shape, dtype=np.int)
        #     bffr     =  np.copy(args[2])
        #     for k in args[0]:
        #         smpl     = ((args[1] == k) * args[2]).reshape((xlen * ylen))
        #
        #         if smpl.sum() == 0:
        #             args[2]  +=  new_idx * (args[1] == k)
        #             new_idx   +=  1
        #
        #         else:
        #             smpl     =   np.delete(smpl, np.where(smpl == 0)[0])
        #             smpl     =   np.median(smpl).astype(int)
        #             bffr     +=  smpl * np.sign(1*(args[1] == k) + 1*(args[2] == smpl))
        #     return bffr
        #
        # cpu_owe   =  multiprocessing.cpu_count()
        # t2  =  time.time()
        # for z in range(zlen):
        #     print z
        #     idxs  =  np.unique(nucs_prj)[1:]
        #     # idxs  =  np.unique(nucs_er[z, :, :])[1:]
        #     tags  =  np.split(idxs, (idxs.size / cpu_owe + 1) * np.arange(1, cpu_owe, dtype=np.int))
        #     args  =  []
        #     for jj in range(len(tags)):
        #         args.append([tags[jj], nucs_er[z, :, :], nucs_prj])
        #
        #     pool     =  multiprocessing.Pool()
        #     results  =  pool.map(f_prj, args)
        #     pool.close()
        #
        #     nucs_prj  =  label(sum(results))
        #     print np.sign(nucs_prj).sum()
        #
        # deltaT2 = time.time() - t2

        # for z in range(1, zlen):
        #     print z
        #     idxs  =  np.unique(nucs_er[z, :, :])[1:]
        #     for k in idxs:
        #         smpl  =  ((nucs_er[z, :, :] == k) * nucs_prj).reshape((xlen * ylen))
        #
        #         if smpl.sum() == 0:
        #             nucs_prj  +=  new_idx * (nucs_er[z, :, :] == k)
        #             new_idx   +=  1
        #
        #         else:
        #             smpl      =   np.delete(smpl, np.where(smpl == 0)[0])
        #             smpl      =   np.median(smpl).astype(int)
        #             nucs_prj  +=  smpl * (nucs_er[z, :, :] == k) * (1 - np.sign(nucs_prj))


        nucs_prj  =  label(nucs_prj)


        def f_dil(inds):
            msk  =  np.zeros(nucs_dil.shape, dtype=np.int)
            for i in inds:
                msk  +=  i * binary_dilation(nucs_dil == i, iterations=4) * (1 - np.sign(msk))
            return msk


        cpu_owe   = multiprocessing.cpu_count()
        # t1        =  time.time()
        nucs_dil  =  np.copy(nucs_prj)
        idxs      =  np.unique(nucs_dil)[1:]
        args      =  np.split(idxs, (idxs.size / cpu_owe + 1) * np.arange(1, cpu_owe, dtype=np.int))
        while (nucs_dil == 0).sum() > 0:
            print (nucs_dil == 0).sum()
            pool     =  multiprocessing.Pool()
            results  =  pool.map(f_dil, args)

            pool.close()

            for j in range(len(results)):
                nucs_dil  +=  results[j] * (1 - np.sign(nucs_dil))

            nucs_dil  =  label(nucs_dil)

        # deltaT1 = time.time() - t1

        # nucs_dil  =  np.copy(nucs_prj)
        # t2 =  time.time()
        # idxs = np.unique(nucs_dil)[1:]
        # while (nucs_dil == 0).sum() > 0:
        #     print (nucs_dil == 0).sum()
        #
        #     for k in idxs:
        #         # print k
        #         nuc_bff    =   binary_dilation(nucs_dil == k, iterations=4)
        #         nucs_dil  *=  (1 - (nucs_dil == k))
        #         nucs_dil  +=  k * nuc_bff * (1 - np.sign(nucs_dil))
        #
        # nucs_dil  =  label(nucs_dil)
        # deltaT2 = time.time() - t2



        ref_spts  =  np.zeros(spts_ctrs.shape[1])
        for k in range(spts_ctrs.shape[1]):
            ref_spts[k]  =  nucs_dil[spts_ctrs[1, k], spts_ctrs[0, k]]

        fls_clrd  =  np.zeros(np.append(nucs_prj.shape, 3))
        for j in range(ref_spts.size):
            fls_clrd[:, :, 0]  +=  (nucs_prj == ref_spts[j]) * 1

        fls_clrd[:, :, 2]  =  np.sign(nucs_prj) * (1 - np.sign(fls_clrd[:, :, 0])) * fls_clrd[:, :, 0].max()
        fls_clrd[:, :, 1]  =  spts.sum(0) * fls_clrd[:, :, 0].max()

        self.nucs_lbls  =  nucs_lbls
        self.nucs_dil       =  nucs_dil
        self.fls_clrd       =  fls_clrd




# class NucleiArea2:
#     def __init__(self, nucs, spts_ctrs, spts):
#
#         val1        =  threshold_otsu(nucs)
#         nucs_t      =  nucs > val1
#         nucs_t_f    =  gaussian(nucs_t, 3)
#         nucs_t_f_s  =  nucs_t_f.sum(0)
#         val2        =  threshold_otsu(nucs_t_f_s)
#         nucs_lbls   =  label(nucs_t_f_s > val2, connectivity=1)
#         nucs_lbls   =  remove_small_objects(nucs_lbls, 80)
#
#         rgp    =  regionprops(nucs_lbls)
#         areas  =  np.zeros(len(rgp))
#         for i in range(areas.size):
#             areas[i]  =  rgp[i]['area']
#
#         areas_av  =  areas.mean()
#         j         =  np.where(areas > 2*areas_av)[0]
#         rad_av    =  np.round(np.sqrt(areas_av)).astype(int)
#
#         if len(j) > 0:
#             j  =  j[0]
#             smpl        =  nucs_lbls == rgp[j]['label']
#             distance    =  distance_transform_edt(smpl)
#             local_maxi  =  peak_local_max(distance, indices=False, footprint=np.ones((rad_av, rad_av)))
#             markers     =  label(local_maxi)
#             lab_water   =  watershed(-distance, markers, mask=smpl)
#
#             nucs_lbls  =  nucs_lbls * (1 - np.sign(smpl*1)) + lab_water + np.sign(lab_water) * 1000
#             nucs_lbls  =  label(nucs_lbls)
#
#
#         idxs  =  np.unique(nucs_lbls)[1:]
#
#         nucs_dil  =  np.copy(nucs_lbls)
#         while (nucs_dil == 0).sum() > 0:
#             # print (nucs_dil == 0).sum()
#             for k in idxs:
#                 nuc_bff   =   binary_dilation(nucs_dil == k)
#                 nucs_dil  *=  (1 - (nucs_dil == k))
#                 nucs_dil  +=  k * nuc_bff * (1 - np.sign(nucs_dil))
#
#         # nucs_dil  =  label(nucs_dil)
#
#         ref_spts  =  np.zeros(spts_ctrs.shape[1])
#         for k in range(spts_ctrs.shape[1]):
#             ref_spts[k]  =  nucs_dil[spts_ctrs[1, k], spts_ctrs[0, k]]
#
#         fls_clrd  =  np.zeros(np.append(nucs_lbls.shape, 3))
#         for j in range(ref_spts.size):
#             fls_clrd[:, :, 0]  +=  (nucs_lbls == ref_spts[j]) * 1
#
#         fls_clrd[:, :, 2]  =  np.sign(nucs_lbls) * (1 - np.sign(fls_clrd[:, :, 0])) * fls_clrd[:, :, 0].max()
#         fls_clrd[:, :, 1]  =  spts.sum(0) * fls_clrd[:, :, 0].max()
#
#         self.nucs_lbls  =  nucs_lbls
#         self.nucs_dil   =  nucs_dil
#         self.fls_clrd   =  fls_clrd
