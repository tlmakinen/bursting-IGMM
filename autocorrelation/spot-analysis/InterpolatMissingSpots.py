"""This function replaces missing spots.

Given a spot that disappears for some frames, this function takes its coordinate
just before and just after disappearing and interpolates its position(s) in the
missing frame(s). Than takes the value of the green raw data in these points,
interpolating the supposed volume too, and replaces the missing information in
this way.
"""


import numpy as np
from skimage.morphology import label
from skimage.measure import regionprops
# import pyqtgraph as pg

import SpotsBuilder


class InterpolatMissingSpots:
    def __init__(self, green4D, spots_tracked_3D, spots_3D, ref_point, numb_zrs):
        
        steps      =  spots_tracked_3D.shape[0]    
        prof_int   =  (spots_3D.spots_ints * (spots_tracked_3D == ref_point)).sum(2).sum(1)                                                     # intensity time series of the spot
        prof_int2  =  (spots_3D.spots_ints * (spots_tracked_3D == ref_point)).sum(2).sum(1)                                                     # intensity time series of the spot
        prof_bin   =  np.sign((spots_tracked_3D == ref_point).sum(2).sum(1))                                                                     # binary time series of the spots presence (1 it is there, 0 it is not)
        prof_lbls  =  label(1 - prof_bin)                                                                                                # labels of the zeros; time series ending with ..., 1, 0, 0 give problems, so we check that the last label is not at the end of the trace. the '-1' compensate for the fact that the first index is 0 

        if np.where(prof_lbls == prof_lbls.max())[0][-1] == steps - 1:                                                                          
            prof_lbls  *=  (1 - (prof_lbls == prof_lbls.max()))


        frms2work  =  []

        for l in range(1, prof_lbls.max() + 1):
            if np.sum(prof_lbls == l) <= numb_zrs:                                                                                               # select zeros segments smaller than a threshold (numb_zrs)
                frms2work.append(np.where(prof_lbls == l))                                                                                         # time coordinates of these points

        for k in range(len(frms2work)):                                                                                                          # segment by segment
            bfr_tzxy     =  np.zeros((4), dtype=np.int)
            aft_tzxy     =  np.zeros((4), dtype=np.int)
            bfr_tzxy[0]  =  (frms2work[k][0].min() - 1).astype(np.int)                                                                                                             # time before disappearing
            aft_tzxy[0]  =  (frms2work[k][0].max() + 1).astype(np.int)                                                                                                             # time after disappearing (when it just reappear)

            bfr_rgp  =  regionprops((spots_tracked_3D[bfr_tzxy[0], :, :] == ref_point).astype(np.int))
            aft_rgp  =  regionprops((spots_tracked_3D[aft_tzxy[0], :, :] == ref_point).astype(np.int))

            bfr_tzxy[2:]  =  np.round(bfr_rgp[0]['centroid'])                                                                                    # txy coordinate of the spot before disappearing
            aft_tzxy[2:]  =  np.round(aft_rgp[0]['centroid'])                                                                                    # txy coordinate of the spot after disappearing

            bfr_loc  =  np.argmin(((spots_3D.spots_tzxy[:, np.array([0, 2, 3])] - bfr_tzxy[np.array([0, 2, 3])])**2).sum(1))                     # locate in tzxy coordinate matrix
            aft_loc  =  np.argmin(((spots_3D.spots_tzxy[:, np.array([0, 2, 3])] - aft_tzxy[np.array([0, 2, 3])])**2).sum(1))                     # for some reason there is not perfect matching between the xy-coordinate of center of mass detected in 3D and the xy-coordinates of the same object projected in xy

            bfr_tzxy[1]  =  spots_3D.spots_tzxy[bfr_loc, 1]                                                                                                                        # z coordinate
            aft_tzxy[1]  =  spots_3D.spots_tzxy[aft_loc, 1]                                                                                                                        # z coordinate

            z_new  =  np.round(np.linspace(bfr_tzxy[1], aft_tzxy[1], 2 + frms2work[k][0].size)).astype(np.int)                                      # interpolation over the z, x, y positions the built spot will have
            z_new  =  z_new[1:-1]                                                                                                                   # interpolations are linear

            x_new  =  np.round(np.linspace(bfr_tzxy[2], aft_tzxy[2], 2 + frms2work[k][0].size)).astype(np.int)
            x_new  =  x_new[1:-1]

            y_new  =  np.round(np.linspace(bfr_tzxy[3], aft_tzxy[3], 2 + frms2work[k][0].size)).astype(np.int)
            y_new  =  y_new[1:-1]

            bfr_vol  =  (spots_3D.spots_vol[bfr_tzxy[0], :, :] * (spots_tracked_3D[bfr_tzxy[0], :, :] == ref_point)).sum()
            aft_vol  =  (spots_3D.spots_vol[aft_tzxy[0], :, :] * (spots_tracked_3D[aft_tzxy[0], :, :] == ref_point)).sum()
            vol_new  =  np.round(np.linspace(bfr_vol, aft_vol, 2 + frms2work[k][0].size)).astype(np.int)                                            # volume interpolation
            vol_new  =  vol_new[1:-1]

            sp_rel_coord  =  SpotsBuilder.SpotsBuilder(vol_new).sp_rel_coord                                                                        # generation of the coordinates of picels of the built spot

            for i in range(len(sp_rel_coord)):
                # ant[frms2work[k][0][i], z_new[i], x_new[i], y_new[i]]  =  2
                prof_int[frms2work[k][0][i]]  +=  green4D[frms2work[k][0][i], z_new[i], x_new[i], y_new[i]].sum()                                   # the time series of the spot is filled in the missing steps with the values of the raw data in the pixels of the built pixels
                for j in range(sp_rel_coord[i].shape[0]):
                    zz  =  z_new[i] + sp_rel_coord[i][j, 0]
                    xx  =  x_new[i] + sp_rel_coord[i][j, 1]
                    yy  =  y_new[i] + sp_rel_coord[i][j, 2]
                    # print(xx, zz, yy, frms2work[0][0][i], i)
                    # ant[frms2work[k][0][i], zz, xx, yy]  =  2
                    try:                                                                                                                            # built spots can go out of the matrix (very rare). In this case we juts skip the pixel that goes out
                        prof_int[frms2work[k][0][i]]  +=  green4D[frms2work[k][0][i], zz, xx, yy].sum()
                        spots_tracked_3D[frms2work[k][0][i], xx, yy]  =  ref_point
                    except IndexError:
                        pass

        self.corrected_trace  =  prof_int  
        self.original_trace   =  prof_int2
        self.spots_tracked_3D  =  spots_tracked_3D

        
        # w = pg.plot(prof_int2, symbol='o')
        # w.plot(prof_int, symbol='x', pen='r')
        # pg.plot(prof_int, symbol='x', pen='r')
