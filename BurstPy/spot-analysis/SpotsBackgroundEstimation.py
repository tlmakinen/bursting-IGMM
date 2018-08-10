"""This function removes the background from spots traces.

Background value is different for each spot. It is extimated
averaging intensity values of pixels around the spots itself.
x, y coordinates are taken from the spots_trackes matrix (T, X, Y)
and the z coordinates is given by the maximum intensity in the row.
A 'cage' is built around the spots and intesities of the cage-pixels
are averaged as background value.
"""


import numpy as np
from skimage.measure import regionprops



class SpotsBackgroundEstimation:
    def __init__(self, spots_tracked, green4d, spt_id):

        z_tot, x_tot, y_tot  =  green4d[0, :, :, :].shape
        spt_mtx              =  (spots_tracked == spt_id)                  # the spots is isolared from the matrix (T, X, Y)
        spt_prof             =  spt_mtx.sum(2).sum(1)                      # its intensity profile             
        # t_prof               =  np.where(spt_prof != 0)[0]                 # list of t-frames in which is on 
        bkg_values           =  np.zeros(spt_prof.size)                     

        for t in range(bkg_values.size):                                                               ###  THIS NEEDS CONTROLS, THE SAMPLING RECTANGLE CAN GO OUT OF THE MATRIX
            if spt_prof[t] != 0:                                                        # check the spot is present in the frame
                x_ctr, y_ctr  =  regionprops(spt_mtx[t].astype(np.int))[0]['centroid']
                x_ctr, y_ctr  =  int(round(x_ctr)), int(round(y_ctr))
                z_ctr         =  np.argmax(green4d[t, :, x_ctr, y_ctr])              # z center coordinate

                bkg_vol    =  np.zeros(green4d[0, :, :, :].shape)
                z_min_ext  =  np.max([z_ctr - 6, 0])                                        # edges of the cage, internal and external. Controls for spots close to the borders
                z_max_ext  =  np.min([z_ctr + 7, z_tot])
                z_min_int  =  np.max([z_ctr - 4, 0])
                z_max_int  =  np.min([z_ctr + 5, z_tot])

                x_min_ext  =  np.max([x_ctr - 6, 0])
                x_max_ext  =  np.min([x_ctr + 7, x_tot])
                x_min_int  =  np.max([x_ctr - 4, 0])
                x_max_int  =  np.min([x_ctr + 5, x_tot])

                y_min_ext  =  np.max([y_ctr - 6, 0])
                y_max_ext  =  np.min([y_ctr + 7, y_tot])
                y_min_int  =  np.max([y_ctr - 4, 0])
                y_max_int  =  np.min([y_ctr + 5, y_tot])
                
                bkg_vol[z_min_ext:z_max_ext, x_min_ext:x_max_ext, y_min_ext:y_max_ext]  =  1            # cage definition (contains more than 1400 pixels)
                bkg_vol[z_min_int:z_max_int, x_min_int:x_max_int, y_min_int:y_max_int]  =  0
                bkg_bffr                                                                =  green4d[t, :, :, :] * bkg_vol
                bkg_values[t]                                                           =  bkg_bffr[bkg_bffr != 0].mean()

        self.bkg_values  =  bkg_values






# ####  HOW TO VISUALLY CHECK  #####

# green4d  =  np.fromfile(foldername + '/raw_data_green4D.bin', 'uint16')
# green4d  =  green4d[4:].reshape((green4d[3], green4d[0], green4d[2], green4d[1]))

# spots_tracked  =  np.fromfile(foldername + '/spots_tracked.bin', 'uint16')
# spots_tracked  =  spots_tracked[3:].reshape((spots_tracked[2], spots_tracked[1], spots_tracked[0]))

# spots_msk4d  =  SpotsDetection3D.SpotsDetection3D([green4d[60:62, :, :, :], 7, 4])    ## only 2 frames to do something less computational demandind, but you can run it on all

# spt_id  =  121              # for example

# z_tot, x_tot, y_tot  =  green4d[0, :, :, :].shape
# spt_mtx              =  (spots_tracked == spt_id)
# spt_prof             =  spt_mtx.sum(2).sum(1)
# t_prof               =  np.where(spt_prof != 0)[0]

# t  =  np.where(t_prof == 60)[0][0]    # beacause we choose 60:62 in spots_msk4d

# x_ctr, y_ctr  =  regionprops(spt_mtx[t_prof[t]].astype(np.int))[0]['centroid']
# x_ctr, y_ctr  =  int(round(x_ctr)), int(round(y_ctr))
# z_ctr         =  np.argmax(green4d[t_prof[t], :, x_ctr, y_ctr])

# bkg_vol                                                                 =  np.zeros(green4d[0, :, :, :].shape)
# z_min_ext  =  np.max([z_ctr - 6, 0])
# z_max_ext  =  np.min([z_ctr + 7, z_tot])
# z_min_int  =  np.max([z_ctr - 4, 0])
# z_max_int  =  np.min([z_ctr + 5, z_tot])

# x_min_ext  =  np.max([x_ctr - 6, 0])
# x_max_ext  =  np.min([x_ctr + 7, x_tot])
# x_min_int  =  np.max([x_ctr - 4, 0])
# x_max_int  =  np.min([x_ctr + 5, x_tot])

# y_min_ext  =  np.max([y_ctr - 6, 0])
# y_max_ext  =  np.min([y_ctr + 7, y_tot])
# y_min_int  =  np.max([y_ctr - 4, 0])
# y_max_int  =  np.min([y_ctr + 5, y_tot])
            
# bkg_vol[z_min_ext:z_max_ext, x_min_ext:x_max_ext, y_min_ext:y_max_ext]  =  1
# bkg_vol[z_min_int:z_max_int, x_min_int:x_max_int, y_min_int:y_max_int]  =  0

# pg.image(np.sign(spots_msk4d.spts_msk4d[0, :, :])*5 + bkg_vol)