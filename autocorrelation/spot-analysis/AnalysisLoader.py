"""This function loads the previously done analysis."""

import numpy as np


class RawData:
    def __init__(self, foldername):

        imarray_red    =  np.fromfile(foldername + '/raw_data_nuclei.bin', 'uint16')
        imarray_red    =  imarray_red[3:].reshape((imarray_red[2], imarray_red[1], imarray_red[0]))

        imarray_green  =  np.fromfile(foldername + '/raw_data_spots.bin', 'uint16')
        imarray_green  =  imarray_green[3:].reshape((imarray_green[2], imarray_green[1], imarray_green[0]))

        green4D  =  np.fromfile(foldername + '/raw_data_green4D.bin', 'uint16')
        green4D  =  green4D[4:].reshape((green4D[3], green4D[0], green4D[2], green4D[1]))

        self.imarray_green  =  imarray_green
        self.imarray_red    =  imarray_red
        self.green4D        =  green4D


class SpotsIntsVol:
    def __init__(self, foldername):

        spots_3D_ints  =  np.fromfile(foldername + '/spots_3D_ints.bin', 'uint16')
        spots_3D_ints  =  spots_3D_ints[3:].reshape((spots_3D_ints[2], spots_3D_ints[1], spots_3D_ints[0]))

        spots_3D_vol  =  np.fromfile(foldername + '/spots_3D_vol.bin', 'uint16')
        spots_3D_vol  =  spots_3D_vol[3:].reshape((spots_3D_vol[2], spots_3D_vol[1], spots_3D_vol[0]))

        spots_3D_tzxy  =  np.fromfile(foldername + '/spots_tzxy.bin', 'uint16')
        spots_3D_tzxy  =  spots_3D_tzxy[2:].reshape((spots_3D_tzxy[1], spots_3D_tzxy[0]))

        self.spots_ints   =  spots_3D_ints
        self.spots_vol    =  spots_3D_vol
        self.spots_tzxy   =  spots_3D_tzxy


class Features:
    def __init__(self, foldername):

        statistics_info  =  np.fromfile(foldername + '/spots_features3D.bin', 'float')
        statistics_info  =  statistics_info[2:].reshape((int(statistics_info[1]), int(statistics_info[0])))

        self.statistics_info  =  statistics_info
