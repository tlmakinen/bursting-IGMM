"""This function builds a random spots in 3D.

Given the coordinate of the center of mass and its volume, this function will
aggregate pixels in 3D at random ti build a fake spot.
"""


import numpy as np
import random
from itertools import product


class SpotsBuilder:
    def __init__(self, vol_new):

        sp_rel_coord  =  []

        for k in range(vol_new.size):

            edge  =  np.fix(vol_new[k]**(1 / 3.0) + .99999)                                         # identify the cubic volum in which the spot must me contained (the 0.99999 drives the aprroximation)

            if edge == 2:                                                                        # p gives the possible random movements to perform from the center for each of the 3 dimension
                p  =  [0, 1]
            if edge == 3:
                p  =  [-1, 0, 1]
            if edge == 4:
                p  =  [-1, 0, 1, 2]

            a  =  set(product(set(p), repeat=3))
            a  =  list(a)
            a.remove(tuple(np.zeros(len(a[0]), dtype=int)))
            random.shuffle(a)                                                                     # list of random movements from the center inside inside the cube

            sp_rel_coord.append(np.array(a[:vol_new[k] - 1]))                                              # random movements needed to reach the given volume

        self.sp_rel_coord  =  sp_rel_coord
