"""This function finds the coordinates of the extrema of a roi segment.

It returns the coordinates of the end points of a roi segment and the
array_region of the same segment over the image.
"""

import numpy as np


class EndPtsArRegFromROI:

    def __init__(self, roi, matrix):  # , imgitm):
        self.roi     =  roi
        self.matrix  =  matrix
        # self.imgitm  =  imgitm

        # ar_reg  =  roi.getArrayRegion(self.matrix, self.imgitm)

        end_pts  =  np.zeros((2, 2))
        coords   =  roi.parentBounds().getCoords()
        delta_x  =  roi.getHandles()[1].x() - roi.getHandles()[0].x()
        delta_y  =  roi.getHandles()[1].y() - roi.getHandles()[0].y()

        if delta_x * delta_y < 0:
            end_pts[0, :]  =  coords[0], coords[3]
            end_pts[1, :]  =  coords[2], coords[1]
        else:
            end_pts[0, :]  =  coords[0], coords[1]
            end_pts[1, :]  =  coords[2], coords[3]

        end_pts  =  np.round(end_pts).astype(np.int)
        end_pts  =  end_pts * (end_pts > 0)                                                                                                     # this puts a zero where it is negative
        constr   =  np.array([[self.matrix.shape[0] - 1, self.matrix.shape[1] - 1], [self.matrix.shape[0] - 1, self.matrix.shape[1] - 1]])      # matrix of the constrain, maximum length in x and y
        end_pts  =  end_pts * (end_pts < constr)  +  constr * (end_pts > constr)                                                                # where a points exceeds the constrain, it is replaced by the maximum value it can have

        # self.ar_reg   =  ar_reg
        self.end_pts  =  end_pts
