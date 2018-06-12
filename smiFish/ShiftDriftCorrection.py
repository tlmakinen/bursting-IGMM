import numpy as np


class ShiftDriftCorrection:
    def __init__(self, stack):

        steps, xlen, ylen  =  stack.shape

        drif_r   =  np.zeros(stack.shape)
        drif_l   =  np.zeros(stack.shape)
        drif_u   =  np.zeros(stack.shape)
        drif_d   =  np.zeros(stack.shape)

        drif_ru  =  np.zeros(stack.shape)
        drif_rd  =  np.zeros(stack.shape)
        drif_lu  =  np.zeros(stack.shape)
        drif_ld  =  np.zeros(stack.shape)

        for t in range(steps):
            drif_r[t, :xlen - t*2, :]  =  stack[t, t*2:, :]
            drif_l[t, :xlen - t*2, :]  =  stack[t, :xlen - t*2, :]
            drif_u[t, :, :xlen - t*2]  =  stack[t, :, t*2:]
            drif_d[t, :, :xlen - t*2]  =  stack[t, :, :ylen - t*2]

            drif_ru[t, :xlen - t*2, :xlen - t*2]  =  stack[t, t*2:, t*2:]
            drif_rd[t, :xlen - t*2, :xlen - t*2]  =  stack[t, t*2:, :ylen - t*2]
            drif_lu[t, :xlen - t*2, :xlen - t*2]  =  stack[t, :xlen - t*2, t*2:]
            drif_ld[t, :xlen - t*2, :xlen - t*2]  =  stack[t, :xlen - t*2, :ylen - t*2]


        corr_r   =  np.zeros((steps-1))
        corr_l   =  np.zeros((steps-1))
        corr_u   =  np.zeros((steps-1))
        corr_d   =  np.zeros((steps-1))

        corr_ru  =  np.zeros((steps-1))
        corr_rd  =  np.zeros((steps-1))
        corr_lu  =  np.zeros((steps-1))
        corr_ld  =  np.zeros((steps-1))

        for t in range(steps-1):
            corr_r[t]  =  np.sum(np.abs(drif_r[t, :, :] - drif_r[t+1, :, :]))
            corr_l[t]  =  np.sum(np.abs(drif_l[t, :, :] - drif_l[t+1, :, :]))
            corr_u[t]  =  np.sum(np.abs(drif_u[t, :, :] - drif_u[t+1, :, :]))
            corr_d[t]  =  np.sum(np.abs(drif_d[t, :, :] - drif_d[t+1, :, :]))

            corr_ru[t]  =  np.sum(np.abs(drif_ru[t, :, :] - drif_ru[t+1, :, :]))
            corr_rd[t]  =  np.sum(np.abs(drif_rd[t, :, :] - drif_rd[t+1, :, :]))
            corr_lu[t]  =  np.sum(np.abs(drif_lu[t, :, :] - drif_lu[t+1, :, :]))
            corr_ld[t]  =  np.sum(np.abs(drif_ld[t, :, :] - drif_ld[t+1, :, :]))


