"""This function provide post-processing to FIshQuant software.

This function takes the excel file generate from output data of FQ-software and gives as output
the matrix of the centre of mass coordinate of the detected mRNA spots.
"""

import numpy as np
import pyqtgraph as pg
from xlrd import open_workbook


class Exc2Spots3D:
    # def __init__(self, xls_data_fname, zlen):
    def __init__(self, xls_data_fname):
        print("Ciao")
        wb     =  open_workbook(xls_data_fname)
        s_wb   =  wb.sheets()[0]
        nrows  =  s_wb.nrows

        ref_start  =  []
        ref_end    =  []

        for i in range(nrows):
            print(i)
            if s_wb.col(0)[i].value == 'SPOTS_START':                           # finds the xls coordinate of the first position value (they can be more than 1 in case you worked with several FQ cells)
                ref_start  =  np.append(ref_start, i + 2)

            if s_wb.col(0)[i].value == 'SPOTS_END':                             # finds the xls coordinate of the last position value (they can be more than 1 in case you worked with several FQ cells)
                ref_end  =  np.append(ref_end, i)

        for j in range(nrows):
            # print(j)
            if s_wb.col(0)[j].value == "Pix-XY":                                # finds pix_xy size position
                break

        pix_size_xy  =  s_wb.col(0)[j + 1].value
        pix_size_z   =  s_wb.col(1)[j + 1].value

        ref_start   =  ref_start.astype(int)
        ref_end     =  ref_end.astype(int)
        n_tot_spts  =  0
        for num_c in range(ref_start.size):                                     # total number of detected spots, in all the cells if there are more than one
            n_tot_spts  +=  ref_end[num_c] - ref_start[num_c]

        spts_ctrs  =  np.zeros((4, n_tot_spts))
        count      =  0
        for num_c in range(ref_start.size):
            # print(num_c)
            for k in range(ref_start[num_c], ref_end[num_c]):
                spts_ctrs[:, count]  =   s_wb.col(1)[k].value, s_wb.col(0)[k].value, s_wb.col(2)[k].value, s_wb.col(3)[k].value          # x, y and z in this order
                count                +=  1

        ampl_av  =  spts_ctrs[3, :].mean()
        j        =  np.where(spts_ctrs[3, :] > 2 * ampl_av)

        w  =  pg.plot(spts_ctrs[3, :], symbol='x')
        w.plot(ampl_av * 2 * np.ones(spts_ctrs[3, :].shape), pen='r')

        # if j[0].size < .1 * spts_ctrs[3, :].size:
        #     spts_ctrs = np.delete(spts_ctrs, j[0], axis=1)

        spts_ctrs         =  spts_ctrs[:3, :]
        spts_ctrs[:2, :]  =  np.round(spts_ctrs[:2, :] / pix_size_xy).astype(int)
        spts_ctrs[2, :]   =  np.round(spts_ctrs[2, :] / (pix_size_z)).astype(int)

        self.spts_ctrs  =  spts_ctrs.astype(np.int)
        self.ampl_av    =  ampl_av


# class Exc2Spots2D:
#     def __init__(self, xls_data_fname):
#
#         wb    =  open_workbook(xls_data_fname)
#         s_wb  =  wb.sheets()[0]
#         nrows  =  s_wb.nrows
#
#         ref_start  =  []
#         ref_end    =  []
#
#         for i in range(nrows):
#             if s_wb.col(0)[i].value == 'SPOTS_START':
#                 ref_start  =  np.append(ref_start, i + 2)
#
#             if s_wb.col(0)[i].value == 'SPOTS_END':
#                 ref_end  =  np.append(ref_end, i)
#
#         for j in range(nrows):
#             if s_wb.col(0)[j].value == "Pix-XY":
#                 break
#
#         pix_size_xy  =  s_wb.col(0)[j + 1].value
#
#         ref_start   =  ref_start.astype(int)
#         ref_end     =  ref_end.astype(int)
#         n_tot_spts  =  0
#         for num_c in range(ref_start.size):
#             n_tot_spts  +=  ref_end[num_c] - ref_start[num_c]
#
#         spts_ctrs  =  np.zeros((2, n_tot_spts))
#         count      =  0
#         for num_c in range(ref_start.size):
#             for k in range(ref_start[num_c], ref_end[num_c]):
#                 spts_ctrs[:, count]  =   s_wb.col(1)[k].value, s_wb.col(0)[k].value                                      # x and y, in this order
#                 count                +=  1
#
#         spts_ctrs  =  np.round(spts_ctrs / pix_size_xy).astype(int)
#
#         self.spts_ctrs  =  spts_ctrs


# class FlagSpotsDetection(QtGui.QDialog):
#     def __init__(self, parent=None):
#         super(FlagSpotsDetection, self).__init__(parent)
#
#         title_lbl  =  QtGui.QLabel("Do you want to delete that spots?", self)
#         title_lbl.setFixedSize(300, 25)
#
#         dothat_btn  =  QtGui.QPushButton("Yes", self)
#         dothat_btn.clicked.connect(self.dothat)
#         # cierra_btn.setToolTip('Write the merged analysis')
#         dothat_btn.setFixedSize(300, 25)
#
#         dont_dothat_btn  =  QtGui.QPushButton("No", self)
#         dont_dothat_btn.clicked.connect(self.dont_dothat)
#         # cierra_btn.setToolTip('Write the merged analysis')
#         dont_dothat_btn.setFixedSize(300, 25)
#
#         cierra_btn  =  QtGui.QPushButton("Ok", self)
#         cierra_btn.clicked.connect(self.cierra)
#         # cierra_btn.setToolTip('Write the merged analysis')
#         cierra_btn.setFixedSize(300, 25)
#
#         yes_and_no  =  QtGui.QHBoxLayout()
#         yes_and_no.addWidget(dothat_btn)
#         yes_and_no.addWidget(dont_dothat_btn)
#
#         layout  =  QtGui.QVBoxLayout()
#         layout.addWidget(title_lbl)
#         layout.addLayout(yes_and_no)
#         layout.addWidget(cierra_btn)
#
#         self.setWindowModality(Qt.ApplicationModal)
#         self.setLayout(layout)
#         self.setGeometry(300, 300, 350, 100)
#         self.setWindowTitle("Choose Algorithm")
#
#
#     def dothat(self):
#         self.dothat_var  =  1
#
#
#     def dont_dothat(self):
#         self.dothat_var  =  0
#
#
#     def cierra(self):
#         self.close()
