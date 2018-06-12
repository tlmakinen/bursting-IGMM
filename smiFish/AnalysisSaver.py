"""This function writes the results of the post-processing on FQ data."""

import numpy as np
import pyqtgraph as pg
import xlwt
from skimage.measure import regionprops
from skimage.color import label2rgb
from PyQt5 import QtGui, QtWidgets

# import ProgressBar
# import RemoveBorderNuclei


class AnalysisSaver:
    def __init__(self, fwritename, nucs_dil, nucs_lbls, nuc_activation, raw_data_fname):

        mycmap       =  np.fromfile('mycmap.bin', 'uint16').reshape((10000, 3))
        nucs_dil_3c  =  label2rgb(nucs_dil, bg_label=0, bg_color=[0, 0, 0], colors=mycmap)

        w        =  pg.image(nucs_dil_3c)
        txt_pos  =  regionprops(nucs_dil)
        for t in range(len(txt_pos)):
            a  =  pg.TextItem(str(txt_pos[t]['label']), 'k')
            w.addItem(a)
            a.setPos(txt_pos[t]['centroid'][0] - 30, txt_pos[t]['centroid'][1] - 30)

        idxs        =  np.unique(nucs_dil)[1:]
        info        =  np.zeros((idxs.size, 4), dtype=np.int)                   # for each nucleus we have tag, number of spots, region volume, nucleus volume
        info[:, 0]  =  idxs
        pbar        =  ProgressBar(total1=idxs.size)
        pbar.show()
        j  =  0
        for k in range(idxs.size):
            pbar.update_progressbar(j)
            j           +=  1
            msk         =  (nucs_dil == idxs[k])
            ref         =  msk * nuc_activation
            ref         =  ref.reshape(ref.size)
            ref         =  np.delete(ref, np.where(ref == 0))
            if ref.size > 0:
                info[k - 1, 1]  =  np.median(ref)
                info[k - 1, 2]  =  msk.sum() * nucs_lbls.shape[0]                   # region volume
                info[k - 1, 3]  =  np.sum(msk * np.sign(nucs_lbls))                 # nucleus volume

        book    =  xlwt.Workbook(encoding='utf-8')
        sheet1  =  book.add_sheet("Sheet 1")

        sheet1.write(0, 0, "Nuc_id")
        sheet1.write(0, 1, "Numb of Spts")
        sheet1.write(0, 2, "Region Volume")
        sheet1.write(0, 3, "Nucleus Volume")
        sheet1.write(0, 4, "X coord")
        sheet1.write(0, 5, "Y coord")

        for k in range(idxs.size):
            sheet1.write(k + 1, 0, info[k - 1, 0])
            sheet1.write(k + 1, 1, info[k - 1, 1])
            sheet1.write(k + 1, 2, info[k - 1, 2])
            sheet1.write(k + 1, 3, info[k - 1, 3])
            sheet1.write(k + 1, 4, txt_pos[k - 1]['centroid'][0])
            sheet1.write(k + 1, 5, txt_pos[k - 1]['centroid'][1])

        sheet1.write(k + 3, 0, raw_data_fname)

        book.save(fwritename + '/journal.xls')

        pbar.close()


class ProgressBar(QtGui.QWidget):
    def __init__(self, parent=None, total1=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line1 = QtGui.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        main_layout  =  QtGui.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()
