"""Main file of GUI to perform post processing on FishQuant data."""

import sys
import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt
from PyQt5 import QtGui, QtWidgets, QtCore
from skimage.color import label2rgb
from importlib import reload

import LoadSegment
import NucleiArea
import AnalysisSaver
import EndPtsArRegFromROI
import LabelsModify


class MainWindow(QtGui.QMainWindow):
    def __init__(self, parent=None):

        QtGui.QMainWindow.__init__(self, parent)

        widget = QtGui.QWidget(self)
        self.setCentralWidget(widget)

        loadDataAction  =  QtGui.QAction(QtGui.QIcon('exit.png'), '&Load data', self)
        loadDataAction.setShortcut('Ctrl+L')
        loadDataAction.setStatusTip('Load lsm files and the xls output of FQ')
        # loadDataAction.triggered.connect(self.load_data)

        save_action  =  QtGui.QAction(QtGui.QIcon('exit.png'), '&Save', self)
        save_action.setShortcut('Ctrl+S')
        save_action.setStatusTip('Save the analysis in a folder')
        # save_action.triggered.connect(self.save_analysis)

        menubar   =  self.menuBar()

        fileMenu  =  menubar.addMenu('&File')
        fileMenu.addAction(loadDataAction)
        fileMenu.addAction(save_action)

        tabs_tot  =  QtGui.QTabWidget()
        tab_nucs  =  QtGui.QWidget()
        tab_spts  =  QtGui.QWidget()

        tabs_tot.addTab(tab_nucs, "Nuclei")
        tabs_tot.addTab(tab_spts, "Spots")

        framepp1  =  pg.ImageView(self)
        framepp1.ui.roiBtn.hide()
        framepp1.ui.menuBtn.hide()

        bt1_btn  =  QtGui.QPushButton("N Segm", self)
        # bt1_btn.clicked.connect(self.nucs_3dsegm)
        bt1_btn.setToolTip('Segmentation of nuclei z by z')
        bt1_btn.setFixedSize(110, 25)


        framepp2  =  pg.ImageView(self)
        framepp2.ui.roiBtn.hide()
        framepp2.ui.menuBtn.hide()

        bt2_btn  =  QtGui.QPushButton("FruFru", self)
        # bt2_btn.clicked.connect(self.nucs_3dsegm)
        bt2_btn.setToolTip('Segmentation of nuclei z by z')
        bt2_btn.setFixedSize(110, 25)

        box1  =  QtGui.QHBoxLayout()
        box1.addWidget(framepp1)
        box1.addWidget(bt1_btn)

        box2  =  QtGui.QHBoxLayout()
        box2.addWidget(framepp2)
        box2.addWidget(bt2_btn)

        tab_nucs.setLayout(box1)
        tab_spts.setLayout(box2)

        layout  =  QtGui.QHBoxLayout(widget)
        layout.addWidget(tabs_tot)

        self.setGeometry(100, 100, 900, 800)
        self.setWindowTitle('smiFIsh')
        self.setWindowIcon(QtGui.QIcon('DrosophilaIcon.png'))
        self.show()


if __name__ == "__main__":

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
