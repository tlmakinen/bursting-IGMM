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
        loadDataAction.triggered.connect(self.load_data)

        save_action  =  QtGui.QAction(QtGui.QIcon('exit.png'), '&Save', self)
        save_action.setShortcut('Ctrl+S')
        save_action.setStatusTip('Save the analysis in a folder')
        save_action.triggered.connect(self.save_analysis)

        # load_analysis_action  =  QtGui.QAction(QtGui.QIcon('exit.png'), '&Load Analysis', self)
        # load_analysis_action.setShortcut('Ctrl+P')
        # load_analysis_action.setStatusTip('Load analysis')
        # load_analysis_action.triggered.connect(self.load_analysis)

        menubar   =  self.menuBar()

        fileMenu  =  menubar.addMenu('&File')
        fileMenu.addAction(loadDataAction)
        fileMenu.addAction(save_action)
        # fileMenu.addAction(load_analysis_action)

        fname_edt = QtGui.QLineEdit(self)
        fname_edt.setToolTip('Name of the file you are working on')

        tabs  =  QtGui.QTabWidget()
        tab1  =  QtGui.QWidget()
        tab2  =  QtGui.QWidget()
        tab3  =  QtGui.QWidget()

        framepp1  =  pg.ImageView(self, name='Frame1')
        framepp1.getImageItem().mouseClickEvent  =  self.click
        framepp1.ui.roiBtn.hide()
        framepp1.ui.menuBtn.hide()
        framepp1.timeLine.sigPositionChanged.connect(self.update_frame2)

        framepp2  =  pg.ImageView(self)
        framepp2.ui.roiBtn.hide()
        framepp2.ui.menuBtn.hide()
        framepp2.view.setXLink('Frame1')
        framepp2.view.setYLink('Frame1')
        framepp2.timeLine.sigPositionChanged.connect(self.update_frame1)

        framepp3  =  pg.ImageView(self)
        framepp3.ui.roiBtn.hide()
        framepp3.ui.menuBtn.hide()
        framepp3.view.setXLink('Frame1')
        framepp3.view.setYLink('Frame1')

        frame1_box  =  QtGui.QHBoxLayout()
        frame1_box.addWidget(framepp1)

        frame2_box  =  QtGui.QHBoxLayout()
        frame2_box.addWidget(framepp2)

        frame3_box  =  QtGui.QHBoxLayout()
        frame3_box.addWidget(framepp3)

        tab1.setLayout(frame1_box)
        tab2.setLayout(frame2_box)
        tab3.setLayout(frame3_box)

        tabs.addTab(tab1, "Segmented")
        tabs.addTab(tab2, "Raw Data")
        tabs.addTab(tab3, "Nuclei Regions")

        nucs_3dsegm_btn  =  QtGui.QPushButton("N Segm", self)
        nucs_3dsegm_btn.clicked.connect(self.nucs_3dsegm)
        nucs_3dsegm_btn.setToolTip('Segmentation of nuclei z by z')
        nucs_3dsegm_btn.setFixedSize(110, 25)

        nucs_2d_dilation_btn  =  QtGui.QPushButton("N Dilation", self)
        nucs_2d_dilation_btn.clicked.connect(self.nucs_2d_dilation)
        nucs_2d_dilation_btn.setToolTip('Nuclei dilation to find nuclei regions')
        nucs_2d_dilation_btn.setFixedSize(110, 25)
        nucs_2d_dilation_btn.setEnabled(False)

        shuffle_clrs_btn  =  QtGui.QPushButton("Shuffle Clrs", self)
        shuffle_clrs_btn.clicked.connect(self.shuffle_clrs)
        shuffle_clrs_btn.setToolTip('Shuffle Colors')
        shuffle_clrs_btn.setFixedSize(110, 25)

        num_iter_lbl  =  QtGui.QLabel('Thr Rate', self)
        num_iter_lbl.setFixedSize(60, 25)

        num_iter_edt  =  QtGui.QLineEdit(self)
        num_iter_edt.textChanged[str].connect(self.num_iter_var)
        num_iter_edt.setToolTip('Sets the threshold in terms of background value; suggested value 5')
        num_iter_edt.setFixedSize(35, 25)
        num_iter_edt.setText('5')

        num_iter_lbl_edt  =  QtGui.QHBoxLayout()
        num_iter_lbl_edt.addWidget(num_iter_lbl)
        num_iter_lbl_edt.addWidget(num_iter_edt)

        man_active_toggle = QtGui.QCheckBox('Hand Cut', self)
        man_active_toggle.setFixedSize(110, 25)
        man_active_toggle.stateChanged.connect(self.man_active)

        manual_cut_btn  =  QtGui.QPushButton("Cut", self)
        manual_cut_btn.clicked.connect(self.manual_cut)
        manual_cut_btn.setToolTip('Manual cutting of nuclei (Ctrl+Suppr)')
        manual_cut_btn.setFixedSize(110, 25)
        manual_cut_btn.setEnabled(True)

        keys  =  QtGui.QVBoxLayout()
        keys.addLayout(num_iter_lbl_edt)
        keys.addWidget(nucs_3dsegm_btn)
        keys.addWidget(man_active_toggle)
        keys.addWidget(manual_cut_btn)
        keys.addWidget(nucs_2d_dilation_btn)
        keys.addStretch()
        keys.addWidget(shuffle_clrs_btn)

        frame_fname  =  QtGui.QVBoxLayout()
        frame_fname.addWidget(fname_edt)
        frame_fname.addWidget(tabs)

        layout  =  QtGui.QHBoxLayout(widget)
        layout.addLayout(frame_fname)
        layout.addLayout(keys)

        self.framepp1         =  framepp1
        self.framepp2         =  framepp2
        self.framepp3         =  framepp3
        self.fname_edt        =  fname_edt
        self.man_active_flag  =  0
        self.c_count          =  0
        self.mycmap           =  np.fromfile('mycmap.bin', 'uint16').reshape((10000, 3)) / 255.0

        self.nucs_2d_dilation_btn  =  nucs_2d_dilation_btn
        self.manual_cut_btn        =  manual_cut_btn

        self.setGeometry(100, 100, 900, 800)
        self.setWindowTitle('FQ_PostAnalysis')
        self.setWindowIcon(QtGui.QIcon('DrosophilaIcon.png'))
        self.show()


    def closeEvent(self, event):
        quit_msg  =  "Are you sure you want to exit the program?"
        reply     =  QtGui.QMessageBox.question(self, 'Message', quit_msg, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()


    def num_iter_var(self, text):
        self.num_iter_var  =  np.int(text)


    def update_frame2(self):
        self.framepp2.setCurrentIndex(self.framepp1.currentIndex)


    def update_frame1(self):
        self.framepp1.setCurrentIndex(self.framepp2.currentIndex)


    def man_active(self, state):
        if state == QtCore.Qt.Checked:
            self.man_active_flag  =  1
            self.nucs_2d_dilation_btn.setEnabled(False)
            self.manual_cut_btn.setEnabled(True)
        else:
            self.man_active_flag  =  0
            self.nucs_2d_dilation_btn.setEnabled(True)
            self.manual_cut_btn.setEnabled(False)


    def click(self, event):
        event.accept()
        pos        =  event.pos()
        modifiers  =  QtGui.QApplication.keyboardModifiers()

        if self.man_active_flag == 1:
            if modifiers  ==  QtCore.Qt.ShiftModifier:
                if self.c_count - 2 * int(self.c_count / 2) == 0:
                    self.pos1          =   pos
                else:
                    self.roi      =  pg.LineSegmentROI([self.pos1, pos], pen='r')
                    self.framepp1.addItem(self.roi)

                self.c_count  +=  1


    def manual_cut(self):
        cif           =  self.framepp1.currentIndex
        hh            =  self.framepp1.view.viewRange()
        ff            =  EndPtsArRegFromROI.EndPtsArRegFromROI(self.roi, self.nucs_3d_det.nucs_lbls[cif, :, :])
        self.end_pts  =  ff.end_pts

        self.bufframe                             =  np.copy(self.nucs_3d_det.nucs_lbls[cif, :, :])
        self.nucs_3d_det.nucs_lbls[cif, :, :]     =  LabelsModify.LabelsModify(self.nucs_3d_det.nucs_lbls[cif, :, :], self.end_pts).labels_fin
        self.nucs_3d_det.nucs_lbls_3c[cif, :, :]  =  label2rgb(self.nucs_3d_det.nucs_lbls[cif, :, :], bg_label=0, bg_color=[0, 0, 0], colors=self.mycmap)
        self.framepp1.setImage(self.nucs_3d_det.nucs_lbls_3c)
        self.framepp1.setCurrentIndex(cif)
        self.framepp1.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
        self.framepp1.view.setYRange(hh[1][0], hh[1][1], padding=.0002)
        self.framepp1.removeItem(self.roi)


    def keyPressEvent(self, event):
        if event.key() == (QtCore.Qt.ControlModifier and Qt.Key_Z):
            cif                                       =  self.framepp1.currentIndex
            hh                                        =  self.framepp1.view.viewRange()
            self.nucs_3d_det.nucs_lbls[cif, :, :]     =  self.bufframe
            self.nucs_3d_det.nucs_lbls_3c[cif, :, :]  =  label2rgb(self.nucs_3d_det.nucs_lbls[cif, :, :], bg_label=0, bg_color=[0, 0, 0], colors=self.mycmap)
            self.framepp1.setImage(self.nucs_3d_det.nucs_lbls_3c)
            self.framepp1.setCurrentIndex(cif)
            self.framepp1.view.setXRange(hh[0][0], hh[0][1], padding=.0002)
            self.framepp1.view.setYRange(hh[1][0], hh[1][1], padding=.0002)

        if event.key() == (QtCore.Qt.ShiftModifier and Qt.Key_Delete):
            self.manual_cut()


    def load_data(self):
        # self.raw_data_fname  =  str(QtGui.QFileDialog.getOpenFileName(None, "Select the lsm file to analyze", filter="*.tif"))
        self.raw_data_fname  =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select the lsm file to analyze", filter='*.lsm *.tif')[0])
        self.fname_edt.setText(self.raw_data_fname)
        self.xls_data_fname  =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select the .xls file of the analyzed lsm file", filter="*.xlsx *.xls")[0])


    def shuffle_clrs(self):
        np.random.shuffle(self.mycmap)
        self.nucs_3d_det.nucs_labls_3c              =   label2rgb(self.nucs_3d_det.nucs_lbls, bg_color=[0, 0, 0], bg_label=0, colors=self.mycmap)
        self.nucs_3d_det.nucs_labls_3c[:, :, :, 0]  *=  (1 - self.nucs_3d_det.spts_ext)
        self.nucs_3d_det.nucs_labls_3c[:, :, :, 1]  *=  (1 - self.nucs_3d_det.spts_ext)
        self.nucs_3d_det.nucs_labls_3c[:, :, :, 2]  *=  (1 - self.nucs_3d_det.spts_ext)
        self.nucs_3d_det.nucs_labls_3c[:, :, :, 1]  +=  self.nucs_3d_det.spts_ext

        cif  =  self.framepp1.currentIndex
        self.framepp1.setImage(self.nucs_3d_det.nucs_labls_3c)
        self.framepp1.setCurrentIndex(cif)

        self.framepp3.setImage(label2rgb(self.nucs_area.nucs_dil, bg_color=[0, 0, 0], bg_label=0, colors=self.mycmap))


    def nucs_3dsegm(self):
        self.nucs_3d_det  =  LoadSegment.LoadSegment3D(self.raw_data_fname, self.xls_data_fname, self.num_iter_var)
        self.framepp1.setImage(self.nucs_3d_det.nucs_lbls_3c)
        self.framepp2.setImage(self.nucs_3d_det.nucs)


    def nucs_2d_dilation(self):
        self.nucs_area  =  NucleiArea.NucleiArea(self.nucs_3d_det.nucs_lbls, self.nucs_3d_det.spts_ext_ctrs, self.nucs_3d_det.spts_ext)
        self.framepp3.setImage(label2rgb(self.nucs_area.nucs_dil, bg_color=[0, 0, 0], bg_label=0, colors=self.mycmap))


    def save_analysis(self):
        folder  =  str(QtGui.QFileDialog.getExistingDirectory(None, "Select Directory"))
        reload(AnalysisSaver)
        AnalysisSaver.AnalysisSaver(folder, self.nucs_area.nucs_dil, self.nucs_3d_det.nucs_lbls, self.nucs_area.fls_clrd[:, :, 0], self.raw_data_fname)


if __name__ == "__main__":

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
