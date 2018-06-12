"""Main file of GUI to perform post processing on FishQuant data."""

import sys
import numpy as np
import pyqtgraph as pg
import tifffile
from PyQt5.QtCore import Qt
from PyQt5 import QtGui, QtWidgets, QtCore
from skimage.color import label2rgb
from importlib import reload

import LoadSegment
import NucleiArea
import AnalysisSaver
# import EndPtsArRegFromROI
import LabelsModify
import SpotsDetection3D_Recursive


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

        exitAction  =  QtGui.QAction(QtGui.QIcon('exit.png'), '&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)

        menubar   =  self.menuBar()

        fileMenu  =  menubar.addMenu('&File')
        fileMenu.addAction(loadDataAction)
        fileMenu.addAction(save_action)
        fileMenu.addAction(exitAction)

        fname_nucs_edt = QtGui.QLineEdit(self)
        fname_nucs_edt.setToolTip('Name of the file you are working on')

        tabs_tot  =  QtGui.QTabWidget()
        tab_nucs  =  QtGui.QWidget()
        tab_spts  =  QtGui.QWidget()

        tabs_tot.addTab(tab_nucs, "Nuclei")
        tabs_tot.addTab(tab_spts, "Spots")

        tabs_nucs  =  QtGui.QTabWidget()
        tab1_nucs  =  QtGui.QWidget()
        tab2_nucs  =  QtGui.QWidget()
        tab3_nucs  =  QtGui.QWidget()

        frame1_nucs  =  pg.ImageView(self, name='Frame1')
        frame1_nucs.getImageItem().mouseClickEvent  =  self.click
        frame1_nucs.ui.roiBtn.hide()
        frame1_nucs.ui.menuBtn.hide()
        frame1_nucs.timeLine.sigPositionChanged.connect(self.update_frame2)

        frame2_nucs  =  pg.ImageView(self)
        frame2_nucs.ui.roiBtn.hide()
        frame2_nucs.ui.menuBtn.hide()
        frame2_nucs.view.setXLink('Frame1')
        frame2_nucs.view.setYLink('Frame1')
        frame2_nucs.timeLine.sigPositionChanged.connect(self.update_frame1)

        frame3_nucs  =  pg.ImageView(self)
        frame3_nucs.ui.roiBtn.hide()
        frame3_nucs.ui.menuBtn.hide()
        frame3_nucs.view.setXLink('Frame1')
        frame3_nucs.view.setYLink('Frame1')

        frame1_box_nucs  =  QtGui.QHBoxLayout()
        frame1_box_nucs.addWidget(frame1_nucs)

        frame2_box_nucs  =  QtGui.QHBoxLayout()
        frame2_box_nucs.addWidget(frame2_nucs)

        frame3_box_nucs  =  QtGui.QHBoxLayout()
        frame3_box_nucs.addWidget(frame3_nucs)

        tab1_nucs.setLayout(frame1_box_nucs)
        tab2_nucs.setLayout(frame2_box_nucs)
        tab3_nucs.setLayout(frame3_box_nucs)

        tabs_nucs.addTab(tab1_nucs, "Segmented")
        tabs_nucs.addTab(tab2_nucs, "Raw Data")
        tabs_nucs.addTab(tab3_nucs, "Nuclei Regions")

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
        frame_fname.addWidget(fname_nucs_edt)
        frame_fname.addWidget(tabs_nucs)

        layout_nucs  =  QtGui.QHBoxLayout()
        layout_nucs.addLayout(frame_fname)
        layout_nucs.addLayout(keys)


        fname_spts_edt = QtGui.QLineEdit(self)
        fname_spts_edt.setToolTip('Name of the file you are working on')

        frame1_spts  =  pg.ImageView(self)
        frame1_spts.ui.roiBtn.hide()
        frame1_spts.ui.menuBtn.hide()

        frame2_spts  =  pg.ImageView(self)
        frame2_spts.ui.roiBtn.hide()
        frame2_spts.ui.menuBtn.hide()

        frame3_spts  =  pg.PlotWidget(self)
        frame3_spts.setFixedSize(140, 100)

        tabs_spts  =  QtGui.QTabWidget()
        tab1_spts  =  QtGui.QWidget()
        tab2_spts  =  QtGui.QWidget()

        frame1_box_spts  =  QtGui.QHBoxLayout()
        frame1_box_spts.addWidget(frame1_spts)

        frame2_box_spts  =  QtGui.QHBoxLayout()
        frame2_box_spts.addWidget(frame2_spts)

        tab1_spts.setLayout(frame1_box_spts)
        tab2_spts.setLayout(frame2_box_spts)

        tabs_spts.addTab(tab1_spts, "Raw Data")
        tabs_spts.addTab(tab2_spts, "Segmented")

        kern_gauss_lbl  =  QtGui.QLabel('Gauss Size', self)
        kern_gauss_lbl.setFixedSize(60, 25)

        kern_gauss_edt  =  QtGui.QLineEdit(self)
        kern_gauss_edt.textChanged[str].connect(self.kern_gauss_var)
        kern_gauss_edt.setToolTip('Sets the size of the Gaussian Filter for the pre-smoothing; suggested value 2')
        kern_gauss_edt.setFixedSize(35, 25)

        kern_gauss_lbl_edt  =  QtGui.QHBoxLayout()
        kern_gauss_lbl_edt.addWidget(kern_gauss_lbl)
        kern_gauss_lbl_edt.addWidget(kern_gauss_edt)

        thresh_lbl  =  QtGui.QLabel('Thresholds', self)
        thresh_lbl.setFixedSize(110, 25)

        min_thr_lbl  =  QtGui.QLabel('Min Thr', self)
        min_thr_lbl.setFixedSize(60, 25)

        min_thr_edt  =  QtGui.QLineEdit(self)
        min_thr_edt.textChanged[str].connect(self.min_thr_var)
        min_thr_edt.setToolTip('Sets the minimum value for the thresgolding; suggested value 5')
        min_thr_edt.setFixedSize(35, 25)

        min_thr_lbl_edt  =  QtGui.QHBoxLayout()
        min_thr_lbl_edt.addWidget(min_thr_lbl)
        min_thr_lbl_edt.addWidget(min_thr_edt)

        max_thr_lbl  =  QtGui.QLabel('Max Thr', self)
        max_thr_lbl.setFixedSize(60, 25)

        max_thr_edt  =  QtGui.QLineEdit(self)
        max_thr_edt.textChanged[str].connect(self.max_thr_var)
        max_thr_edt.setToolTip('Sets the maximum value for the thresholding; suggested value 7')
        max_thr_edt.setFixedSize(35, 25)

        max_thr_lbl_edt  =  QtGui.QHBoxLayout()
        max_thr_lbl_edt.addWidget(max_thr_lbl)
        max_thr_lbl_edt.addWidget(max_thr_edt)

        steps_thr_lbl  =  QtGui.QLabel('Steps', self)
        steps_thr_lbl.setFixedSize(60, 25)

        steps_thr_edt  =  QtGui.QLineEdit(self)
        steps_thr_edt.textChanged[str].connect(self.steps_thr_var)
        steps_thr_edt.setToolTip('Sets the number of different thresholds between min and max to test')
        steps_thr_edt.setFixedSize(35, 25)

        steps_thr_lbl_edt  =  QtGui.QHBoxLayout()
        steps_thr_lbl_edt.addWidget(steps_thr_lbl)
        steps_thr_lbl_edt.addWidget(steps_thr_edt)

        spts_3d_segm_btn  =  QtGui.QPushButton("Thr Study", self)
        spts_3d_segm_btn.clicked.connect(self.spts_3d_segm)
        spts_3d_segm_btn.setToolTip('Segmentation study over several thresholds')
        spts_3d_segm_btn.setFixedSize(110, 25)

        sld_thr  =  QtGui.QScrollBar(QtCore.Qt.Horizontal, self)
        sld_thr.valueChanged.connect(self.sld_thr_update)

        keys  =  QtGui.QVBoxLayout()
        keys.addLayout(kern_gauss_lbl_edt)
        keys.addStretch()
        keys.addWidget(thresh_lbl)
        keys.addLayout(min_thr_lbl_edt)
        keys.addLayout(max_thr_lbl_edt)
        keys.addLayout(steps_thr_lbl_edt)
        keys.addWidget(spts_3d_segm_btn)
        keys.addStretch()
        keys.addWidget(frame3_spts)
        keys.addWidget(sld_thr)
        keys.addStretch()

        tabs_spts_fname  =  QtGui.QVBoxLayout()
        tabs_spts_fname.addWidget(fname_spts_edt)
        tabs_spts_fname.addWidget(tabs_spts)

        layout_spts  =  QtGui.QHBoxLayout()
        layout_spts.addLayout(tabs_spts_fname)
        layout_spts.addLayout(keys)

        tab_nucs.setLayout(layout_nucs)
        tab_spts.setLayout(layout_spts)

        layout  =  QtGui.QHBoxLayout(widget)
        layout.addWidget(tabs_tot)

        self.frame1_nucs      =  frame1_nucs
        self.frame2_nucs      =  frame2_nucs
        self.frame3_nucs      =  frame3_nucs
        self.fname_nucs_edt   =  fname_nucs_edt
        self.fname_spts_edt   =  fname_spts_edt
        self.frame1_spts      =  frame1_spts
        self.frame2_spts      =  frame2_spts
        self.frame3_spts      =  frame3_spts
        self.sld_thr          =  sld_thr
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
        self.num_iter_value  =  np.int(text)


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
        # ff            =  EndPtsArRegFromROI.EndPtsArRegFromROI(self.roi, self.nucs_3d_det.nucs_lbls[cif, :, :])
        # self.end_pts  =  ff.end_pts
        pp       =  self.roi.getHandles()
        pp       =  [self.roi.mapToItem(self.framepp1.imageItem, p.pos()) for p in pp]
        end_pts  =  np.array([[int(pp[0].x()), int(pp[0].y())], [int(pp[1].x()), int(pp[1].y())]])

        self.bufframe                             =  np.copy(self.nucs_3d_det.nucs_lbls[cif, :, :])
        # self.nucs_3d_det.nucs_lbls[cif, :, :]     =  LabelsModify.LabelsModify(self.nucs_3d_det.nucs_lbls[cif, :, :], self.end_pts).labels_fin
        self.nucs_3d_det.nucs_lbls[cif, :, :]     =  LabelsModify.LabelsModify(self.nucs_3d_det.nucs_lbls[cif, :, :], end_pts).labels_fin
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
        self.raw_nucs_fname  =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select nuclei file to analyze", filter='*.tif')[0])
        self.raw_spts_fname  =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select spots file to analyze", filter='*.tif')[0])
        self.fname_nucs_edt.setText(self.raw_nucs_fname)
        self.fname_spts_edt.setText(self.raw_spts_fname)
        # self.xls_data_fname  =  str(QtWidgets.QFileDialog.getOpenFileName(None, "Select the .xls file of the analyzed lsm file", filter="*.xlsx *.xls")[0])
        self.spts  =  tifffile.imread(self.raw_spts_fname)
        self.spts  =  self.spts[:, 1000:1512, 1000:1512]
        self.frame1_spts.setImage(self.spts)


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
        self.nucs_3d_det  =  LoadSegment.LoadSegment3D(self.raw_nucs_fname, self.xls_data_fname, self.num_iter_value)
        self.framepp1.setImage(self.nucs_3d_det.nucs_lbls_3c)
        self.framepp2.setImage(self.nucs_3d_det.nucs)


    def nucs_2d_dilation(self):
        self.nucs_area  =  NucleiArea.NucleiArea(self.nucs_3d_det.nucs_lbls, self.nucs_3d_det.spts_ext_ctrs, self.nucs_3d_det.spts_ext)
        self.framepp3.setImage(label2rgb(self.nucs_area.nucs_dil, bg_color=[0, 0, 0], bg_label=0, colors=self.mycmap))


    def save_analysis(self):
        folder  =  str(QtGui.QFileDialog.getExistingDirectory(None, "Select Directory"))
        reload(AnalysisSaver)
        AnalysisSaver.AnalysisSaver(folder, self.nucs_area.nucs_dil, self.nucs_3d_det.nucs_lbls, self.nucs_area.fls_clrd[:, :, 0], self.raw_data_fname)


    def kern_gauss_var(self, text):
        self.kern_gauss_value  =  np.float(text)


    def min_thr_var(self, text):
        self.min_thr_value  =  np.float(text)


    def max_thr_var(self, text):
        self.max_thr_value  =  np.float(text)


    def steps_thr_var(self, text):
        self.steps_thr_var  =  np.float(text)
        self.thr_vals   =  np.linspace(self.min_thr_value, self.max_thr_value, self.steps_thr_var)
        self.sld_thr.setMaximum(self.steps_thr_var)


    def spts_3d_segm(self):
        # print(np.linspace(self.min_thr_value, self.max_thr_value, self.steps_thr_var))
        self.spts_numb  =  SpotsDetection3D_Recursive.SpotsDetection3D_Recursive(self.spts, self.kern_gauss_value, self.thr_vals).spts_num
        self.frame3_spts.plot(self.thr_vals, self.spts_numb)
        self.roi  =  pg.LineSegmentROI([[self.thr_vals[0], 0], [self.thr_vals[0], self.spts_numb.max()]], pen='b')
        self.frame3_spts.addItem(self.roi)


    def sld_thr_update(self):
        print(self.sld_thr.value())
        print(self.thr_vals[self.sld_thr.value()])
        # self.roi.setPos([self.thr_vals[self.sld_thr.value(), 0])
        self.roi.setPos([self.thr_vals[self.sld_thr.value()], 0], update=True)



if __name__ == "__main__":

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
