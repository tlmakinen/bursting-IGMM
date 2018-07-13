import sys
import numpy as np
import pyqtgraph as pg
# from PyQt4.QtCore import Qt
from PyQt4 import QtCore, QtGui


class MainWindow(QtGui.QMainWindow):

    def __init__(self, parent=None):

        QtGui.QMainWindow.__init__(self, parent)

        widget = QtGui.QWidget(self)
        self.setCentralWidget(widget)

        self.frame1  =  pg.ImageView(self)
        self.frame1.ui.roiBtn.hide()
        self.frame1.ui.menuBtn.hide()

        chop_btn  =  QtGui.QPushButton("Chop Stack", self)
        chop_btn.clicked.connect(self.chop)
        chop_btn.setToolTip('Chop the stack in 3 parts: before, during and after mitosis')
        chop_btn.setFixedSize(110, 17)

        layout  =  QtGui.QVBoxLayout(widget)
        layout.addWidget(self.frame1)
        layout.addWidget(chop_btn)

        self.setGeometry(100, 100, 600, 400)
        self.setWindowTitle('SegmentTrack4MemoryGUI')
        # self.setWindowIcon(QtGui.QIcon('DrosophilaIcon.png'))
        self.show()

    def chop(self):
        self.mpp = PopUp()
        self.mpp.show()
        self.mpp.procStart.connect(self.PopUp)

    @QtCore.pyqtSlot(str)
    def PopUp(self, message):
        print message
        # self.lineEdit.setText("From B: " + message)
        self.frame1.setImage(np.random.rand(50, 50))
        self.mpp.close()


class PopUp(QtGui.QWidget):
    procStart = QtCore.pyqtSignal(int)

    def __init__(self):
        QtGui.QWidget.__init__(self)

        frame1  =  pg.ImageView(self)
        frame1.ui.roiBtn.hide()
        frame1.ui.menuBtn.hide()

        chop_btn  =  QtGui.QPushButton("Chop Stack", self)
        chop_btn.clicked.connect(self.chop)
        chop_btn.setToolTip('Chop the stack in 3 parts: before, during and after mitosis')
        chop_btn.setFixedSize(110, 17)

        layout  =  QtGui.QVBoxLayout()
        layout.addWidget(frame1)
        layout.addWidget(chop_btn)

        self.setLayout(layout)
        self.setGeometry(300, 300, 600, 400)
        self.setWindowTitle("Modifier Tool")

    @QtCore.pyqtSlot()
    def chop(self):
        val  =  9
        self.procStart.emit(val)


if __name__ == "__main__":

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
