"""This is a ProgressBar widget used by other functions."""
from PyQt4 import QtGui



class ProgressBar(QtGui.QWidget):
    def __init__(self, parent=None, total1=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line1 = QtGui.QLineEdit()

        self.progressbar1 = QtGui.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        main_layout = QtGui.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        self.progressbar1.setValue(val1)
