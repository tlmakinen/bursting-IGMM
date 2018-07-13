import numpy as np
import tifffile
from skimage.filters import gaussian
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl


raw_data_fname  =  '/home/atrullo/Dropbox/FQ_PostAnalysis/border.lsm'
raw_data        =  tifffile.imread(raw_data_fname)
# data            =  raw_data[0, :, 1, :, :]
raw_data        =  raw_data[0, :, 1, :, :]
data            =  gaussian(raw_data, 4)


positive        =  np.log(np.clip(data, 0, data.max())**2)
negative        =  np.log(np.clip(-data, 0, -data.min())**2)
d2              =  np.empty(data.shape + (4,), dtype=np.ubyte)
d2[..., 0]      =  positive * (255./positive.max())
d2[..., 1]      =  negative * (255./negative.max())
d2[..., 2]      =  d2[..., 1]
d2[..., 3]      =  d2[..., 0]*0.3 + d2[..., 1]*0.3
d2[..., 3]      =  (d2[..., 3].astype(float) / 255.) **2 * 255
d2[:, 0, 0]     =  [255, 0, 0, 100]
d2[0, :, 0]     =  [0, 255, 0, 100]
d2[0, 0, :]     =  [0, 0, 255, 100]


app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.opts['distance'] = 200
w.show()
w.setWindowTitle('pyqtgraph example: GLVolumeItem')
v = gl.GLVolumeItem(d2)
v.translate(-50, -50, -100)
w.addItem(v)

ax = gl.GLAxisItem()
w.addItem(ax)
