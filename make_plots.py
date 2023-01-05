import os
import pandas as pd
import numpy as np
from itertools import chain
#from sklearn.preprocessing import StandardScaler
#from matplotlib import pyplot as plt, cm as cm, mlab as mlab
from matplotlib import pyplot as plt
import sys
from bighist.bighist import StratifiedTimeSeries
#from bighist.bighist import createGridForPC12
import math

# Load the Seshat dataset from the 2017 PNAS article and confirm that PC1
# explains 77.2% of the variance
ts = StratifiedTimeSeries.loadSeshatPNAS2017Data()
ts.addFlowAnalysis()
Dsqr = [v**2 for v in ts.D]
PC1_var = Dsqr[0] / np.sum(Dsqr)
if not math.isclose(PC1_var, .772, abs_tol=.001):
    raise Exception('PC1 variance is not .772 to necessary tolerance')
else:
    print('The variance explained by the first PC is ' + str(PC1_var))

NGAs = list(set(ts.df[ts.subseriesColumn]))
arrowAlpha = 1
newWorldCol = (1,0,0,arrowAlpha)
oldWorldCol = (0,0,1,arrowAlpha)
subseriesColors = StratifiedTimeSeries.\
    buildNewOldWorldSubseriesColors(newWorldCol, oldWorldCol, NGAs)

ts.plotPC1PC2Movements(subseriesColors=subseriesColors)

r0 = 1.5
minPoints = 20
dGrid = .2
#u0Vect,v0Vect = createGridForPC12(dGrid, ts.velArrayOut)
velScaling = 100

# Set the figure size
plt.rcParams["figure.figsize"] = (15,8)

#fig = plt.figure(figsize=(15,8))
#plt.axis()
#plt.xlim(-4,6)
#plt.ylim(-2.5,2.5)
plt.xticks(size=15)
plt.yticks(size=15)
plt.xlabel("PC1", size=15)
plt.ylabel("PC2", size=15)
plt.savefig(os.path.join('figures',\
    "pc12_movement_plot_original_data_original_NGAs.png"),
    dpi=600)
    #"pc12_movement_plot_original_data_original_NGAs.png"),, transparent=True)

# Set the figure size back to the default
plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]