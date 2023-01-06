import os
import pandas as pd
import numpy as np
from itertools import chain
#from sklearn.preprocessing import StandardScaler
#from matplotlib import pyplot as plt, cm as cm, mlab as mlab
from matplotlib import pyplot as plt
import sys
from bighist.bighist import StratifiedTimeSeries, getNGAs
#from bighist.bighist import createGridForPC12
import math

# Two datasets Seshat datasets are used:
# Original: The 2017 PNAS article dataset
# Equinox:  The more recent Equinox dataset
#
# In both cases, what is used are the complexity characteristics (CCs). There
# are 9 CCs for the Original dataset and 8 CCs for the Equinox dataset. In
# addition to the data themseleves, the Equinox dataset is filtered to contain
# only original NGAs (there are 5 new NGAs in the Equinox dataset). There are
# thus three distinct StratifiedTimeSeries (see the bighistory Python package):
#
# (a) tsOrigData
# (b) tsEquinoxDataOrigNGAs
# (b) tsEquinoxDataAllNGAs
#
# For each StratifiedTimeSeries, the following three plots are made:
#
# (1) PC1 histogram
# (2) PC1-PC2 movemont plot
# (3) PC2 versus PC1, "windowing" by ranges of PC1 values

# Load the Seshat dataset from the 2017 PNAS article and confirm that PC1
# explains 77.2% of the variance
tsOrigData = StratifiedTimeSeries.loadSeshatPNAS2017Data()
tsOrigData.addFlowAnalysis()
Dsqr = [v**2 for v in tsOrigData.D]
PC1_var = Dsqr[0] / np.sum(Dsqr)
if not math.isclose(PC1_var, .772, abs_tol=.001):
    raise Exception('PC1 variance is not .772 to necessary tolerance')
else:
    print('The variance explained by the first PC is ' + str(PC1_var))

origNGAs = getNGAs('PNAS2017')
allNGAs = getNGAs('Equinox')

tsEquinoxDataOrigNGAs = StratifiedTimeSeries.\
    loadSeshatEquinoxData(subseriesToKeep=origNGAs)
tsEquinoxDataOrigNGAs.addFlowAnalysis()
tsEquinoxDataAllNGAs = StratifiedTimeSeries.\
    loadSeshatEquinoxData(subseriesToKeep=allNGAs)
tsEquinoxDataAllNGAs.addFlowAnalysis()

def makeMovementPlot(ts, file_stem):
    NGAs = list(set(ts.df[ts.subseriesColumn]))
    arrowAlpha = 1
    newWorldCol = (1,0,0,arrowAlpha)
    oldWorldCol = (0,0,1,arrowAlpha)
    subseriesColors = StratifiedTimeSeries.\
        buildNewOldWorldSubseriesColors(newWorldCol, oldWorldCol, NGAs)
    ts.plotPC1PC2Movements(subseriesColors=subseriesColors)
    plt.xticks(size=15)
    plt.yticks(size=15)
    plt.xlabel("PC1", size=15)
    plt.ylabel("PC2", size=15)

    plt.savefig(os.path.join('figures',\
        'PC1_PC2_movement_plot_' + file_stem + '.png' ),
        dpi=600)
# Set the figure size (set back to default below)
plt.rcParams["figure.figsize"] = (15,8)

# Call makeMovementPlot to make the three movement plots
makeMovementPlot(tsOrigData, 'original_data')
makeMovementPlot(tsEquinoxDataOrigNGAs, 'equinox_data_orig_NGAs')
makeMovementPlot(tsEquinoxDataAllNGAs, 'equinox_data_all_NGAs')

# Set the figure size back to the default
plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]