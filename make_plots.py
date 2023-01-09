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
    print('The variance explained by the first PC is, correctly, ' +\
        str(PC1_var) + ' for the Seshat2017 dataset')

origNGAs = getNGAs('PNAS2017')
allNGAs = getNGAs('Equinox')

tsEquinoxDataOrigNGAs = StratifiedTimeSeries.\
    loadSeshatEquinoxData(subseriesToKeep=origNGAs)
tsEquinoxDataOrigNGAs.addFlowAnalysis()
tsEquinoxDataAllNGAs = StratifiedTimeSeries.\
    loadSeshatEquinoxData(subseriesToKeep=allNGAs)
tsEquinoxDataAllNGAs.addFlowAnalysis()

# Use the same PC1 and PC2 limits for all plots
pc1Min = -7
pc1Max = 7
pc2Min = -2.5
pc2Max = 2.5

def makeSlidingWindowPlot(ts, file_stem, pc1Min, pc1Max, pc2Min, pc2Max):
    PC1 = ts.velArrayOut[:,0,0]
    PC2 = ts.velArrayOut[:,1,0]
    PC1_vel = ts.velArrayOut[:,0,1]*100
    PC2_vel = ts.velArrayOut[:,1,1]*100
    
    window_width =1.0
    overlap = .5
    
    score_list = []
    vel_list = []
    score_std_list = []
    vel_std_list = []
    
    score_error_list = []
    vel_error_list = []
    
    center_list = []
    
    PC1_min = np.min(PC1)
    PC1_max = np.max(PC1)
    
    n_window = np.ceil((PC1_max - PC1_min - window_width)/
                       (window_width-overlap)).astype(int)
    
    for i in range(n_window):
        window = np.array([PC1_min+i*(window_width-overlap),
                           PC1_min+i*(window_width-overlap)+window_width])
        center = np.mean(window)
        loc = (window[0]<=PC1) * (PC1<window[1])
        
        PC2_in_window = PC2[loc]
        PC2_vel_in_window = PC2_vel[loc]
        PC2_vel_in_window = PC2_vel_in_window[~np.isnan(PC2_vel_in_window)]
        
        score = np.mean(PC2_in_window)
        vel = np.mean(PC2_vel_in_window)
        score_std = np.std(PC2_in_window)
        vel_std = np.std(PC2_vel_in_window)
    
        score_error = score_std/np.sqrt(len(PC2_in_window) )
        vel_error = vel_std/np.sqrt( len(PC2_vel_in_window) )
    
        
        center_list.append(center)
        score_list.append(score)
        vel_list.append(vel)
        score_std_list.append(score_std)
        vel_std_list.append(vel_std)
        
        score_error_list.append(score_error)
        vel_error_list.append(vel_error)
       
    plt.axis()
    plt.xlim(pc1Min, pc1Max)
    plt.ylim(pc2Min, pc2Max)
    plt.plot(center_list, score_list, 'b-o')
    plt.errorbar(center_list, score_list, yerr=score_error_list, capthick=2, capsize=3)
    plt.xlabel("PC1 (center of window)")
    plt.ylabel("Average PC2")

    plt.savefig(os.path.join('figures',\
        'sliding_window_plot_' + file_stem + '.png' ),
        dpi=600)
    plt.clf()

# Call makeMovementPlot to make the three movement plots. Due to the averaging,
# different plot limits are needed for PC2.
makeSlidingWindowPlot(tsOrigData, 'original_data',
                 pc1Min, pc1Max, -1.0, 1.5)
makeSlidingWindowPlot(tsEquinoxDataOrigNGAs, 'equinox_data_orig_NGAs',
                 pc1Min, pc1Max, -1.0, 1.5)
makeSlidingWindowPlot(tsEquinoxDataAllNGAs, 'equinox_data_all_NGAs',
                 pc1Min, pc1Max, -1.0, 1.5)

def makeMovementPlot(ts, file_stem, pc1Min, pc1Max, pc2Min, pc2Max):
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
    plt.xlim(pc1Min, pc1Max)
    plt.ylim(pc2Min, pc2Max)

    plt.savefig(os.path.join('figures',\
        'PC1_PC2_movement_plot_' + file_stem + '.png' ),
        dpi=600)
    plt.clf()

# Set the figure size (set back to default below)
plt.rcParams["figure.figsize"] = (15,8)

# Call makeMovementPlot to make the three movement plots
makeMovementPlot(tsOrigData, 'original_data',
                 pc1Min, pc1Max, pc2Min, pc2Max)
makeMovementPlot(tsEquinoxDataOrigNGAs, 'equinox_data_orig_NGAs',
                 pc1Min, pc1Max, pc2Min, pc2Max)
makeMovementPlot(tsEquinoxDataAllNGAs, 'equinox_data_all_NGAs',
                 pc1Min, pc1Max, pc2Min, pc2Max)

# Set the figure size back to the defaultd woplt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]

def makeHistogram(ts, fileStem, pc1Min, pc1Max):
    num_bins = 50
    n, bins, patches = plt.hist(ts.pcMatrix[:,0],
                                num_bins,
                                density=0,
                                facecolor='blue',
                                alpha=0.5)
    plt.xlabel("Projection onto first Principal Component")
    plt.ylabel("Counts")
    plt.xlim(pc1Min, pc1Max)
    plt.savefig(os.path.join('figures',\
        'PC1_histogram_' + fileStem + '.png' ),
        dpi=600)
    plt.clf()

# Call makeHistogram to make the three histograms
makeHistogram(tsOrigData, 'original_data',
              pc1Min, pc1Max)
makeHistogram(tsEquinoxDataOrigNGAs, 'equinox_data_orig_NGAs',
              pc1Min, pc1Max)
makeHistogram(tsEquinoxDataAllNGAs, 'equinox_data_all_NGAs',
              pc1Min, pc1Max)

