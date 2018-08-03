import genutils as gu 
from numpy import matlib        
contents=gu.prContents(media=['Glu 2%'], strains=['77.WT', 'blabla', 'B', 'C', 'D', 'E', 'F', 'G'], swapRowCol=True)

#cd to the accesspr directory

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as pch
from accesspr import *
import scipy.interpolate as scint
from prPlottingFunctions import *
import colors
import re
import platereader as pr

mfolder='/Users/s1259407/Dropbox/PhD/phd_peter_swain/data/plate_reader_data/20180601minsheng'


p=pr.platereader(mfolder+'/20180601minsheng.xlsx', mfolder+'/contents.xls')

#Generating a plot with all strains in one media.
cols={'229.WT': 'red',
 '252.Hxt1ko': '#587a67',
 '253.Hxt2ko': '#292573',
 '254.Hxt3ko': '#309545',
 '255.Hxt4ko': '#d620a0',
 '256.Hxt5ko': '#2a0d62'}

plotRawStatPerMediaRobust(p, dtype='OD', addLegend=True, strainColors=cols)

#draw plate layout of all wells
ppf.experimentOverview(p, dtype='OD')   


p.correctmedia()
#get growth statistics
p.getstats(cvfn='nn')

#store the pr file with processed data
pickle.dump(p, open('minsheng.pkl', 'wb'))


#p.makedataframe