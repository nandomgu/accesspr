import platereader as pr
import numpy as np
import _pickle as pickle 
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import os.path as path
import os 
#import pandas as pd
import modin.pandas as pd
import scipy
import scipy.interpolate as scint
from functools import partial
import prPlottingFunctions as ppf
from random import choice, randint, sample
import statsmodels.stats.api as sms
from matplotlib.colors import hex2color
from datetime import datetime
from xlrd import XLRDError
import re
import colors
from sklearn.preprocessing import normalize, StandardScaler
from sklearn.decomposition import PCA  
import time
import pdb
#from matplotlib.mlab import find
from matplotlib.widgets import RectangleSelector
from platform import system, release 
#from decimal import *
#getcontext().prec = 3
import modin.pandas as pd
#useful small functions   
#find the column number that corresponds to a specific column name
findloc=lambda loc, platelocs: np.nonzero([j==loc for j in platelocs])[0][0]      
grepcolumns= lambda name, df: np.nonzero([not not re.findall( f , name )[0][0] for f in df.columns.values])     
formatlist=lambda t: [t] if np.size(t)==1 else t
formatstring=lambda t: t if np.size(t)==1 and len(t)>1 else t[0]
formatval=lambda t: t[0] if isinstance(t, list) and np.size(t)==1 and isinstance(t, int) else np.array(t)
formatnum=lambda t: t[0] if np.size(t)==1 and type(t)==np.ndarray  else t 
formatfloat=lambda t: np.array([t]) if np.size(t)==1 and isinstance(t, float) else np.array(t)
treatmulti=lambda x: x[0] if np.size(x)>1 else x
fixempty=lambda x: np.nan if not x else x
removenulls= lambda x: x[x['strain']!='null']   #function for hassle free removal of nulls
def rotatelabels(ax, angle=90):
    [item.set_rotation(angle) for item in ax.get_xticklabels()]   
def internalStd(rawdata, stat, wtindices, fr=0.4, to=0.7, len=30):
    '''
    regression, transformedData=internalStd(rawdata, stat, wtindices, fr=0.4, to=0.7, len=30)
    
    internalStd standardises fluorescence measurements of any gain and machine in units of WT fluorescence.
    By doing so, fluorescence is expressed in units of the change of fluorescence stat in the wt.
    Why this is important:
    The relative fluorescence measurements can change dramatically among machine and gain. This standardisation
    brings all measurements to the same scale (the rate of change of fl in the wt).
    It does so by:
        a) calculating the rate of change of fluorescence channel stat between OD=fr and OD=to in the wells wtindices,
        which should correspond to an untagged strain (if possible in a standard control medium). 
        Higher values of 'to' ensure a more precise, convergent regression.
        b) the standard fluorescence is obtained by subtracting the intercept and dividing by the slope of the regression
        in all wells. 
    rawdata- a rawdata dictionary where rawdata[stat] is a m wells x n timepoints dataframe for the given stat.
    stat- the fluorescence measurement label of interest, such as 'GFP80'
    wtindices- the indices in the rawdata[stat] matrix that correspond to the reference condition, normally the untagged strain 
    in a standard medium.
    fr, to- the initial and final OD in the wt for which to calculate the rate of change in stat. can be altered o fit the 
    actual growth values of the experiment. 0.4 and 0.7 work great for BY4741 in 2% glucose and sucrose.
    '''
    vals=np.zeros([np.size(wtindices),len]) #prepopulating a matrix
    xvals=np.matlib.repmat(np.linspace(fr, to, len), np.size(vals, 0),1)  # prepopulating x values for such matrix   
    for d in range(0, size(wtindices)):
        func= interp1d(rawdata['OD'].iloc[wtindices[d],: ].values, rawdata[stat].iloc[wtindices[d],: ].values)
        vals[d, :]=func(linspace(fr, to, len))
    reg=scipy.stats.linregress(xvals.flatten(),vals.flatten())
    transformeddata= (rawdata[stat]-reg.intercept)/reg.slope
    return reg, transformeddata

def mergedicts(dictlist):
    final={}
    for dic in dictlist:
        for key in dic.keys():
            final[key]=dic[key]
    return final

def DFTransform(df, targetVariable, groupVariable, group, scale=1, shift=0):
    df2=df
    a=df
    a[targetVariable][a[groupVariable]==group]=df[targetVariable][df[groupVariable]==group]*scale+shift
    df=df2
    print(df[targetVariable][df[groupVariable]==group])
    print(a[targetVariable][a[groupVariable]==group])
    return a
    
def createPatches(colorDict, loc='upper right'):
	patches= [[pch.Patch(color=c) for c in colorDict[j]] for j in colorDict.keys() ]
	lgnds= colorDict.keys()
	plt.figlegend(patches, lgnds, loc)


def halfAUCIndex(vector):
    auc= np.trapz(vector)
    if np.isnan(auc):
        return np.nan
    areavector=np.zeros(np.size(vector))
    for j in range(0, np.size(vector)):
        areavector[j]= np.trapz(vector[range(0, j)]) ###each element is the area till that point.
    areavector=areavector/auc
    halfAreaIndex=np.where(abs(areavector-0.5)==np.nanmin(abs(areavector-0.5)))
    return(halfAreaIndex)

def halfAUCTime(p, media, strain, stat, aligned=True):
    try:
        hai=halfAUCIndex(p.d[media][strain][stat])
        if np.size(hai)>1:
            print('there is more than one half AUC index for '+strain+' in '+media+'. Choosing the first value.')
            hai=hai[0][0]
            print(hai)
        if np.isnan(hai):
            return np.nan
        else:
            if aligned== False:
                return(p.t[hai])
            else:
                return(p.d[media][strain]['Time centered at gr peak'][hai])
    except:
        return(np.nan)

def extractMagnitude(inputstr, munits='%'):
    d=inputstr.split(' ')[1].split(munits)[0]
    #print('number is ', d)
    return(float(d))
    
def filterDataFrame(df, filterkeys, filtervalues):###filter keys and values must be the same size!
    for x in range(0, np.size(filterkeys)):
        #try:
        df=df[df[filterkeys[x]]==filtervalues[x]]
        #except:
        #print('filtering of variable', filterkeys[x], 'failed. Carrying on...')
    return(df)

def DFsubset(df, variable, values):
    '''
    this function filters a dataframe such that dv[variable] only contains the values provided. useful when there is too many elementsfor variable.
    '''
    #values=formatlist(values) ##making sure that the value is a list even though it might be only one string.
    df2=pd.DataFrame(columns= df.columns)
    for x in values:
        df2=pd.concat([df2, filterDataFrame(df, [variable], [x])])
    df2=df2.reset_index(drop=True)
    return(df2)

def DFexclude(df, variable, values):
    '''
    this function filters a dataframe such that dv[variable] DOES NOT contain the entries with the values provided.    '''
    values=formatlist(values)
    newvalues=list(set(df[variable].values) - set(values))  #we obtain all the values in df[variable] discarding the ones in values
    df2=DFsubset(df, variable, newvalues)
    return(df2)

def digitizeCentered(vector, nbins=10, bins=False, extendBy=2):
    vector=np.around(vector,3)
    bins=symmetricMap(np.max(abs(vector)), nbins=nbins, extendBy=extendBy)
    digitized= np.digitize(vector, bins=bins)
    return(digitized, bins) 
    
def flatten(listoflists):
    ''' 
    flatten(listoflists)
    takes a list of lists and return a concatenated list of the individual elements of each list.
    '''
    flattened = [val for sublist in listoflists for val in sublist]
    return flattened
    
def colorScatter(p, media, strain, xstat=False, ystat='OD mean', colorBy='d/dtgr', symmetric=True, cmap='bwr',nbins=100, extendBy=2, alpha=1, markersize=12, marker='o', addLegend=True, vmin=0, vmax=0, xlabel=1,ylabel=1):
    if xstat==False:
        xstat='time'
    x=p.d[media][strain][xstat]
    y=p.d[media][strain][ystat]
    if symmetric==True:
        vmin=0
        vmax=nbins+extendBy-1
        cols, bins=digitizeCentered(p.d[media][strain][colorBy], nbins=nbins, extendBy=extendBy)
        #cols=cols.reshape(-1,1)
        #print('cols=', cols)
        plt.scatter(x, y,c=cols, cmap=cmap, edgecolor='none', vmin=vmin, vmax=vmax, alpha=alpha, s=markersize)
    else:
        if vmin!=0 or vmax!=0: ### if the limits have actually been specified, if no symmetry is demanded then  we do simple routine
            bins= np.linspace(np.min(np.around(p.d[media][strain][colorBy],2)), np.max(np.around(p.d[media][strain][colorBy],2)), nbins)
            #colvec=np.digitize(np.around(p.d[media][strain][colorBy],2), bins) #leaving unbinned to see if there is any inconvenient
            colvec=np.around(p.d[media][strain][colorBy],2) ##rounding the colors and the just plotting
            plt.scatter(x, y,c=colvec, cmap=cmap, edgecolor='none', vmin=vmin, vmax= vmax, alpha=alpha, s=markersize)
        else: ####uses the real limits of the colouring variable, binned
            bins= np.linspace(np.min(np.around(p.d[media][strain][colorBy],2)), np.max(np.around(p.d[media][strain][colorBy],2)), nbins)
            colvec=np.digitize(np.around(p.d[media][strain][colorBy],2), bins)
            vmin=np.min(colvec)
            vmax=np.max(colvec)
            plt.scatter(x, y,c=colvec, cmap=cmap, edgecolor='none', vmin=vmin, vmax= vmax, alpha=alpha, s=markersize)
    if xlabel==1:
        plt.xlabel(xstat)
    if ylabel==1:
        plt.ylabel(ystat)
    if addLegend==True:
        f=plt.gcf()
        refaxis = f.add_axes([0.90, .1, 0.05, 0.7])
        gradient = np.linspace(vmin, vmax, 256)
        gradient = np.flipud(np.vstack((gradient, gradient)).T)
        refaxis.imshow(gradient, aspect='auto', cmap=cmap)
        refaxis.set_title(colorBy)
        refaxis.get_xaxis().set_visible(False)
        refaxis.get_yaxis().set_ticks( np.linspace(256,0, 2) )
        refaxis.get_yaxis().tick_right()
        refaxis.set_yticklabels([np.round(bins[0],2), np.round(bins[-1],2)])

def symmetricMap(poslimit, nbins=10, extendBy=False):
    ''' 
    This function creates a symmetric bin set where the edges are -poslimit and poslimit, equally spaced in 10 bins.
    Extend adds extra bins at the edges.
    '''
    if extendBy != False:
        x, step= np.linspace(-poslimit, poslimit, nbins, retstep=True)
        bins= np.arange(-poslimit-step*extendBy, poslimit+step*(extendBy+1), step)
    else:
        bins=np.linspace(-poslimit, poslimit, nbins)
    return(bins)
    
normalflperod='c-GFP80'
def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))
def fadeColor(color, fadeIntensity=0.5):
    flag=0
    while(flag<1):
        colst= color.split('#')[1]
        chunkstr=list(chunkstring(colst,2))
        el1=hex(int(chunkstr[0],16)+int( fadeIntensity*float(int('aa',16)))).split('0x')[1]
        el2=hex(int(chunkstr[1],16)+int( fadeIntensity*float(int('aa',16)))).split('0x')[1]
        el3=hex(int(chunkstr[2],16)+int( fadeIntensity*float(int('aa',16)))).split('0x')[1]
        fadedColor= '#'+el1+el2+el3
        if int(el1+el2+el3,16)> int('ffffff',16):
            fadeIntensity= fadeIntensity-0.1 
            print('Warning: Fade intensity overflow. attempting weaker fading')
        else:
            flag=1
            print('color chosen:'+fadedColor)
    return fadedColor

def fitsBounds(x, bounds):
    if x>= bounds[0] and x<= bounds[1]:
        return True
    else:
        return False
        
def cycleThroughVector(vector, ind):
    '''
    function for generalized indexing. if index goes longer than the length of vector,
    then search carries on from the beginning.
    '''
     
def rewriteExperiments(xpr, experiments='all'):
    if experiments=='all':
        experiments= xpr.data.keys()
    for key in experiments:
        if path.isdir(xpr.source):
            pickle.dump(xpr.data[key], open(xpr.source + '/' +key, 'wb'))
            print('Experiment pickle file'+ xpr.source+'/'+key+' has been replaced.')

def normalizePerCategory(df, groupingVar, valueVar):
    cats= np.unique(df[groupingVar])
    
    for cat in cats:
        vec=df[df[groupingVar]==cat][valueVar]
        df[df[groupingVar]==cat][valueVar]= (vec- np.min(vec))/np.max(vec)
    return df

def alignTime(p, media, strain, FL='c-GFPperod', centerAtPeakFL=0):
	"centers time in the plate reader experiment p at the time of maximum growth rate (default) or time of peak fluorescence. this function is to be used by accesspr"
	# time t where growth rate is max for given curve
	if strain== 'null':
	    return np.nan, np.nan
	if not ppf.hasKey(p.d[media][strain], 'gr'):
	    out=np.nan;
	    centeredTime=np.nan;
	else:
	    if centerAtPeakFL==0:
	        centeredTime=p.t[np.where(p.d[media][strain]['gr']==max(p.d[media][strain]['gr']))]
	        alignedTimeVector=p.t-centeredTime
	        out=alignedTimeVector
	    else:
	        centeredTime=p.t[np.where(p.d[media][strain][normalflperod]==max(p.d[media][strain][FL]))]
	        #out= {'alignedTime': alignedTimeVector , 'rawFL':rawFLVector ,'normalizedFL':normalizedFLVector , 'peakFL':peak,  'peakTime':peakTime, 'gr':p.d[media][strain]['gr']}
	        alignedTimeVector=p.t-centeredTime
	        out=alignedTimeVector
	return out, centeredTime
	
def normalizeOverTime(p, media, strain, dtype=normalflperod, subtractBg=0):
    rawFLVector=p.d[media][strain][dtype]
    noBgFLVector=p.d[media][strain][dtype]-min(p.d[media][strain][dtype])
    normalizedFLVector=(noBgFLVector-np.mean(noBgFLVector))/np.std(noBgFLVector)
    return normalizedFLVector
	
def alignStats(p, media, strain, dtype, subtractBg=0):
    if ppf.hasKey(p.d[media][strain], dtype) and strain!='null':
        rawFLVector=p.d[media][strain][dtype]
        shp=np.shape(p.d[media][strain][dtype])
        #if we've got more than one column we average them.
        if len(shp)>1:
            rawFLVector=np.nanmean(p.d[media][strain][dtype], 1)
        noBgFLVector=rawFLVector-np.nanmin(rawFLVector)
        normalizedFLVector=(noBgFLVector-np.nanmean(noBgFLVector))/np.nanstd(noBgFLVector)
        normalizedFLPeak=np.nanmax(normalizedFLVector)
        alignedTimeVector, centeredTime=alignTime(p,media,strain)
        try:
            alignedPeakTime=alignedTimeVector[np.where(normalizedFLVector==normalizedFLPeak)]
        except:
            alignedPeakTime=np.nan
        FLPeak=np.nanmax(rawFLVector)
        halfFL=FLPeak/2
        #print('halffl='+str(halfFL))
        absolutePeakTime=p.t[np.where(rawFLVector==np.nanmax(rawFLVector))]
        ##find the first value at which half the expression is reached, finding the timepoint where the difference is minimal.
        #FUTURE WORK: could be improved to be found by interpolation to maximise precision. and to make sure that only 
        #points before the peak time are recovered.
        y=rawFLVector-halfFL
        g= scint.interp1d(p.t, y)
        tresampled=np.linspace(p.t[0], p.t[-1], 10000)
        yresampled=g(tresampled)
        f = scint.UnivariateSpline(tresampled, yresampled, s=0)
        try:
            roots=f.roots() #out of all the roots we get the earliest
            ##out of all the roots, we get the one closest to the absolute peak time
            subtraction= roots-absolutePeakTime
            ##we want negative and closest
            negatives= subtraction<0 ##those who are before the peak time
            if sum(negatives)==0:
                negatives=False(negatives) ##if there are no negatives well, whateever, try with anything available.
            closest= np.where(abs(subtraction[negatives])== np.min(abs(subtraction[negatives])))
            timeToHalfFL=roots[closest]
            #whichishalf=np.where(abs(rawFLVector-halfFL)==np.nanmin(abs(rawFLVector-halfFL)))
            #recvering the timepoint at which the half is reached.
        except:
            print('problem finding half fluorescence time for '+strain+' in '+media+'. obtaining time stats without FL background.' )
            try:
            #removing background fluorescence to improve the calculation
                noBGPeak= np.nanmax(noBgFLVector)
                y=noBgFLVector-noBGPeak/2
                g= scint.interp1d(p.t, y)
                tresampled=np.linspace(p.t[0], p.t[-1], 10000)
                yresampled=g(tresampled)
                f = scint.UnivariateSpline(tresampled, yresampled, s=0)
                roots=f.roots() #out of all the roots we get the earliest
                ##out of all the roots, we get the one closest to the absolute peak time
                subtraction= roots-absolutePeakTime #eak time
                ##we want negative and closest
                negatives= subtraction<0 ##those who are before the peak time
                if sum(negatives)==0:
                    negatives=False(negatives) ##if there are no negatives well, whateever, try with anything available.
                closest= np.where(abs(subtraction[negatives])== np.min(abs(subtraction[negatives]))) ### time from half fl to full fl. the sign of the steepness determines whether this was found before (negative) or after the peak.
                timeToHalfFL=roots[closest][0]
            except Exception as err:
                print('problem finding half fluorescence time for'+strain+' in '+media+': '+ str(err) )
                halfFLToPeakTime=np.nan
                slope=np.nan
                timeToHalfFL=np.nan 
        responseTimeAligned=timeToHalfFL-centeredTime
        halfFLToPeakTime= timeToHalfFL-absolutePeakTime ### time from half fl to full fl. the sign of the steepness determines whether this was found before (negative) or after the peak.
        steepness=1/treatmulti(halfFLToPeakTime)
        slope= (FLPeak-halfFL)/(absolutePeakTime-timeToHalfFL)
        #timeToHalfFL=p.t[np.where(rawFLVector==minDiffHalf)]
    else:
        rawFLVector=np.nan
        normalizedFLVector=np.nan
        normalizedFLPeak=np.nan
        FLPeak=np.nan
        halfFL=np.nan
        alignedPeakTime=np.nan
        absolutePeakTime=np.array([np.nan])
        slope=np.array([np.nan])
        responseTime=np.nan
        timeToHalfFL=np.nan
        steepness=np.nan
        halfFLToPeakTime=np.array([np.nan])
        responseTimeAligned=np.nan
    out= {'rawFL': rawFLVector, 'normalizedFL': normalizedFLVector, 'FLPeak': FLPeak, 'normalizedFLPeak':normalizedFLPeak, 'alignedPeakTime':alignedPeakTime, 'absolutePeakTime':treatmulti(absolutePeakTime), 'slope': treatmulti(slope), 'responseTime': timeToHalfFL, 'steepness': steepness, 'half2PeakTime': treatmulti(halfFLToPeakTime), 'responseTimeAligned': responseTimeAligned, 'halfFL': halfFL }
    return out
def str2numTS(tsString):
    tsArray= np.array(tsString.split())
    return tsArray

def extractTS(df, media, strain, tstype='FLperODTS'):
    a=df['strain']==strain 
    b=df['media']==media
    out=df[a & b][tstype].values
    #out=df[a & b][tstype].values[0].split()
    return( out)
def getRawData(expt, sheet=[]):
    if isinstance(expt, str): #process the string if it is one, otherwise assemble it.
        xl=pd.ExcelFile(expt)
        if not sheet:
            df=xl.parse(0)
        else:
            df=xl.parse(sheet)
    else:
        xl= pd.ExcelFile(expt.wdir + expt.name+'.xlsx')
        df= xl.parse(expt.dsheetnumber)
    dlabels= df[df.columns[0]].values ##get column names
    itime= np.nonzero(dlabels == 'Time [s]')[0] ##indices of the time vectors for each datatype, which are a good reference for upcoming data entries
    #itime + 0 is time. -1 is cycle number. +1 is temperature. +2 is the first well whose label is found at df.ix[itime[i]+2][0]
    #there is a one row blank between tables of different measurements.
    #last well from previous measurement type  df.ix[itime[i]-4][0]
    #blank df.ix[itime[i]-3][0]
    #measurement type df.ix[itime[i]-2][0]
    #cycle number df.ix[itime[i]-1][0]
    #time df.ix[itime[i]][0]
    #temperature df.ix[itime[i]+1][0]
    #first well df.ix[itime[i]+2][0]
    timepoints=df.ix[itime[0]][1:].dropna().values.astype('float')/3600 #
    datastarts= itime+2
    dataends= itime -4
    if np.size(datastarts)>1:
        matrixrows= dataends[1]-datastarts[0]+1 ###the end of the first matrix is 4 indices above the next time vector
    else:
        matrixrows=[i for i, x in enumerate([pd.isnull(j)  for j in dlabels[datastarts[0]:]]) if x ][0] 
    rawdata={} 
    for i in itime:
        d= i+2
        #rawdata[df.ix[i-2][0]]=pd.DataFrame(df.ix[d:d+matrixrows-1, 1:1+np.size(timepoints)], index=df.ix[d:d+matrixrows-1,0], columns=timepoints) #getting rid of incomplete timepoints
        dt=df.ix[i-2][0] ##the datatype
        actualdata=df.ix[d:d+matrixrows-1, 1:1+np.size(timepoints)].values #measurements
        welllabels=df.ix[d:d+matrixrows-1, 0] #first column is invariably the well labels.
        timelabels=[str(d) for d in range(0,np.size(timepoints))]
        #we choose the timepoints vector because it is recorded on the go so it tells you how many timepoints you have actually recorded.
        #number of cycles gets fully set from the beginning of the expt so it is not trustworthy.
        #we replace the overflow value with nans.
        rawdata[dt]=pd.DataFrame(data=actualdata, index=welllabels,columns=timelabels).replace('OVER', value=np.nan) 
    return rawdata, df, itime

def preprocessExpt(expt, fillnans=True, standardize=False, normOD=0.8, extension='preprocessed.xlsx', plot=False, main='GFP80', supporting='GFP60', machine=None, correctmachine=True, finalpath=''):
    '''
    preprocessExpt(expt, fillnans=True, standardize=True, mapFL=True, normOD=0.8, extension='preprocessed.xlsx' )
    preprocess data in excel sheet and output a filenamepreprocesed.xls
    '''
    rawdata, df, itime=getRawData(expt)
    #indices of the main channel which are not nans
    if fillnans==True:
        notnans=np.isnan(rawdata[main])==False
        #fitting a regression between the the supporting (lower) measurement and no nan values of the main measurement
        reg=scipy.stats.linregress(rawdata[supporting].values[notnans].flatten(),rawdata[main].values[notnans].flatten())
        ln=np.linspace(0, np.nanmax(rawdata[supporting]), 300)
        ln2=np.linspace(0, np.nanmax(rawdata[supporting]), 300)*reg.slope+reg.intercept
        if plot==True:
            plt.scatter(rawdata[supporting],rawdata[main])
            plt.xlabel(supporting)
            plt.ylabel(main)
            plt.plot(ln, ln2, color='red')
            #plt.legend([ 'y='+str(reg.slope)+'x'+stringSign(reg.intercept)+str(reg.intercept)])
        #take the supporting fluorescence and transform it into the  main
        rawdata[main]=rawdata[supporting]*reg.slope+reg.intercept  
        if correctmachine and machine=='Plate Reader 1':
            rawdata[main]= rawdata[main]*2.26+2.48
            print('PR1 experiment'+expt.name+'mapped to PR2\n')
        #take the values of supporting at the indices of the missing values of main and then apply the transformation
        #substitution=rawdata[supporting][np.isnan(rawdata[main])]*reg.slope+reg.intercept     
        #rawdata[main][np.isnan(substitution)==False]=substitution[np.isnan(substitution)==False]
    #performing the transformation based on WT fluorescence
    tdata=rawdata
    #reg1, tdata[stat1]=internalStd(rawdata, stat1, wtindices, fr=0.4, to=normOD)
    #reg2, tdata[stat2]=internalStd(rawdata, stat2, wtindices, fr=0.4, to=normOD)
    #### trying to create a new excel file with the new preprocessed data
    df2=df #create a copy of the original data
    datastarts= itime+2
    dataends= itime -4
    matrixrows= dataends[1]-datastarts[0]+1 ###the end of the first matrix is 4 indices above the next time vector
    timepoints=df.ix[itime[0]][1:].dropna().values.astype('float')/3600 #
    for i in itime:
        d= i+2
        #rawdata[df.ix[i-2][0]]=pd.DataFrame(df.ix[d:d+matrixrows-1, 1:1+np.size(timepoints)], index=df.ix[d:d+matrixrows-1,0], columns=timepoints) #getting rid of incomplete timepoints
        dt=df.ix[i-2][0] ##the datatype
        df2.ix[d:d+matrixrows-1, 1:1+np.size(timepoints)]=tdata[dt].values #measurements
    ##now that we have got the raw data we can apply functions to the matrices directly
    #preprocesspath=finalpath+datetime.now().strftime("%Y-%m-%d")+'/'+expt.name+extension;
    preprocesspath=expt.name+extension;
    df2.to_excel(preprocesspath, header=True, index=False)
    return tdata, preprocesspath

def bootstrapInterps(obj, centeringVariable, dtype, bootstrap=0, plot=True, col='black', alpha=0.15, marker=None, markevery=10):
    '''
    bootstrapInterps(obj, centeringVariable, dtype, bootstrap=0, plot=True, col='black', alpha=0.15)
        creates bootstrap replicates from timeseriesxreplicate matrix obj[dtype], and plots the mean and std deviation of these replicates.
        it samples random numbers of columns n=bootstrap times, calculates the mean of those replicates and piles up such 
        means into a bootstrap mean matrix B of len(obj[centeringVariable]) timepoints x n bootstrapmeans
        if bootstrap==0, then it plots the mean and sd of the replicates.
        output: 
            dictionary final with fields
                dtype= B, matrix of timepoints x bootstrapmeans
                centeringVariable= x variable (of interpolated x values for all replciates) used in plotting, gotten from obj.
    '''
    mn=np.nanmean(obj[dtype],1)
    sd=np.nanstd(obj[dtype],1)
    #calculating the coefficient of variation
    ##generate bootstrap replicates
    ncols=np.size(obj[dtype], 1)
    nrows=np.size(obj[dtype], 0)
    ###generating new traces  from combining other traces
    final={}
    final[centeringVariable]=obj[centeringVariable]
    for j in range(0, bootstrap): ###till the number of bootstrap replicates is reached
        ###sample a number of replicates
        if j==0:
            reps=np.nanmean(obj[dtype][:, sample(range(0,ncols),randint(1,ncols))], 1)  ###sample from one to ncols -1 , and get those specific replicates from the set
        if j>0:
            addmn=np.nanmean(obj[dtype][:, sample(range(0,ncols),randint(1,ncols))], 1)  ###sample from one to ncols -1 , and get those specific replicates from the set
            reps=np.column_stack([reps, addmn])
    if bootstrap>0:
        totalmn=np.nanmean(reps,1)
        totalsd=np.nanstd(reps,1)
        plt.plot(obj[centeringVariable], totalmn, color=col, marker=marker, markevery=markevery)
        plt.fill_between(obj[centeringVariable], totalmn-totalsd, totalmn+totalsd, color=col, alpha=alpha)
        plt.xlabel(centeringVariable)
        plt.ylabel(dtype)
        final[dtype]=reps
        final[dtype+' mean']=totalmn
        final[dtype+'bootstrapsd']=totalsd
    else:
        final[dtype]=interpolated[dtype]
        final[dtype+' mean']=mn
        final[dtype+'sd']=sd
        plt.plot(obj[centeringVariable], mn, color=col)
        plt.fill_between(obj[centeringVariable], mn-sd, mn+sd, color=col, alpha=alpha )
        plt.xlabel(centeringVariable)
        plt.ylabel(dtype)
    return final

class accesspr:	
    '''
    accesspr version 4.89(the number matches with that of the compatible platereader software version)
    
    accesspr is a class that allows to integrate and organize the information
    from many plate reader experiments to produce publication-grade plots and tables. 
    To exploit the functionality of accesspr, it is ideal
    to learn about pandas dataframes and the seaborn graphics package. 
    
    ***Newest FEATURES***
    pcaclick- get a bird's eye view of all your dataset's curves by using PCA, then browse and select curves of interest.
    -----------------------------------------------------------------------------------------------------------
    ATTRIBUTES: 
    
    The attributes of accesspr allow you to know what is in your experiments and what stage of the processing they are in. 
    
    source: 
        the path of the folder that contains the pickle files
    data:
        a dictionary of platereader objects, corresponding to each experiment in each pickle file. 
        The keys are each experiment's pickle file name. e.g. xpr.data[picklefilename].d[media][strain]
    allreplicates:
        a pandas dataframe of the experiment, media, strain and plate location of each replicate in the experiment ensemble.
    allcontents:
        a pandas dataframe containing all unique combinations of strains and media found in all experiments. 
    statcontents:
        a pandas dataframe describing what statistics are in which experiment, and the fraction (0-1) of conditions 
        in that experiment that have such statistic. By default it looks for tme varying statistics, but you can search for other statistics using
        self.containsstat(statname)
    extractionFieldsScalar: 
        Point statistics to be extracted when making dataframes through self.makedataframe('notime'). By default this includes all point statistics
        found in the experiments, but more can be added.
    extractionFieldsTime: 
       Time varying statistics/measurements that are present in at least one experiment. Suc statistics can be extracted 
       using xpr.makedataframe('timerow') or xpr.makedataframe('timecol', times=[x1,x2,x3...])
    allstrains: 
        contains a nonredundant list of all strains found in the experiments.
    allmedia:
        contains a nonredundant list of all media found in the experiments.
    exptInfo:
        dictionary of paths for all the experiments loaded. 
        self.exptInfo[exptname]
                               ['datapath']
                               ['contentspath']
    allexperiments:
        contains a nonredundant list of all experiments in the accesspr object. 
    FL: 
        a dictionary that contains the default fluorescence channels to be used for every experiment. can be changed manually or by calling
        xpr.assignFL()
    machines:
        a dictionary (with experiment keys) that contains the machine in which this experiment was run. currently identified machines are 'Plate Reader 1' and 'Plate Reader 2'
    serialnumbers:
        a dictionary (with experiment keys) that contains the machine serial number. even if the machine has no identified human-readable name (see machines), the serialnumber will be stored here.
    -------------------------------------------------------------------
    METHODS: (*:runs during initialization)
    To create an accesspr object, 
        xpr=accesspr(PICKLEFOLDERPATH, encoding='latin1', ignoreFiles=False, FL='noFL', FLperod='noFLperod', onlyFiles=False, analyseFL=True, preprocess=True):
            initializes an accesspr instance by using the files in path or list of paths PICKLEFOLDERPATH. 
            ignoreFiles is a list of the file names in this folder that should be ignored.
    The paths of any experiments that fail to import will be added to the xpr.failedfiles list.
    
    To interrogate about experiment contents and keep track of experiment processing,
    
        * getallcontents(self):
             produces self.allconditions attribute (see attributes)
             
        * containsstat(self, STAT):
            produces the statcontents attribute (see attributes)
            
        * listcontents(self, verbose=True):
            fills the allstrains and allmedia attributes. If verbose==True, produces a human-readable display of the contents of each experiment.
            
        *assignFL(mainFL='noFL', mainFLperod='noFLperod', experiments=False)
            if arguments are specified, sets the given channels as the  default fluorescence channels for each experiment. 
            if arguments are unspecified, sets GFP as default.
            You can specify a list of fluorescence channels to be found in an experiment in order of priority. (if one is not found, then the next one is searched)
        
        containssetup(self, MEDIA,STRAIN):
            retrieves a list of experiments that contain media MEDIA and strain STRAIN. 
            Internally sets the temporary property containslist to this experiment list. 
        
        plotrawreplicates(self, MEDIA, STRAIN, dtype='OD', exptColors=False):
            plots the raw repicates of dtype VS time of a specific strain or media, for all the experiments that contain that specific condition.
            if exptColors is not specified, then random colors are assigned to each experiment and returned as a function output.
        plotrepmean(self, MEDIA, STRAIN, dtype='')
            uses interpTimes to align and interpolate replicates across experiments and creates a mean± std line plot for each media in MEDIA
        findMachines()
            finds the machine wher this experiment was run by looking into several places: first, in pr.machine, then in the experiment file.
            Returns dictionary self.serialnumbers[EXPERIMENT] ad self.machines[EXPERIMENT], and if no machine is found, self[expt]=np.nan.
    To do platereader-related and other data processing, experiment by experiment,
    
        * alignAll(self): 
        for each media, strain and experiment, produces 2 alternative time vectors,
        found along the other stats at xpr.data[EXPT][MEDIA][STRAIN]
            'Time centered at gr peak' for which t=0 corresponds to the max growth rate time.
            'Time centered at flperod peak' for which t=0 corresponds to the peak time of the default fluorescence channel.
            
        correctauto(self, f=['GFP', 'AutoFL'], experiments='all',media='all', strains='all', refstrain=['WT'], figs=True, correctOD=True, noruns=2, bd=False, no1samples=100, rewrite=False, rerun=False):
            corrects autofluorescence for each condition in every experiment. 
            if rewrite ==True, then it rewrites the pickle file of the relevant experiments *** CAUTION: PORTABILITY ISSUES REPORTED FOR REWRITING
            if rerun==True, it reprocesses all experiments regardless of whether they have been previously processed.
            
        getstats(self,experiments='all', media='all', strain='all', dtype=['OD'], rewrite=False, bd=False, cvfn='sqexp', esterrs=False, stats=True, plotodgr=False, rerun=False):
            calculates the growth rate over time, as well as the lag time, final OD and other related statistics per experiment. 
            if dtype is not 'OD', then it calculates the d/dt of the particular statistic.
            if rewrite ==True, then it rewrites the pickle file of the relevant experiments *** CAUTION: PORTABILITY ISSUES REPORTED FOR REWRITING. 
            if rerun==True, it reprocesses all experiments regardless of whether they have been previously processed.
        
    Data processing and extraction for a specific media and strain
        
        interptimesnew(conditionsDF, dtype='FLperod', centeringVariable='gr', upperLim=16, exceptionShift=0.01, ignoreExps=False, factor=False)
            Creates a sub dataset such that time-varying property dtype is sampled at the same times across replicates.
            This is important because the exact sampling times might vary across replicates.
            returns a dictionary with two keys:
                ['time']- the general sampling times for all replicates
                [dtype]- the  dtype statistic sampled at times ['time']for all replicates found. 
            the function 
                a) aligns the replicates to their max growth rate time
                b)trims experiments so that all have data in a same shared time interval
                c)performs interpolation using an arbitrary experiment's timepoints. 

                *** useful dataFrame-specific functions:
                    
                    DFsubset(dataFrame, variable, values)
                        Allows to create a smaller data frame out of a bigger one such that
                        a variable (a key/column) only has specific values. 
                        This is useful when an entire dataframe is too big or messy to work with in its entirety.
                        
                        #example. create a small dataframe such that the new one only contains 3 specific strains. 
                        #first we check the total number of strains in the original dataset
                        np.size(np.unique(df['strain']))
                        #Out[]: 46
                        #then we subset this dataframe
                        df2= DFsubset(df, 'strain', ['409.Hxt4', '506.Hxt4rgt2', '403.Hxt4std1']) 
                        #then we check whether the dataframe's strain number was reduced.
                        np.size(np.unique(df2['strain']))
                        #Out[1366]: 3
    Data processing and extraction of all accesspr information
        
        extractallinfonew(self, excludeNull=False)
            extracts the information of all replicates for all conditions. 
            returns a dataframe most single point statistics and string-formatted- time series of each condition's replicate info. 
            the statistics deployed and the order (currently fixed) can be found in the self.extractionFields attribute.
            The function concatenates the output of extractRepInfo iterating over all conditions in accesspr.
            It is highly recommended to exclude the null especially when plotting with pandas built-in functions or seaborn, as they cannot handle 
            NaNs very well.
    
    Plotting information contained in time series: 
    
        timeStat(self, media, strain, times=[4], dtype='FLperod', includeMediaStrain=True)
            Extracts the value of statistic dtype at Absolute times times for all replicates of strain in media.
            the values in the times vector that exceed the duration of the shortest experiment (self.interpLimit) are excluded.
            returns a dataframe where each time of the form
                    experiment  media   strain  t0 t1 ...
                
                
                
            if includeMediaStrain=False, the media/strain columns are excluded (useful when no grouping 
            is needed.
        timeStatAligned(self, media, strain, times=[0], dtype='FLperod', includeMediaStrain=True)
            Extracts the value of statistic dtype at ALIGNED times times (see interpTimes) for all replicates of strain in media. 
            returns a dataframe where each time of the form
                        media   strain  t0 t1 ...
                exp1
                exp2
                
            if includeMediaStrain=False, the media/strain columns are excluded (useful when no grouping 
            is needed.
            
        timeStatAll(times, aligned=False)
            Performs timeStat (if aligned==False) or timeStatAligned (if timeStat==True) 
            on all conditions in the xpr object. 
            Returns a dataframe with all conditions.
            It is suggested to reorganize the dataset this way:  
                
                tidydf=pd.melt(df, id_vars=['media', 'strain'], value_vars=times, var_name='time', value_name=dtype)
                
                Doing this tells seaborn that all the time columns are actually the same kind of information:
                It collapses all time columns into one, such that the resulting  dataframe look like this
                (e.g., for times -1 and 1):
                    
                index   media   strain  time    FLperod
                    0   Glu 2%      Hxt4   -1  1446.09
                    1   Glu 2%      Hxt4   -1  1446.09
                    2   Glu 2%      Hxt4    1  10058.7
                will allow sns.factorplot and sns.swarmplot  to plot the time vs FLperod relationship
                appropriately. 
            
        colorScatter(self, media, strain, experiments=False, xstat=False, ystat='FLperod', colorBy='d/dtgr', symmetric=True, cmap='bwr',nbins=40, extendBy=2, alpha=1, markersize=12, marker='o', addLegend=False)
            creates a scatter plot of array xstat (defaults to time) vs array ystat, where the color of each point is determined by the magnitude array colorBy
            all three xstat, ystat and colorBy must have the same number of arguments.
            If experiments ==False (default) all experiments of with media/strain are used.
            color map properties:
                cmap- determines the default color map used (defaults to blue-white-red 'bwr')
                nbins- the number of bins by which colorBy is to be discretized.
                symmetric- if True(default), the bins are generated in symmetrical fashion around 0 (e.g. very useful for z score binning).
                    otherwise it just spans from the minimum to the maximum value of colorBy.
                extendBy- adds extra bins not contemplated by colorBy. useful to distinguish outliers.
            addLegend- adds a colormap scale reference. must be set to False until the last scatter plot has been added to the figure.
                
        plotrepmean(self, media, strain, dtype='FLperod', col='Black', alpha=0.2, exceptionShift=0.01):
            for a given strain and condition, runs interpTimes and plots the mean ± standard deviation (shaded area) over time for the combined replicates. 
            It operates on the current figure, so multiple plots can be overlayed.
            
        plotrawreplicates(self, MEDIA, STRAIN, dtype='OD', exptColors=False):
            plots the raw within-experiment replicates of dtype VS time of a specific strain or media, for all the experiments that contain that specific condition.
            if exptColors is not specified, then random colors are assigned to each experiment. 
            Not suitable for FLperod or gr, as it needs the statistic
            to be present per timepoint per well.
            returns a dictionary of colors assigned to each experiment.
        
        plotConditionStd(self, media, strain, dtype='OD', normalize=0, centerAtFLPeak=0,col='black', range=0):
            plots the individual experiments of strain in media, without time processing. 
            It is ideal for statistics that collapse multiple within-experiment replicates, like FLperod or gr.
            This function is quite basic and doesn't interpolate times nor it trims edges.
            The advantage of this function is that it allows to normalize each replicate individually.
            
        plotConditionAligned(self, media, strain, dtype='OD', normalize=0, centerAtFLPeak=0,col='black', range=0)
            plots the individual experiments of strain in media, centering at the time of max growth rate. 
            It is ideal for statistics that collapse multiple within-experiment replicates, like FLperod or gr.
            This function is quite basic and doesn't interpolate times nor it trims edges.
            The advantage of this function is that it allows to normalize each replicate individually.
    
    Plotting single-point statistics
        
        text2d(self, media=False, strains=False, xstat='FinalOD', ystat='maxGR', strainColors=False, xlim=False, ylim=False, markersize=500, newFig=True):
            plots single point statistics xstat versus ystat, using the media string as a plotting marker. 
            If no strains or media are included, then it plots all media, colouring each strain in a different way.
            returns a strainColors dictionary, where each key is a strain and each value is a hexadecimal color string ('#rrggbb'). 
        
        Plotting single point statistics through seaborn
            One of the greatest advantages of contstructing a pandas dataframe is  the possibility to 
            use the seaborn graphics package as seaborn was developed with dataframes in mind. 
            I recommend skimming through all the plotting possibilities with seaborn in these 2 links
                https://stanford.edu/~mwaskom/software/seaborn/examples/index.html
                https://stanford.edu/~mwaskom/software/seaborn/generated/seaborn.barplot.html
            Examples.
            #create an accesspr object
            xpr=accesspr(PATH)
            
            #create a dataframe with all conditions in the experiment.
            prDataFrame= extractAllInfo(excludeNull=True)
            
            #Create a barplot of the fluorescence peak for each strain
            plt.figure()
            sns.barplot('strain', 'FLPeak', data=prDataFrame )
            
            #Create a barplot of the fluorescence peak by strain, this time grouping by different kinds of media.
            plt.figure()
            sns.barplot('strain', 'FLPeak', data=prDataFrame, hue='media' )

    '''

    def __init__(self, source=[], params={}, encoding='latin1', ignoreFiles=False,analyseFL=True, onlyFiles=False, FL='noFL', FLperod='noFLperod',  preprocess=False, ignorepickles=False, ignorexls=False, filecontains=[]):
        '''
        obj=acesspr( source, params={}; encoding='latin1', ignoreFiles=False, FL='noFL', FLperod='noFLperod', onlyFiles=False, analyseFL=True, preprocess=True)
        initializes an accesspr instance from a path specifying either a directory or a pickle file. if this fails, try changing the encoding to ascii, utf8, utf16 or latin1 (so far tried).
        source- the path or paths where the plate reader data is located. 
        ignoreFiles= list of file names to ignore 
        onlyFiles= list of experiments to exculsively add (i.e. ignore all experiments except for this list). 
            Useful when a folder has too many irreevant experiments
        FL= indicates the main fluorescence to be indexed in the experiments. defaults to noFL
        FLperod= name of the variable that contains the fluorescence per od. only edit if your FL per od has a non standard name. 
        analyseFL= if False, fluorescence routines are ignored.
        preprocess- if working with FL experiments run in different plate readers and different gains, preprocess corrects for all these differences.
            For details, see self.preprocessAll2().
        '''
        self.version='4.94'
        self.ignorexls=ignorexls
        self.ignorepickles=ignorepickles
        self.date=datetime.today().strftime('%Y-%m-%d')
        self.savedir='./accessprdata'
        self.objectpath=[]
        self.stagepath=[]
        self.stagecount=[]
        self.prtype='Tecan'
        self.ignoredwells=[]
        self.filecontains=filecontains
        #source is a directory with several pickle files. if no pickles are found, then it tries to find experiments.
        self.source = source
        self.failedfiles=[]
        self.transformpr1pr2=lambda x: 2.26*x+2.48
        self.ignoreFiles=ignoreFiles
        self.onlyFiles=onlyFiles
        self.analyseFL=analyseFL
        self.encoding=encoding
        self.interpLimit=20
        self.preprocess=preprocess
        #exptInfo contains the source of each experiment, whether it is a pickle or just raw pr data
        self.exptInfo={}
        self.data = {}
        self.experimentDuration={}
        #activeExperiments is an array the list of currently relevant pickle files meant to be loaded onto memory.
        #this list is modified by contains functions
        self.activeExpts=[];
        self.releaseNotes='xpr.FL contains default fluorescence per experiment'
        self.defaultStats=['gr','GFP','c-GFPperod', 'GFP100','c-GFP100perod','GFP90','c-GFP90perod','GFP80','c-GFP80perod','GFP70', 'c-GFP70perod','GFP60', 'c-GFP60perod','GFP50', 'c-GFP50perod']
        self.extractionFieldsScalar=[]
        self.extractionFieldsTime=[]
        self.extractionFields=['experiment', 'machine','media','strain','InitialOD','FinalOD','lagTime','realTime','alignedTime','maxGR','maxGRvar','maxGRTime','grAUC','InitialRawFL','InitialFLperOD','FinalFLperOD','FLPeak','FLAbsPeakTime','FLAlignedPeakTime','FLperodAUC','slope','responseTime','steepness','responseTimeAligned']
        self.extractionFieldsOD=['experiment', 'machine','media','strain','InitialOD','FinalOD','lagTime','realTime','alignedTime','maxGR','maxGRvar','maxGRTime','grAUC']
        self.extractionFieldsGrowth=['d/dtgr var','time of max gr','time of max gr var','local max gr','OD gp','doubling time var', 'doubling time', 'max OD', 'lag time', 'max gr', 'max gr var']
        self.extractionFieldsFL=['FLPeak','absolutePeakTime','alignedPeakTime','responseTime','responseTimeAligned','halfFL','normalizedFLPeak','slope','steepness']
        self.mediaValue={'Glu 0.05%': 0.05, 'Glu 0.1%': 0.1, 'Glu 0.2%': 0.2,'Glu 0.4%': 0.4,'Glu 0.6%': 0.6,'Glu 0.8%': 0.8,'Glu 1%': 1,'Glu 1.5%': 1.5,'Glu 2%': 2,'2% Glu': 2,'0.2% Glu': 0.2, 'SucGlu 1%': 0.5, 'SucGlu 1.8% 0.2%': 0.2/2, 'SucGlu 0.2% 1.8%': 1.8/2}
        self.strainAlias={'YST_498': 'Hxt1', 'YST_499': 'Hxt1', 'Hxt4': 'Hxt4', 'Hxt2': 'Hxt2','Hxt3': 'Hxt3','Hxt5': 'Hxt5','Hxt6': 'Hxt6','Hxt7n': 'Hxt7' }
        self.allconditions=pd.DataFrame()
        if analyseFL==True: #if fl analysis is required, we first fill up the FL fields
            self.analyseFL=True
            self.FL={}
            self.consensusFLs=[]
            self.consensusFL=FL
            self.consensusFLperod=FLperod
            self.supportingFL=['GFP60', 'GFP70', 'GFP75']
        if source:
            self.loadInitial() #load experiments, potentially doing preprocessing if this is indicated during initialisation.
        if self.data:#once experiments have been loaded, we assign fluorescence to each of the individual experiments.
            self.listcontents(verbose=False) ## indexes experiments, runs getvariables(), checkallstats() and getallcontents
            #self.getallcontents()
            self.strainColors=False ##### list of colours to be assigned to each strain in plotting routines. currently not in use.
            self.mediaColors=False ##### list of colours to be assigned to each media during plotting routines. currently not in use.
            self.exptColors=False  ##### list of colours to be assigned to each media during plotting routines. currently not in use.
            self.aligned=False
            self.processRecord=pd.DataFrame(columns=['experiment', 'correctauto', 'getstatsOD', 'getstatsFLperod'], index= list(self.data.keys())) ### stamps are incomplete, error, 1 if processing is complete but date unknown  and date if processing date is known.
            self.statcontents=pd.DataFrame( index= list(self.data.keys())) ### stamps are incomplete, error, 1 if processing is complete but date unknown  and date if processing date is known.
            self.wildTypeList=['WT', '77.WT', '229.WT']
            self.machines={}
            self.serialnumbers={}
            ##This stores the normalisation factor for each experiment.
            self.normalizationFactor={}
            #this provides the conditions to obtain the normaliation factor. the channel, and the od for which the normalizing factor is obtained
            self.normStandard={}
            self.normStandard['channel']='AutoFL'
            self.normStandard['OD']=0.4 
            #self.getNormFactor()
            self.refstrains=dict()
            self.getvariables(verbose=False)
            self.replicateLocations()
            self.checkallstats()
            self.findMachines()
            if analyseFL==True and not self.FL:
                for key in self.data.keys():
                    self.FL[key]={}  
                self.assignFL(mainFL=FL, mainFLperod=FLperod)
                self.assignsupportingFL()
            else:
                self.analyseFL=False
            #             try:
            #                 self.alignAll(rerun=True)
            #             except Exception as err:
            #                 print('failed to align: '+str(err))
            #self.interpLimit= np.min(np.array(list(self.experimentDuration.values())))
            #this stores the absolute duration (in hrs) of the shortest experiment
            #local lambdas
            self.getprops= lambda num: (self.allreplicates.iloc[num, 0], self.allreplicates.iloc[num, 1], self.allreplicates.iloc[num, 2], self.allreplicates.iloc[num, 3])
            self.getindcurve= lambda expt, media, strain, plateloc, dtype:  self.data[expt].d[media][strain][dtype][:, findloc(plateloc, self.data[expt].d[media][strain]['plateloc'])]
            self.getenscurve= lambda expt, media, strain, dtype:  self.data[expt].d[media][strain][dtype]
            self.findnull= lambda e, m, s, dtype: self.getenscurve(e,m,'null',dtype).mean(1)
            self.findwt= lambda e, m, s, dtype: self.getenscurve(e,m,'WT',dtype).mean(1)
            self.findref= lambda e, m, s, dtype, ref: self.getenscurve(e,m,ref,dtype).mean(1)
            self.nullsubtract= lambda e, m, s, dtype: self.getenscurve(e,m,s,dtype).mean(2)-self.getenscurve(e,m,'null',dtype).mean(1) 
            self.wtdivide= lambda e, m, s, dtype: self.getenscurve(e,m,s,dtype).mean(2)/self.getenscurve(e,m,'WT',dtype).mean(1)  
            self.wtfsubtract= lambda e, m, s, dtype: self.getenscurve(e,m,s,dtype).mean(2)-self.getenscurve(e,m,'WT',dtype).mean(1) 
            self.ensrefdivide= lambda e, m, s, dtype, ref: self.getenscurve(e,m,s,dtype).mean(2)/self.getenscurve(e,m,ref,dtype).mean(1)  
            self.ensrefsubtract= lambda e, m, s, dtype, ref: self.getenscurve(e,m,s,dtype).mean(2)-self.getenscurve(e,m,ref,dtype).mean(1) 
            self.refdivide= lambda e, m, s, pl, dtype, ref: self.getindcurve(e,m,s,pl,dtype).mean(2)/self.getenscurve(e,m,ref,dtype).mean(1)  
            self.refsubtract= lambda e, m, s, pl, dtype, ref: self.getindcurve(e,m,s,pl,dtype).mean(2)-self.getenscurve(e,m,ref,dtype).mean(1) 
    def getvariables(self, verbose=True):
        '''
        getvariables(self, verbose=True)
        looks into al the experiments and catalogues all available variables for each condition.
        Names of all variables found are stored in:
            self.extractionFieldsScalar- scalar variables, which will be brought to the dataframe.
            self.sextractionFieldsTime- time variables, which can serve as a reference for time plotting
        '''
        scalarvars=[]
        timevars=[]
        self.extractionFieldsScalar=[]
        self.extractionFieldsTime=[]
        for j in self.allexperiments:
            localvars=self.data[j].datavariables(out=True)
            [scalarvars.append(j) for j in localvars[2]]
            [timevars.append(j) for j in localvars[1]]
        self.extractionFieldsScalar=list(np.unique(scalarvars))
        self.extractionFieldsTime=list(np.unique(timevars))
        if verbose==True:
            print('Scalar variables registered (found in self.extractionFieldsScalar:\n'+'\n'.join(self.extractionFieldsScalar)+'\n')
            print('Time-varying variables registered (found in self.extractionFieldsTime:\n'+'\n'.join(self.extractionFieldsTime)+'\n')
    def checkallstats(self):
        self.statcontents=pd.DataFrame( index= list(self.data.keys()))
        for k in self.defaultStats+self.extractionFieldsTime:
            self.containsstat(k, printstats=False)
            #self.containsstat('GFP', printstats=False)
            #self.containsstat('c-GFPperod', printstats=False)
            #self.containsstat('GFP100', printstats=False)
            #self.containsstat('c-GFP100perod', printstats=False)
            #self.containsstat('GFP90', printstats=False)
            #self.containsstat('c-GFP90perod', printstats=False)
            #self.containsstat('GFP80', printstats=False)
            #self.containsstat('c-GFP80perod', printstats=False)
            #self.containsstat('GFP70', printstats=False)
            #self.containsstat('c-GFP70perod', printstats=False)
            #self.containsstat('GFP60', printstats=False)
            #self.containsstat('c-GFP60perod', printstats=False)
            #self.containsstat('GFP50', printstats=False)
            #self.containsstat('c-GFP50perod', printstats=False)
            #self.containsstat('FLperod', printstats=False)
    def getNormFactor(self):
        #if self.allexperiments != self.activeExpts:
        #    self.loadFresh()
        for expt in self.allexperiments:
            for m in self.data[expt].allconditions:
                #this is the wildtype in this experiment
                s=[j in self.wildTypeList for j in self.data[expt].allstrains]
                wt=str(np.array(self.data[expt].allstrains)[np.where(s)][0])
                print(wt)
                print(str(np.array(self.data[expt].allstrains)[np.where(s)][0]))
                try:
                    mmin=np.min(self.data[expt].d[m][wt]['OD mean']-self.normStandard['OD'])
                    ind=np.where((self.data[expt].d[m][wt]['OD mean']-self.normStandard['OD'])==mmin)
                except:
                    print('couldn''t get to normalise in '+str(m)+' '+str(wt))
                self.normalizationFactor[expt]=np.mean(self.data[expt].d[m][wt][self.normStandard['channel']][ind, :])
    def loadInitial(self):
        '''
        Search for files in the source path. the pr files can be either:
        1. pickles created with the pickle library
        2. folders, each containing two excel datasheets: a) data file and b)contents file (characterised by ending in contents.xls or contents.xlsx
        '''
        importerrors='The following files could not be imported:\n'
        print('accesspr version'+self.version)
        if isinstance(self.source, str):
            self.source=[self.source]
        for src in self.source:
            if path.exists(src) and path.isdir(src):
                dirs = os.listdir(src)
                for entry in dirs:
                    if entry.endswith('.pkl') and (self.ignoreFiles == False or ((entry in self.ignoreFiles)==False)) and (self.onlyFiles==False or entry in self.onlyFiles) and (not self.data or not (entry in list(self.data.keys())) ) and self.ignorepickles==False:
                        print('trying to import pickle '+entry)
                        pkl = open(src + '/' + entry, 'rb')
                        p=pickle.load(pkl,encoding= self.encoding)
                        exptname=re.split(r'[/\\\\]', p.name)[-1]
                        #print('processed name of'+ p.name+':- '+exptname+'\n')
                        self.data[exptname] = p
                        self.activeExpts.append(exptname)
                        pkl.close()
                        self.experimentDuration[exptname]=self.data[exptname].t[-1]
            ##if the path its                
                    elif path.exists(src) and os.path.isdir(src+'/'+entry) and (self.ignoreFiles == False or ((entry in self.ignoreFiles)==False)) and (self.onlyFiles==False or entry in self.onlyFiles):
                        #this generally means that the path is a directory. we look inside to see if there are any excel files that look like pr files.
                        datafile=[]
                        contentsfile=[]
                        pth=src+'/'+entry; #store the full path of the directory.
                        #look for 
                        files=os.listdir(pth)
                        hasExcel=0;
                        for f in files:
                            if f.endswith('.xls') or f.endswith('.xlsx')  :
                                hasExcel+=1;
                                if f.endswith('contents.xls') or f.endswith('contents.xlsx') or np.array([isinstance(j, str) for j in re.findall("contents",f, flags=re.IGNORECASE)]).any():
                                    contentsfile=src+'/'+entry+'/'+f
                                    #print('contentsfile:'+contentsfile)
                                else: #it is not a contents file so we check for, or do preprocessing.
                                    if not f.endswith('preprocessed.xlsx'): #if no preprocessing required import file wo preprocessing
                                        datafile=src+'/'+entry+'/'+f
                                            #print('datafile:'+datafile+'\n')
                        if hasExcel>=2 and not self.ignorexls :
                            print('directory with excel files will be incorporated.')
                            #the folder may contain more than one experiments and the name of the experiment may differ from
                            #that of the folder. therefore we use the experiment filename as a unifying reference
                            #print('data and contents files: '+ formatstring(datafile)+formatstring(contentsfile))
                            if contentsfile and datafile and ((self.onlyFiles==False and not self.filecontains) or  np.array([j in datafile for j in self.filecontains]).any()):
                                #print(contentsfile, datafile)
                                #creating pr file from scratch
                                try:
                                    p=pr.platereader(datafile, contentsfile, info=False)
                                    exptname=re.split(r'[/\\\\]', p.name)[-1]
                                except XLRDError:
                                    print('failed to import experiment: excel file appears to be corrupt.')
                                    importerrors=importerrors+'\n'+exptname+'\t file appears corrupt'
                                    self.failedfiles.append([datafile, contentsfile])
                                    continue
                                except Exception as err:
                                    print('problem while preprocessing:'+str(err))
                                    importerrors=importerrors+'\n'+exptname+'\t seems like a wrong sheet'
                                    self.failedfiles.append([datafile, contentsfile])
                                if ppf.hasKey(self.data, exptname):
                                    continue
                                else:
                                    self.data[exptname]=p
                                    self.activeExpts.append(exptname)
                                    print('platereader experiment '+exptname+':\nsuccessfully loaded\n')
                                    if self.preprocess==True: #we want to preprocess the 
                                        try:
                                            self.assignFL(experiments=formatlist(exptname), mainFL=self.consensusFL)
                                            self.assignsupportingFL(experiments=formatlist(exptname))
                                            self.preprocessAll2(experiments=formatlist(exptname), fillnans=True, normOD=0.8, extension='preprocessed.xlsx',  reload=True)
                                            print('successfully preprocessed\n')
                                        except Exception as err:
                                            print('failed to preprocess experiment: '+str(err))
        print(importerrors)
#                                 self.failedfiles.append([datafile, contentsfile, err])
#                                 continue
            #pickle.dump(self.data, open('xpr_startingdata.pkl', 'wb'))
            #pickle.dump(self, open('xprBackup.pkl', 'wb'))
    def getexpt(self, num):
        return self.data[self.allexperiments[num]]
    def preprocessAll(self, fillnans=True, normOD=0.8, extension='preprocessed.xlsx',  reload=True):
        '''
        preprocessExpt(expt, fillnans=True, standardize=True, mapFL=True, normOD=0.8, extension='preprocessed.xlsx' )
        preprocess all experiments to fix problems in measurement and normalize data
        '''
        preprocessed=[];
        for exp in  self.allexperiments:
            print('experiment'+exp)
            for candidate in self.supportingFL:
                print('candidate is '+candidate)
                if candidate in self.statcontents.columns and round(self.statcontents.loc[exp, candidate])==1:
                    supporting= candidate
                    print(candidate+' found for '+exp )
                    break
                else:
                    print(candidate+' not present, trying another')
                    continue
            try:
                preprocessExpt(self.data[exp], fillnans=fillnans, extension=extension, supporting=supporting, main=self.consensusFL, machine=xpr.machines[expt])
                preprocessed.append(1)
                print('experiment '+exp+'.- '+ supporting+' used for linear correction of nans in '+self.consensusFL) 
                supporting=[] #clearing up supporting to avoid runons
            except: 
                preprocessed.append(0)
                print('experiment '+exp+'.- not processed due to error') 
        print('processed experiments list:')
        print([self.allexperiments[j] for j in np.where(preprocessed)[0]])
        self.loadFresh(exptList= [self.allexperiments[j] for j in np.where(preprocessed)[0]], extension='preprocessed.xlsx') 
    def assignsupportingFL(self, experiments=[], candidates=[]):
        '''
        assignsupportingFL(self, experiments, candidates)
        Finds the supporting fluorescences to be used in each experiment
        '''
        if not candidates:
            candidates=self.supportingFL
        if not experiments:
            experiments=self.allexperiments
        else:
            experiments=formatlist(experiments)
        for exp in  experiments:
            exp=formatstring(exp)
            #print('experiment'+exp)
            candidatecounter=0
            for candidate in self.supportingFL:
                #print('candidate is '+candidate)
                if candidate in self.data[exp].datatypes: #self.statcontents.columns and round(self.statcontents.loc[exp, candidate])==1:
                    supporting= candidate
                    print(candidate+' found for '+exp )
                    self.FL[exp]['supportingFL']=supporting
                    break
                else:
                    print(candidate+' not present, trying another')
                    candidatecounter+=1
                    if candidatecounter>=np.size(self.supportingFL): #if we went through all candidates and none of them were found, then supporting fluorescence becomes unassigned.  
                        self.FL[exp]['supportingFL']=['unassigned']
                    continue
    def preprocessAll2(self, experiments=[], fillnans=True, normOD=0.8, extension='preprocessed.xlsx'):
        '''
        preprocessExpt(expt, fillnans=True, standardize=True, mapFL=True, normOD=0.8, extension='preprocessed.xlsx' )
        preprocess all experiments to fix problems in measurement and correct for plate reader differences. this version directly replaces the data arrays experiment by experiment.
        fillnans.- fills overflow values (nans) in the main FL by scaling a supporting FL (lower gain) to the scale of the main FL.
            e.g. GFP80 gives high signal to noise but some points overflow. therefore GFP60 is linearly transformed to replace GFP80
        experiments.- list of experiments to preprocess.
        extension- file with preprocessed data is created with the extension 'preprocessed.xslx' as default.    
            *WARNING* preprocessed xlsx file is not reliable as of version 4.89. 
        '''
        #pdb.set_trace()
        main=self.consensusFL
        preprocessed=[];
        if not experiments:
            experiments=self.allexperiments
        else:
            experiments=formatlist(experiments)
        for exp in  experiments:
            exp=formatstring(exp)
            if not self.FL:
                self.FL={}
            if not ppf.hasKey(self.FL, exp):
                self.FL[exp]={}
                self.FL[exp]['mainFL']=self.consensusFL
            if not ppf.hasKey(self.FL[exp], 'supportingFL'):
                self.assignsupportingFL(exp)
            if self.FL[exp]['supportingFL']=='unassigned':
                print('experiment does not have a supporting fluorescence')
                continue
            else:
                supporting=self.FL[exp]['supportingFL']
            tdata, p=preprocessExpt(self.data[exp], fillnans=fillnans, extension=extension, supporting=supporting, main=self.consensusFL, machine=self.data[exp].machine)
            print('experiment '+exp+'.- '+ supporting+' used for linear correction of nans in '+self.consensusFL+'. Processed experiment file created')
            cl=ppf.conditionList(self.data[exp])
            for j in range(0, len(cl)):
                cond=cl['media'].values[j]
                strn=cl['strain'].values[j]
                self.data[exp].d[cond][strn][main]= np.array(tdata[main].loc[self.data[exp].d[cond][strn]['plateloc'],:].T)
                #    except Exception as err:
                #        print(err)    
            print('Preprocessed data of '+exp+' updated in structure.\n')
            supporting=[] #clearing up supporting to avoid runons
            preprocessed.append(1)
            #except: 
            #    preprocessed.append(0)
            #    print('experiment '+exp+'.- not processed due to error') 
            #    supporting=[]
        print('processed experiments list:')
        print([experiments[j] for j in np.where(preprocessed)[0]])
        #this version of the function does not load fresh.
        #self.loadFresh(exptList= [self.allexperiments[j] for j in np.where(preprocessed)[0]], extension='preprocessed.xlsx') 
    def loadPickles(self, exptList=[]):
        '''
        loadExpts(exptList=self.activeExpts)
        bring out experiments in exptList to the workspace. this function is to become medular in order to not saturate
        the memory of the computer. whenever not in use, the experiments will be removed from memory
        '''
        if not exptList:
            exptList=self.activeExpts
        #clearing up previously loaded data.
        self.data={}
        #from pickle.
        for entry in exptList:
            if entry.endswith('.pkl') and (self.ignoreFiles == False or ((entry in self.ignoreFiles)==False)) and (self.onlyFiles==False or entry in self.onlyFiles) :
                print('trying to open ', entry)
                try:
                    pkl = open(self.source + '/' + entry, 'rb')
                    self.data[entry] = pickle.load(pkl,encoding= encoding)
                except:
                    pkl.close()
    def loadFresh(self, exptList=[], extension=None):
        '''
        loadFresh(exptList=self.activeExpts)
        load experiments fresh from the original files.
        '''
        if not exptList:
            exptList=self.allexperiments
        #clearing up previously loaded data.
        #self.data={}
        #self.checkallstats()
        #from pickle.
        for entry in exptList:
            #making sure the entry is not in ignorefiles or is in onlyfiles
            #if (self.ignoreFiles == False or ((entry in self.ignoreFiles)==False)) and (self.onlyFiles==False or entry in self.onlyFiles) :
            print('trying to open ', entry, 'from the information in exptInfo')
            if extension: # if there is an extension like preprocessed.xslx, then we reload with such extension.
                try:
                    self.data[entry] = pr.platereader(self.exptInfo[entry]['datapath'].split('.')[0]+extension,self.exptInfo[entry]['contentspath'])
                except Exception as err:
                    print('file with extension '+extension+' could not be added. adding to failedfiles and importing unprocessed file.')
                    self.failedfiles.append([self.exptInfo[entry]['datapath'].split('.')[0]+extension, self.exptInfo[entry]['contentspath'], err])
                    self.data[entry] = pr.platereader(self.exptInfo[entry]['datapath'],self.exptInfo[entry]['contentspath'])
            else:
                self.data[entry] = pr.platereader(self.exptInfo[entry]['datapath'],self.exptInfo[entry]['contentspath'])
    def findMachines(self):
        snref={}
        snref['1006006292']='Plate Reader 1'
        snref['1411011275']='Plate Reader 2'
        snref['906006717']='Plate Reader 3'
        snref['906006718']='Plate Reader 4'
        snref['906006719']='Plate Reader 5'
        snref['']=np.nan
        for expt in self.allexperiments:
            try:
                self.machines[expt]= self.data[expt].machine
                print('success in finding machine attribute')
            except:
                print('experiment '+expt+' lacking ''machine'' attribute. Attempting retrieval from file...')
                pth=self.exptInfo[expt]['datapath']
                #print(pth)
                xl= pd.ExcelFile(pth)
                df= xl.parse()#fixed to sheet 0 just as a temporary patch
                #print(df.head(10))
                if self.prtype == 'Tecan':
                    ##getting the serial number of the machine to know whether it is PR1 or PR2
                    snField=df.keys()[['Tecan i-control 'in j for j in df.keys()]][0] #this gives us the name of the 'infinite' field
                    try:
                        self.serialnumbers[expt]=df[snField][0].split('Serial number: ')[1] 
                    except:
                        print('problem finding serial number for experiment '+expt+'.')
                        self.serialnumbers[expt]=''
                    try:
                        self.machines[expt]= snref[self.serialnumbers[expt]]
                    except Exception as err:
                        print('Couldn''t find machine for '+expt+' :'+str(err))
                        self.machines[expt]= np.nan
    def assignFL(self, mainFL='noFL', mainFLperod='noFLperod', experiments=False):
        '''
        assignFL(self, mainFL='noFL', mainFLperod='noFLperod', experiments=False)
        
        Assigns fluorescence channel(s) mainFL as the main fluorescence channel(s) to be extracted.
        This allows to deal with heterogeneous channel names that are equivalent data.
        
        if no mainFL is declared, GFP is assigned as the main fluorescence channel. 
        such main FL channel is stored under the xpr.FL dictionary, and is called by 
        xpr.makedataframe('notime') to create fields 'FL'...  and 'FLOD'...
        You can specify a list of fluorescence channels to be found in an experiment in order of priority.
        (if one is not found, then the next one is searched)
        if experiments is given, then the fluorescence is assigned only for the given list of experiment, otherwise to all experiments.
        This routine runs by default during initial loading of the accesspr object.
        '''
        if experiments== False:
            experiments=list(self.data.keys())
        for j in experiments:
            self.FL[j]={}
        GFPFields=['GFP' in key for key in self.statcontents.keys()]
        ###Getting fields that contain GFP and are present in all experiments		
        self.consensusFLs=self.statcontents.keys()[GFPFields][np.where([np.around(sum(self.statcontents[field]))== np.size(self.statcontents,0) for field in  self.statcontents.keys()[GFPFields]])]
        if mainFL=='noFL' and mainFLperod=='noFLperod' and self.analyseFL==True:
            print('No consensus fluorescence was found. please specify fluorescence channels to work with. Defaulting to GFP.')
            fl='GFP'
            for expt in self.data.keys():
                self.FL[expt]['mainFL']=fl
                self.FL[expt]['mainFLperod']='c-'+fl+'perod'
                self.FL[expt]['mainFLperodvar']='c-'+fl+'perod var'
        else:
            if isinstance( mainFL, str):
                mainFL=[mainFL]
            if isinstance(mainFLperod, str):
                mainFLperod=[mainFLperod]
            assigned={};
            for expt in experiments:
                assigned[expt]=0
            for j in range(0, np.size(mainFL)):
                if mainFL[j]!= 'noFL':
                    for expt in experiments:
                        #self.containsstat(mainFL[j], printstats=False) #making sure it is there if at all , if the field is not a string or is unassigned
                        if ( mainFL[j] in self.statcontents.columns and self.statcontents.loc[expt, mainFL[j]]>0): # ((ppf.hasKey(self.FL[expt], 'mainFL') and ((not isinstance(self.FL[expt]['mainFL'], str)) or self.FL[expt]['mainFL']=='unassigned'))) :
                            #print('about to substitute'+self.FL[expt]['mainFL'] )
                            if assigned[expt]==0:
                                self.FL[expt]['mainFL']=mainFL[j]
                                self.FL[expt]['mainFLperod']='c-'+mainFL[j]+'perod'
                                self.FL[expt]['mainFLperodvar']='c-'+mainFL[j]+ 'perod var'
                                assigned[expt]=1
                        else:
                            self.FL[expt]['mainFL']='unassigned'
            for j in range(0, np.size(mainFLperod)):
                if mainFLperod[j]!= 'noFLperod':  ###priority to specified fluorescence per OD
                    for expt in experiments:
                        self.containsstat(mainFLperod[j], printstats=False) #making sure it is there if at all  
                        self.FL[expt]['mainFLperod']=mainFLperod[j]
                        self.FL[expt]['mainFLperodvar']='c-'+mainFL[j]+'perod var'
                        #print(mainFLperod[j]+' does not seem to exist for '+expt)
    #         if mainFL!= 'noFL':  ###priority to specified fluorescence
#             for expt in experiments:
#                 self.FL[expt]['mainFL']=mainFL
#                 if mainFL in list(self.consensusFLs)==False:
#                     print('Warning: '+mainFL+' is not a consensus Fluorescence. Consider changing or leaving unspecified')
#         if mainFLperod!= 'noFLperod':  ###priority to specified fluorescence per OD
#             for expt in experiments:
#                 self.FL[expt]['mainFLperod']=mainFLperod
#                 self.FL[expt]['mainFLperodvar']=mainFLperod+' var'
#                 if mainFLperod in list(self.consensusFLs)==False:
#                     print('Warning: '+mainFLperod+' is not a consensus Fluorescence per od. Consider changing or leaving unspecified')
# 
#         if mainFLperod in list(self.consensusFLs)==False and mainFL in list(self.consensusFLs): ##if the fl channel is there but there is no cfl for that:
#             mainFLperod='c-'+mainFL+'perod' ##we force the FLperod to match FL
#             for expt in experiments:
#                 self.FL[expt]['mainFLperod']=mainFLperod
#                 self.FL[expt]['mainFLperodvar']=mainFLperod+' var'
#         
#         #if there is no specification, we look for consensus
#         if mainFL=='noFL' and mainFLperod=='noFLperod' and np.size(self.consensusFLs)>0: # there is a consensus but was not specified in the arguments
#             flag=0;
#             allFLODs=self.consensusFLs[['perod' in j for j in list(self.consensusFLs)]] 
#             allFLs=self.consensusFLs[['perod' in j==False for j in list(self.consensusFLs)]] 
#             if np.size(allFLODs)>0: ###if there was any perod in the consensus
#                 if np.size(allFLODs)>1: ## if there is more than one, we select the first one as a default
#                     print('More than one consensus FLperod detected. Assigning '+str(allFLODs[0])+' as the default.')
#                 mainFLperod=str(allFLODs[0])
#                 for expt in self.data.keys():
#                     self.FL[expt]['mainFLperod']=mainFLperod
#                     self.FL[expt]['mainFLperodvar']=mainFLperod+'var'
#             if np.size(allFLs)>0: ###if there was any FL in the consensus
#                 if np.size(allFLs)>1: ## if there is more than one, we select the first one as a default
#                     print('More than one consensus FL detected. Assigning '+str(allFLs[0])+' as the default.')
#                 mainFL=str(allFLs[0])
#                 for expt in self.data.keys():
#                     self.FL[expt]['mainFL']=mainFL
#         else:
#             if np.size(self.consensusFLs)==0:
#                 print('No consensus fluorescence was found. please specify fluorescence channels to work with. Defaulting to GFP.')
#                 fl='GFP'
#                 for expt in self.data.keys():
#                     self.FL[expt]['mainFL']=fl
#                     self.FL[expt]['mainFLperod']='c-'+fl+'perod'
#                     self.FL[expt]['mainFLperodvar']='c-'+fl+'perodvar'
#             if mainFL!='noFL':
#                 print('\n Adding '+mainFL+' as a consensus. Please correct autofluorescence using this channel.') ### if no consensus but there was mainFL inserted
#                 for expt in self.data.keys():
#                     self.FL[expt]['mainFL']=mainFL
#                     self.FL[expt]['mainFLperod']='c-'+mainFL+'perod'
#                     self.FL[expt]['mainFLperodvar']='c-'+mainFL+'perodvar'
#             else:
#                 if mainFLperod!='noFLperod':
#                     fl= mainFLperod.split('perod')[0]
#                     try:
#                         fl=fl.split('c-')[1]
#                     except:
#                         print('channel name does not include leading c-') 
#                     print('adding '+mainFLperod+' as a consensus corrected fluorescence. assuming raw Fluorescence to be '+fl+'. Please correct if necessary') ### if no consensus but there was mainFL inserted
#                     for expt in self.data.keys():
#                         self.FL[expt]['mainFL']=fl
#                         self.FL[expt]['mainFLperod']='c-'+mainFL+'perod'
#                         self.FL[expt]['mainFLperodvar']='c-'+mainFL+'perodvar'
#                 else:
#                     print(' \nplease add fluorescence channels manually. \ncheck the contents using xpr.containsstat() and xpr.statcontents.')
                
                #                 
                #                 
                #                 
                #                 for fl in allFLs:
                #                     if 'perod' in fl:
                #                         flag=flag+1
                #                         for expt in self.data.keys():
                #                             self.FL[expt]['mainFLperod']=fl
                #                         self.consensusFLperod=fl
                #                         if fl != 'FLperod': #if it is flperod then we don't adjust the names
                #                             self.FL[expt]['mainFL']=fl.split('perod')[0].split('c-')[1] 
                #                             self.consensusFL=fl.split('perod')[0].split('c-')[1] 
                #                         self.FL[expt]['mainFLperodvar']=fl+'var'
                #                 ###second priority goes to the mainFL if it is part of the consensus.
                #                     else:                     
                #                         if mainFL in self.consensusFLs: #if there is a mainFL specified
                #                             for expt in self.data.keys():
                #                                 self.FL[expt]['mainFL']=fl #make this FL be the main fluorescence 
                #                                 if mainFLperod in consensusFLs==False: #if mainFLperod is not in the consensus but the
                #                                     self.FL[expt]['mainFLperod']='c-'+fl+'perod'
                #                                     self.FL[expt]['mainFLperodvar']='c-'+fl+'perodvar'
                #                                 else:
                #                                     self.FL[expt]['mainFLperod']=mainFLperod
                #                                     self.FL[expt]['mainFLperodvar']=mainFLperod+perodvar'
                #else: 
                    #print('no consensus Fluorescence was found in all experiments. please check the contents using xpr.containsstat() and xpr.statcontents.')
        #else:
        ###Block for when there is no same fluorescence in all experiments To be filled later
    def listcontents(self, verbose=True):
        '''
        Lists all the different conditions and strains that are contained 
        in the different experiments. Assembles xpr.allstrains, xpr.allmedia and self.allexperiments
        Arguments:
        verbose: prints the contents of every experiment to screen.
        '''
        self.allmedia=set()
        self.allstrains=set()
        self.allexperiments=[]
        self.getvariables()
        self.checkallstats()
        self.getallcontents()
        for exp in self.data.keys():
            self.allexperiments.append(exp)
            self.exptInfo[exp]= {}
            self.exptInfo[exp]['type']='datasheet'
            self.exptInfo[exp]['datapath']=self.data[exp].name+'.xlsx'
            self.exptInfo[exp]['contentspath']=self.data[exp].aname
            self.experimentDuration[exp]=self.data[exp].t[-1]
            try:
                self.exptInfo[exp]['contentspath']=self.data[exp].aname
            except:
                print('failed to find a contents file name')
            if verbose==True:
                print('\n\nExperiment ('+exp+') contains:\n------------------ ')
            for media in self.data[exp].d.keys():
                if verbose==True:
                    print('\nMedium ('+media+'): ',)
                self.allmedia.add(media)
                for strain in sorted(self.data[exp].d[media].keys()):
                    if verbose==True:
                        print( strain +',')
                    self.allstrains.add(strain)
        self.getvariables()
        self.checkallstats()
        self.getallcontents()
    def containsstat(self, stat, printstats=True):
        '''checks whether the experiments contain the given statistic stat.
        '''
        for key in self.data.keys():
            #print('condition list of ', key, ':') 
            cl= ppf.conditionList(self.data[key], excludeNull=True)
            #print( cl)
            percentage=0
            for x in range(0, np.size(cl, 0)):
                media=cl.values[x,1]
                strain=cl.values[x,2]
                if ppf.hasKey(self.data[key].d[media][strain], stat):
                    percentage=percentage+ 1/np.size(cl, 0)
                    #print(percentage)
                    self.statcontents.loc[key, stat]=percentage
                else:
                    if stat=='gr': #in the special case that we are dealing with gr, we really want to know if it is 0.
                        self.statcontents.loc[key, stat]=0
        if printstats==True:
            print(self.statcontents)
    def resetFL(self):
        '''resets the processing of fluorescence in order to start from scratch in case something went wrong'''
        for expt in self.allexperiments:
            self.data[expt].reset()
    def containssetup(self, media, strain, verbose=False, strict=True, musthave='gr'):
        '''
        Creates a list of experiments that contain the specified conditions.
        This list is used when plotting or statistics methods are called.
        Arguments:
        media: the media to be searched
        strain: the strain to be searched
        verbose: prints a little message whether a given experiment contains the setup (both media and strain)
        musthave: a statistic that must be present in the experiments in order to be returned by the function
        strict: retrieves just the experiments that contain media, strain and musthave. if False, experiments are filtered on the basis of media and strain only.
         This is relevant for some routines that require growth rate.
        '''
        self.media = media
        self.strain = strain
        self.nofod = []
        self.containslist = []
        for key in self.data.keys():
            if strict==True:
                if ppf.hasKey(self.data[key].d, media) and ppf.hasKey(self.data[key].d[media], strain) and ppf.hasKey(self.data[key].d[media][strain], musthave):
                    self.containslist.append(key)
            else:
                if ppf.hasKey(self.data[key].d, media) and ppf.hasKey(self.data[key].d[media], strain) :
                    self.containslist.append(key)
            if verbose != False: 
                    print('The experiment %s contains the requested setup.' % (key))
            else:
                if verbose != False: 
                    print('The experiment %s does not contain the requested setup.' % (key))
                
    def strainContains(self, strng):
        '''retrieves all strains whose name contains a specific identifier strng
        '''
        sd=list(self.allstrains)
        return [sd[c] for c in np.where([strng in j for j in  self.allstrains])[0]]

    def mediaContains(self, strng):
        '''retrieves all strains whose name contains a specific identifier strng
        '''
        sd=list(self.allmedia)
        return [sd[c] for c in np.where([strng in j for j in  self.allmedia])[0]]
    def getallcontents(self):
        '''
        getallcontents(self)
        creates two dataframes in the accesspr object:
            self.allconditions.- a dataframe of unique media/strain conditions in all experiments (useful for looping through all conditions once)
            self.allcontents.- a dataframe containing reating experiment, media, strain and plate locations of each condition.
            This object is useful when you want to only consider or merge certain wells or biological replicates 
        '''
        cl= pd.DataFrame()
        rl=pd.DataFrame() #dataframe of the replicates
        for key in self.data.keys():
            #print(cl)
            cl= pd.concat([cl, ppf.conditionList(self.data[key])], ignore_index=True)
            rl= pd.concat([rl, ppf.replicateList(self.data[key])], ignore_index=True)
        self.allreplicates=rl
        self.allcontents=cl
        self.allreplicates['experiment']= [re.split(r'[/\\\\]', exptname)[-1]  for exptname in self.allreplicates['experiment']]
        self.allcontents['experiment']=[re.split(r'[/\\\\]', exptname)[-1]   for exptname in self.allcontents['experiment']]
        self.allconditions= pd.DataFrame(np.column_stack([cl.values[:,1], cl.values[:,2]]), columns=['media', 'strain'])
        self.allconditions=self.allconditions.drop_duplicates()
    def correctauto(self, f=[], experiments='all',media='all', strains='all', refstrain=['WT', '229.WT', '77.WT', 'REF'], figs=False, correctOD=True, noruns=2, bd=False, no1samples=100, rewrite=False, rerun=False, correctmedia=True, mediausemean=False, ignoreneg=True):
        ''' function designed to run correctauto on all experiments, or combinations of media and strains. 
        '''
        if experiments=='all':
            experiments=self.allexperiments
        #self.containsstat('c-'+f[0]+'perod',printstats=False)
        for key in experiments: 
            #if self.statcontents.loc[key, 'c-'+f[0]+'perod']==1 and rerun==False:
            #    print('Experiment ', key, ' already contains FLperod')
            #    continue
            whichref= np.where(np.array([c in self.data[key].allstrains for c in refstrain])) #which ref of the ones added is in this experiment
            localref= np.array(refstrain)[whichref] #hopefully there is never 2 different refs
            allfl=[];
            if not localref:
                print('experiment', key, ' does not have a reference strain. going to the next.')
                continue
            else:
                localref=formatstring(localref)
            try:
                print('experiment', key, ':')
                if not f:
                    fl=[self.FL[key]['mainFL'], 'AutoFL']
                else:
                    fl=f
                    allfl.append(self.FL[key]['mainFL'])
                print('fluorescences used:')
                print(fl)
                print('\n')
                self.data[key].correctauto(f=fl, conditions=media, strains=strains, refstrain=localref, figs=figs, correctOD=correctOD, noruns=noruns, bd=bd, no1samples=no1samples, correctmedia=correctmedia , ignoreneg=ignoreneg, mediausemean=mediausemean)
                if figs==False:
                    plt.close('all')
            except Exception as err: #LinAlgErr:
                print(err)
                for e in range(2,5):
                    try:
                        print('try number '+str(e))
                        self.data[key].correctauto(f=fl, conditions=media, strains=strains, refstrain=localref, figs=figs, correctOD=correctOD, noruns=noruns, bd=bd, no1samples=no1samples, correctmedia=correctmedia, mediausemean=mediausemean, ignoreneg=ignoreneg )
                    except Exception as err2:
                        print('something went wrong. '+str(err2))
                plt.close('all')
        if rewrite==True:
            if path.isdir(self.source):
                pickle.dump(self.data[key], open(self.source + '/' +key, 'wb'))
                print('Experiment pickle file'+ self.source+key+' has been replaced.')
        [self.containsstat('c-'+j+'perod',printstats=False) for j in np.unique(allfl)]

    def statDerivative(self, dtype='FLperod', experiments='all',media='all', strains='all', rewrite=False, rerun=False):
        ''' function designed to obtain a simple gradient (derivative) of any stat inside the experiments
        '''
        if experiments=='all':
            experiments=list(self.data.keys())
        for key in experiments:
            cl=self.allcontents[self.allcontents['experiment']==key][['media','strain']]
            for j in range(0, np.size(cl,0)):
                #try:
                #print(j)
                print('experiment', key,'media', cl.values[j,0], 'strain', cl.values[j,1],  ':')
                self.data[key].d[cl.values[j,0]][cl.values[j,1]]['d/dt'+dtype]= ppf.statDerivative(self.data[key], cl.values[j,0],cl.values[j,1], dtype=dtype)
                maxd=np.nanmax(self.data[key].d[cl.values[j,0]][cl.values[j,1]]['d/dt'+dtype])
                self.data[key].d[cl.values[j,0]][cl.values[j,1]]['maxd/dt'+dtype]= maxd
                miin=np.nanmin(self.data[key].d[cl.values[j,0]][cl.values[j,1]]['d/dt'+dtype]-maxd)
                pktime=np.where(self.data[key].d[cl.values[j,0]][cl.values[j,1]]['d/dt'+dtype]==miin)
                self.data[key].d[cl.values[j,0]][cl.values[j,1]]['maxd/dt'+dtype+'time']= pktime
                #print(self.data[key].d[cl.values[j,0]][cl.values[j,1]]['d/dt'+dtype])
                #except:
                #    print('something went wrong. carrying on.')
                #    continue
        
        
    def statByStatScatter(self, media, strain,xstat='OD', ystat='GFP',experiments='all', color='black', indexrange=False):
        for m in media:
            self.containssetup(m,strain)
            experiments = self.containslist
            print(experiments)
            plt.title(xstat+' vs '+ystat)
            plt.xlabel(xstat)
            plt.ylabel(ystat)
            for exp in experiments:
                if indexrange==True:
                    plt.scatter(self.data[exp].d[m][strain][xstat][indexrange],self.data[exp].d[m][strain][ystat][indexrange], color=color)
                else:
                    plt.scatter(self.data[exp].d[m][strain][xstat],self.data[exp].d[m][strain][ystat], color=color)

    def statHist(self, media, strain, stat='GFP',experiments='all', color='black', bins=10):
        self.containssetup(media,strain)
        if experiments=='all':
            experiments = self.containslist
        plt.title('Frequency of '+ stat)
        for exp in experiments:
                plt.hist(self.data[exp].d[media][strain][stat], color=color, bins=bins)
            
    def rewriteExperiments(self, experiments='all'):
        if experiments=='all':
            experiments= self.data.keys()
        for key in experiments:
            if path.isdir(self.source):
                pickle.dump(self.data[key], open(self.source + '/' +key, 'wb'))
                print('Experiment pickle file'+ self.source+key+' has been replaced.')
                
    def nullDistribution(self, experiments='all', exptColors=0):
        ''' creates a scatterplot of the OD and Fluorescence of each experiment. the fluorecence
        chosen depends on the main fluorescence assigned to that experiment.
        '''
        if exptColors==0:
            exptColors= colorDict(keys=self.allexperiments)
        cl= ac.DFsubset(self.allconditions, 'strain', ['null'])
        cl=cl.reset_index(drop=True)
        for j in range(0, np.size(cl, 0)):
            md=cl.loc[j, 'media']
            st=cl.loc[j, 'strain']
            self.containssetup(cl.loc[j, 'media'], cl.loc[j, 'strain'], strict=False)
            expts=self.containslist
        for expt in expts:
            plt.scatter(self.data[expt].d[md][st]['OD'], self.data[expt].d[md][st][self.FL[expt]['mainFL']], color=exptColors[expt])
        plt.xlabel('OD of null')
        plt.ylabel('FL of null')
        createFigLegend(dic=exptColors)
    def getstats(self, experiments='all', conditionsDF=False, media='all', strain='all', dtype='OD', savestate=False, bd=False, cvfn='matern', esterrs=False, noruns=5, stats=True, plotodgr=False, rerun=False, figs=False, iskip=False):
        '''
        This method ensures that all experiments that are used for plotting and subsequent 
        statistical analysis have run odstats().
        It is recommended to run checkodstats() before performing further analysis to 
        access the full power of the plate reader data.
        Overwrites existing pickle files.
        '''
        #if isinstance(replicateDF, pd.DataFrame)
        #    try:
        #        expt, media, strain, plateloc= replicateDF.values[j]
        #    except:
        #        expt, media, strain= replicateDF.values[j]
        #print "Fitting requested ODs for experiments"
        #if isinstance(conditionsDF, pd.DataFrame):
        #    
        if experiments=='all':
            experiments = self.data.keys()
        if type(experiments)==str and experiments != 'all':
            experiments=[experiments]
        print('Performing non-parametric fitting of growth statistics using gaussian processes.\n')
        for key in experiments: 
            print('...\nExperiment', key, ':\n')
            if self.statcontents.loc[key, 'gr']>.99 and dtype=='OD' and rerun==False:
                print('Experiment ', key, ' already contains gr')
                continue
            else:
                try:
                    self.data[key].getstats(dtype=dtype, esterrs=esterrs, bd=bd, cvfn=cvfn, stats=stats, plotodgr=plotodgr, noruns=noruns, figs=figs, iskip=iskip)
                    if savestate==True:
                        picklefilename='./xprdata/'+self.date+'/'+'getstats/'+key+'.pkl'
                        pickle.dump(self.data[key], open(picklefilename, 'rb'))
                        print('processing state saved locally in'+picklefilename)
                except Exception as err:
                    print('something went wrong while obtaining stats for experiment'+key+': '+str(err))
            #if rewrite==True:
            #    if path.isdir(self.source):
            #        pickle.dump(self.data[key], open(self.source + '/' +key, 'wb'))
            #        print('Experiment pickle file'+ self.source+key+' has been replaced.')
        self.getvariables()
        self.checkallstats()
    def align(self, media, strain, alignFL=False):
        '''
        In the current version 'plotalign(self,media,strain,type)' sets the time vector <t0,...,tn> 
        to zero at t.where(data == max).
        In other words, the time vector is shifted to the left.
        A new key 'taligned' containing the aligned time vector is created for each experiment 
        that contains the requested setup.
        '''
        self.containssetup(media, strain, strict=True)
        #print "Aligning times for each strain from data[]"
        for key in self.containslist:
            flag=0
            if strain== 'null':
                alignedTimeGR=np.nan
                alignedTimeFL=np.nan
            else:
                centeredTimeGR=self.data[key].t[np.where(self.data[key].d[media][strain]['gr']==max(self.data[key].d[media][strain]['gr']))]
                alignedTimeGR=self.data[key].t-centeredTimeGR
                if alignFL==True:
                    mainFL=self.FL[key]['mainFL']
                    mainFLperod=self.FL[key]['mainFLperod']
                    try:
                        centeredTimeFL=self.data[key].t[np.where(self.data[key].d[media][strain][mainFLperod]==max(self.data[key].d[media][strain][mainFLperod]))]
                    except KeyError:
                        try:
                            print('Warning: experiment '+key+' does not contain corrected '+mainFLperod+'. Attempting to use raw '+mainFL+' peak') 
                            centeredTimeFL=self.data[key].t[np.where(self.data[key].d[media][strain][mainFL+' mean']==max(self.data[key].d[media][strain][mainFL+' mean']))]
                        except:
                            print('Fluorescence peak failed to be found. setting centered FL time to NaN')
                            #centeredTimeFL=np.nan #np.matlib.repmat(np.nan,np.size(xpr.data[key].t,0),1).reshape((np.size(xpr.data[key].t,0),))
                            flag=1
                    if flag==1:
                        alignedTimeFL=np.matlib.repmat(np.nan,np.size(self.data[key].t,0),1).reshape((np.size(self.data[key].t,0),))
                        self.data[key].d[media][strain]['Time centered at FL peak'] = alignedTimeFL
                    else:
                        alignedTimeFL=self.data[key].t-centeredTimeFL
                        self.data[key].d[media][strain]['Time centered at FL peak'] = alignedTimeFL
            self.data[key].d[media][strain]['Time centered at gr peak'] = alignedTimeGR

    def alignNew(self, conditions=[], alignFL=False):
        '''
        In the current version 'plotalign(self,media,strain,type)' sets the time vector <t0,...,tn> 
        to zero at t.where(data == max).
        In other words, the time vector is shifted to the left.
        A new key 'taligned' containing the aligned time vector is created for each experiment 
        that contains the requested setup.
        '''
        if not(conditions):
            conditions=self.allcontents
        for j in range(0, np.size(self.allcontents,0)+1):
            expt, media, strain, plateloc= self.allcontents.values[j] 
            flag=0
            if strain== 'null':
                alignedTimeGR=np.nan
                alignedTimeFL=np.nan
            else:
                centeredTimeGR=self.data[key].t[np.where(self.data[key].d[media][strain]['gr']==max(self.data[key].d[media][strain]['gr']))]
                alignedTimeGR=self.data[key].t-centeredTimeGR
                if alignFL==True:
                    mainFL=self.FL[key]['mainFL']
                    mainFLperod=self.FL[key]['mainFLperod']
                    try:
                        centeredTimeFL=self.data[key].t[np.where(self.data[key].d[media][strain][mainFLperod]==max(self.data[key].d[media][strain][mainFLperod]))]
                    except KeyError:
                        try:
                            print('Warning: experiment '+key+' does not contain corrected '+mainFLperod+'. Attempting to use raw '+mainFL+' peak') 
                            centeredTimeFL=self.data[key].t[np.where(self.data[key].d[media][strain][mainFL+' mean']==max(self.data[key].d[media][strain][mainFL+' mean']))]
                        except:
                            print('Fluorescence peak failed to be found. setting centered FL time to NaN')
                            #centeredTimeFL=np.nan #np.matlib.repmat(np.nan,np.size(xpr.data[key].t,0),1).reshape((np.size(xpr.data[key].t,0),))
                            flag=1
                    if flag==1:
                        alignedTimeFL=np.matlib.repmat(np.nan,np.size(self.data[key].t,0),1).reshape((np.size(self.data[key].t,0),))
                        self.data[key].d[media][strain]['Time centered at FL peak'] = alignedTimeFL
                    else:
                        alignedTimeFL=self.data[key].t-centeredTimeFL
                        self.data[key].d[media][strain]['Time centered at FL peak'] = alignedTimeFL
            self.data[key].d[media][strain]['Time centered at gr peak'] = alignedTimeGR

    def alignAll(self, rerun=False, alignFL=False):
        '''
        In the current version 'plotalign(self,media,strain,type)' sets the time vector <t0,...,tn> 
        to zero at t.where(data == max).
        In other words, the time vector is shifted to the left.
        A new key 'taligned' containing the aligned time vector is created for each experiment 
        that contains the requested setup.
        '''
        if self.aligned==False or rerun==True:
            for expt in list(self.data.keys()):
                print('aligning experiment ', expt)
                cl=ppf.conditionList(self.data[expt])[['media', 'strain']]
                #print(cl)
                media= cl.values[:,0]
                strains= cl.values[:,1]
                for x in range(0, np.size(cl.values, 0)):
                    if strains[x]=='null':
                        continue ###this is to be replaced by the use of raw fluorescence
                    else:
                        self.align(media[x], strains[x], alignFL=alignFL)
                        #print('aligning media '+media[x]+';strain: '+strains[x])
                        #plt.plot(media[x], strains[x], alignFL=alignFL)
            self.containsstat('gr',printstats=False)
            self.containsstat('Time centered at gr peak',printstats=False)
            self.containsstat('Time centered at FL peak',printstats=False)
            if 'Time centered at gr peak' in list(self.statcontents.columns.values):
                grAlignedFraction=self.statcontents['Time centered at gr peak'].sum()/ np.size(self.statcontents['Time centered at gr peak'])
                print( str(self.statcontents['Time centered at gr peak'].sum())+" out of "+str(np.size(self.statcontents['Time centered at gr peak']))+"experiments aligned. \n Tips: \n .getstats() calculates growth statistics from all experiments. \n.statcontents lets you see which experiments need to get have gr or FLperod. \n.alignAll(rerun=True) to try aligning again")
            else:
                grAlignedFraction=0
            if grAlignedFraction>.95:
                self.aligned=True
                print('Experiments aligned successfully.')
            else:
                self.aligned=False
        else:
            print('Experiments have already been aligned. to realign, try rerun=True')
            self.getvariables()
            self.checkallstats()
    def plotrepmean(self, media=False, strain=False, conditionsDF=False, experiments='all', ignoreExps=False, dtype='OD', col='Black', alpha=0.2, normalise=False, excludeFirst=0, excludeLast=-1, bootstrap=5, centeringVariable='time', factor=False, factorColors=False, factorMarkers=False, loc='upper left', markevery=30, addLegend=True):
        markers=['.', 'o', 's', '^', '+', 'v', '*', 'p', '<', '>', 'h', 'x', 'D']*100
        '''
        intepdata, factorcolors, factormarkers=plotrepmean(self, media=False, strain=False, conditionsDF=False, experiments='all', ignoreExps=False, dtype='OD', col='Black', alpha=0.2, normalise=False, excludeFirst=0, excludeLast=-1, bootstrap=5, centeringVariable='time', factor=False, factorColors=False, factorMarkers=False, loc='upper left')
        Interpolates and bootstraps replicates of the same condition, and plots the mean and standard deviation of such bootstraps. 
        **TIP** 
        You can provide a conditions Dataframe (conditionsDF=DF)  in the format of self.allreplicates, and provide a factor=FACTOR, 
        to subdivide all those replicates by their category in FACTOR. e.g. factor='strain' will group the replicates by strain
        '''
        if centeringVariable=='Time centered at gr peak':
            print('removing null and expts without estimated growth rate')
            if isinstance(strain, list):
                strain.remove('null')
        if isinstance(conditionsDF, pd.DataFrame):
            conditionsDF= conditionsDF[conditionsDF['strain'].values != 'null']
            conditionsDF=conditionsDF.reset_index(drop=True)
            if centeringVariable=='Time centered at gr peak' or dtype=='gr':
                print('removing expts without estimated growth rate')
                self.checkallstats()
                contents=self.statcontents
                contents=contents.reset_index()
                expswithgr=contents[contents['gr']>0.99]['index'].values
                conditionsDF=conditionsDF[[j  in expswithgr for j in conditionsDF['experiment']]]
                conditionsDF=conditionsDF.reset_index(drop=True)
            #print(conditionsDF)
        if experiments=='all':
            experiments=self.allexperiments
        else:
            if isinstance(experiments, str):
                experiments=[experiments]
        if isinstance(ignoreExps, list):
            [experiments.remove(j) for j in ignoreExps]
        if isinstance(media, str):
            media=[media]
        if isinstance(strain, str):
            strain=[strain]
        if dtype=='':
            try:
                dtype=self.consensusFLperod
            except:
                dtype='OD'
        if np.size(media)==1:
            mediaName=media
            media=[media]
        cv=[]
        flag=0
        if media and strain: #we 
            for m in media:
                for s in strain:
                    conditionsDF=DFsubset(DFsubset(DFsubset(self.allreplicates, 'media', [m]), 'strain', s), 'experiment', experiments)
                    #print('processing '+m)
                    interpolated= self.interptimesnew(replicateMatrix=conditionsDF, dtype=dtype, centeringVariable=centeringVariable)
                    interpolated[dtype]=interpolated[dtype][excludeFirst:excludeLast]
                    interpolated[centeringVariable]=interpolated[centeringVariable][excludeFirst:excludeLast]
                    if normalise==True:
                        interpolated[dtype]=interpolated[dtype]/np.nanmax(flatten(interpolated[dtype]))
                    interpolatedFinal=bootstrapInterps(interpolated, centeringVariable, dtype, bootstrap=bootstrap, col=col, alpha=alpha)
                    return(interpolatedFinal) 
        if isinstance(conditionsDF, pd.DataFrame): #if the conditions matrix is 
            interpolated= self.interptimesnew( replicateMatrix=conditionsDF, dtype=dtype, centeringVariable=centeringVariable)
            interpolated[dtype]=interpolated[dtype][excludeFirst:excludeLast]
            interpolated[centeringVariable]=interpolated[centeringVariable][excludeFirst:excludeLast]
            if normalise==True:
                interpolated[dtype]=interpolated[dtype]/np.nanmax(flatten(interpolated[dtype])) 
            if isinstance(factor, str): 
            ### This if block interpolates of full set of conditions together, but are plotted by group based on variable factor. Lucia Bandiera 20180501
                createdColorsFlag=0
                createdMarkersFlag=0
                if not(isinstance(factorColors, dict)):  #if the color dictionary never existed we create it but remember we did with the flag
                    factorColors=ppf.colorDict(keys=np.unique(conditionsDF[factor].values), colors=colors.strongColors*100)
                    createdColorsFlag=1
                if not(isinstance(factorMarkers, dict)): #if the marker dictionary never existed we create it but remember we did with the flag
                    factorMarkers=ppf.colorDict(keys=np.unique(conditionsDF[factor].values), colors=markers) #this dict assembly works well for markers so it is great
                    createdMarkersFlag=1
                for i_fact,fact in enumerate(conditionsDF[factor].unique()):# iterate over the factor levels
                    print(fact)
                    i_DF = conditionsDF.index[conditionsDF[factor]== fact].get_values() # find the indexes of the dataframe in which the factor level occurs
                    i_interpolatedCols = [conditionsDF.index.get_loc(i_e) for i_e in i_DF] # find the corresponding position in the data frame (equals the number of columns in interpolate)
                    print(i_interpolatedCols)
                    interpolatedToUse={}
                    interpolatedToUse[dtype] = interpolated[dtype][:,np.where(conditionsDF[factor]== fact)[0]]
                    interpolatedToUse[centeringVariable]=interpolated[centeringVariable]
                    #final object here is a dictionary for each factor in order to obtain separate means and sds for each
                    #if not(isinstance(factorColors, dict)): #if the dictionary never existed we create it but remember we did with the flag
                    #    factorColors={}
                    #    factorMarkers={}
                    #    createdColorsFlag=0
                    #else:
                    #    createdColorsFlag=1
                    if ppf.hasKey(factorColors,fact)==False: 
                        #if there is no factorColors or if the argument does not contain a color for the given factor we make one
                        factorColors[fact]=(colors.strongColors*100)[i_fact] 
                    if ppf.hasKey(factorMarkers,fact)==False: 
                        factorMarkers[fact]=(markers*100)[i_fact]
                    interpolatedFinal=bootstrapInterps(interpolatedToUse, centeringVariable, dtype, bootstrap=bootstrap, col=factorColors[fact], alpha=alpha, marker=factorMarkers[fact], markevery=markevery)
                    interpolated['Group_'+fact+'_'+dtype+' mean']=interpolatedFinal[dtype+' mean']
                    if bootstrap>0:
                        st=dtype+'bootstrapsd'
                    else:
                        st=dtype+'sd'
                    interpolated['Group_'+fact+'_'+st]=interpolatedFinal[st]
                if addLegend==True:
                    ppf.createFigLegend(keys=list(factorColors.keys()), cols= list(factorColors.values()), markers=list(factorMarkers.values()), loc=loc ) #creating a legend for the figure
                return(interpolated, factorColors, factorMarkers) 
            else:
                    #interpolated final contains the mean and sd that were plotted
                    interpolatedFinal=bootstrapInterps(interpolated, centeringVariable, dtype, bootstrap=bootstrap, col=col, alpha=alpha)
                    return(interpolatedFinal)

    def getMediaValues(self, df=None):
        mediaValue=[]
        if df is None:
            df=self.extractAllInfo(excludeNull=False)
        for j in range(0, np.size(df['media'])):
            try:
                #print( df['media'][j] + 'is'+ str(self.mediaValue[df['media'][j]]))
                mediaValue.append(self.mediaValue[df['media'][j]])
            except:
                mediaValue.append(np.nan)
        df['mediaValue']=mediaValue
        return df
    def colorScatter(self, media, strain, experiments=False, xstat=False, ystat='FLperod', colorBy='d/dtgr', symmetric=True, cmap='bwr',nbins=40, extendBy=2, alpha=1, markersize=12, marker='o', addLegend=False, vmin=0, vmax=0, xlabel=1, ylabel=1):
        if experiments==False:
            self.containssetup(media,strain)
            experiments=self.containslist
        lastexpt=experiments[-1]
        for n in range(0, np.size(experiments)-1):
            colorScatter(self.data[experiments[n]], media, strain, xstat=xstat, ystat=ystat, colorBy=colorBy, symmetric=symmetric, cmap=cmap,nbins=nbins, extendBy=extendBy, alpha=alpha, markersize=markersize, marker=marker, addLegend=False, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel)
        colorScatter(self.data[lastexpt], media, strain, xstat=xstat, ystat=ystat, colorBy=colorBy, symmetric=symmetric, cmap=cmap,nbins=nbins, extendBy=extendBy, alpha=alpha, markersize=markersize, marker=marker, addLegend=addLegend, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel)

    def plotDefault(self, media, strain, typeofdata, aligned=False):
        D = self.data
        self.media = media
        self.strain = strain
        self.typeofdata = typeofdata
        plt.figure()
        plt.title(self.typeofdata + ' of' + self.strain + ' in ' + \
                      self.media + ' over time', color = 'b')
        if aligned:
            plt.xlabel('Aligned time in h')
        else:
            plt.ylabel(self.typeofdata + ' in a.u.')
        legends = []
        if aligned:
            self.align(self.media, self.strain)
            for key in self.containslist: 
                if ppf.hasKey(D[key].d, media) and ppf.hasKey(D[key].d[media], strain)\
                    and ppf.hasKey(D[key].d[media][strain], typeofdata):
                    plt.plot(self.aligned[key], D[key].d[self.media][self.strain][self.typeofdata])
                    legends.append([key])
                else: continue
        else:
            for key in self.containslist: 
                if ppf.hasKey(D[key].d, media) and ppf.hasKey(D[key].d[media], strain)\
                    and ppf.hasKey(D[key].d[media][strain], typeofdata):
                    plt.plot(D[key].t, D[key].d[self.media][self.strain][self.typeofdata])
                    legends.append([key])
                else: continue
        plt.legend(legends, loc = 'upper right')
        plt.show(block=False)

    def plotrawreplicates(self, media, strain, dtype='OD',xlim=False, ylim=False, experiments='all', exptColors=False, addLegend=True, xstat=[], color=False):
        markers=['o', 's', '^', 'x', '*', 'h','+', 'v', 'p', '>','o', 's', '^', 'x', '*', 'h','+', 'v', 'p', '>','o', 's', '^', 'x', '*', 'h','+', 'v', 'p', '>']
        self.containssetup(media, strain, strict=False)
        if exptColors==False:
            if color:
                exptColors= ppf.colorDict(keys= self.allexperiments, colors= color* np.size(self.allexperiments))
            else:
                exptColors=dict(zip(self.containslist, ppf.randomColors(np.size(self.containslist))))
        patches=[]
        for x in self.containslist:
            if experiments!='all' and (x in experiments)==False:
                continue
            patches.append(pch.Patch(color=exptColors[x]))
            arr=self.data[x].d[media][strain][dtype]
            plateloc=self.data[x].d[media][strain]['plateloc']
            if not xstat:
                arrx=self.data[x].t
            else:
                arrx=self.data[x].d[media][strain][xstat]
            [plt.plot(arrx, arr[:,k], label= x[0:8]+'...'+x[-4:]+' '+plateloc[k], marker=markers[k], markevery=10, color= exptColors[x]) for k in range(0, len(plateloc))]
        if addLegend==True:
            #plt.figlegend(patches, self.containslist, 'upper right')
            plt.legend()
        plt.title('Replicates of '+strain+ ' in '+media)
        if ylim:
            plt.ylim(ylim)
        if xlim:
            plt.xlim(xlim)
        return exptColors
    def timeStatAll(self, times, media='all', strains='all', xstat='time', dtype='OD', scale=False, subtractBackground=False):
        df=pd.DataFrame()
        if scale != False:
            maxdf=self.extractallinfonew(excludeNull=True) #extract all info to get maxima.
            if media !='all': #if media subset is entered then filter dataframe by that subset
                maxdf=DFsubset(maxdf, 'media', media)
            if strains !='all':#if strain subset is entered then filter dataframe by that subset
                maxdf=DFsubset(maxdf, 'strain', strains)
            if dtype=='FLperod': ####depending on the type, the max will be given by different properties.
                mx= np.max(maxdf['FLPeak'])
            if dtype=='gr':
                mx= np.max(maxdf['maxGR'])
            if dtype=='OD':
                mx= np.max(maxdf['FinalOD'])
        else:
            mx=False
        for j in range(0, np.size(self.allreplicates,0)):
            e,m,s,pl=self.allreplicates.values[j,0],self.allreplicates.values[j,1],self.allreplicates.values[j,2],self.allreplicates.values[j,3]
            #print('e '+e+' '+'m '+m+' s '+s+'\n')
            if self.allreplicates.values[j,2]=='null':
                continue
            else:
                df=pd.concat([df, self.timeStat(experiment=e,  media=m,strain=s, plateloc=pl, times=times, dtype=dtype, scale=scale, xstat=xstat, max=mx)])
            #except:
            #    print('condition ', self.allconditions.values[j,1], ' in ' ,self.allconditions.values[j,0], 'showed extraction problems.')
        df.index= range(0,np.size(df,0))
        if media !='all': #if media subset is entered then filter dataframe by that subset
            df=DFsubset(df, 'media', media)
        if strains !='all':#if strain subset is entered then filter dataframe by that subset
            df=DFsubset(df, 'strain', strains)
        return df
    def timeStat(self,  media, strain, experiment=[], plateloc=[], times= [4], dtype='OD', includedescriptors=True, scale=False, max=False, subtractBackground=False, background=False, xstat='time' ):
        '''
        df=timeStat(self, media, strain, times=[4], dtype='FLperod', includeMediaStrain=True):
        Generates a dataframe where columns are media, strain and the values of dtype at times (one column per element in times).
        and and each row is the value of dtype at those times for one replicate. 
        The times are absolute from the beginning of each experiment. Default is at 4 hrs.
        '''
        self.containssetup(media,strain, strict=False) #finding the experiments with a given media and strain
        if not experiment:
            expts=self.containslist ##here are the compliant experiments.
        else:
            expts= formatlist(experiment)
        cols= ['experiment','media', 'strain', 'plateloc']+ formatlist(times)
        isWithinBound= np.array([j<=self.interpLimit for j in formatval(times)]) ###true or false array saying 
        #print(type(isWithinBound))
        #whether the times fall within the length of the shortest experiment
        if np.size(np.where(isWithinBound==False))>0: #warn that times out of the bound will be excluded.
            print('Warning: times above ', self.interpLimit, ' are being excluded')
        ###we immediately get rid of the times outside the duration of the shortest experiment.
        times=formatnum(times)
        if includedescriptors==True: ##whether the user wants to iinclude the media and strain columns
            fin=pd.DataFrame( columns=cols) #preparing the output dataframe with media and strain
        else:
            fin=pd.DataFrame( columns= times) #preparing the output dataframe w/o media and strain
        for j in range(0, np.size(expts)): ##loop through experiment names
            #print(expts[j])
            #try:
            if not ppf.hasKey(self.data[expts[j]].d[media][strain], dtype):
                datavals= np.zeros([np.size(self.data[expts[j]].t, 0), 1])*np.nan
                continue
            if np.ndim( self.data[expts[j]].d[media][strain][dtype])>1:
                #print('enter1. multiple columns')
                if not plateloc or isinstance(plateloc, list): #if plateloc is a list then average all locs
                    #print('enter2 there is no plateloc or plateloc is a string')
                    datavals= self.data[expts[j]].d[media][strain][dtype].mean(1)
                else:
                    if isinstance(plateloc, str): #if it is a string then it is only one well. 
                        #print('enter3 there is plateloc string')
                        plateloclist=self.data[expts[j]].d[media][strain]['plateloc']
                        if plateloc in plateloclist: #if that well is in the list (has not been removed)
                            #print('enter4 plateloc is in the plateloc list. datavals assigned')
                            datavals=  self.data[expts[j]].d[media][strain][dtype][:, findloc(plateloc, plateloclist)]
                        else:
                            if not plateloc in plateloclist or plateloc in self.data[expts[j]].ignoredwells or self.ignoredwells[data[expts[j]]]:
                                #print('enter5 well seems ignored. filling nans')
                                datavals= np.zeros([np.size(self.data[expts[j]].t, 0), 1])*np.nan
            else:
                #print('enter6 one column. datavals assigned')
                datavals=self.data[expts[j]].d[media][strain][dtype]
            if subtractBackground==True: ##if the user wants to subtract a minimum value of the signal
                if background == False: ##if no general minimum is given we use the minimum level of the signal
                    background= np.min(datavals)
            else:#if subtraction of the minimum is not needed we make background 0.
                background=0
            datavals=datavals-background #whatever the value of background is at this point, we subtract it from the signal datavals.
            if scale==True: ##if the user wants to scale the response relative to a maximum
                if max != False: ##if no general maximum is given we use the maximum level in the signal
                    f = scint.interp1d(self.data[expts[j]].d[media][strain][xstat],  datavals/max, bounds_error=False, fill_value=np.nan) #create interpolation function
                else:#if a max is given we divide the signal by it.
                    f = scint.interp1d(self.data[expts[j]].d[media][strain][xstat],  datavals/np.max(datavals) , bounds_error=False, fill_value=np.nan) #create interpolation function
            else: #if scaling is not required then we just proceed with the original time series.
                f = scint.interp1d(self.data[expts[j]].d[media][strain][xstat],  datavals, bounds_error=False, fill_value=np.nan) #create interpolation function
            fin.loc[j,times]= [j for j in f(times)] #the dataframe at the given expt's index and  and at the times columns requested will be obtained by interpolating from the data
            if includedescriptors==True: # add media and strain if required
                fin.loc[j,'experiment']= expts[j]
                fin.loc[j,'plateloc']=plateloc
                fin.loc[j,'media']= media
                fin.loc[j,'strain']= strain
        fin[times]=fin[times].astype('float')
        return fin        
    def replicateLocations(self):
        if not self.allconditions.empty:
            daf=pd.DataFrame(index=[self.allconditions['strain'], self.allconditions['media']], columns=self.allexperiments)
            daf=daf.fillna(0)
            for j in daf.index:
                self.containssetup(media=j[1], strain=j[0], strict=False)
                for c in self.containslist:
                    daf.loc[j, c]=1 
            self.conditionLocTable=daf
            self.numReplicates=pd.DataFrame(columns=['strain', 'media', 'numReplicates'])
            self.numReplicates['numReplicates']=daf.sum(1).values
            self.numReplicates['strain']=self.allconditions['strain']
            self.numReplicates['media']=self.allconditions['media']
    def makedataframe(self, type='notime', times=False, dtype='OD', xstat='time', conditionsDF=False):
        ''' makedataframe(self, type='notime', times=False, dtype='OD', xstat='time')
        Exports a dataframe with data in 3 different formats:
        type= the format of the dataframe to be exported. three possible types:
            'notime'.-  point statistics for each replicate in the experiments. anything that can be exported will be exported.
            Interpolated time statistics: statistics sampled at equivalent times across experiments and replicates...
            'timecol'.- specify timpepoints times and get statistic dtype evaluated at points times of  variable xstat (defaults to times)
            times vector must be specified. each  is a replicate and columns are descriptor variables of such replicates an the values at times times.
        'timerow'.-  retrieve a dataframe of time traces (in rows) of data type dtype, each column being a different time series. times can be specified to retrieve specific points but it is not necessary.
        times.- the times at which extract the value of xstat for each replicate when using type timerow and timecol
        xstat= the time varying statistic that serves as a reference for interpolation. defaults to time.
            other potential options for it are 'Time centered at gr peak' for aligned values, or other time driven variables like OD for example, as long as there is one per
            condition.
        '''
        if not isinstance(conditionsDF, (pd.DataFrame)):
            conditionsDF=self.allreplicates
        if type=='notime':
            df=self.extractallinfonew(replicateDF=conditionsDF)
        if type=='timecol':
            if not isinstance(times, (int, list, np.ndarray)):
                raise ValueError('For timerow please provide a list of times at which to extract data: times=[x1, x2, x3].')
            else:
                df= self.timeStatAll(times=times, dtype=dtype, xstat=xstat)
        if type=='timerow':
            if not isinstance(times, (int, list, np.ndarray)):
                raise ValueError('For timerow please provide a list of times at which to extract data: times=[x1, x2, x3].')
            alldtypes= formatlist(dtype)
            df=pd.DataFrame();
            c=1;
            for j in alldtypes:
                if c==1 and not isinstance(times, (int, list, np.ndarray)):
                    interpolated=self.interptimesnew(replicateMatrix=conditionsDF, dtype=j, centeringVariable=xstat)
                    times= interpolated['time']
                    df=pd.DataFrame(interpolated[dtype],index=times, columns=inerpolated['names'])
                else:
                    interpolated=self.interptimesnew(replicateMatrix=conditionsDF, dtype=j, centeringVariable=xstat, xvals=times)
                    df=pd.concat([df, pd.DataFrame(interpolated[j],index=interpolated[xstat], columns=interpolated['names'])], axis=1)
                
                c+=1
        return df
    def interptimesnew(self, replicateMatrix=False, dtype='OD', centeringVariable='time', descriptors=False, xvals=False):    
        '''
        interpTimes(self, media, strain, dtype='OD', centeringVariable='time')   
        interpolate all replicates across the same time scale.
        replicateMatrix is a matrix of the format of self.allreplicates with the specific replicates to interpolate
        dtype: data type to extract. it must be present in all wells.
        centeringVariable: the variable used for the x axis
        descriptors: the final object will contain full descriptor factors for each condition:  expt, media, strain, plateloc
        '''
        if not isinstance(replicateMatrix, pd.DataFrame):
            replicateMatrix=self.allreplicates
        #replicateMatrix=replicateMatrix[replicateMatrix['strain']!='null'] ###removing annoying nulls!!! which should have nans in all values yet they don't because.
        #replicateMatrix=replicateMatrix.reset_index(drop=True)
        #print(replicateMatrix)
        interpRange=[]
        maxLengths=[]
        startPoints=[]
        removelist=[]
        adjustedTimes=dict()
        adjustedTimes[dtype]=[]
        adjustedTimes[centeringVariable]=[] #working with ordered lists
        filter=[];
        for j in range(0, np.size(replicateMatrix, 0)): ###retireving the limiting values for the interpolation amongst all experiments.
            try:
                expt, media, strain, plateloc=replicateMatrix.values[j, [0,1,2,3]]
                maxLengths.append(np.around(self.data[expt].d[media][strain][centeringVariable][-1],2))
                startPoints.append(np.around(self.data[expt].d[media][strain][centeringVariable][0],2))
            except Exception as err:
                print('problem measuring length of data trace. '+str(err)+'\n'+centeringVariable+' may not exist for '+expt+' '+media+' '+strain+'\n excluding expriment '+expt)
                removelist.append(expt)
                continue
        removelist=np.unique(removelist)
        #for j in removelist:
        #    replicateMatrix=replicateMatrix[replicateMatrix['experiment']!=j]
        #    replicateMatrix=replicateMatrix.reset_index(drop=True)
        interpRange=[np.max(np.array(startPoints)), np.min(np.array(maxLengths))]
        print('setting interpolation limit with new information')
        self.interpLimit=np.min(np.array(maxLengths))
        func= lambda x: fitsBounds(x, interpRange) ### this lambda creates function func which will evaluate whether x falls within the interpRange
        namesarray=[]
        exptsarray=[]
        mediaarray=[]
        strainsarray=[]
        platelocarray=[]
        dtypearray=[]
        flag=[]
        for j in range(0, np.size(replicateMatrix, 0)):
            #print('processing replicate number', j)
            #print(j)
            expt, media, strain, plateloc=replicateMatrix.values[j, [0,1,2,3]]
            #finding the points that fit the range
            fitpoints=np.where([func(x) for x in self.data[expt].d[media][strain][centeringVariable]])
            if j==0:
                firstx=np.around(self.data[expt].d[media][strain][centeringVariable][fitpoints],2)
                finalDict={}
                finalDict[dtype]=np.zeros([len(firstx), np.size(replicateMatrix,0)]) #working with ordered lists
            namesarray.append(expt+' '+media+' '+strain+' '+' '+plateloc+ ' '+ dtype)
            if descriptors==True:
                exptsarray.append(expt)
                mediaarray.append(media)
                strainsarray.append(strain)
                platelocarray.append(plateloc)
                dtypearray.append(dtype)
            #print(np.shape(self.data[expt]. d[media][strain][dtype])[1])
            #Index in the data matrix that corresponds to this specific well
            if not(plateloc in self.data[expt].ignoredwells):
                platelocIndex=np.where([plateloc==k for k in self.data[expt].d[media][strain]['plateloc']])[0][0]
            else:
                adjustedTimes[dtype].append(np.zeros([np.size(fitpoints),1])*np.nan)
                adjustedTimes[centeringVariable].append(np.zeros([np.size(fitpoints),1])*np.nan)
                continue
            if not(dtype.endswith(' mean')) and not(dtype.endswith('gr')) and not(dtype.endswith('perod')) and not(dtype.endswith('var')):
                adjustedTimes[dtype].append(self.data[expt].d[media][strain][dtype][fitpoints, platelocIndex])
            else:
                try:
                    adjustedTimes[dtype].append(self.data[expt].d[media][strain][dtype][fitpoints])
                except Exception as err:
                    print('Error occured during processing of'+expt+' '+media+' '+strain+' '+' '+plateloc+ ' '+ dtype+':\n'+str(err)+'\n filling vector with nans')
                    adjustedTimes[dtype].append(np.zeros([np.size(fitpoints),1])*np.nan)
                    #adjustedTimes[centeringVariable].append(np.zeros([np.size(fitpoints),1])*np.nan)
            adjustedTimes[centeringVariable].append(np.around(self.data[expt].d[media][strain][centeringVariable][fitpoints],2))
        #finalDict={};
        if isinstance(xvals, (bool)): #if it is not a boolean
            xvals=np.around(adjustedTimes[centeringVariable][0],2)  #arbitrarily taking the times of the first condition as a reference
        finalDict[dtype]=np.empty([np.size(xvals), len(adjustedTimes[dtype])], dtype=None)
        finalDict[centeringVariable]= np.array(xvals)
        for j in range(0, len(adjustedTimes[dtype])): 
            try:
                fint=scint.interp1d(adjustedTimes[centeringVariable][j],adjustedTimes[dtype][j], bounds_error=False, fill_value=np.nan) #interpolate times j and response j
                finalDict[dtype][:, j]=[fint(q) for q in xvals]
            except:
                finalDict[dtype][:, j]=np.empty([np.size(xvals)])*np.nan
        finalDict['names']=namesarray
        if descriptors==True:
            finalDict['experiment']=exptsarray
            finalDict['media']=mediaarray
            finalDict['strain']=strainsarray
            finalDict['plateloc']=platelocarray
            finalDict['dtype']=dtypearray
        return finalDict

    def addLegend(self, strains=False, media=False, strainColors=False, mediaColors=False):
        patches=[]
        legends=[]
        if strains==True:
            legends= self.allstrains[self.allstrains != 'null']
            legends= legends[ ['WT' in x for x in legends]]
        patches= [pch.Patch(color=hxt4col),pch.Patch(color=std1col),pch.Patch(color=rgt2col)]
        legends=['409.Hxt4', '403.std1', '506.rgt2']
        plt.figlegend(patches, legends, 'upper right')
        #             adjustedTimes[expt]= self.data[expt][np.where(self.data[expt].t<self.data[shortestExperiment].t[-1])]
        #maxLengths[m]=self.data[experimentList[m]].t[-1]
        # for m in range(0, np.size(experimentList)-1): ##we get the max time duration of all the experiments
        #             maxLengths[m]=self.data[experimentList[m]].t[-1]
        #         #now we look for the experiment with the minimum length
        #         shortestExperiment=experimentList[maxLengths.index(min(maxLengths)[0])]
        #         #now using the shortest experiment as a reference, we get the timepoints that are within the hour-wise length of that experiment.
        #         sharedTimes=dict()
        #         for expt in experimentList:
        #             adjustedTimes[expt]= self.data[expt][np.where(self.data[expt].t<self.data[shortestExperiment].t[-1])]
        #             
            
        #self.data[shortestExperiment].t[-1] 
        #first we take the real time point indices that do not exceed the upper limit, which is in hours. say, if indices 161-200 are above, then you consider 1-160
        #then we interpolate the timepoints
    def extractallinfonew(self, replicateDF=[], growthstats=[], flstats=[]):
        self.getvariables() #making sure all the latest added variables are available.
        if not isinstance(replicateDF, pd.DataFrame):
            replicateDF= self.allcontents
        df=pd.DataFrame(index= replicateDF.index, columns=self.extractionFields)
        #for j in range(0, np.size(replicateDF,0)):
        #    expt, media, strain, plateloc=  formatstring(replicateDF.values[j][grepcolumns('experiment', df)]),formatstring(replicateDF.values[j][grepcolumns('media', df)]),formatstring(replicateDF.values[j][grepcolumns( 'strain', df)]),formatstring(replicateDF.values[j][grepcolumns('plateloc', df)])             ##Pending: generate modular extraction functions depending on expt, media, strain, plateloc
        filler=[np.nan for j in range(0, len(replicateDF))]
        if not growthstats:
            growthstats=self.extractionFieldsScalar
        if not flstats:
            flstats= self.extractionFieldsFL
        growthdata={};
        for stat in  growthstats:   
            try:
                growthdata[stat]=[self.extractStat(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], stat) for j in range(0, len(replicateDF))]
            except:
                growthdata[stat]=[self.extractStat(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], stat) for j in range(0, len(replicateDF))]
        growthdata['gr auc']=[self.extractGRAUC(replicateDF.loc[j,:].values[0],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], replicateDF.loc[j,:].values[3]) for j in range(0, len(replicateDF))]
        growthdata['initial OD']=[self.extractInitialODFull(  replicateDF.loc[j,:].values[0],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], replicateDF.loc[j,:].values[3]) for j in range(0, len(replicateDF))]
        growthdata['final OD raw']=[self.extractFinalODFull(  replicateDF.loc[j,:].values[0],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], replicateDF.loc[j,:].values[3]) for j in range(0, len(replicateDF))]
        growthdata['max OD raw']=[self.extractMaxODFull(  replicateDF.loc[j,:].values[0],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], replicateDF.loc[j,:].values[3]) for j in range(0, len(replicateDF))]
        nfl=self.consensusFL
        FLData={};
        getalstFL=lambda expt, m,s,fl: alignStats(expt, m, s, dtype=fl)
        treatmulti=lambda x: x[0] if np.size(x)>1 else x
        #print('entering problematic part')
        alloutsFL=[getalstFL(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1],replicateDF.loc[j,:].values[2], self.FL[replicateDF.loc[j,:].values[0]]['mainFL']) for j in range(0, len(replicateDF))]
        #print('transforming into floats')
        for field in formatlist(flstats):
            #print('now for '+field)
            #pdb.set_trace()
            FLData[field]=[np.float(fixempty(treatmulti(alloutsFL[j][field]))) for j in range(0, len(replicateDF))]
        nflod=self.consensusFLperod
        FLODData={};
        getalstFLOD=lambda expt, m,s,fl: alignStats(expt, m, s, dtype=fl)
        #aligned fluorescence statistics calculation
        alloutsFLOD=[getalstFLOD(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1],replicateDF.loc[j,:].values[2],self.FL[replicateDF.loc[j,:].values[0]]['mainFLperod']) for j in range(0, len(replicateDF))]
        alloutsFLODvar=[getalstFLOD(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1],replicateDF.loc[j,:].values[2],self.FL[replicateDF.loc[j,:].values[0]]['mainFLperodvar']) for j in range(0, len(replicateDF))]
        for field in flstats:
            try:
                FLODData['FLOD'+field]=[np.float(treatmulti(alloutsFLOD[j][field])) for j in range(0, len(replicateDF))]   
            except:
                FLODData['FLOD'+field]=filler
        return pd.concat([replicateDF, pd.DataFrame.from_dict(mergedicts([growthdata, FLData, FLODData])).drop('plateloc', axis=1)], axis=1 )
    def fixNaNs(self, dtype='GFP80', repairWith='GFP70', experiments=False):
        self.statcontents
        if experiments==False:
            experiments=self.allexperiments
        for expt in experiments:
            plt.figure()
            if self.machines[expt]=='Plate Reader 2' :
                if np.isnan(self.statcontents.loc[expt, repairWith])==False:
                    ppf.transformEach(self.data[expt], stat1=repairWith, stat2=dtype, plot=True)
                    ppf.replaceStat(self.data[expt], dtype, replaceWith='transformed'+repairWith)
                    self.containsstat('transformed'+repairWith, printstats=False)
                else:
                    try:
                        if np.isnan(self.statcontents.loc[expt, 'GFP75'])==False:
                            ppf.transformEach(self.data[expt], stat1='GFP75', stat2=dtype, plot=True)
                            ppf.replaceStat(self.data[expt], 'GFP80', replaceWith='transformedGFP75')
                            self.containsstat('transformedGFP75',printstats=False)
                    except:
                        print('no GFP75')
                        if np.isnan(self.statcontents.loc[expt, 'GFP60'])==False:
                            ppf.transformEach(self.data[expt], stat1='GFP60', stat2=dtype, plot=True)
                            ppf.replaceStat(self.data[expt], 'GFP80', replaceWith='transformedGFP60')
                            self.containsstat('transformedGFP60',printstats=False)
    #Data extraction minifunctions
    def extractMaxGR(self, expt, media, strain, plateloc):
        return self.data[expt].d[media][strain]['max growth rate']
    def extractMaxGRTime(self, expt, media, strain, plateloc):
        return self.data[expt].d[media][strain]['time of max growth rate']
    def extractMaxGRVar(self, expt, media, strain, plateloc):
        return self.data[expt].d[media][strain]['max growth rate var']
    def extractLagTime(self, expt, media, strain, plateloc):
        return self.data[expt].d[media][strain]['lag time']
    def extractGRAUC(self, expt, media, strain, plateloc):
        try:
            return ppf.statArea(self.data[expt], media, strain, 'gr')
        except:
            return np.nan
    def extractInitialODplateloc(self, expt, media, strain, plateloc):
        whichisplateloc=np.where([plateloc==j for j in self.data[expt].d[media][strain]['plateloc']])[0]
        return self.data[expt].d[media][strain]['OD'][:,whichisplateloc][0]
    def extractFinalODPlateloc(self, expt, media, strain, plateloc):
        whichisplateloc=np.where([plateloc==j for j in self.data[expt].d[media][strain]['plateloc']])[0]
    def extractMaxODPlateloc(self, expt, media, strain, plateloc):
        whichisplateloc=np.where([plateloc==j for j in self.data[expt].d[media][strain]['plateloc']])[0]
        return np.nanmax(self.data[expt].d[media][strain]['OD'][:,whichisplateloc])
    def extractGRFull(self, expt, media, strain, plateloc=False):
        if isinstance(conditionsDF['plateloc'], list):
            whichisplateloc=np.where([plateloc==j for j in self.data[expt].d[media][strain]['plateloc']])[0]
            return self.data[expt].d[media][strain]['OD'][:,whichisplateloc][0]
        else:
            return np.nanmean(self.data[expt].d[media][strain]['OD'][0,:])
    def extractInitialODFull(self, expt, media, strain, plateloc=False):
        if not isinstance(plateloc, list):
            whichisplateloc=np.where([plateloc==j for j in self.data[expt].d[media][strain]['plateloc']])[0]
            return np.nanmean(self.data[expt].d[media][strain]['OD'][:,whichisplateloc][0])
        else:
            return np.nanmean(self.data[expt].d[media][strain]['OD'][0,:])
    def extractFinalODFull(self, expt, media, strain, plateloc=False):
        return np.nanmean(self.data[expt].d[media][strain]['OD'][-1,:])
    def extractMaxODFull(self, expt, media, strain, plateloc=False):
        return np.nanmax(self.data[expt].d[media][strain]['OD'])
    def extractStat(self, expt, m, s, stat):
        if ppf.hasKey(expt.d[m][s], stat):
            return expt.d[m][s][stat]
        else:
            return np.nan
    def pcaclick(self, reps=None, dtype='OD', clicknumber=0, components=[0,1,2], times=[0,1,2,3,4,5,6,7,8,9,10], rownorm=True, colorby=[], color=colors.nicePastels+colors.strongColors, dotsize=.6, alpha=0.5, label='pcaclick_0'):
        '''
        pcaclick(self, reps=None, , dtype='OD', clicknumber=0,components=[0,1,2], times=[0,1,2,3,4,5,6,7,8,9,10], rownorm=True, colorby=[], color=colors.nicePastels+colors.strongColors, dotsize=.6, alpha=0.5)
        
        Explore all curves of dtype using Principal Components Analysis. Then click and select observations of interest.
        reps= pandas dataframe with replicate information (defaults to self.allreplicates)
        dtype= the type of curves to project. defaults to 'OD'
        clicknumber=the number of clicks intended. assumes no clicks are intended so make sure to change to >0
        components= components of interest in the PCA. recommended to leave as [0,1,2]
        times.- the times at which to sample the curves of dtype. make sure this spans the length of all your curves
        rownorm.- (boolean) whether to normalise each curve by centering its values to their mean and scaling them  by the standard deviation.
        colorby.- *RECOMMENDED* variable in reps to group and colour the observations. e.g. if 'media', each media type will be coloured differently.
        color.- the colours to use when using colorby.
        dotsiz.- the size of the points in the scatterplot
        alpha.- (0..1) the transparency of the points in the scatterplot. more curves call for a smaller number.
        '''
        if not isinstance(reps, pd.DataFrame):
            reps=self.allreplicates ##dataframe of conditions
        fig, ax=plt.subplots(2,2)
        V=self.makedataframe('timerow', dtype=dtype, times=times, conditionsDF=reps).T.fillna(method='ffill').fillna(0)#get the values of the data frame and transpose them to have the time series.
        ##we try to fill the nans by extending the last value in the column whenever possible
        if np.isnan(V).any().any():
            print('removing columns that contain nans...\n')
            nancols=sum(np.isnan(V.values),0)
            print(V.columns.values[nancols>0]) ##these are the columns that have nans
            print('\n')
            V=V[V.columns[nancols==0]].values
            reps=reps.iloc[np.nonzero(nancols==0)[0][0], :]
        ##we fill the trailing NAs with the last value inserted using the last value of the series. then make all remaining nans a 0 in case some curves are completely nan.
        if rownorm==True:
            normrow= lambda row: (row-np.nanmean(row))/np.nanstd(row)
            xn=V.apply(normrow, axis=1).values
        else:
            xn=V.values
        pca = PCA(n_components=5)
        pca.fit(xn) # run the pca
        xpca=pca.transform(xn) # projecting xn onto the pca components
        if colorby:
            if np.size(colorby)>1:
                colorby=colorby[0]
            if not colorby in reps.columns:
                print('variable '+colorby+'not in the available features.')
                #return 0
            else:
                classvector=reps[colorby].values
                varindex=np.nonzero(reps.columns.values==colorby)[0][0]; 
                d=ppf.colorDict(keys=np.unique(classvector), colors=color)
                plt.sca(ax[0][0])
                plt.scatter(xpca[:, components[0]], xpca[:, components[1]], c=[d[j] for j in classvector], lw=dotsize, alpha=alpha)
                plt.sca(ax[1][0])
                plt.scatter(xpca[:, components[0]], xpca[:, components[1]], lw=dotsize, alpha=alpha, color='grey')
                #plt.sca(ax[1][1])
                #plt.scatter(xpca[:, components[1]], xpca[:, components[2]], lw=dotsize)
                ppf.createFigLegend(dic=d)  
        else:
            cl='red'
            plt.scatter(xpca[:, components[0]], xpca[:, components[1]]) 
        time.sleep(0.5)
        plt.sca(ax[1][0])
        plt.xlabel('Component '+str(components[0]))  
        plt.ylabel('Component '+str(components[1]))
        plt.title('PCA of '+dtype)
        plt.sca(ax[1][1])
        plt.xlabel('Component '+str(components[0]))  
        plt.ylabel('Component '+str(components[2]))
        plt.title('PCA of '+dtype)
        pointvector=[]
        incrementlabel= lambda label, j: label if not label in reps.columns else label.split('_')[0]+'_'+str(j)
        for j in range(1, 100):
            newlabel= incrementlabel(label,j)
            if not newlabel in reps.columns:
                break
        if clicknumber >0:
            g=0
            while(g<clicknumber):
                plt.figure(plt.gcf().number)
                plt.suptitle('Click on the scatterplots to explore curves.\n '+str(g)+'/'+str(clicknumber)+' clicks')
                plt.sca(ax[0][0]);
                a=plt.ginput(1);
                subst=[ scipy.spatial.distance.euclidean(np.array([xpca[j, 0], xpca[j, 1]]), a[0]) for j in range(0, np.size(xpca, 0))];
                point=np.argmin(subst)
                pointvector.append(point)
                exptname=reps.iloc[point, 0]
                plateloc=reps.iloc[point, 3]
                print(reps.iloc[point, :]);
                if colorby:
                    cl=d[reps.iloc[point, formatnum(varindex)]];
                #plt.sca(ax[0][0]);
                #plt.scatter(xpca[point, 0], xpca[point, 1], marker='x', c='red', lw=3);
                plt.sca(ax[1][0]);
                plt.scatter(xpca[point, 0], xpca[point, 1], marker='x', c=cl, lw=3);
                plt.sca(ax[0][1]); plt.plot(xn[point, :], color=cl);
                #adding a schematic of in which experiment this is
                ax[1][1].clear()
                axplate=ax[1][1]
                plt.sca(axplate);
                #plt.scatter(xpca[point, 1], xpca[point, 2], marker='x', c=cl, lw=3);
                plt.xlim([0, 12])
                plt.ylim([0,8])
                plt.title(exptname+' '+plateloc)  
                axplate.set_yticks([.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]) 
                xticks=[1,2,3,4,5,6,7,8,9,10,11,12]
                axplate.set_xticks(np.array(xticks)-.5) 
                axplate.set_xticklabels(['%d' % (j) for j in xticks]  )
                letters=['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']
                position=np.where([plateloc[0]== j for j in letters])[0][0]   
                axplate.set_yticklabels(letters)   
                [plt.axvline(x=j) for j in xticks]
                [plt.axhline(y=j) for j in np.linspace(0, 8, 9)]
                #print(letters)
                #print('\nplateloc '+plateloc+' position '+str(position))
                rct=pch.Rectangle(xy=(np.float(plateloc[1:])-1, position), width=1, height=1, color='red') 
                axplate.add_patch(rct)  
                fig.canvas.draw();
                fig.canvas.flush_events();
                g+=1
                factor= np.full((len(reps), 1), False, dtype=bool)
                factor[np.unique(pointvector)]=True
                reps[newlabel]=factor
            plt.suptitle('Click on the scatterplots to explore curves.\n '+str(g)+'/'+str(clicknumber)+' clicks')
        print('clicked curves marked at self.allreplicates['+newlabel+']\n')
        return reps, reps.iloc[np.unique(pointvector),:  ], pca
    def excludereps(self, reps=[], byindex=False):
        ''' excludereps(self, reps=[], byindex=False)
            excludes replicates based on their index in the self.allreplicates dataframe
            reps.- dataframe of replicates to exclude (in the format of self.allreplicates
            byindex(internal use mainly).- exclude replicates assuming the indices of self.allreplicates and reps match. more efficient, but not always guaranteed to work fine due to python indexing.'''
        if byindex:
            self.allreplicates=self.allreplicates[[x not in reps.index.values for x in self.allreplicates.index.values]]                                                                               
        else:
            findrep= lambda a: np.nonzero(np.array(df['experiment'].values==a[0])[0][0] & np.array(df['media'].values==a[1]) & np.array(df['strain'].values==a[2]) & np.array(df['plateloc'].values==a[3]  ) )[0]  
            inds=[findrep(self.allreplicates.iloc[j, :].values) for j in reps.index]  
            self.allreplicates=self.allreplicates.iloc[[x not in inds for x in range(0, len(self.allreplicates))], :]
            self.allreplicates=self.allreplicates.reset_index(drop=True)
    def openfailed(self, number):
        '''try to create a plate reader experiment in the failedfiles dictionary'''
        p=pr.platereader(self.failedfiles[number][0], self.failedfiles[number][1])
        return p 
    def reloadrobust(self, picklefile, only=[], exclude=[]):
        '''load experiments one by one from an xpr.data pickle into an existing xpr data object'''
        dat=pickle.load(open(picklefile, 'rb'))
        #self.allexperiments= np.unique(self.allexperiments)
        for j in list(dat.keys()):
            if not only or j in only or not(j in exclude):
                self.data[j]= dat[j]
        self.listcontents()
        self.getvariables()
        self.getallcontents()
        self.checkallstats()
    def getresiduals(self, dtype='OD'):
        for j in self.allreplicates:
            ex,m,s,pl= df.iloc[j][0],df.iloc[j][1],df.iloc[j][2],df.iloc[j][3]
            self.data[ex].d[m][s][dtype+'residual']=self.data[ex].d[m][s][dtype]-self.data[ex].d[m][s][f+dtype]
            self.data[ex].d[m][s]['sqdif-'+dtype]=np.sum(self.data[ex].d[m][s][dtype+'residual']**2)
            self.data[ex].d[m][s]['normsqdif-'+dtype]=np.sum(self.data[ex].d[m][s][dtype+'residual']**2)/np.len(self.data[ex].d[m][s][dtype+'residual'])
    def overview(self, exptnum, experiment=[], dtype= 'OD', colormap='cool', colorMapRange=False, timeRange=False, addFL=False):
        if experiment:
            p=self.data[experiment]
            st=experiment
        else:
            p=self.getexpt(exptnum)
            st=self.allexperiments[exptnum]
        axarray=ppf.experimentOverview(p,  dtype= dtype, colormap=colormap, colorMapRange=colorMapRange, timeRange=timeRange, addFL=addFL)
        plt.suptitle(st)
        return axarray
    def ignorewells(self, experiments=[]):
        if not experiments:
            experiments=self.allexperiments
        if not not self.ignoredwells:
            for expt in list(self.ignoredwells.keys()):
                self.data[expt].ignorewells(self.ignoredwells[expt])
    def createdirs(self):
        if not os.path.isdir(self.savedir):
            os.mkdir(self.savedir)
        datefolder=self.savedir+'/'+self.date
        if not os.path.isdir(datefolder):
            os.mkdir(datefolder)
        #         if not os.path.isdir(self.savedir+'/'+self.date+'/pickles'):    
        #             self.picklepath=self.savedir+'/'+self.date+'/pickles'
        #             os.mkdir(self.savedir+'/'+self.date+'/pickles')
        #         if not os.path.isdir(self.savedir+'/'+self.date+'/objectdata'):    
        #             self.objectpath=self.savedir+'/'+self.date+'/objectdata'
        #             os.mkdir(self.savedir+'/'+self.date+'/objectdata')
        #         if not os.path.isdir(self.savedir+'/'+self.date+'/staged'):    
        #             self.stagepath=self.savedir+'/'+self.date+'/staged'
        #             os.mkdir(self.savedir+'/'+self.date+'/staged')
    def stagedata(self, label='', experiments=False):
        try:
            self.stagecount+=1
            newstagepath=self.savedir+'/'+self.date+'/accessprstage_'+str(self.stagecount)+'_'+label+'_'+datetime.now().strftime('%Y-%m-%d-%H%Mhrs')
            print('Creating path for the new stage...\n')
            if not os.path.isdir(newstagepath):
                os.mkdir(newstagepath)  
                pickle.dump(self, open(newstagepath+'/accessprobject_'+label+'_'+self.stagecount+'.pkl', 'wb'))  
                pickle.dump(self.data, open(newstagepath+'/datastage_'+str(self.stagecount)+label+datetime.now().strftime('%Y-%m-%d-%H%Mhrs')+'.pkl', 'wb'))
                if experiments:
                    for j in self.allexperiments:
                        pickle.dump(self.data[j], open(newstagepath+'/experiments/'+j+'.pkl', 'wb'))
            self.currentstagepath=newstagepath
        except Exception as err:
            print('data could not be staged: '+err)
            self.stagecount-=1
    def changestage(self, number=[], label=[]):
        if not number and not label:
            return('Please specify either stage number or stage label')
        else:
            files= os.listdir(stagepath)
            if number:
                self.reloadRobust(files[find(['datastage_'+str(number) in j for j in files])[0]])
    def relativetimepoints(self, df=[], dtype='OD', xstat='time',reftime=[], max=[], maxbefore=[], maxafter=[], nsteps=11):
        '''For each condition in the dataset, extract data at each side of reference timepoint reftime, which may be different for every replicate in the experiment. 
        Example:
            Extract OD over time within a window of 5 hours to each side of the fluorescence peak time, exactly 11 timepoints
                relativetimepoints(self, dtype='OD', reftime='FLAbsolutePeakTime', max=5, nsteps=11)
            Extract OD over time 3 hours before and 2 hours before max growth rate
                relativetimepoints(self, dtype='OD', reftime='time of max gr', maxbefore=3, maxafter=2, nsteps=11)
        '''
        if nsteps%2==0:
            raise Warning('Warning: nsteps must be odd in order to retrieve values at query point')
        intervalsize=1
        if not isinstance(df, pd.DataFrame):
            df=self.makedataframe()##retrieving only the hxt4 lines, let's ignore that we have different strains.
        #hxt4df=DFsubset(hxt4df, 'media', )
        mag=intervalsize#magnitude of interval
        nbacksteps=np.floor(nsteps/2)-1
        nfrontsteps=np.ceil(nsteps/2)
        if not max and not maxbefore and not maxafter:
            print('please indicate span length for both sides (max), span before the reftime (maxbefore) or span after the reftime (maxafter)') 
        if max:
            maxbefore=np.floor(max)
            maxafter=np.ceil(max)
        else:
            if maxbefore and not maxafter:
                maxafter=maxbefore
            else:
                if maxafter and not maxbefore:
                    maxbefore=maxafter
                else:
                    if maxbefore and maxafter:
                        nbacksteps= np.floor(nsteps* (maxbefore/(maxbefore+maxafter)))-1  
                        nfrontsteps= np.ceil(nsteps* (maxafter/(maxbefore+maxafter))) 
        steps= np.linspace(-maxbefore*mag, maxafter*mag, nsteps ) #points surrounding the peak time
        #stepsbefore=np.linspace(((-max)*mag), -1*mag, nbacksteps) #points before the peak time. negative training set
        #stepsafter=np.linspace( 0,(max)*mag, nfrontsteps)
        intervals=[np.around(steps+j,2) for j in df[reftime]]
        #nintervalsb=[np.around(stepsbefore+j,2) for j in df[variable]]
        #nintervalsa=[np.around(stepsafter+j,2) for j in df[variable]]
        timevalues=np.unique(intervals) #we find all the unique times obtained from this expansion
        #timevaluesb=np.unique(nintervalsb)
        #timevaluesa=np.unique(nintervalsa)
        #the timevalues array has to be within the interpLimit to work
        #to simplify the process for now we extract  datapoints
        timedf=self.makedataframe(conditionsDF= df, type='timerow', dtype=dtype, times=list(timevalues)) #list(timevalues[0:85])); 
        #timedfb=self.makedataframe(conditionsDF= df, type='timerow', dtype='OD', times=list(timevaluesb)) #list(timevalues[0:85])); 
        #timedfa=self.makedataframe(conditionsDF= df, type='timerow', dtype='OD', times=list(timevaluesa)) #list(timevalues[0:85])); 
        dmOD=np.zeros([np.size(intervals,0), np.size(intervals[0])])  #initialising the design matrix for our training.
        #dmODb=np.zeros([np.size(intervals,0), np.size(nintervalsb[0])])  #initialising the design matrix for our training.
        #dmODa=np.zeros([np.size(intervals,0), np.size(nintervalsa[0])])  #initialising the design matrix for our training.
        for j,k in enumerate(df.index):
            #a,b,c, d=hxt4df.iloc[k, (0, 1,2,3)]
            #try:
            e,m,s,pl= df.iloc[k, 0],df.iloc[k, 1],df.iloc[k, 2],df.iloc[k, 3]
        
            try:
                vec=timedf.transpose().iloc[k, ([np.where(timedf.transpose().columns==j)[0][0] for j in intervals[k]])].values
            except Exception as err:
                print(' '.join([e, m, s, pl])+' error:'+ str(err) )
                vec=[np.nan for j in intervals[k]]
            #vecb=timedfb.transpose().iloc[k, ([np.where(timedfb.transpose().columns==j)[0][0] for j in nintervalsb[k]])].values
            #veca=timedfa.transpose().iloc[k, ([np.where(timedfa.transpose().columns==j)[0][0] for j in nintervalsa[k]])].values
            #except:
                #vec=np.nan* np.size(intervals[0])
                #vecb=np.nan* np.size(intervals[0])
                #veca=np.nan* np.size(intervals[0])
                #vec0=np.nan* np.size(intervals[0])
            dmOD[k, :]=  vec
            #dmODb[k, :]=  vecb
            #dmODa[k, :]=  veca
        beforeindices=[j for j in range(0, int(nbacksteps+1))]  
        afterindices=[j for j in range(int(nbacksteps+1), int(nbacksteps+1+nfrontsteps))]
        return dmOD,steps, beforeindices, afterindices
    def plotraw(self, df=None, dtype='OD', x='time', cmap=colors.strongColors+colors.nicePastels, include='all', exclude=[], row='strain', rowvalues=None, col='media', colvalues=None, hue='experiment', huevalues=None, marker=None, markervalues=None, markers=['.', 'o', 's', '^', '+', 'v', '*', 'p', '<', '>', 'h', 'x', 'D']*10, xlim=[], ylim=[], separatelabel=True, centerat=None, ignore='clicked'):
        '''plotraw(self, df=None, dtype='OD', x='time', cmap=colors.strongColors+colors.nicePastels, include='all', exclude=[], row='strain', rowvalues=None, col='media', colvalues=None, hue='experiment', huevalues=None, marker=None, markervalues=None, markers=['.', 'o', 's', '^', '+', 'v', '*', 'p', '<', '>', 'h', 'x', 'D']*10, xlim=[], ylim=[], separatelabel=True):
        Plot categorised raw curves of x vs dtype. 
        row, col, hue and marker are factors to categorise by. They must be columns in self.allreplicates.
        rowvalues, colvalues etc. are the categories themselves (defaults to all existent categories) 
            Examples:
            self.plotraw(dtype='OD', row='strain', hue='experiment', col='media', colvalues=['m1', 'm2'] )
            #create a plot of m rows x 2 columns where each row corresponds
            # to m different strains and columns are curves grown in 2 media: m1 (first) and m2 (second)
            #the colour of each curve  denotes what experiment where the curve is found. 
        '''
        colnumber=None
        rownumber=None
        if not isinstance( df, pd.DataFrame):
            df=self.allreplicates

        if hue is None:
            hue='strain'
            huevalues=self.allstrains
            huenumber= 1
            huedict= ppf.colorDict(keys=self.allstrains, colors=['blue']*100)
        else:
            if huevalues is None:
                huevalues=np.unique(df[hue])
            huedict= ppf.colorDict(keys=huevalues, colors=cmap)
        if marker is None:
            markervalues=self.allstrains
            marker='strain'
            markerdict= ppf.colorDict(keys=self.allstrains, colors=['']*100)
        else:
            if markervalues is None:
                markervalues=np.unique(df[marker])
            markerdict= ppf.colorDict(keys=markervalues, colors=markers)
        if row is None:
            rowvalues=list(self.allstrains)
            row='strain'
            rownumber=1
            rowdict= ppf.colorDict(keys=rowvalues, colors=[0]*100)
            print(rowdict)
        else:
            if rowvalues is None:
                rowvalues=np.unique(df[row])
            rownumber= np.size(rowvalues)
            rowdict= ppf.colorDict(keys=rowvalues, colors=[j for j in range(0, rownumber)])
            print(rowdict)
        if col is None:
            col='strain'
            colvalues=list(self.allstrains)
            colnumber=1
            coldict= ppf.colorDict(keys=colvalues, colors=[0]*100)
            print(coldict)
        else:
            if colvalues is None:
                colvalues=np.unique(df[col])
            colnumber= np.size(colvalues)
            coldict= ppf.colorDict(keys=colvalues, colors=[j for j in range(0, colnumber)])
            print(coldict)
        #subsetting
        #which=[j.all() for j in pd.concat([df[rowvar]==rowvalue, df[colvar]==colvalue, df[huevar]==huevalue,df[markervar]== markervalue], axis=1).values] 
        #test if the dataframe row has the given features.
        hasfeatures= lambda e, m, s, pl: [all([e in j, m in j, s in j, pl in j]) for j in df.values]  
        #check if the individual line has an element of the subsets
        isinsubset= lambda j: all([df[row][j] in rowvalues, df[col][j] in colvalues, df[hue][j] in huevalues, df[marker][j] in markervalues])
        isincluded= lambda j: any( [x in include for x in df[[row, col, hue, marker]].iloc[j].values])
        isexcluded= lambda j: any( [x in exclude for x in df[[row, col, hue, marker]].iloc[j].values])
        fig, axarray=plt.subplots(rownumber, colnumber, figsize= [10,10], squeeze=False) ## stupid squeeze!!
        print('DIMENSIONS \n'+ str(rownumber)+'  by  '+str(colnumber))
        getprops= lambda num: (df.iloc[num, 0], df.iloc[num, 1], df.iloc[num, 2], df.iloc[num, 3])
        #getcurve = lambda e, m, s, pl, dtype:  getindcurve(e,m,s,pl) if np.size(self.data[e].d[m][s][dtype],1)>1 and len(self.data[e].d[m][s]['plateloc']) >1 else getenscurve(e,m,s)
        #example line of how to  use getcurve and getprops.
        #plt.figure(); plt.plot(np.array([getcurve(*getprops(j)) for j in [1,2,3,4,5,6]]).T)    
        for j in range(0, len(df)):
            a={}
            a['experiment'], a['media'], a['strain'], a['plateloc']= getprops(j)
            a[hue]=df.iloc[j][hue]
            a[col]= df.iloc[j][col]
            a[row]=df.iloc[j][row]
            a[marker]=df.iloc[j][marker] 
            #if a['plateloc'] in self.data[a['experiment']].ignoredwells or a['plateloc'] in self.ignoredwells[a['experiment']] or df[ignore].values[j]:
            #    continue
            if (isinsubset(j) or isincluded(j)) and not(isexcluded(j)):
                try:
                    curve=self.getindcurve(a['experiment'], a['media'], a['strain'], a['plateloc'], dtype) 
                except:
                    curve=self.getenscurve(a['experiment'], a['media'], a['strain'], dtype) 
                #try:
                ax=axarray[rowdict[df[row][j]]][coldict[df[col][j]]]
                #except:
                #    ax=axarray
                xx=self.getenscurve(a['experiment'], a['media'], a['strain'], x)
                if centerat:
                    xx=xx-df.iloc[j][centerat]
                ax.plot(xx, curve, color=huedict[a[hue]], marker=markerdict[a[marker]], markevery=20, label=a[marker])
                if xlim:
                    ax.set_xlim(xlim)
                if ylim:
                    ax.set_ylim(ylim)
                #axarray[rowdict[df[row][j]]][coldict[df[col][j]]].title.set_text(df[row][j]+' '+df[col][j])
                if coldict[df.iloc[j][col]]==0:
                    ax.set_ylabel(a[row])
                if rowdict[df.iloc[j][row]]==0:
                    ax.set_title(a[col])
        fig.text(0.5, 0.04, x, ha='center')
        fig.text(0.04, 0.5,  dtype, va='center', rotation='vertical')
        if separatelabel:
            plt.figure()
            ppf.createFigLegend(dic=huedict, markers=markers)
        try:
            return axarray
        except:
            return ax
def drawplatelayout():
    axplate=plt.axes()
    plt.xlim([0, 12])
    plt.ylim([0,8])
    #plt.title(exptname+' '+plateloc)  
    axplate.set_yticks([.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]) 
    xticks=[1,2,3,4,5,6,7,8,9,10,11,12]
    axplate.set_xticks(np.array(xticks)-.5) 
    axplate.set_xticklabels(['%d' % (j) for j in xticks]  )
    letters=['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']  
    axplate.set_yticklabels(letters)   
    [plt.axvline(x=j) for j in xticks]
    [plt.axhline(y=j) for j in np.linspace(0, 8, 9)]
    return axplate
def selectwells(n=1, axplate=False, fig=False, color='red'):
    if not axplate:
        axplate=plt.gca()
    if not fig:
        fig= plt.gcf()
    plateloc=[];
    coords=plt.ginput(n, mouse_stop=3, mouse_pop=2, mouse_add=1)
    for j in range(0, len(coords)):
        x=int(np.floor(coords[j][0]))
        y=int(np.floor(coords[j][1]))
        rct=pch.Rectangle(xy=(x, y), width=1, height=1, color=color, alpha=.3)
        xticks=[1,2,3,4,5,6,7,8,9,10,11,12]
        letters=['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']  
        axplate.add_patch(rct)
        plateloc.append(letters[y]+str(xticks[x]))  
    fig.canvas.draw();
    fig.canvas.flush_events();
    return plateloc
def clickwells(self, expts=[], dtype= 'OD', colormap='cool', colorMapRange=False, size=4, timeRange=False, addFL=False, label='clicked', color='red'):
    if not expts:
        expts=self.allexperiments
    self.ignoredwells={};
    self.allreplicates[label]=False
    for j in expts:
        clickfig=plt.figure()
        ax=showplate(self.data[j], dtype= dtype, colormap=colormap, colorMapRange=colorMapRange, size=size, timeRange=timeRange, addFL=addFL) 
        plt.suptitle(j+'\nleftclick: ignore well \n rclick/enter: next expt. delete: unclick last well')
        wells=selectwells(n=96, axplate=ax, color=color)
        self.ignoredwells[j]=wells;
        self.allreplicates[label]=np.logical_or(self.allreplicates[label], np.logical_and( self.allreplicates['experiment']==j, [x in wells for x in self.allreplicates['plateloc']]))      
def showplate(p, dtype= 'OD', colormap='cool', colorMapRange=False, size=4, timeRange=False, addFL=False):
    ax=plt.axes()
    defaultRange={'OD': [0,1.5], 'GFP':[0,8000], 'AutoFL': [0,1000], 'mCherry': [0, 5000],'GFP80':[0,8000], 'GFP60':[0,8000], 'GFP100':[0,8000]}
    if addFL!= False:
        defaultRange.update(addFL)
    if colorMapRange==False:
        colorMapRange= defaultRange[dtype]
    xstat=p.t
    drawplatelayout()
    ylim= [0,  1]
    xlim=[0, 1];
    xdata= (xstat-xlim[0])/ np.max(xstat-xlim[0])
    cl=ppf.conditionList(p)[['media', 'strain']]
    problematicWells=[]
    for x in range(0, np.size(cl,0)):
        media=cl.values[x,0]
        strain=cl.values[x,1]
        conditionCoordinates=[ppf.plateCoordinates(a) for a in p.d[media][strain]['plateloc']]
        xpositions=np.array(conditionCoordinates)[:,1] #A1, 0 0, is in fact located at 11, 0 visually
        ypositions=7-np.array(conditionCoordinates)[:, 0]          
        #print('plateloc: ', p.d[media][strain]['plateloc'])
        #print('conditionCoordinates: ', conditionCoordinates)
        for cc in range(0,np.size(conditionCoordinates,0)):
            ccdata= p.d[media][strain][dtype][:,cc]
            ccdata= (ccdata-np.min(ccdata))/np.max(ccdata-np.min(ccdata)) #scaling the data of this well btn 0 and 1
            plt.ylim([ylim[0], ylim[1]*8])
            plt.xlim([xlim[0], xlim[1]*12])
            if p.d[media][strain]['plateloc'][cc] in p.ignoredwells:
                continue
            else:
                try:
                    plt.scatter(xdata+xpositions[cc], ccdata+ypositions[cc], c=p.d[media][strain][dtype][:,cc], s= size, cmap=colormap, vmin=colorMapRange[0], vmax=colorMapRange[1], edgecolor='None')
                except IndexError:
                    #print('time length: ', np.size(p.t))
                    #print('vector length: ', np.size(p.d[media][strain][dtype],0))
                    print('Index error:  well ', p.d[media][strain]['plateloc'][cc], ', media ', media, ', strain ', strain)
    plt.colorbar()
    return ax

def line_select_callback(eclick, erelease):
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
    print(" The buttons you used were: %s %s" % (eclick.button, erelease.button))


def toggle_selector(event):
    print(' Key pressed.')
    if event.key in ['Q', 'q'] and toggle_selector.RS.active:
        print(' RectangleSelector deactivated.')
        toggle_selector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print(' RectangleSelector activated.')
        toggle_selector.RS.set_active(True)


def selectarea(current_ax):
    # drawtype is 'box' or 'line' or 'none'
    toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1, 3],  # don't use middle button
                                           minspanx=5, minspany=5,
                                           spancoords='pixels',
                                           interactive=True)
    plt.connect('key_press_event', toggle_selector)
    plt.show()

def selectpoint(ax=False, fig=False):
    if not ax:
        ax=plt.gca()
    if not fig:
        fig=plt.gcf()
    x, y= plt.ginput()
    return x,y

def clickpoints(self, conditionsDF=False, dtype='OD', x='time', xlab='xpoint', ylab='ypoint'):
    if not conditionsDF:
        conditionsDF=self.allreplicates
    xvec=[]
    yvec=[]
    print(conditionsDF)
    for j in range(0, len(conditionsDF)):
        expt, media, strain,plateloc= conditionsDF.values[j, [0,1,2,3]]
        plt.figure()
        plt.suptitle(str(np.round(j/len(conditionsDF)))+expt)
        plt.title(media+strain+plateloc)
        plt.plot(self.data[expt].d[media][strain][x], self.data[expt].d[media][strain][dtype][:, findloc(plateloc, self.data[expt].d[media][strain]['plateloc'])])
        xy=plt.ginput(timeout=60)
        xvec.append(xy[0][0])
        yvec.append(xy[0][1])
        plt.close(plt.gcf())
    conditionsDF[xlab]=xvec
    conditionsDF[ylab]=yvec
    def storeproperties(self):
        propertylist=self.__dict__.keys() 
        finaldict={}
        for j in propertylist:
            finaldict[j]= self.__dict__[j]
        return finaldict
def mergedicts(dictlist):
    final={}
    for dic in dictlist:
        for key in dic.keys():
            final[key]=dic[key]
    return final

# mag=1#magnitude of interval
# max=5
# nsteps=21 #number of values in the interval. must be odd to pass by zero
# steps= np.linspace(-max*mag, max*mag, nsteps ) #points surrounding the peak time
# stepsbefore=np.linspace((2*(-max)*mag)-1, -1*mag, nsteps ) #points before the peak time. negative training set
# stepsafter=np.linspace( 1*mag,2*(max)*mag+1, nsteps )
# stepstart=np.linspace( 0*mag,2*(max)*mag, nsteps ) #points at the beginning of some time series.
# intervals=[np.around(steps+j,2) for j in hxt4df['absolutePeakTime']]
# nintervalsb=[np.around(stepsbefore+j,2) for j in hxt4df['absolutePeakTime']]
# nintervalsa=[np.around(stepsafter+j,2) for j in hxt4df['absolutePeakTime']]
# nintervals0=[np.around(stepstart,2) for j in hxt4df['absolutePeakTime']]
# timevalues=np.unique(intervals) #we find all the unique times obtained from this expansion
# timevaluesb=np.unique(nintervalsb)
# timevaluesa=np.unique(nintervalsa)
# timevalues0=np.unique(nintervals0)
# #the timevalues array has to be within the interpLimit to work
# timedf=xpr.makedataframe(conditionsDF= hxt4df, type='timerow', dtype='OD', times=list(timevalues)) #list(timevalues[0:85])); 
# timedfb=xpr.makedataframe(conditionsDF= hxt4df, type='timerow', dtype='OD', times=list(timevaluesb)) #list(timevalues[0:85])); 
# timedfa=xpr.makedataframe(conditionsDF= hxt4df, type='timerow', dtype='OD', times=list(timevaluesa)) #list(timevalues[0:85])); 
# timedf0=xpr.makedataframe(conditionsDF= hxt4df, type='timerow', dtype='OD', times=list(timevalues0)) #list(timevalues[0:85])); 


# expdic={}
# attributelist=[
# p.version
# p.dsheetnumber
# p.asheetnumber
# p.wdir
# p.name
# p.aname
# p.standardgain
# p.ignoredwells
# p.ignored
# p.negativevalues
# p.serialnumber
# p.machine
# p.expdate
# p.datatypes
# p.gains
# p.mediaGP
# p.platelabels
# p.allstrains
# p.allconditions
# p.alldata
# p.importtime
# p.gc
# 
# 
# #migrating a p structure into a dictionary
# expdic={}
# for expt in self.allexperiments:
#     p=self.data[expt]
#     cl=ppf.conditionList(p)
#     expdic[expt]={}
#     expdic[expt]['version']=p.version
#     expdic[expt]['dsheetnumber']=p.dsheetnumber
#     expdic[expt]['asheetnumber']=p.asheetnumber
#     expdic[expt]['wdir']=p.wdir
#     expdic[expt]['name']=p.name
#     expdic[expt]['aname']=p.aname
#     expdic[expt]['standardgain']=p.standardgain
#     expdic[expt]['ignoredwells']=p.ignoredwells
#     expdic[expt]['ignored']=p.ignored
#     expdic[expt]['negativevalues']=p.negativevalues
#     expdic[expt]['serialnumber']=p.serialnumber
#     expdic[expt]['machine']=p.machine
#     expdic[expt]['expdate']=p.expdate
#     expdic[expt]['datatypes']=p.datatypes
#     expdic[expt]['gains']=p.gains
#     expdic[expt]['platelabels']=p.platelabels
#     expdic[expt]['allstrains']=p.allstrains
#     expdic[expt]['allconditions']=p.allconditions
#     expdic[expt]['alldata']=p.alldata
#     expdic[expt]['importtime']=p.importtime
#     expdic[expt]['gc']=p.gc
#     for j in cl.values:
#         m= j[1]
#         s=j[2]
#         expdic[expt]['data'][m][s]


### wizard


# if not self.data:
#     print('There is no data in your experiment. Try:\na) Check whether the path you inserted is correct\nb) try loading data from a previously stored accesspr data file with self.reloadrobust(filename)
# if not FL:
#     print('No fluorescence has been assigned. try running self.assignFL(mainFL=YOUR_FLUORESCENCE_CHANNEL))


def odtime(self, od='0.25', default='first'):
    from scipy.interpolate import interp1d
    arr=np.zeros(len(self.allreplicates))*np.nan
    for j in range(0, len(self.allreplicates)):
        e,m,s,pl= self.getprops(j)
        x=self.getenscurve(e,m,s, 'time')
        y=self.getindcurve(e,m,s,pl, 'OD')
        try:
            func=interp1d(y,x, bounds_error=False, fill_value=np.nan)
            if default=='first':
                arr[j]=func(od)
        except Exception as err:
            print(str(j)+':'+e+' '+m+' '+s+' '+pl+':'+str(err))
            arr[j]= np.nan
    self.allreplicates['time to OD of '+str(od)]=arr

def makeinterpolant(self, emsp=None,  x='time', y='OD', xcenter=0, ycenter=0, ymultiply=1, ydivide=1):
    from scipy.interpolate import interp1d
    if not not emsp:
        e,m,s,pl=emsp[0], emsp[1], emsp[2], emsp[3]
    else:
        e,m,s,pl= self.getprops(num)
    x=self.getenscurve(e,m,s, x)-xcenter
    y=(self.getindcurve(e,m,s,pl, y)-ycenter)*ymultiply/ydivide
    func=interp1d(x, y, bounds_error=False, fill_value=np.nan)
    return func

#ensemble interpolant
def makeinterpolantens(self, ems=None,  x='time', y='OD', xcenter=0, ycenter=0, ymultiply=1, ydivide=1):
    from scipy.interpolate import interp1d
    if not not ems:
        e,m,s=ems[0], ems[1], ems[2]
    else:
        e,m,s,pl= self.getprops(num)
    x=self.getenscurve(e,m,s, x)-xcenter
    y=(self.getenscurve(e,m,s, y).mean(1)-ycenter)*ymultiply/ydivide
    func=interp1d(x, y, bounds_error=False, fill_value=np.nan)
    return func

def makeinterpolantens2(self, ems=None,  x='time', y='OD', xcenter=None, ycenter=0, ymultiply=1, ydivide=1):
    from scipy.interpolate import interp1d
    if not not ems:
        e,m,s=ems[0], ems[1], ems[2]
    else:
        e,m,s,pl= self.getprops(num)
    if not not xcenter: # if it is a number
        x=self.getenscurve(e,m,s, x)-xcenter
    if isinstance(xcenter, str):
        x=self.getenscurve(e,m,s, x)-self.xcenter
    y=(self.getenscurve(e,m,s, y).mean(1)-ycenter)*ymultiply/ydivide
    func=interp1d(x, y, bounds_error=False, fill_value=np.nan)
    return func



def processts(self, x, y, nums=None, times=[0], xcenter=None, ycenter=None, ymultiply=None, ydivide=None):
    if not nums:
        nums=range(0, len(self.allreplicates))
    inn=[]
    for j in nums:
        if not not xcenter:
            xc=self.allreplicates[xcenter][j]
        else:
            xc=0
        if not not ycenter:
            yc=self.allreplicates[ycenter][j]
        else:
            yc=0
        if not not ymultiply:
            ym=self.allreplicates[ymultiply][j]
        else:
            ym=1
        if not not ydivide:
            yd=self.allreplicates[ydivide][j]
        else:
            yd=1
        inn.append(makeinterpolant(self, j, x=x, y=y, xcenter= xc, ycenter=yc, ymultiply=ym)(times))
    return inn
    


#####writing all curves to a table!!!

# def save(filename, datatypes
# 
# writecurve_time_OD= lambda fil, j: fil.write(','.join(list(self.getprops(j))) +','+ printarray(self.getenscurve(*self.getprops(j)[0:3], 'time'))+','+printarray( self.getindcurve(*self.getprops(j), 'OD'))+'\n')  
# 
# [writecurve_time_OD(fil, j) for j in rng] 


###contents reformatting. goal: attach experiment media strain plateloc dtype to data
#platecodes=pd.DataFrame([[j+k for j in list(testcnt.index)] for k in [str(x) for x in testcnt.columns]]) 





#def contentcolumns(contentsfile):

#[[testcnt.loc[j,int(k)].split(' in ')[1],testcnt.loc[j,int(k)].split(' in ')[0],j+k ]  for j in list(testcnt.index) for k in [str(x) for x in testcnt.columns]]  


