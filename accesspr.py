import platereader as pr
import numpy as np
import _pickle as pickle 
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import os.path as path
import os 
import pandas as pd
import scipy.interpolate as scint
from functools import partial
import prPlottingFunctions as ppf
from random import choice, randint, sample
import statsmodels.stats.api as sms
from matplotlib.colors import hex2color
import colors
#from decimal import *
#getcontext().prec = 3

### RAW DATA FUNCTIONS. functions that receive a rawdata dictionary where rawdata[stat] is a matrix of m wells x n timepoints,
###coming directly from the plate reader excel sheet.

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
    df2=pd.DataFrame(columns= df.columns)
    for x in values:
        df2=pd.concat([df2, filterDataFrame(df, [variable], [x])])
    df2=df2.reset_index(drop=True)
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
    
def colorScatter(p, media, strain, xstat=False, ystat='ODmn', colorBy='d/dtgr', symmetric=True, cmap='bwr',nbins=100, extendBy=2, alpha=1, markersize=12, marker='o', addLegend=True, vmin=0, vmax=0, xlabel=1,ylabel=1):
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
            except:
                print('problem finding half fluorescence time for'+strain+' in '+media+'\n' )
                halfFLToPeakTime=np.nan
                slope=np.nan
                timeToHalfFL=np.nan 
        responseTimeAligned=timeToHalfFL-centeredTime
        halfFLToPeakTime= timeToHalfFL-absolutePeakTime ### time from half fl to full fl. the sign of the steepness determines whether this was found before (negative) or after the peak.
        steepness=1/halfFLToPeakTime[0]
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
    out= {'rawFL': rawFLVector, 'normalizedFL': normalizedFLVector, 'FLPeak': FLPeak, 'normalizedFLPeak':normalizedFLPeak, 'alignedPeakTime':alignedPeakTime, 'absolutePeakTime':absolutePeakTime[0], 'slope': slope[0], 'responseTime': timeToHalfFL, 'steepness': steepness, 'half2PeakTime': halfFLToPeakTime[0], 'responseTimeAligned': responseTimeAligned, 'halfFL': halfFL }
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
def getRawData(expt):
    if isinstance(expt, str): #process the string if it is one, otherwise assemble it.
        xl=pd.ExcelFile(expt)
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
    matrixrows= dataends[1]-datastarts[0]+1 ###the end of the first matrix is 4 indices above the next time vector
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
    return rawdata

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
        final[dtype+'mn']=totalmn
        final[dtype+'bootstrapsd']=totalsd
    else:
        final[dtype]=interpolated[dtype]
        final[dtype+'mn']=mn
        final[dtype+'sd']=sd
        plt.plot(obj[centeringVariable], mn, color=col)
        plt.fill_between(obj[centeringVariable], mn-sd, mn+sd, color=col, alpha=alpha )
        plt.xlabel(centeringVariable)
        plt.ylabel(dtype)
    return final

def preprocessExpt(expt, fillnans=True, standardize=False, normOD=0.8, extension='preprocessed.xlsx', plot=False, main='GFP80', supporting='GFP60'):
    '''
    preprocessExpt(expt, fillnans=True, standardize=True, mapFL=True, normOD=0.8, extension='preprocessed.xlsx' )
    preprocess data in excel sheet and output a filenamepreprocesed.xls
    '''
    rawdata=getRawData(expt)
    #indices of the main channel which are not nans
    if fillnans==True:
        notnans=np.isnan(rawdata[main])==False
        #fitting a regression between the the supporting (lower) measurement and no nan values of the main measurement
        reg=scipy.stats.linregress(rawdata[supporting].values[notnans].flatten(),rawdata[main].values[notnans].flatten())
        ln=np.linspace(0, np.nanmax(rawdata[supporting]), 300)
        ln2=np.linspace(0, np.nanmax(rawdata[supporting]), 300)*reg.slope+reg.intercept
        if plot==True:
            plt.scatter(rawdata[supporting],rawdata[main])
            plt.xlabel(main)
            plt.ylabel(supporting)
            plt.plot(ln, ln2, color='red')
            plt.legend([ 'y='+str(reg.slope)+'x'+stringSign(reg.intercept)+str(reg.intercept)])
        #fixing the nans in stat2 by applying a regression from values in stat 1
        substitution=rawdata[main][np.isnan(rawdata[supporting])]*reg.slope+reg.intercept     
        rawdata[main][np.isnan(substitution)==False]=substitution[np.isnan(substitution)==False]
    #performing the transformation based on WT fluorescence
    tdata={}
    tdata['OD']=rawdata['OD']
    #reg1, tdata[stat1]=internalStd(rawdata, stat1, wtindices, fr=0.4, to=normOD)
    #reg2, tdata[stat2]=internalStd(rawdata, stat2, wtindices, fr=0.4, to=normOD)
    #### trying to create a new excel file with the new preprocessed data
    df2=df #create a copy of the original data
    for i in itime:
        d= i+2
        #rawdata[df.ix[i-2][0]]=pd.DataFrame(df.ix[d:d+matrixrows-1, 1:1+np.size(timepoints)], index=df.ix[d:d+matrixrows-1,0], columns=timepoints) #getting rid of incomplete timepoints
        dt=df.ix[i-2][0] ##the datatype
        df2.ix[d:d+matrixrows-1, 1:1+np.size(timepoints)]=tdata[dt].values #measurements
    ##now that we have got the raw data we can apply functions to the matrices directly
    df2.to_excel(expt.name+extension, header=True, index=False)
    return tdata


class accesspr:	
    '''
    accesspr version 4.81(the number matches with that of the compatible platereader software version)
    
    accesspr is a class that allows to integrate and organize the information
    from many plate reader experiments to produce publication-grade plots and tables. 
    To fully understand the functionality of accesspr, it is ideal
    to learn about pandas (especially dataframes) and the seaborn graphics package. 
    
    ***Newest FEATURES***
    self.allReplicates contains the condition in every well in every experiment
    self.plotReplicateMeanNew() and self.interptimesnew now work with replicate dataframes as an input to include specific replicates
    self.getRawData() creates a raw data dictionary to handle the full pr output altogether.
    ATTRIBUTES:
    source: 
        the path of the folder that contains the pickle files
    data:
        a dictionary of platereader objects, corresponding to each experiment in each pickle file. 
        The keys are each experiment's pickle file name. e.g. xpr.data[picklefilename].d[media][strain]
    allContents:
        a pandas dataframe containing all unique combinations of strains and media found in all experiments. 
    statcontents:
        a pandas dataframe describing whether each expe
    extractionFields: 
        name of the columns of all the features to be extracted when calling self.extractRepInfo() and self.extractAllInfo()
    allstrains: 
        contains a nonredundant list of all strains found in the experiments.
    allmedia:
        contains a nonredundant list of all media found in the experiments.
    allExperiments:
        contains a nonredundant list of all experiments in the accesspr object.
    aligned:
        indicates whether the experiments have been aligned in time by the max growth rate and the max Fluorescence (see method alignAll) 
    FL: 
        a dictionary that contains the default fluorescence channels to be used for every experiment. can be changed manually or by calling
        xpr.assignFL()
    machines:
        a dictionary (with experiment keys) that contains the machine in which this experiment was run. currently identified machines are 'Plate Reader 1' and 'Plate Reader 2'
    serialnumbers:
        a dictionary (with experiment keys) that contains the machine serial number. even if the machine has no identified human-readable name (see machines), the serialnumber will be stored here.
    METHODS: (*:runs during initialization)
    To create an accesspr object, 
    
        xpr=accesspr(PICKLEFOLDERPATH, ignoreFiles=False):
            initializes an accesspr instance by using the files in path PICKLEFOLDERPATH. 
            ignoreFiles is a list of the '*.pkl' file names in this folder that wish to be ignored.
            
    To interrogate about experiment contents and keep track of experiment processing,
    
        * getAllContents(self):
             produces self.allConditions attribute (see attributes)
             
        * containsstat(self, STAT):
            produces the statcontents attribute (see attributes)
            
        * listcontents(self, verbose=True):
            fills the allstrains and allmedia attributes. If verbose==True, produces a human-readable display of the contents of each experiment.
            
        *assignFL(mainFL='noFL', mainFLperod='noFLperod')
            if arguments are specified, sets the given channels as the  default fluorescence channels for each experiment. 
            if arguments are unspecified, looks for consensus fluorescence amongst experiments. Otherwise sets GFP as default.
        
        containssetup(self, MEDIA,STRAIN):
            retrieves a list of experiments that contain media MEDIA and strain STRAIN. 
            Internally sets the temporary property containslist to this experiment list. 
        
        plotRawReplicates(self, MEDIA, STRAIN, dtype='OD', exptColors=False):
            plots the raw repicates of dtype VS time of a specific strain or media, for all the experiments that contain that specific condition.
            if exptColors is not specified, then random colors are assigned to each experiment and returned as a function output.
        plotReplicateMean(self, MEDIA, STRAIN, dtype='')
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
        
        interpTimes(self, media, strain, dtype='FLperod', centeringVariable='gr', upperLim=16, exceptionShift=0.01, ignoreExps=False)
            Creates a sub dataset such that time-varying property dtype is sampled at the same times across replicates.
            This is important because the exact sampling times might vary across replicates.
            returns a dictionary with two keys:
                ['time']- the general sampling times for all replicates
                [dtype]- the  dtype statistic sampled at times ['time']for all replicates found. 
            the function 
                a) aligns the replicates to their max growth rate time
                b)trims experiments so that all have data in a same shared time interval
                c)performs interpolation using an arbitrary experiment's timepoints. 

        extractRepInfo(self, media, strain)
            One of the most useful functions in accesspr. 
            Finds all experiments with replicates for strain in media, and then 
            returns a pandas dataframe with most single point statistics and string-formatted- time series of each replicate's info. 
            the statistics deployed and the order (currently fixed) can be found in the self.extractionFields attribute. 
            To understand the power of having the data in a pandas dataframe, I recommend the following websites:
            10 minutes to pandas (check out the data frame section)
                http://pandas.pydata.org/pandas-docs/stable/10min.html
            plotting directly with pandas
                http://pandas.pydata.org/pandas-docs/stable/visualization.html
            
            for extraction of all experimental conditions present in your accesspr structure, see method extractAllInfo() below.
            output dataframe Field descriptions: (p- one point statistic. s-array stored in a string)
                'experiment'(p)-        name of the pickle file containing this replicate
                'media'(p)-             media identifier string for this replicate
                'strain'(p)-            strain identifier string for this replicate
                'InitialOD'(p)-         the raw OD at the first timepoint of the experiment
                'FinalOD'-(p)           the raw  OD reached at the last timepoint of the experiment
                'InitialRawFL'(p)-      the mean raw GFP  measure across all within experiment replicates
                'InitialFLperOD'(p)-    the corrected FLperod at the beginning of the experiment. 
                'FinalFLperOD'(p), -    the corrected FLperod at the last timepoint of the experiment
                'FLPeak'(p)       -     the max value of the FLperod reached
                'FLAbsPeakTime'(p)-     the time (in hrs after the beginning of the experiment) at which the max value of the FLperod is reached
                'FLAlignedPeakTime'(p)- the time (in hrs relative to the max growth rate time) at which the max Value of the FLperod is reached.
                'lagTime'(p)-           the time (in hrs after the beginning of the experiment) to reach the max growth rate.
                'realTime'(s)-          string containing the data sampling times (in hours after the beginning of the experiment)
                'alignedTime'(s)-       string containing the data sampling times (in hours relative to the max growth rate time)
                'FLperODTS'(s)-         string containing the fluorescence per od data samples for the whole experiment.
                'maxGR' (p)-            max growth rate reached by the culture
                'maxGRvar' (p)-         max growt rate variance for all within-experiment replicates
                'maxGRTime'(p)-         time (in hours after the beginning of the experiment) at which the max growth rate is reached (possibly redundant with lag Time but may be small differences in the way it is calculated)
                'FLperodAUC'(p)-        total area under the curve (AUC)for the FLperod. this measure does not include any background subtraction
                'grAUC'(p)-             total area under the curve (AUC) of the growth rate. may come in handy
                'halfFLAreaTime'(p)-    time at which half the total FLperod AUC in the experiment has been accumulated.
                'halfGRAreaTime'(p)-    time at which half the total growth rate AUC in the experiment has been accumulated.
    
                Future versions will allow to customize the pandas dataframe exported but extracting all features will always be recommended.
                
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
        
        extractAllInfo(self, excludeNull=False)
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
                
        plotReplicateMean(self, media, strain, dtype='FLperod', col='Black', alpha=0.2, exceptionShift=0.01):
            for a given strain and condition, runs interpTimes and plots the mean ± standard deviation (shaded area) over time for the combined replicates. 
            It operates on the current figure, so multiple plots can be overlayed.
            
        plotRawReplicates(self, MEDIA, STRAIN, dtype='OD', exptColors=False):
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

    def __init__(self, source, encoding='latin1', ignoreFiles=False, FL='noFL', FLperod='noFLperod', onlyFiles=False, analyseFL=True, preprocess=True):
        '''
        initializes an accesspr instance from a path specifying either a directory or a pickle file. if this fails, try changing the encoding to ascii, utf8, utf16 or latin1 (so far tried).
        '''
        self.prtype='Tecan'
        #source is a directory with several pickle files. if no pickles are found, then it tries to find experiments.
        self.source = source
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
        self.version='4.81'
        self.releaseNotes='xpr.FL contains default fluorescence per experiment'
        self.allConditions=pd.DataFrame()
        self.loadInitial()
        #self.interpLimit= np.min(np.array(list(self.experimentDuration.values())))
        #this stores the absolute duration (in hrs) of the shortest experiment
        self.defaultStats=['gr','GFP','c-GFPperod', 'GFP100','c-GFP100perod','GFP90','c-GFP90perod','GFP80','c-GFP80perod','GFP70', 'c-GFP70perod','GFP60', 'c-GFP60perod','GFP50', 'c-GFP50perod']
        self.extractionFields=['experiment', 'machine','media','strain','InitialOD','FinalOD','lagTime','realTime','alignedTime','maxGR','maxGRvar','maxGRTime','grAUC','InitialRawFL','InitialFLperOD','FinalFLperOD','FLPeak','FLAbsPeakTime','FLAlignedPeakTime','FLperodAUC','slope','responseTime','steepness','responseTimeAligned']
        self.extractionFieldsOD=['experiment', 'machine','media','strain','InitialOD','FinalOD','lagTime','realTime','alignedTime','maxGR','maxGRvar','maxGRTime','grAUC']
        self.extractionFieldsGrowth=['d/dtgrvar','time of max growth rate','time of max growth rate var','local max growth rate','OD gp','doubling time var', 'doubling time', 'max OD', 'lag time', 'max growth rate', 'max growth rate var']
        self.extractionFieldsFL=['FLPeak','absolutePeakTime','alignedPeakTime','responseTime','responseTimeAligned','halfFL','normalizedFLPeak','slope','steepness']
        self.mediaValue={'Glu 0.2%': 0.2,'Glu 0.4%': 0.4,'Glu 0.6%': 0.6,'Glu 0.8%': 0.8,'Glu 1%': 1,'Glu 1.5%': 1.5,'Glu 2%': 2,'2% Glu': 2,'0.2% Glu': 0.2, 'SucGlu 1%': 0.5, 'SucGlu 1.8% 0.2%': 0.2/2, 'SucGlu 0.2% 1.8%': 1.8/2}
        self.strainAlias={'YST_498': 'Hxt1', 'YST_499': 'Hxt1', 'Hxt4': 'Hxt4', 'Hxt2': 'Hxt2','Hxt3': 'Hxt3','Hxt5': 'Hxt5','Hxt6': 'Hxt6','Hxt7n': 'Hxt7' }
        self.listcontents(verbose=False)
        #self.getAllContents()
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
        try:
            self.getAllContents()
        except:
            print('Impossible to get  allContents list for all experiments.')
        self.refstrains=dict()
        self.checkallstats()
        if analyseFL==True:
            self.analyseFL=True
            self.FL={}
            self.consensusFLs=[]
            self.consensusFL=FL
            self.consensusFLperod=FLperod
            self.supportingFL='GFP60'
            for key in self.data.keys():
                self.FL[key]={}
            try:    
                self.assignFL(mainFL=FL, mainFLperod=FLperod)
            except:
                print('problems finding consensus fluorescence')
                self.assignFL(mainFL='GFP')
        else:
            self.analyseFL=False
        try:
            self.alignAll(rerun=True)
        except:
            print('Problem aligning experiments, possibly owing to missing growth rate. Consider running .getstats()')
        self.findMachines()
        self.replicateLocations()
    def checkallstats(self):
        self.statcontents=pd.DataFrame( index= list(self.data.keys()))
        for k in self.defaultStats:
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
                    mmin=np.min(self.data[expt].d[m][wt]['ODmn']-self.normStandard['OD'])
                    ind=np.where((self.data[expt].d[m][wt]['ODmn']-self.normStandard['OD'])==mmin)
                except:
                    print('couldn''t get to normalise in '+str(m)+' '+str(wt))
                self.normalizationFactor[expt]=np.mean(self.data[expt].d[m][wt][self.normStandard['channel']][ind, :])
    def loadInitial(self):
        '''
        Search for files in the source path. the pr files can be either:
        1. pickles created with the pickle library
        2. folders, each containing two excel datasheets: a) data file and b)contents file (characterised by ending in contents.xls or contents.xlsx
        '''
        if path.exists(self.source) and path.isdir(self.source):
            dirs = os.listdir(self.source)
            for entry in dirs:
                print('trying to import '+entry)
                if entry.endswith('.pkl') and (self.ignoreFiles == False or ((entry in self.ignoreFiles)==False)) and (self.onlyFiles==False or entry in self.onlyFiles) :
                    pkl = open(self.source + '/' + entry, 'rb')
                    p=pickle.load(pkl,encoding= self.encoding)
                    exptname=p.name.split('/')[-1]
                    self.data[exptname] = p
                    self.exptInfo[exptname]= {}
                    self.exptInfo[exptname]['type']='pickle'
                    self.exptInfo[exptname]['datapath']=self.data[exptname].name+'.xlsx'
                    self.exptInfo[exptname]['dsheetnumber']=self.data[exptname].dsheetnumber
                    try:
                        self.exptInfo[exptname]['contentspath']=self.data[exptname].aname
                    except:
                        print('failed to find a contents file name')
                    self.activeExpts.append(exptname)
                    pkl.close()
                    self.experimentDuration[exptname]=self.data[exptname].t[-1]
        ##if the path its                
                elif path.exists(self.source) and os.path.isdir(self.source+'/'+entry) and (self.ignoreFiles == False or ((entry in self.ignoreFiles)==False)) and (self.onlyFiles==False or entry in self.onlyFiles):
                    #this generally means that the path is a directory. we look inside to see if there are any excel files that look like pr files.
                    datafile=[]
                    contentsfile=[]
                    pth=self.source+'/'+entry; #store the full path of the directory.
                    #look for 
                    files=os.listdir(pth)
                    hasExcel=0;
                    for f in files:
                        if f.endswith('.xls') or f.endswith('.xlsx'):
                            hasExcel+=1;
                            if f.endswith('contents.xls') or f.endswith('contents.xlsx'):
                                contentsfile=self.source+'/'+entry+'/'+f
                            else: #it is not a contents file so we check for, or do preprocessing.
                                if self.preprocess==True:
                                    #                                     if f.endswith('.xslx') and not f.endswith('preprocessed.xlsx'): #freshly preprocess everytime.
                                    #                                         preprocessExpt(expt, main=self.consensusFL, supporting=self.supportingFL)
                                    #                                         datafile=self.source+'/'+entry+'/'+f.split('.xlsx')[0]+'preprocessed.xslx'
                                    #                                 else:
                                    if not f.endswith('preprocessed.xlsx'): #if no preprocessing required import file wo preprocessing
                                        datafile=self.source+'/'+entry+'/'+f
                        if hasExcel>=2:
                            print('directory with excel files will be incorporated.')
                            #the folder may contain more than one experiments and the name of the experiment may differ from
                            #that of the folder. therefore we use the experiment filename as a unifying reference
                            if contentsfile and datafile:
                                #creating pr file from scratch
                                p=pr.platereader(datafile, contentsfile)
                                exptname=p.name.split('/')[-1]
                                self.data[exptname]=p 
                                print('successfully loaded a platereader experiment '+exptname)
                                self.activeExpts.append(exptname)
                                self.exptInfo[exptname]= {}
                                self.exptInfo[exptname]['type']='datasheet'
                                self.exptInfo[exptname]['datapath']=self.data[exptname].name+'.xlsx'
                                self.exptInfo[exptname]['contentspath']=self.data[exptname].aname
                                self.experimentDuration[exptname]=self.data[exptname].t[-1]
        #pickle.dump(self.data, open('xpr_startingdata.pkl', 'wb'))
        #pickle.dump(self, open('xprBackup.pkl', 'wb'))
    def preprocessAll(self, fillnans=True, normOD=0.8, extension='preprocessed.xlsx'):
        '''
        preprocessExpt(expt, fillnans=True, standardize=True, mapFL=True, normOD=0.8, extension='preprocessed.xlsx' )
        preprocess all experiments to fix problems in measurement and normalize data
        '''
        for exp in self.allexperiments:
            preprocessExpt(self.data[exp], fillnans=fillnans, extension=extension, supporting='GFP60', main='GFP80')
    def stageData():
        ##future: add way to verify that data changes have been made
        pickle.dump(self.data, open('xpr_datastage'+str(self.currentStage+1).zfill(2)+'.pkl', 'wb'))
        self.currentStage=self.currentStage+1
        #used to save the data processing at different stages so that processing can be reversed at any point
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
    def loadFresh(self, exptList=[]):
        '''
        loadFresh(exptList=self.activeExpts)
        load experiments fresh from the original files.
        '''
        if not exptList:
            exptList=self.activeExpts
        #clearing up previously loaded data.
        self.data={}
        self.checkallstats()
        #from pickle.
        for entry in exptList:
            #making sure the entry is not in ignorefiles or is in onlyfiles
            #if (self.ignoreFiles == False or ((entry in self.ignoreFiles)==False)) and (self.onlyFiles==False or entry in self.onlyFiles) :
            print('trying to open ', entry, 'from the information in exptInfo')
            self.data[entry] = pr.platereader(self.exptInfo[entry]['datapath'],self.exptInfo[entry]['contentspath'])
    def findMachines(self):
        for expt in self.allexperiments:
            try:
                self.machines[expt]= self.data[expt].machine
            except:
                print('experiment '+expt+' lacking ''machine'' attribute. Attempting retrieval from file...')
                try:
                    pth=self.data[expt].wdir + self.data[expt].dname
                    #print(pth)
                    xl= pd.ExcelFile(self.data[expt].wdir + self.data[expt].dname)
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
                    if self.serialnumbers[expt]=='1006006292':
                        self.machines[expt]='Plate Reader 1'
                    else:
                        if self.serialnumbers[expt]=='1411011275':
                            self.machines[expt]='Plate Reader 2'
                        else:
                            self.machines[expt]='unknown'
                except:
                    print('could not extract machine info from source file. Machine left unknown.')					
                    self.machines[expt]= np.nan
    def assignFL(self, mainFL='noFL', mainFLperod='noFLperod'):
        '''
        every experiment can have their own fluorescence. 
        we search for the consensus, and select the consensus as the main fluorescence
        '''
        
        GFPFields=['GFP' in key for key in self.statcontents.keys()]
        ###Getting fields that contain GFP and are present in all experiments		
        self.consensusFLs=self.statcontents.keys()[GFPFields][np.where([np.around(sum(self.statcontents[field]))== np.size(self.statcontents,0) for field in  self.statcontents.keys()[GFPFields]])]
        
        if mainFL!= 'noFL':  ###priority to specified fluorescence
            for expt in self.data.keys():
                self.FL[expt]['mainFL']=mainFL
                if mainFL in list(self.consensusFLs)==False:
                    print('Warning: '+mainFL+' is not a consensus Fluorescence. Consider changing or leaving unspecified')
        if mainFLperod!= 'noFLperod':  ###priority to specified fluorescence per OD
            for expt in self.data.keys():
                self.FL[expt]['mainFLperod']=mainFLperod
                if mainFLperod in list(self.consensusFLs)==False:
                    print('Warning: '+mainFLperod+' is not a consensus Fluorescence per od. Consider changing or leaving unspecified')

        if mainFLperod in list(self.consensusFLs)==False and mainFL in list(self.consensusFLs): ##if the fl channel is there but there is no cfl for that:
            mainFLperod='c-'+mainFL+'perod' ##we force the FLperod to match FL
            for expt in self.data.keys():
                self.FL[expt]['mainFLperod']=mainFLperod
                self.FL[expt]['mainFLperodvar']=mainFLperod+'var'
        
        #if there is no specification, we look for consensus
        if mainFL=='noFL' and mainFLperod=='noFLperod' and np.size(self.consensusFLs)>0: # there is a consensus but was not specified in the arguments
            flag=0;
            allFLODs=self.consensusFLs[['perod' in j for j in list(self.consensusFLs)]] 
            allFLs=self.consensusFLs[['perod' in j==False for j in list(self.consensusFLs)]] 
            if np.size(allFLODs)>0: ###if there was any perod in the consensus
                if np.size(allFLODs)>1: ## if there is more than one, we select the first one as a default
                    print('More than one consensus FLperod detected. Assigning '+str(allFLODs[0])+' as the default.')
                mainFLperod=str(allFLODs[0])
                for expt in self.data.keys():
                    self.FL[expt]['mainFLperod']=mainFLperod
                    self.FL[expt]['mainFLperodvar']=mainFLperod+'var'
                    
            if np.size(allFLs)>0: ###if there was any FL in the consensus
                if np.size(allFLs)>1: ## if there is more than one, we select the first one as a default
                    print('More than one consensus FL detected. Assigning '+str(allFLs[0])+' as the default.')
                mainFL=str(allFLs[0])
                for expt in self.data.keys():
                    self.FL[expt]['mainFL']=mainFL
        else:
            if np.size(self.consensusFLs)==0:
                print('No consensus fluorescence was found. please specify fluorescence channels to work with. Defaulting to GFP.')
                fl='GFP'
                for expt in self.data.keys():
                    self.FL[expt]['mainFL']=fl
                    self.FL[expt]['mainFLperod']='c-'+fl+'perod'
                    self.FL[expt]['mainFLperodvar']='c-'+fl+'perodvar'
            if mainFL!='noFL':
                print('\n Adding '+mainFL+' as a consensus. Please correct autofluorescence using this channel.') ### if no consensus but there was mainFL inserted
                for expt in self.data.keys():
                    self.FL[expt]['mainFL']=mainFL
                    self.FL[expt]['mainFLperod']='c-'+mainFL+'perod'
                    self.FL[expt]['mainFLperodvar']='c-'+mainFL+'perodvar'
            else:
                if mainFLperod!='noFLperod':
                    fl= mainFLperod.split('perod')[0]
                    try:
                        fl=fl.split('c-')[1]
                    except:
                        print('channel name does not include leading c-') 
                    print('adding '+mainFLperod+' as a consensus corrected fluorescence. assuming raw Fluorescence to be '+fl+'. Please correct if necessary') ### if no consensus but there was mainFL inserted
                    for expt in self.data.keys():
                        self.FL[expt]['mainFL']=fl
                        self.FL[expt]['mainFLperod']='c-'+mainFL+'perod'
                        self.FL[expt]['mainFLperodvar']='c-'+mainFL+'perodvar'
                else:
                    print(' \nplease add fluorescence channels manually. \ncheck the contents using xpr.containsstat() and xpr.statcontents.')
                
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
        D = self.data
        self.allmedia=set()
        self.allstrains=set()
        self.allexperiments=[]
        for exp in D.keys():
            self.allexperiments.append(exp)
            if verbose==True:
                print('\n\nExperiment ('+exp+') contains:\n------------------ ')
            for media in D[exp].d.keys():
                if verbose==True:
                    print('\nMedium ('+media+'): ',)
                self.allmedia.add(media)
                for strain in sorted(D[exp].d[media].keys()):
                    if verbose==True:
                        print( strain +',')
                    self.allstrains.add(strain)
                
    def containsstat(self, stat, printstats=False):
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
    #def backup(self):
    #    '''generates numbered backup pickles to backtrack processing in case something fails'''
    def resetFL(self):
        '''resets the processing of fluorescence in order to start from scratch in case something went wrong'''
        for expt in self.allexperiments:
            xpr.data[expt].reset()
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
        D = self.data
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
    
    def getAllContents(self):
        '''
        getAllContents(self)
        creates two dataframes in the accesspr object:
            self.allConditions.- a dataframe of unique media/strain conditions in all experiments (useful for looping through all conditions once)
            self.allcontents.- a dataframe containing reating experiment, media, strain and plate locations of each condition.
            This object is useful when you want to only consider or merge certain wells or biological replicates 
        '''
        cl= pd.DataFrame()
        rl=pd.DataFrame() #dataframe of the replicates
        for key in self.data.keys():
            #print(cl)
            cl= pd.concat([cl, ppf.conditionList(self.data[key])], ignore_index=True)
            rl= pd.concat([rl, ppf.replicateList(self.data[key])], ignore_index=True)
        self.allReplicates=rl
        self.allcontents=cl
        self.allConditions= pd.DataFrame(np.column_stack([cl.values[:,1], cl.values[:,2]]), columns=['media', 'strain'])
        self.allConditions=self.allConditions.drop_duplicates()

    def correctauto(self, f=['GFP', 'AutoFL'], experiments='all',media='all', strains='all', refstrain=['WT'], figs=True, correctOD=True, noruns=2, bd=False, no1samples=100, rewrite=False, rerun=False, correctmedia=True, mediausemean=False, ignoreneg=True):
        ''' function designed to run correctauto on all experiments, or combinations of media and strains. 
        '''
        if experiments=='all':
            experiments=list(self.data.keys())
        self.containsstat('c-'+f[0]+'perod')
        for key in experiments: 
            #if self.statcontents.loc[key, 'c-'+f[0]+'perod']==1 and rerun==False:
            #    print('Experiment ', key, ' already contains FLperod')
            #    continue
            whichref= np.where(np.array([c in self.data[key].allstrains for c in refstrain])) #which ref of the ones added is in this experiment
            localref= np.array(refstrain)[whichref] #hopefully there is never 2 different refs
            try:
                print('experiment', key, ':')
                self.data[key].correctauto(f=f, conditions=media, strains=strains, refstrain=localref, figs=figs, correctOD=correctOD, noruns=noruns, bd=bd, no1samples=no1samples, correctmedia=correctmedia , ignoreneg=ignoreneg, mediausemean=mediausemean)
                plt.close('all')
                break
            except: #LinAlgErr:
                for e in range(2,20):
                    try:
                        print('try number '+str(e))
                        self.data[key].correctauto(f=f, conditions=media, strains=strains, refstrain=localref, figs=figs, correctOD=correctOD, noruns=noruns, bd=bd, no1samples=no1samples, correctmedia=correctmedia, mediausemean=mediausemean, ignoreneg=ignoreneg )
                    except:
                        print('something went wrong.')
                plt.close('all')
        if rewrite==True:
            if path.isdir(self.source):
                pickle.dump(self.data[key], open(self.source + '/' +key, 'wb'))
                print('Experiment pickle file'+ self.source+key+' has been replaced.')
        self.containsstat('c-'+f[0]+'perod')

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
        cl= ac.DFsubset(self.allConditions, 'strain', ['null'])
        cl=cl.reset_index(drop=True)
        for j in range(0, np.size(cl, 0)):
            md=cl.loc[j, 'media']
            st=cl.loc[j, 'strain']
            xpr.containssetup(cl.loc[j, 'media'], cl.loc[j, 'strain'], strict=False)
            expts=xpr.containslist
        for expt in expts:
            plt.scatter(self.data[expt].d[md][st]['OD'], self.data[expt].d[md][st][self.FL[expt]['mainFL']], color=exptColors[expt])
        plt.xlabel('OD of null')
        plt.ylabel('FL of null')
        createFigLegend(dic=exptColors)
    def getstats(self, experiments='all', conditionsDF=False, media='all', strain='all', dtype='OD', rewrite=False, bd=False, cvfn='sqexp', esterrs=False, stats=True, plotodgr=False, rerun=False, exitearly= False, linalgmax= 5):
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
        if experiments=='all':
            experiments = self.data.keys()
        if type(experiments)==str and experiments != 'all':
            experiments=[experiments]
        print('Performing non-parametric fitting of growth statistics using gaussian processes.\n')
        for key in experiments: 
            print('...\nExperiment', key, ':\n')
            if self.statcontents.loc[key, 'gr']>.95 and dtype=='OD' and rerun==False:
                print('Experiment ', key, ' already contains gr')
                continue
            #if self.statcontents.loc[key, 'd/dt FLperod']==1 and dtype==['FLperod'] and rerun==False:
            #   print('Experiment ', key, ' already contains d/dt FLperod')
            #   continue
            cl=ppf.conditionList(self.data[key])[['media', 'strain']].values
            for x in range(0, np.size(cl,0)):
                if dtype=='FLperod' and ('WT' in cl[x,1] or 'null' in cl[x,1]):
                    continue
                else:
                        self.data[key].getstats(dtype, conditions=cl[x,0], strains=cl[x,1], esterrs=esterrs, bd=bd, cvfn=cvfn, stats=stats, plotodgr=plotodgr)
                        plt.close('all')
                        #print('something went wrong. carrying on.')
            if rewrite==True:
                if path.isdir(self.source):
                    pickle.dump(self.data[key], open(self.source + '/' +key, 'wb'))
                    print('Experiment pickle file'+ self.source+key+' has been replaced.')
        self.containsstat('gr')

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
                            centeredTimeFL=self.data[key].t[np.where(self.data[key].d[media][strain][mainFL+'mn']==max(self.data[key].d[media][strain][mainFL+'mn']))]
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
                            centeredTimeFL=self.data[key].t[np.where(self.data[key].d[media][strain][mainFL+'mn']==max(self.data[key].d[media][strain][mainFL+'mn']))]
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
            self.containsstat('gr')
            self.containsstat('Time centered at gr peak')
            self.containsstat('Time centered at FL peak')
            if 'Time centered at gr peak' in list(self.statcontents.columns.values):
                grAlignedFraction=self.statcontents['Time centered at gr peak'].sum()/ np.size(self.statcontents['Time centered at gr peak'])
            else:
                grAlignedFraction=0
            if grAlignedFraction==1:
                self.aligned=True
                print('Experiments aligned successfully.')
            else:
                self.aligned=False
                print( str(self.statcontents['Time centered at gr peak'].sum())+" out of "+str(np.size(self.statcontents['Time centered at gr peak']))+"experiments aligned. \n Tips: \n .getstats() calculates growth statistics from all experiments. \n.statcontents lets you see which experiments need to get have gr or FLperod. \n.alignAll(rerun=True) to try aligning again")
            print('Experiments have already been aligned. to realign, try rerun=True')
    def plotReplicateMeanNew(self, media=False, strain=False, conditionsDF=False, experiments='all', ignoreExps=False, dtype='OD', col='Black', alpha=0.2, exceptionShift=0.01, normalise=False, excludeFirst=0, excludeLast=-1, bootstrap=0, centeringVariable='time', factor=False, factorColors=False, factorMarkers=False, loc='upper left'):
        markers=['.', 'o', 's', '^', '+', 'v', '*', 'p', '<', '>', 'h', 'x', 'D']*100
        '''plots mean plus shaded area across all replicates. returns the mean coefficient of variation across replicates.'''
        if centeringVariable=='Time centered at gr peak':
            print('removing null from strain list as it does not have growth rate')
            if isinstance(strain, list):
                strain.remove('null')
            if isinstance(conditionsDF, pd.DataFrame):
                conditionsDF= conditionsDF[conditionsDF['strain'].values != 'null']
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
                    conditionsDF=DFsubset(DFsubset(DFsubset(self.allReplicates, 'media', [m]), 'strain', s), 'experiment', experiments)
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
                    interpolatedFinal=bootstrapInterps(interpolatedToUse, centeringVariable, dtype, bootstrap=bootstrap, col=factorColors[fact], alpha=alpha, marker=factorMarkers[fact], markevery=10)
                    interpolated['Group_'+fact+'_'+dtype+'mn']=interpolatedFinal[dtype+'mn']
                    if bootstrap>0:
                        st=dtype+'bootstrapsd'
                    else:
                        st=dtype+'sd'
                    interpolated['Group_'+fact+'_'+st]=interpolatedFinal[st]
                ppf.createFigLegend(keys=list(factorColors.keys()), cols= list(factorColors.values()), markers=list(factorMarkers.values()), loc=loc ) #creating a legend for the figure
                return(interpolated, factorColors, factorMarkers) 
            else:
                    #interpolated final contains the mean and sd that were plotted
                    interpolatedFinal=bootstrapInterps(interpolated, centeringVariable, dtype, bootstrap=bootstrap, col=col, alpha=alpha)
                    return(interpolatedFinal)
    def plotReplicateMean(self, media, strain, experiments='all', ignoreExps=False, dtype='', col='Black', alpha=0.2, exceptionShift=0.01, normalise=False, excludeFirst=0, excludeLast=-1, bootstrap=0, centeringVariable='Time centered at gr peak'):
        '''plots mean plus shaded area across all replicates. returns the mean coefficient of variation across replicates.'''
        if dtype=='':
            try:
                dtype=self.consensusFLperod
            except:
                dtype='OD'
        if np.size(media)==1:
            mediaName=media
            media=[media]
        cv=[]
        for m in media:
            #print('processing '+m)
            interpolated= self.interpTimes(m, strain, dtype=dtype, experiments=experiments, ignoreExps=ignoreExps, centeringVariable=centeringVariable)
            interpolated[dtype]=interpolated[dtype][excludeFirst:excludeLast]
            interpolated[centeringVariable]=interpolated[centeringVariable][excludeFirst:excludeLast]
            if normalise==True:
                interpolated[dtype]=interpolated[dtype]/np.nanmax(flatten(interpolated[dtype])) 
            mn=np.nanmean(interpolated[dtype],1)
            sd=np.nanstd(interpolated[dtype],1)
            #calculating the coefficient of variation
            cv.append(np.nanmean(sd/mn))
            ##generate bootstrap replicates
            ncols=np.size(interpolated[dtype], 1)
            nrows=np.size(interpolated[dtype], 0)
            reps=mn#np.nans(nrows)#create a dummy column for the bootstraps
            ###generating new traces  from combining other traces
            for j in range(0, bootstrap): ###till the number of bootstrap replicates is reached
                ###sample a number of replicates
                addmn=interpolated[dtype][:, sample(range(0,ncols),randint(1,ncols))].mean(1)  ###sample from one to ncols -1 , and get those specific replicates from the set
                reps=np.column_stack([reps, addmn])
            if bootstrap>0:
                totalmn=np.nanmean(reps,1)
                totalsd=np.nanstd(reps,1)
                plt.plot(interpolated[centeringVariable], totalmn, color=col)
                plt.fill_between(interpolated[centeringVariable], mn-totalsd, mn+totalsd, color=col, alpha=alpha, label=strain+' in '+m )
                plt.xlabel(centeringVariable)
                plt.ylabel(dtype)
            else:
                plt.plot(interpolated[centeringVariable], mn, color=col)
                plt.fill_between(interpolated[centeringVariable], mn-sd, mn+sd, color=col, alpha=alpha, label=strain+' in '+m )
                plt.xlabel(centeringVariable)
                plt.ylabel(dtype)
                ###then get their mean
                ###then add it to the column         
        return cv

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

    def plotConditionAligned(self, media, strain, dtype='OD', experiments='all', normalize=0, centerAtFLPeak=0,col='black', range=0):
        '''
        Plots the specified datatype for specified condition for all suitable
        experiments in one plot, colouring all of them the same way, allowing to normalize or cnter at fluorescence
        If aligned=True, the data are plotted against the aligned time vector.
        '''
        self.alignAll()
        self.containssetup(media, strain)
        self.media = media
        self.strain = strain
        expts=self.containslist
        plt.title('aligned '+ dtype + ' over time across replicates', color = 'b')
        timestr=['Time relative to growth rate peak time (hrs)', 'Time relative to Fluorescence peak time (hrs)']
        plt.xlabel(timestr[centerAtFLPeak])
        plt.ylabel(dtype + ' in a.u.')
        legends = []
        if normalize==1:
            for key in expts:
                if experiments!='all': 
                    if key in experiments==False:
                        continue 
                if ppf.hasKey(self.data[key].d, media) and ppf.hasKey(self.data[key].d[media], strain)\
                    and ppf.hasKey(self.data[key].d[media][strain], dtype):
                    plt.plot(self.data[key].d[media][strain]['Time centered at gr peak'], normalizeOverTime(self.data[key], self.media, self.strain, dtype=dtype),color=col)
                    plt.ylabel('normalized '+ dtype + ' in a.u.')
                else: continue
        else:
            for key in expts:
                if experiments!='all': 
                    if key in experiments==False:
                        continue  
                if ppf.hasKey(self.data[key].d, media) and ppf.hasKey(self.data[key].d[media], strain)\
                    and ppf.hasKey(self.data[key].d[media][strain], dtype):
                    plt.plot(self.data[key].d[media][strain]['Time centered at gr peak'], self.data[key].d[self.media][self.strain][dtype],color=col)
                else: continue 
        #plt.show(block=False)

    def plotConditionStd(self, media, strain, dtype='OD', normalize=0, centerAtFLPeak=0,col='black', range=0):
        '''
        Purpose: compare conditions from the results of different experiments
        Plots the specified datatype for specified condition for all suitable
        experiments in one plot, colouring all of them the same way, allowing to normalize or cnter at fluorescence
        If aligned=True, the data are plotted against the aligned time vector.
        '''
        self.containssetup(media, strain)
        plt.title('aligned '+ dtype + ' over time across replicates', color = 'b')
        timestr=['Time relative to growth rate peak time (hrs)', 'Time relative to Fluorescence peak time (hrs)']
        plt.xlabel(timestr[centerAtFLPeak])
        plt.ylabel(dtype + ' in a.u.')
        legends = []
        self.align(self.media, self.strain)
        if normalize==1:
            for key in self.containslist: 
                if ppf.hasKey(self.data[key].d, media) and ppf.hasKey(self.data[key].d[media], strain)\
                    and ppf.hasKey(self.data[key].d[media][strain], dtype):
                    plt.plot(self.data[key].t, normalizeOverTime(self.data[key], self.media, self.strain), color=col)
                else: continue
        else:
            for key in self.containslist: 
                if ppf.hasKey(self.data[key].d, media) and ppf.hasKey(self.data[key].d[media], strain)\
                    and ppf.hasKey(self.data[key].d[media][strain], dtype):
                    plt.plot(self.data[key].t, self.data[key].d[self.media][self.strain][dtype], color=col)
                else: continue
        plt.show(block=False)

    def plotRawReplicates(self, media, strain, dtype='OD',xlim=False, ylim=False, experiments='all', exptColors=False, addLegend=True, xstat=[], color=False):
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
            [plt.plot(arrx, arr[:,k], label= x+' '+plateloc[k], marker=markers[k], markevery=10, color= exptColors[x]) for k in range(0, len(plateloc))]
        if addLegend==True:
            #plt.figlegend(patches, self.containslist, 'upper right')
            plt.legend()
        plt.title('Replicates of '+strain+ ' in '+media)
        if ylim:
            plt.ylim(ylim)
        if xlim:
            plt.xlim(xlim)
        return exptColors

    def timeStatAligned(self, media, strain, times=[0], dtype='OD', includeMediaStrain=True):
        '''
        df=timeStatAligned(self, media, strain, times=[0], dtype='FLperod', includeMediaStrain=True):
        Generates a dataframe where columns are media, strain and the values of dtype at times (one column per element in times).
         and and each row is the value of dtype at those times for one replicate.
        '''
        it= self.interpTimesNew(replicateMatrix=DFsubset(DFsubset(self.allReplicates, 'media', [media]), 'strain', [strain]),dtype=dtype)
        if includeMediaStrain==True:
            fin=pd.DataFrame( columns= ['experiment', 'machine','media', 'strain']+ times)
        else:
            fin=pd.DataFrame(columns= times)
        for j in range(0, np.size(it[dtype],1)):
            f = scint.interp1d(it['time'], it[dtype][:,j])
            fin.loc[j, 'experiment']=self.containslist[j]
            fin.loc[j,'machine']=self.machines[self.containslist[j]]
            fin.loc[j,'media']= media
            fin.loc[j,'strain']= strain
            fin.loc[j,times]= f(times)
        return fin
    def timeStatAll(self, times, media='all', strains='all', aligned=False, dtype='OD', scale=False, subtractBackground=False):
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
        for j in range(0, np.size(self.allConditions,0)):
            if self.allConditions.values[j,1]=='null':
                continue
            try:
                if aligned==True:
                    df=pd.concat([df, self.timeStatAligned(self.allConditions.values[j,0],self.allConditions.values[j,1], times, dtype=dtype)])
                if aligned==False:
                        df=pd.concat([df, self.timeStat(self.allConditions.values[j,0],self.allConditions.values[j,1], times, dtype=dtype, scale=scale, max=mx)])
            except:
                print('condition ', self.allConditions.values[j,1], ' in ' ,self.allConditions.values[j,0], 'showed extraction problems.')
        df.index= range(0,np.size(df,0))
        if media !='all': #if media subset is entered then filter dataframe by that subset
            df=DFsubset(df, 'media', media)
        if strains !='all':#if strain subset is entered then filter dataframe by that subset
            df=DFsubset(df, 'strain', strains)
        return df
        
    def timeStat(self, media, strain, times= [4], dtype='OD', includeMediaStrain=True, scale=False, max=False, subtractBackground=False, background=False ):
        '''
        df=timeStat(self, media, strain, times=[4], dtype='FLperod', includeMediaStrain=True):
        Generates a dataframe where columns are media, strain and the values of dtype at times (one column per element in times).
        and and each row is the value of dtype at those times for one replicate. 
        The times are absolute from the beginning of each experiment. Default is at 4 hrs.
        '''
        self.containssetup(media,strain, strict=False) #finding the experiments with a given media and strain
        expts=self.containslist ##here are the compliant experiments.
        isWithinBound= [j<=xpr.interpLimit for j in times] ###true or false array saying 
        #whether the times fall within the length of the shortest experiment
        if np.size(np.where(isWithinBound==False))>0: #warn that times out of the bound will be excluded.
            print('Warning: times above ', self.interpLimit, ' are being excluded')
        times= list(np.array(times)[isWithinBound]) ###we immediately get rid of the times outside the duration of the shortest experiment.
        
        if includeMediaStrain==True: ##whether the user wants to iinclude the media and strain columns
            fin=pd.DataFrame( columns= ['experiment','machine','media', 'strain']+ times) #preparing the output dataframe with media and strain
        else:
            fin=pd.DataFrame( columns= times) #preparing the output dataframe w/o media and strain
        for j in range(0, np.size(expts)): ##loop through experiment names
            if np.ndim( self.data[expts[j]].d[media][strain][dtype])>1:
                datavals= self.data[expts[j]].d[media][strain][dtype].mean(1)
            else:
                datavals=self.data[expts[j]].d[media][strain][dtype]
            if subtractBackground==True: ##if the user wants to subtract a minimum value of the signal
                if background == False: ##if no general minimum is given we use the minimum level of the signal
                    background= np.min(datavals)
            else:#if subtraction of the minimum is not needed we make background 0.
                background=0
            datavals=datavals-background #whatever the value of background is at this point, we subtract it from the signal datavals.
            if scale==True: ##if the user wants to scale the response relative to a maximum
                if max != False: ##if no general maximum is given we use the maximum level in the signal
                    f = scint.interp1d(self.data[expts[j]].d[media][strain]['time'],  datavals/max) #create interpolation function
                else:#if a max is given we divide the signal by it.
                    f = scint.interp1d(self.data[expts[j]].d[media][strain]['time'],  datavals/np.max(datavals) ) #create interpolation function
            else: #if scaling is not required then we just proceed with the original time series.
                f = scint.interp1d(self.data[expts[j]].d[media][strain]['time'],  datavals) #create interpolation function
            if includeMediaStrain==True: # add media and strain if required
                fin.loc[j,'experiment']= expts[j]
                fin.loc[j,'machine']=self.machines[self.containslist[j]]
                fin.loc[j,'media']= media
                fin.loc[j,'strain']= strain
            fin.loc[j,times]= f(times) #the dataframe at the given expt's index and  and at the times columns requested will be obtained by interpolating from the data
        return fin
        
    def replicateLocations(self):
        if not self.allConditions.empty:
            daf=pd.DataFrame(index=[self.allConditions['strain'], self.allConditions['media']], columns=self.allexperiments)
            daf=daf.fillna(0)
            for j in daf.index:
                self.containssetup(media=j[1], strain=j[0], strict=False)
                for c in self.containslist:
                    daf.loc[j, c]=1 
            self.conditionLocTable=daf
            self.numReplicates=pd.DataFrame(columns=['strain', 'media', 'numReplicates'])
            self.numReplicates['numReplicates']=daf.sum(1).values
            self.numReplicates['strain']=self.allConditions['strain']
            self.numReplicates['media']=self.allConditions['media']
            
        
    def text2d(self, media=False, strains=False, dimx='FinalOD', dimy='maxGR', strainColors=False, xlim=False, ylim=False, markersize=500, newFig=True):
        if newFig==True:
            plt.figure()
        patches=[]
        legs=[]
        plt.xlabel(dimx); plt.ylabel(dimy)
        if strains==False:
            strains=np.unique(list(self.allConditions['strain'].values))
        if strainColors==False:
            strainColors= dict(zip( self.allConditions['strain'].values, ppf.randomColors(np.size(self.allConditions['strain'].values))))
        print(strainColors)
        for strain in  strains:
            print('current strain is'+strain) 
            patches.append(pch.Patch(color= strainColors[strain]))
            legs.append(strain)
            if np.all(media==False):
                media=self.allConditions[self.allConditions['strain']==strain]['media'].values
            for m in media:
                 df=self.extractRepInfo(m, strain) 
                 plt.scatter(df[dimx], df[dimy], marker=r"${}$".format(m, markersize,markersize), s= markersize, color= strainColors[strain])
        plt.figlegend(patches, legs, 'upper right')
        return strainColors
    def interpTimes(self, media, strain, dtype='OD', centeringVariable='time', upperLim=16, exceptionShift=0.01, ignoreExps=False, experiments='all'):    
        '''
        interpTimes(self, media, strain, dtype='FLperod', centeringVariable='time', upperLim=16, exceptionShift=0.01, ignoreExps=False, experiments='all')   
        interpolate all replicates across the same time scale.
        
        '''
        self.containssetup(media, strain, strict=True, musthave=dtype)
        if experiments!='all':
            experimentList=experiments
        else:
            experimentList=self.containslist
        print(experimentList)
        if ignoreExps!=False:
            try:
                [experimentList.remove(j) for j in ignoreExps]
            except:
                print(str(ignoreExps)+'is not on the list')
        interpRange=[]
        maxLengths=[]
        startPoints=[]
        adjustedTimes=dict()
        for expt in experimentList: ###retireving the limiting values for interpolation amongst all experiments.
            #print('processing experiment'+expt+'...')
            adjustedTimes[expt]=dict()
            maxLengths.append(np.around(self.data[expt].d[media][strain][centeringVariable][-1],2))
            startPoints.append(np.around(self.data[expt].d[media][strain][centeringVariable][0],2))
            #startpoints=np.array(startPoints)
            #maxLengths=np.array(maxLengths)
        #print('startpoints: ', startPoints)
        #print('maxlengths: ', maxLengths)
        interpRange=[np.max(np.array(startPoints)), np.min(np.array(maxLengths))]
        #print('Interprange: ', interpRange)
        func= lambda x: fitsBounds(x, interpRange) ### this lambda creates function func which will evaluate whether x falls within the interpRange
        for expt in experimentList:
            #print('entered storage loop')
            fitPoints=np.where([func(x) for x in self.data[expt].d[media][strain][centeringVariable]])
            #print(np.shape(self.data[expt].d[media][strain][dtype])[1])
            #print('fitpoints of '+expt+': ', fitPoints)
            try:
                adjustedTimes[expt][dtype]= self.data[expt].d[media][strain][dtype][fitPoints] 
            except:
                try:
                    print(dtype+' cannot be processed as single dimension. attempting column-averaging...')
                    adjustedTimes[expt][dtype]= self.data[expt].d[media][strain][dtype][fitPoints,:].mean(1)
                except:
                    print('impossible to process '+dtype+'. Please make sure the statistic exists and is time x replicate array.')
            try:
                adjustedTimes[expt][centeringVariable]=np.around(self.data[expt].d[media][strain][centeringVariable][fitPoints],2)
            except:
                print('problems retrieving '+centeringVariable+'. attempting to align specific condition...')
                self.align(media, strain)
                adjustedTimes[expt][centeringVariable]=np.around(self.data[expt].d[media][strain][centeringVariable][fitPoints],2)
            #print('adjustedTimes of '+expt+': ', adjustedTimes[expt]['time']) 
            #print(np.column_stack([adjustedTimes[expt]['time'],adjustedTimes[expt][dtype]]))
        finalDict={};
        finalDict['experiments']=experimentList
        finalDict[dtype]=np.empty([np.size(adjustedTimes[experimentList[0]][centeringVariable]), np.size(experimentList)], dtype=None)

        try:
            for j in range(0, np.size(experimentList)): #arbitrarily taking the first experiment in the list as reference
                fint=scint.interp1d(np.around(adjustedTimes[experimentList[j]][centeringVariable],2),adjustedTimes[experimentList[j]][dtype])
                finalDict[centeringVariable]=np.around(adjustedTimes[experimentList[0]][centeringVariable],2)
                finalDict[dtype][:, j]=fint(np.around(adjustedTimes[experimentList[0]][centeringVariable],2))

        except ValueError:
            print('Warning: time out of range. trying interpolation from full time vectors...')
            #try:
            for j in range(1, np.size(experimentList)): #arbitrarily taking the first experiment in the list as reference
                try:
                    fint=scint.interp1d(np.around(self.data[experimentList[j]].d[media][strain][centeringVariable],2),self.data[experimentList[j]].d[media][strain][dtype])
                except:
                    print(dtype+' cannot be interpolated as single dimension. attempting column-averaging...')
                    fint=scint.interp1d(np.around(self.data[experimentList[j]].d[media][strain][centeringVariable],2),self.data[experimentList[j]].d[media][strain][dtype].mean(1))
                temptime= adjustedTimes[experimentList[0]][centeringVariable] ##[0:-1]
                #temptime[0]=temptime[0]+exceptionShift
                #temptime[-1]=temptime[-1]- exceptionShift
                #print('expt time: ' ,np.around(adjustedTimes[experimentList[j]]['time'],2))
                #print('interpolation time: ' ,temptime)
                finalDict[dtype][:, j]=fint(temptime)
                finalDict[centeringVariable]=temptime
            #except ValueError:
            #print('Error: interpolation time out of bounds. please screen for time discrepancies')
        return finalDict
    def interptimesnew(self, replicateMatrix=False, dtype='OD', centeringVariable='time', descriptors=False):    
        '''
        interpTimes(self, media, strain, dtype='OD', centeringVariable='time')   
        interpolate all replicates across the same time scale.
        replicateMatrix is a matrix of the format of self.allReplicates with the specific replicates to interpolate
        dtype: data type to extract. it must be present in all wells.
        centeringVariable: the variable used for the x axis
        descriptors: the final object will contain full descriptor factors for each condition:  expt, media, strain, plateloc
        '''
        interpRange=[]
        maxLengths=[]
        startPoints=[]
        adjustedTimes=dict()
        adjustedTimes[dtype]=[] #working with ordered lists
        adjustedTimes[centeringVariable]=[] #working with ordered lists
        for j in range(0, np.size(replicateMatrix, 0)): ###retireving the limiting values for the interpolation amongst all experiments.
            expt, media, strain=replicateMatrix.values[j, [0,1,2]]
            maxLengths.append(np.around(self.data[expt].d[media][strain][centeringVariable][-1],2))
            startPoints.append(np.around(self.data[expt].d[media][strain][centeringVariable][0],2))
        interpRange=[np.max(np.array(startPoints)), np.min(np.array(maxLengths))]
        func= lambda x: fitsBounds(x, interpRange) ### this lambda creates function func which will evaluate whether x falls within the interpRange
        namesarray=[]
        exptsarray=[]
        mediaarray=[]
        strainsarray=[]
        platelocarray=[]
        for j in range(0, np.size(replicateMatrix, 0)):
            #print('processing replicate number', j)
            expt, media, strain, plateloc=replicateMatrix.values[j, [0,1,2,3]]
            namesarray.append(expt+' '+media+' '+strain+' '+' '+plateloc)
            if descriptors==True:
                exptsarray.append(expt)
                mediaarray.append(media)
                strainsarray.append(strain)
                platelocarray.append(plateloc)
            #finding the points that fit the range
            fitPoints=np.where([func(x) for x in self.data[expt].d[media][strain][centeringVariable]])
            #print(np.shape(self.data[expt]. d[media][strain][dtype])[1])
            #Index in the data matrix that corresponds to this specific well
            if not(plateloc in self.data[expt].ignoredwells):
                platelocIndex=np.where([plateloc==k for k in self.data[expt].d[media][strain]['plateloc']])[0][0]
            else:
                numpoints=np.size(self.data[expt].d[media][strain][dtype][fitPoints, :],0)
                adjustedTimes[dtype].append(np.zeros([numpoints,1])*np.nan)
                continue
            if not(dtype.endswith('mn')) and not(dtype.endswith('gr')) and not(dtype.endswith('perod')) and not(dtype.endswith('var')):
                adjustedTimes[dtype].append(self.data[expt].d[media][strain][dtype][fitPoints, platelocIndex])
            else:
                adjustedTimes[dtype].append(self.data[expt].d[media][strain][dtype][fitPoints])
            adjustedTimes[centeringVariable].append(np.around(self.data[expt].d[media][strain][centeringVariable][fitPoints],2))
        finalDict={};
        finalDict[dtype]=np.empty([np.size(adjustedTimes[centeringVariable][0]), len(adjustedTimes[dtype])], dtype=None)
        finalDict[centeringVariable]=np.around(adjustedTimes[centeringVariable][0],2) #arbitrarily taking the times of the first condition as a reference
        for j in range(0, len(adjustedTimes[dtype])): 
            try:
                fint=scint.interp1d(adjustedTimes[centeringVariable][j],adjustedTimes[dtype][j]) #interpolate times j and response j
                finalDict[dtype][:, j]=fint(finalDict[centeringVariable])
            except:
                finalDict[dtype][:, j]=np.empty([np.size(finalDict[centeringVariable],0)])*np.nan
        finalDict['names']=namesarray
        if descriptors==True:
            finalDict['experiment']=exptsarray
            finalDict['media']=mediaarray
            finalDict['strain']=strainsarray
            finalDict['plateloc']=platelocarray
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
    def mergedicts(dictlist):
        final={}
        for dic in dictlist:
            for key in dic.keys():
                final[key]=dic[key]
        return final

    def extractallinfonew(self, replicateDF=[], growthstats=[], flstats=[]):
        if not isinstance(replicateDF, pd.DataFrame):
            replicateDF= self.allcontents
        df=pd.DataFrame(index= replicateDF.index, columns=self.extractionFields)
        for j in range(0, np.size(replicateDF,0)):
            expt, media, strain, plateloc= replicateDF.values[j]
            ##Pending: generate modular extraction functions depending on expt, media, strain, plateloc
        filler=[np.nan for j in range(0, len(replicateDF))]
        if not growthstats:
            growthstats=self.extractionFieldsGrowth
        if not flstats:
            flstats= self.extractionFieldsFL
        growthdata={};
        for stat in  growthstats:   
            try:
                growthdata[stat]=[self.extractStat(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], stat) for j in range(0, len(replicateDF))]
            except:
                growthdata[stat]=[self.extractStat(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], stat) for j in range(0, len(replicateDF))]
        growthdata['growth rate auc']=[self.extractGRAUC(replicateDF.loc[j,:].values[0],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], replicateDF.loc[j,:].values[3]) for j in range(0, len(replicateDF))]
        growthdata['initial OD']=[self.extractInitialODFull(  replicateDF.loc[j,:].values[0],replicateDF.loc[j,:].values[1], replicateDF.loc[j,:].values[2], replicateDF.loc[j,:].values[3]) for j in range(0, len(replicateDF))]
        nfl=self.consensusFL
        FLData={};
        getalstFL=lambda expt, m,s,fl: alignStats(expt, m, s, dtype=nfl)
        alloutsFL=[getalstFL(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1],replicateDF.loc[j,:].values[2],'GFP80') for j in range(0, len(replicateDF))]
        for field in flstats:
            FLData[field]=[np.float(alloutsFL[j][field]) for j in range(0, len(replicateDF))]
        nflod=self.consensusFLperod
        FLODData={};
        getalstFLOD=lambda expt, m,s,fl: alignStats(expt, m, s, dtype=nflod)
        #aligned fluorescence statistics calculation
        alloutsFLOD=[getalstFLOD(self.data[replicateDF.loc[j,:].values[0]],replicateDF.loc[j,:].values[1],replicateDF.loc[j,:].values[2],'c-GFP80perod') for j in range(0, len(replicateDF))]
        for field in flstats:
            try:
                FLODData['FLOD'+field]=[np.float(alloutsFLOD[j][field]) for j in range(0, len(replicateDF))]   
            except:
                FLODData['FLOD'+field]=filler
        return pd.concat([replicateDF, pd.DataFrame.from_dict(mergedicts([growthdata, FLData, FLODData]))], axis=1 )
    def extractRepInfo(self, media, strain, strict=True, onlyGrowth=False):
        '''
        generates a table of basic info for all replicates
        '''
        print('extracting '+strain+' in '+media)
        self.containssetup(media, strain, strict=strict)
        self.align(media, strain)
        containslist=self.containslist
        timepointRange=[]
        FLVector={}
        initialOD=[]
        FinalOD=[]
        InitialRawFL=[]
        InitialFLperOD=[]
        FinalFLperOD=[]
        FLPeak=[]
        FLAbsPeakTime=[]
        FLAlignedPeakTime=[]
        FLperODTS=[]
        FLRawTS=[]
        lagTime=[]
        medialist=[]
        strainlist=[]
        dataset=[]
        realTime=[]
        alignedTime=[]
        experiment=[]
        machine=[]
        maxGR=[]
        maxGRvar=[]
        maxGRTime=[]
        grAUC=[]
        FLperodAUC=[]
        halfFLAreaTime=[]
        halfGRAreaTime=[]
        responseTime=[]
        responseTimeAligned=[]
        steepness=[]
        slope=[]
        FLLag=[] #defined as the time to reach half the total fl from the beginning
        for i in range(0, np.shape(self.containslist)[0]): 
            #if the mock assignment below crashes, there is no FLperod
            #checkingFL= self.data[containslist[i]].d[media][strain][nflod]
            #except:
            #    InitialFLperOD.append(np.nan)
            #    FinalFLperOD.append(np.nan)
            #    FLPeak.append(np.nan)
            #    FLAbsPeakTime.append(np.nan)
            #    FLAlignedPeakTime.append(np.nan)
            #    FLperodAUC.append(np.nan)
            #print('nfl= '+nfl)
            if ppf.hasKey(self.data[containslist[i]].d, media) and ppf.hasKey(self.data[containslist[i]].d[media], strain):
                #plt.plot(self.aligned[key], normalizeOverTime(self.data[key], self.media, self.strain), color=col)
                #out= {'normalizedFL': normalizedFLVector, 'FLPeak': max(rawFLVector), 'normalizedFLPeak':normalizedFLPeak, 'alignedPeakTime':alignedPeakTime, 'absolutePeakTime':absolutePeakTime
                #print 'strain is ', strain, 'strain is not WT: ', strain != "WT"
                if strain != 'null':
                    ods=self.data[containslist[i]].d[media][strain]['OD']
                    maxGR.append(self.data[containslist[i]].d[media][strain]['max growth rate'])
                    maxGRvar.append(self.data[containslist[i]].d[media][strain]['max growth rate var'])
                    maxGRTime.append(self.data[containslist[i]].d[media][strain]['time of max growth rate'])
                    tpr=str(0)+'-'+str(np.size(self.data[containslist[i]].t))
                    timepointRange.append(tpr)
                    iod=np.nanmean(ods[0,:])
                    initialOD.append(iod)
                    lagTime.append(self.data[containslist[i]].d[media][strain]['lag time'])
                    FinalOD.append(np.nanmean(ods[np.size(ods,0)-1,:]))
                    medialist.append(media)
                    strainlist.append(strain)
                    grAUC.append(ppf.statArea(self.data[containslist[i]], media, strain, 'gr'))
                    if self.analyseFL:
                        nflod=self.FL[containslist[i]]['mainFLperod']
                        nfl=self.FL[containslist[i]]['mainFL']
                        #alignStats should give nans in case the fluorescence isn't there
                        out=alignStats(self.data[containslist[i]], media, strain, dtype=nflod)
                        FLAbsPeakTime.append(float(out['absolutePeakTime']))
                        FLAlignedPeakTime.append(float(out['alignedPeakTime']))
                        FLPeak.append(out['FLPeak'])
                        try:
                            responseTime.append(out['responseTime'][0])
                        except:
                            responseTime.append(out['responseTime'])
                        try:
                            responseTimeAligned.append(out['responseTimeAligned'][0])
                        except:
                            responseTimeAligned.append(out['responseTimeAligned'])
                        steepness.append(out['steepness'])
                        slope.append(out['slope'])
                        if hasKey(self.data[containslist[i]].d[media][strain], nflod):
                            InitialFLperOD.append(self.data[containslist[i]].d[media][strain][nflod][0])
                            FinalFLperOD.append(self.data[containslist[i]].d[media][strain][nflod][np.size(self.data[containslist[i]].d[media][strain][nflod])-1])
                            FLperodAUC.append(ppf.statArea(self.data[containslist[i]], media, strain, nflod))
                        else:
                            InitialFLperOD.append(np.nan)
                            FinalFLperOD.append(np.nan)
                            FLperodAUC.append(np.nan)
                        if hasKey(self.data[containslist[i]].d[media][strain], nfl):
                            FLs=self.data[containslist[i]].d[media][strain][nfl]
                            InitialRawFL.append(np.nanmean(FLs[0, :]))
                        else:
                            InitialRawFL.append(np.nan)
                if self.analyseFL==False or strain=='null':
                    FLperODTS.append(np.nan)
                    InitialFLperOD.append(np.nan)
                    FinalFLperOD.append(np.nan)
                    FLPeak.append(np.nan)
                    FLAbsPeakTime.append(np.nan)
                    FLAlignedPeakTime.append(np.nan)
                    FLperodAUC.append(np.nan)
                    halfFLAreaTime.append(np.nan)
                    responseTime.append(np.nan)
                    responseTimeAligned.append(np.nan)
                    slope.append(np.nan)
                    steepness.append(np.nan)
                machine.append(self.machines[containslist[i]])
                experiment.append(containslist[i])
                realTime.append(np.str(self.data[containslist[i]].t.T).replace(']', '').replace('[', '').replace('\n', ''))
                at, ct=alignTime(self.data[containslist[i]], media, strain)
                alignedTime.append(str(at).replace(']', '').replace('[', '').replace('\n', '') )
                #print('FLs is', FLs)
                #halfGRAreaTime.append(float(halfAUCTime(self.data[containslist[i]], media, strain, 'gr')))	
        if self.analyseFL:
            dataset= list(zip(containslist,machine,medialist,strainlist, initialOD,FinalOD, lagTime, realTime, alignedTime, maxGR, maxGRvar,maxGRTime,  grAUC,InitialRawFL, InitialFLperOD, FinalFLperOD, FLPeak, FLAbsPeakTime, FLAlignedPeakTime,FLperodAUC, slope, responseTime, steepness, responseTimeAligned))
            df= pd.DataFrame(data=dataset, columns=self.extractionFields)
        else:
            dataset= list(zip(containslist,machine,medialist,strainlist, initialOD,FinalOD, lagTime, realTime, alignedTime, maxGR, maxGRvar,maxGRTime,  grAUC))
            df= pd.DataFrame(data=dataset, columns=self.extractionFieldsOD)
        return df
    def extractAllInfo(self, excludeNull= False, strict=True, onlyGrowth=False):
        #making a decision on wether to use only growth data or not
        if self.analyseFL or onlyGrowth==False:
            extflds=self.extractionFields
        else:
            extflds=self.extractionFieldsOD
            
        df=pd.DataFrame(columns=extflds)
        self.getAllContents()
        for x in range(0, np.size(self.allConditions,0)):
            if self.allConditions.values[x,1]=='null' and excludeNull==True:
                continue
            #try:
            out=self.extractRepInfo(self.allConditions.values[x,0], self.allConditions.values[x,1], strict=strict)
            df=df.append(out, ignore_index=True)
            #except:
            #    print('Problem extracting '+ self.allConditions.values[x,1]+ ' in '+self.allConditions.values[x,0])
            #    continue
            #print ' begin df:', df, "end df\n"
        df= self.getMediaValues(df)
        return df
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
                    self.containsstat('transformed'+repairWith)
                else:
                    try:
                        if np.isnan(self.statcontents.loc[expt, 'GFP75'])==False:
                            ppf.transformEach(self.data[expt], stat1='GFP75', stat2=dtype, plot=True)
                            ppf.replaceStat(self.data[expt], 'GFP80', replaceWith='transformedGFP75')
                            self.containsstat('transformedGFP75')
                    except:
                        print('no GFP75')
                        if np.isnan(self.statcontents.loc[expt, 'GFP60'])==False:
                            ppf.transformEach(self.data[expt], stat1='GFP60', stat2=dtype, plot=True)
                            ppf.replaceStat(self.data[expt], 'GFP80', replaceWith='transformedGFP60')
                            self.containsstat('transformedGFP60')
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
        return self.data[expt].d[media][strain]['OD'][:,whichisplateloc][-1]
    def extractInitialODFull(self, expt, media, strain, plateloc=False):
        return np.nanmean(self.data[expt].d[media][strain]['OD'][0,:])
    def extractFinalODFull(self, expt, media, strain, plateloc=False):
        return np.nanmean(self.data[expt].d[media][strain]['OD'][-1,:])
    def extractStat(self, expt, m, s, stat):
        if ppf.hasKey(expt.d[m][s], stat):
            return expt.d[m][s][stat]
        else:
            return np.nan
            
    #                     if self.analyseFL:
    #                         nflod=self.FL[containslist[i]]['mainFLperod']
    #                         nfl=self.FL[containslist[i]]['mainFL']
    # 
    #                         print('nfl= '+nflod)
    #                         InitialFLperOD.append(self.data[containslist[i]].d[media][strain][nflod][0])
    #                         FinalFLperOD.append(self.data[containslist[i]].d[media][strain][nflod][np.size(self.data[containslist[i]].d[media][strain][nflod])-1])
    #                         FLPeak.append(out['FLPeak'])
    #                         FLAbsPeakTime.append(float(out['absolutePeakTime']))
    #                         #print('apt = '+str(FLAbsPeakTime))
    #                         FLAlignedPeakTime.append(float(out['alignedPeakTime']))
    #                         FLperODTS.append(np.str(self.data[containslist[i]].d[media][strain][nflod]).replace(']', '').replace('[', ''))
    #                         FLperodAUC.append(ppf.statArea(self.data[containslist[i]], media, strain, nflod))
    #                         FLs=self.data[containslist[i]].d[media][strain][nfl]
    #                         InitialRawFL.append(np.nanmean(FLs[0, :]))
    #                         autoFLs=self.data[containslist[i]].d[media][strain]['AutoFL']
    #                         try:
    #                             responseTime.append(out['responseTime'][0])
    #                         except:
    #                             responseTime.append(out['responseTime'])
    #                         try:
    #                             responseTimeAligned.append(out['responseTimeAligned'][0])
    #                         except:
    #                             responseTimeAligned.append(out['responseTimeAligned'])
    #                         steepness.append(out['steepness'])
    #                         slope.append(out['slope'])
    #                         #hat=halfAUCTime(self.data[containslist[i]], media, strain, nflod)
    #                         #print(hat)
    #                         #try:
    #                         #    halfFLAreaTime.append(float(hat))
    #                         #except:
    #                         #    print('problem calculating the half area time for '+strain+' in '+media)
    #                         #    halfFLAreaTime.append(np.nan)
    #                 if self.analyseFL==False or strain=='null':
    #                     FLperODTS.append(np.nan)
    #                     InitialFLperOD.append(np.nan)
    #                     FinalFLperOD.append(np.nan)
    #                     FLPeak.append(np.nan)
    #                     FLAbsPeakTime.append(np.nan)
    #                     FLAlignedPeakTime.append(np.nan)
    #                     FLperodAUC.append(np.nan)
    #                     halfFLAreaTime.append(np.nan)
    #                     responseTime.append(np.nan)
    #                     responseTimeAligned.append(np.nan)
    #                     slope.append(np.nan)
    #                     steepness.append(np.nan)
    #                 machine.append(self.machines[containslist[i]])
    #                 experiment.append(containslist[i])
    #                 realTime.append(np.str(self.data[containslist[i]].t.T).replace(']', '').replace('[', '').replace('\n', ''))
    #                 at, ct=alignTime(self.data[containslist[i]], media, strain)
    #                 alignedTime.append(str(at).replace(']', '').replace('[', '').replace('\n', '') )