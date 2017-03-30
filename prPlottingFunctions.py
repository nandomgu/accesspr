import platereader as pr
import matplotlib.pyplot as plt
import matplotlib.patches as pch
import numpy as np
import pickle
import pandas as pd
from scipy import integrate as integrate
from scipy.interpolate import interp1d

###PLATE READER PLOTTING FUNCTIONS VERSION 3.1.0.
###By Luis Fernando Montaño. Swain lab, 2016.

normalflperodName='FLperod'
cflperodName='cFLperod'
normalodName='OD'
fitodName='fOD'
grName='gr'

def transformStat(p, stat, scale, shift, genericLabel=True):
    '''Generates a linear transformation of expression values'''
    cl=conditionList(p, excludeNull=True).values #condition list array
    for media in cl[:,0]:
        for strain in cl[:,1]:
            if genericLabel==False:
                p.d[media][strain][stat+'x'+scale+'+'+shift]=p.d[media][strain][stat]*scale+shift
            else:
                p.d[media][strain]['transformed'+stat]=p.d[media][strain][stat]*scale+shift

#### P MEDIA STRAIN modular functions. functions that run given these 3 parameters
	
def containsFL(p, media, strain):
    if strain=='null':
        return 0
    else:
        return 1
    
def growthInflectionPoint(p, media, strain):
    return( p.t[    np.where(np.diff(np.sign( np.gradient(p.d[media][strain]['gr']))))[0][0] ])
    

def halfAUCIndex(vector):
    auc= np.trapz(vector)
    areavector=np.zeros(np.size(vector))
    for j in range(0, np.size(vector)):
        areavector[j]= np.trapz(vector[range(0, j)]) ###each element is the area till that point.
    areavector=areavector/auc
    halfAreaIndex=np.where(abs(areavector-0.5)==np.min(abs(areavector-0.5)))
    return(halfAreaIndex)
    
def halfAUCTime(p, media, strain, stat, aligned=True):
    if aligned== False:
        return(p.t[halfAUCIndex(p.d[media][strain][stat])])
    else:
        return(p.d[media][strain][stat]['Time centered at gr peak'][int(halfAUCIndex(p.d[media][strain][stat]))])

def alignByGrinfp(p, media, strain):

    return(p.t-growthInflectionPoint(p, media, strain))

def normalizeOverTime(p, media, strain, subtractBg=0):
    rawFLVector=p.d[media][strain][normalflperodName]
    noBgFLVector=p.d[media][strain][normalflperodName]-min(p.d[media][strain][normalflperodName])
    normalizedFLVector=(noBgFLVector-np.mean(noBgFLVector))/np.std(noBgFLVector)
    return normalizedFLVector
	
def alignStats(p, media, strain, subtractBg=0):
    rawFLVector=p.d[media][strain][normalflperodName]
    noBgFLVector=p.d[media][strain][normalflperodName]-min(p.d[media][strain][normalflperodName])      
    normalizedFLVector=(noBgFLVector-np.mean(noBgFLVector))/np.std(noBgFLVector)
    normalizedFLPeak=max(normalizedFLVector)
    alignedTimeVector=alignTime(p,media,strain)
    alignedPeakTime=alignedTimeVector[np.where(normalizedFLVector==max(normalizedFLVector))]
    absolutePeakTime=p.t[np.where(rawFLVector==max(rawFLVector))]
    out= {'FLPeak': max(rawFLVector), 'normalizedFLPeak':normalizedFLPeak, 'alignedPeakTime':alignedPeakTime, 'absolutePeakTime':absolutePeakTime}
    return out
    
def plotalign(p, media, strain, col='black', normalize=1, centerAtPeakFL=0):
	"aligns fitted curves for given strains respective to normalized peak or maximum growth rate. this is a standalone function to plot."
	# time t where growth rate is max for given curve
	if centerAtPeakFL==0:
	    centeredTime=p.t[np.where(p.d[media][strain]['gr']==max(p.d[media][strain]['gr']))]
	else:
	    centeredTime=p.t[np.where(p.d[media][strain][normalflperodName]==max(p.d[media][strain][normalflperodName]))]
	# calculates difference betwee
	alignedTimeVector=p.t-centeredTime
	rawFLVector=p.d[media][strain][normalflperodName]
	noBgFLVector=p.d[media][strain][normalflperodName]-min(p.d[media][strain][normalflperodName])
	normalizedFLVector=(noBgFLVector-np.mean(noBgFLVector))/np.std(noBgFLVector)
	peak=np.where(rawFLVector==max(rawFLVector))
	normalizedPeak=np.where(normalizedFLVector==max(normalizedFLVector))
	peakTime=alignedTimeVector[np.where(rawFLVector==max(rawFLVector))]
	#plt.figure()
	#plt.plot(alignedTimeVector, rawFLVector, color=col  )
	#plt.figure()
	if normalize==1:
	    plt.plot(alignedTimeVector, normalizedFLVector, color=col  )
	else:
	    plt.plot(alignedTimeVector, rawFLVector, color=col  )
	out= {'alignedTime': alignedTimeVector , 'rawFL':rawFLVector ,'normalizedFL':normalizedFLVector , 'peakFL':peak,  'peakTime':peakTime, 'gr':p.d[media][strain]['gr']}
	return out

def maxGrowthRate(p, media, strain):
	mgr=max(p.d[media][strain]['gr'])
	return mgr
	
def expressionDerivative(p, strain, media, plotResult=0, col="black"):
	der=np.gradient(p.d[media][strain][normalflperodName])
	if plotResult!=0:
		plt.plot(p.t, der[10:size(der)], color=col)
	return der	

def statDerivative(p, media, strain, dtype='FLperod', plotResult=0, col="black"):
	der=np.gradient(p.d[media][strain][dtype])
	if plotResult!=0:
		plt.plot(p.t, der[10:size(der)], color=col)
	return der	
	
def maxExpRateTime(p, strain, media, plotResult=0, col="black"):
	#make sure to load smooth function
	der=np.gradient(smooth(p.d[media][strain][normalflperodName]))
	###the derivative is not defined for the first timepoint in the series, therefore we have to push the indices one to the right
	tmer=p.t[np.where(der==max(der))-array([10])] 
	return tmer		
		
def maxGrowthRateTime(p,strain, media):
    tmgr=timeOfMaxGrowthRate=p.t[np.where(p.d[media][strain]['gr']==max(p.d[media][strain]['gr']))]
    return tmgr
        
def expLevelAtMaxGrowthRateTime(p, strain, media):
	
	timeIndex= min(p.t-maxGrowthRateTime(p, strain, media))
	##print( timeIndex
	elmgrt=p.d[media][strain][normalflperodName][timeIndex]
	
	return elmgrt
	

def maxExpressionLevelTime(p, strain, media):
    tmel=timeOfMaxGrowthRate=p.t[np.where(p.d[media][strain][normalflperodName]==max(p.d[media][strain][normalflperodName]))]
    return tmel	
	
def maxExpressionLevel(p, strain, media):
	mexp=max(p.d[media][strain][normalflperodName])
	return mexp

def integrateExpression(p, strain, media):	
	intgrl=np.trapz(p.d[media][strain]['FLperod'], p.t)	
	return intgrl

def lagTime(p, strain, media):
	lt=p.d[media][strain]['odstats']['lag time']
	return lt

####### P FUNC MEDIA STRAIN FUNCTIONS. functions that receive an arbitrary function as an input. the output of func is usually a single value. 


def MeltMinusMgrt(p, func, strain, conc, media, col, plotResult=0, supress=1):
	
	melt=funcByConcentration(p, func, strain, conc, media,col=col, plotResult=0, supress=1)

	mgrt=funcByConcentration(p, func, strain,conc, media,col=col, plotResult=0, supress=1)
	
	if plotResult != 0:
		plt.plot(conc, melt-mgrt, col=col)
		
	return melt-mgrt	


def funcByInitialFL(p, func, media, strain, col='black'):
    mexp=func(p,strain, media)
    plt.plot(p.d[media][strain][normalflperodName][0], mexp , 'o', color=col)
    plt.xlabel('Initial Fluorescence per OD')
    plt.ylabel(func.__name__)
    return mexp
	
	
def funcByConcentration(p, func, strain, concentrations, col, plotResult=0, supress=1):
    ###this function assumes that there is only one concentration per media.
    ###first get the media for which concentrations is not nan
    noNaNKeys=[list(concentrations.keys())[y] for y in np.where([np.isnan(x)==False for x in c.values()])[0]]  
    concentrations=dict(zip(noNaNKeys,[concentrations[x] for x in noNaNKeys])) ### assemble a new concentrations object without nans
    ###now order the concentrations from lowest to highest
    newOrder =[np.where(list(concentrations.values())==x)[0][0] for x in np.sort(list(concentrations.values()))]
    concentrations= dict(zip( [list(concentrations.keys())[x] for x in newOrder], [list(concentrations.values())[x] for x in newOrder]))
    ###calls a function that is dependent on p, strain, and different kinds of media
    funcs=np.zeros(np.size(concentrations.keys()))
    for i in list(concentrations.keys()):
        funcs[i]= func(p, i,strain)
    if plotResult==1:	
        plt.plot(conc, funcs, color=col, label= strain)
        if supress != 1:
            plt.xlabel('Concentration')
            plt.ylabel( func.__name__)
            plt.legend(loc=4)
            return funcs
		

def funcByConcentrationByStrain(p, func, strains, concentrations, orderedMedia, strainColors, markerlist=['o', 'x', 's', '^', '+', 'v', 'D','<', '*'], normalizeToExt=0, subtractBg=0, supp=0):
    plt.figure()
    for i in range(0, np.size(strains)):
        plt.title(func.__name__+' by Concentration')
        plt.ylabel(func.__name__)
        plt.xlabel('Concentration')
        values=funcByConcentration(p, func, strains[i], concentrations, orderedMedia, strainColors[i])
        ##print( values
        plt.plot(concentrations, values, markerfacecolor=strainColors[i], label=strains[i], color=strainColors[i], linewidth=2, markersize=10, markeredgewidth=2, marker=markerlist[i], markeredgecolor='black')
        plt.legend(loc=2) 


def monodCurve(p, strain, conc, media, col, supress=1):
	mgrs=np.zeros(np.size(media))
	for i in range(0, (np.size(media))):
		mgrs[i]= maxGrowthRate(p, strain, media[i])
	plt.plot(conc, mgrs, color=col, label= strain)
	if supress != 1:
		plt.xlabel('Concentration')
		plt.ylabel('Growth rate')
		plt.legend(loc=4)
	return mgrs
	
	
	
#### P STAT MEDIA STRAIN FUNCTIONS. these functions operate on a statistic that can be found in the plate reader object in p.d[MEDIA][STRAIN][STAT]. examples are normalflperodName, and 'gr'

def statByConcentration(p, stat, strain, concentrations, color=0, supress=1):
    if color==0:
        color=randomColor()
    media= list(concentrations.keys())
    for i in range(0, (size(stats))):
        stats[i]= p.d[media[i]][strain][stat]
    plt.plot(np.array(list(concentrations.values())), stats, color=color, label= strain)
    if supress != 1:
        plt.xlabel('Concentration')
        plt.ylabel(stat)
        plt.legend(loc=4)
    return stats

def plotByOD(p, stat, strain, media, col, normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1):
	if normalizeToSelf==1:
		plt.plot(p.d[media][strain]['odata'][:, 0], (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/max(p.d[media][strain][stat]), color=col, label=media, linewidth=2)
	elif normalizeToExt != 0:
		plt.plot(p.d[media][strain]['odata'][:, 0], (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/normalizeToExt, color=col, label=media, linewidth=2)
	elif subtractBg!=0:		
		plt.plot(p.d[media][strain]['odata'][:, 0], p.d[media][strain][stat]-min(p.d[media][strain][stat]), color=col, label=media, linewidth=2)
	else:	
		plt.plot(p.d[media][strain]['odata'][:, 0], p.d[media][strain][stat], color=col, label=media, linewidth=2)
	
	if supress != 1:
		plt.xlabel('OD')
		plt.ylabel(stat)
		plt.legend(loc=4)
		
		

def plotByGrowthRate(p, stat, strain, media, col='blue', normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1, colorSlope=1):
	if colorSlope==1:
	    cols= sign2color(np.gradient(p.d[media][strain]['gr']))
	else:
	    cols=col
	if normalizeToSelf==1:
		plt.scatter(p.d[media][strain]['gr'], (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/max(p.d[media][strain][stat]), color=cols, label=media, linewidth=0.1)
	elif normalizeToExt != 0:
		plt.scatter(p.d[media][strain]['gr'], (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/normalizeToExt, color=cols, label=media, linewidth=0.1)
	elif subtractBg!=0:		
		plt.scatter(p.d[media][strain]['gr'], p.d[media][strain][stat]-min(p.d[media][strain][stat]), color=cols, label=media, linewidth=0.1)
	else:	
		plt.scatter(p.d[media][strain]['gr'], p.d[media][strain][stat], color=cols, label=media, linewidth=0.1)
	
	if supress != 1:
		plt.legend(loc=4)
	plt.xlabel('Growth Rate')
	plt.ylabel(stat)	
	


def plotByGrowthAccel(p, stat, strain, media, col='blue', normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1, colorSlope=1):
	if colorSlope==1:
	    cols= sign2color(np.gradient(p.d[media][strain]['gr']))
	else:
	    cols=col
	if normalizeToSelf==1:
		plt.scatter(np.gradient(p.d[media][strain]['gr']), (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/max(p.d[media][strain][stat]), color=cols, label=media, linewidth=0.1)
	elif normalizeToExt != 0:
		plt.scatter(np.gradient(p.d[media][strain]['gr']), (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/normalizeToExt, color=cols, label=media, linewidth=0.1)
	elif subtractBg!=0:		
		plt.scatter(np.gradient(p.d[media][strain]['gr']), p.d[media][strain][stat]-min(p.d[media][strain][stat]), color=cols, label=media, linewidth=0.1)
	else:	
		plt.scatter(np.gradient(p.d[media][strain]['gr']), p.d[media][strain][stat], color=cols, label=media, linewidth=0.1)
	
	if supress != 1:
		plt.legend(loc=4)
	plt.xlabel('Growth acceleration')
	plt.ylabel(stat)	


def plotRateByGrowthAccel(p, stat, strain, media, col='blue', normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1, colorSlope=1):
	if colorSlope==1:
	    cols= sign2color(np.gradient(p.d[media][strain]['gr']))
	else:
	    cols=col
	if normalizeToSelf==1:
		plt.scatter(np.gradient(p.d[media][strain]['gr']), np.gradient((p.d[media][strain][stat]-min(p.d[media][strain][stat]))/max(p.d[media][strain][stat])), color=cols, label=media, linewidth=0.1)
	elif normalizeToExt != 0:
		plt.scatter(np.gradient(p.d[media][strain]['gr']), np.gradient((p.d[media][strain][stat]-min(p.d[media][strain][stat]))/normalizeToExt), color=cols, label=media, linewidth=0.1)
	elif subtractBg!=0:		
		plt.scatter(np.gradient(p.d[media][strain]['gr']), np.gradient(p.d[media][strain][stat]-min(p.d[media][strain][stat])), color=cols, label=media, linewidth=0.1)
	else:	
		plt.scatter(np.gradient(p.d[media][strain]['gr']), np.gradient(p.d[media][strain][stat]), color=cols, label=media, linewidth=0.1)
	
	if supress != 1:
		plt.legend(loc=4)
	plt.xlabel('Growth acceleration')
	plt.ylabel(stat+' rate')	
	
	
	
	def plotAccelByGrowthAccel(p, stat, strain, media, col='blue', normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1, colorSlope=1):
	    if colorSlope==1:
		    cols= sign2color(np.gradient(p.d[media][strain]['gr']))
	    else:
		    cols=col
	    if normalizeToSelf==1:
		    plt.scatter(np.gradient(p.d[media][strain]['gr']), np.gradient(np.gradient((p.d[media][strain][stat]-min(p.d[media][strain][stat]))/max(p.d[media][strain][stat]))), color=cols, label=media, linewidth=0.1)
	    elif normalizeToExt != 0:
		    plt.scatter(np.gradient(p.d[media][strain]['gr']), np.gradient(np.gradient((p.d[media][strain][stat]-min(p.d[media][strain][stat]))/normalizeToExt)), color=cols, label=media, linewidth=0.1)
	    elif subtractBg!=0:		
		    plt.scatter(np.gradient(p.d[media][strain]['gr']), np.gradient(np.gradient(p.d[media][strain][stat]-min(p.d[media][strain][stat]))), color=cols, label=media, linewidth=0.1)
	    else:	
		    plt.scatter(np.gradient(p.d[media][strain]['gr']), np.gradient(np.gradient(p.d[media][strain][stat])), color=cols, label=media, linewidth=0.1)
	
	    if supress != 1:
		    plt.legend(loc=4)
	    plt.xlabel('Growth acceleration')
	    plt.ylabel(stat+' rate')	


	
def plotScatterOverTime(p, stat, strain, media, col='blue', normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1, colorSlope=1):
	if colorSlope==1:
	    cols= sign2color(np.gradient(p.d[media][strain]['gr']))
	else:
	    cols=col
	if normalizeToSelf==1:
		plt.scatter(p.t, (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/max(p.d[media][strain][stat]), color=cols, label=media, linewidth=0.1)
	elif normalizeToExt != 0:
		plt.scatter(p.t, (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/normalizeToExt, color=cols, label=media, linewidth=0.1)
	elif subtractBg!=0:		
		plt.scatter(p.t, p.d[media][strain][stat]-min(p.d[media][strain][stat]), color=cols, label=media, linewidth=0.1)
	else:	
		plt.scatter(p.t, p.d[media][strain][stat], color=cols, label=media, linewidth=0.1)
	
	if supress != 1:
		plt.legend(loc=4)
	plt.xlabel('Time (Hrs)')
	plt.ylabel(stat)	



def plotByTotalGrowth(p, stat, strain, media, col, normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1):
    #this function treats the max OD reached in a particular condition as 1, so the changes in expression can be seen relative to the max growth rate achieved 
    totalGrowth=p.d[media][strain]['odata'][:, 0]/ max(p.d[media][strain]['odata'][:, 0]) ### we dont subtract by the starting od as it may cause differences between 2 concentrations that reach the same eventual od but started at different points
    if normalizeToSelf==1:
    	plt.plot(totalGrowth, (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/max(p.d[media][strain][stat]), color=col, label=media, linewidth=2)
    elif normalizeToExt != 0:
    	plt.plot(totalGrowth, (p.d[media][strain][stat]-min(p.d[media][strain][stat]))/normalizeToExt, color=col, label=media, linewidth=2)
    elif subtractBg!=0:		
    	plt.plot(totalGrowth, p.d[media][strain][stat]-min(p.d[media][strain][stat]), color=col, label=media, linewidth=2)
    else:	
    	plt.plot(totalGrowth, p.d[media][strain][stat], color=col, label=media, linewidth=2)
    
    if supress != 1:
    	plt.xlabel('OD')
    	plt.ylabel(stat)
    	plt.legend(loc=4)
	
	
def plotByODAllMedia(p, stat, strain, medialist, mediaColors, normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1):
    plt.figure()
    for i in range(0, np.size(medialist)):
        plotByOD(p, stat, strain, medialist[i], mediaColors[i], normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1)
        plt.xlim([0.2, 1.5])
        plt.title(stat+' of '+strain)
        plt.ylabel(stat)
        plt.xlabel('OD')
        plt.legend(loc=2)
        
        
def statArea(p,media, strain, stat):##make sure it is a time series
    return np.trapz(p.d[media][strain][stat])
    
def plotByODAllMediaByStrain(p, stat, strains, medialist, mediaColors, normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1):
    for i in range(0, np.size(strains)):
        plotByODAllMedia(p, stat, strains[i], medialist, mediaColors, normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1)


def statByConcentrationTimeSlice(p, stat, strain, concentrations, time, color=0, supress=1, lw=2,labelStrain=1):
    noNaNKeys=[list(concentrations.keys())[y] for y in np.where([np.isnan(x)==False for x in concentrations.values()])[0]]  
    concentrations=dict(zip(noNaNKeys,[concentrations[x] for x in noNaNKeys])) ### assemble a new concentrations object without nans
    ###now order the concentrations from lowest to highest
    print('unordered keys:', list(concentrations.keys()))
    print('unordered values:', list(concentrations.values()))
    newOrder =[np.where(list(concentrations.values())==x)[0][0] for x in np.sort(list(concentrations.values()))]
    print('new order should be', newOrder)
    media= [list(concentrations.keys())[x] for x in newOrder] ###ordered media
    print('ordered keys:', media)
    conc= [list(concentrations.values())[x] for x in newOrder] ###ordered concentrations
    print('ordered values:', conc)
    ##time is in hours. it gets converted into the index of the array
    timeIndex=np.round(np.where(abs(p.t-time)==min(abs(p.t-time))))
    ##print((timeIndex)
    stats=np.zeros(np.size(media))
    errors=np.zeros(np.size(media))
    for i in range(0, (np.size(media))):
        print(p.d[media[i]][strain][stat][timeIndex])
        stats[i]= p.d[media[i]][strain][stat][timeIndex].flatten().mean()
        errors[i]=np.sqrt(p.d[media[i]][strain][stat+'var'][timeIndex])
        if labelStrain==1:
            lab=strain
        else:
            lab=str(time)+'h'
    if color==0:
        color=randomColor()
    plt.errorbar(conc, stats, yerr= errors, color=color, label=lab , linewidth=lw)
    if supress != 1:
        plt.title('Dose response of '+strain+' at '+ str(time)+ ' hours' )
        plt.xlabel('Concentration')
        plt.ylabel(stat)
        plt.legend(loc=2)
    return stats


def statByConcentrationManyTimeSlices(p, stat, strain, conc, times, cols=0, lw=2):
    if cols==0:
        cols=randomColors(np.size(times))
    plt.figure()
    plt.title('Dose response of '+strain+' at several times' )
    plt.xlabel('Concentration')
    plt.ylabel(stat)
    plt.legend(loc=2)
    for i in range(0, np.size(times)):
        statByConcentrationTimeSlice(p, stat, strain, conc, times[i], cols[i], lw=lw, labelStrain=0)
        
def statByConcentrationManyStrains(p, stat, strains, conc, time, cols, supress=1, lw=2):
    plt.figure()
    plt.title('Dose response of at '+ str(time)+ ' hours' )
    plt.xlabel('Concentration')
    plt.ylabel(stat)
    for i in range(0, np.size(strains)):
        statByConcentrationTimeSlice(p, stat, strains[i], conc, time, cols[i], supress=1, lw=lw, labelStrain=1)
    plt.legend(loc=2)


def statCorrelation(p, media, strain, stat1, stat2, degree=0):
    '''this function tries to explain one stat in relation to another. there should be the same points
        for both stats, and often these will be raw data. 
        degree is whether we are correlating the normaldata or a particular derivative of it.
    '''
    plt.scatter(p.d[media][strain][stat1],p.d[media][strain][stat2])
    plt.xlabel(stat1)
    plt.ylabel(stat2)
    reg=scipy.stats.linregress(p.d[media][strain][stat1].flatten(),p.d[media][strain][stat2].flatten())
    ln=np.linspace(0, np.max(p.d[media][strain][stat1]), 300)
    ln2=np.linspace(0, np.max(p.d[media][strain][stat1]), 300)*reg.slope+reg.intercept
    plt.plot(ln, ln2, color='red')
    plt.legend([ 'y='+str(reg.slope)+'x'+stringSign(reg.intercept)+str(reg.intercept),strain+' in '+media])
    return reg

def stringSign(x):
    if sign(x)==1:
        return '+'
    else: 
        return ''

### REPORT FUNCTIONS. these functions are used to plot experiment reports, usually plotting a particular variable for all strains or for one strains in all media.
########concentration report functions: these functions assume there are many concentrations of the same variable across conditions, and thus can be ordered. 
######## they also assume all combinations of media and concentrations exist, but this could improve in the future. the concentrations argument usually corresponts to a magnitude scalar per medium in medialist, 
######## such that each medium is associated with each scalar. 
#These functions generally include mediacolors and straincolors as optional inputs. these are ordered lists of the colors used, each color corresponding to each  strain or media in the list of the function.
### it will soon be upgraded to a dictionary as in the robust functions.
### KEY ADVANTAGE TO USING THIS: knowledge about the concentrations of stuff allows you to generate dose response curves and dependency of stats on concentration. You can also use different timepoints as reference. 
# if no colors are provided, a random set of colors will be generated. 
# the functions also output the specific color dictionary used for plotting, such that it can be input to other functions in the future (for example, if you liked a color combination a lot'. 

def plotRawFLPerStrain(p, concentrations, strainColors, mediaColors, strains, medialist):
	for x in range(0,np.size(strains)):
		plt.figure()
		plt.title('Raw FL of '+strains[x]+' in all media')
		plt.ylabel('Raw FL')
		plt.xlabel('Time (Hrs)')
		for y in range(0, np.size(medialist)):
			arr=p.d[medialist[y]][strains[x]]['odata'][:,range(1,np.size(p.d[medialist[y]][strains[x]]['odata'],1),3)]
			plt.plot(p.t, arr, color=mediaColors[y])
			
			
def plotStatOverTimePerStrain(p,stat,concentrations,strainColors,mediaColors,strains,medialist, plotvar=0):
	for x in range(0,np.size(strains)):
# 		#print( x
# 		#print( strains[x]
# 		#print( medialist
		plt.figure()
		plt.title(stat+' of '+strains[x])
		plt.ylabel(stat)
		plt.xlabel('Time (Hrs)')
		##print( mediaColors
		for y in range(0, np.size(medialist)):
# 		    #print( "x= %d, y= %d", x,y
# 		    ##print( "strains[x]",strains[x] , "medialist[y]", medialist[y], "mediaColors[y]", mediaColors[y]
		    arr=p.d[medialist[y]][strains[x]][stat]
		    arrvar=p.d[medialist[y]][strains[x]][stat+'var']
		    #if size(arr)>1:
		    plt.plot(p.t, arr, color= mediaColors[y], label=medialist[y], linewidth=2)
		    if plotvar!=0:
		        plt.plot(p.t, arr+np.sqrt(arrvar), color= mediaColors[y], label=medialist[y])
		        plt.plot(p.t, arr-np.sqrt(arrvar), color= mediaColors[y], label=medialist[y])
		    plt.legend(loc=1)



def funcByMaxOD(p, func, strain, media, col, normalizeToExt=0, subtractBg=0, supressPlot=1, supress=1):
	value=func(p, strain, media)
	mOD=p.d[media][strain]['odstats']['max OD']
	if supressPlot != 1:
	    if normalizeToExt != 0:
		    plt.plot(p.d[media][strain]['odstats']['max OD'], value/normalizetoExt, 'o',color=col, label=strain)
	    elif subtractBg!=0:		
		    plt.plot(p.d[media][strain]['odstats']['max OD'], value-subtractBg, 'o',color=col, label=strain)
	    else:	
		    plt.plot(p.d[media][strain]['odstats']['max OD'], value, 'o', color=col, label=strain)
	
	    if supress != 1:
		    plt.xlabel('OD')
		    plt.ylabel(func.__name__)
		    plt.legend(loc=4)
	
	valuepair= [mOD, value]	
	return valuepair
	
def funcByMaxODAllMedia(p, func, strain, medialist, col, normalizeToExt=0, subtractBg=0, supress=1):
	values=np.zeros((np.size(medialist),2))
	for i in range(0, (np.size(medialist))):
	    if i==np.size(medialist):
	        sup=0
	    else:
	        sup=1
	    va=funcByMaxOD(p, func, strain, medialist[i], col, normalizeToExt, subtractBg, supress=sup)
	    values[i,:]=funcByMaxOD(p, func, strain, medialist[i], col, normalizeToExt, subtractBg, supress=sup)
	return values

def funcByMaxODAllMediaByStrain(p, func, strains, medialist, strainColors, markerlist=['o', 'x', 's', '^', '+', 'v', 'D','<', '*'], normalizeToExt=0, subtractBg=0, supp=0):
    plt.figure()
    for i in range(0, np.size(strains)):
        plt.title(strains[i]+': '+func.__name__+' explained by final OD')
        plt.ylabel(func.__name__)
        plt.xlabel('Max OD')
        values=funcByMaxODAllMedia(p, func, strains[i], medialist, strainColors[i], normalizeToExt=0, subtractBg=0, supress=supp)
        ##print( values
        plt.plot(values[:,0], values[:,1], markerfacecolor=strainColors[i], label=strains[i], color=strainColors[i], linewidth=2, markersize=10, markeredgewidth=2, marker=markerlist[i], markeredgecolor='black')
        plt.legend(loc=2) 


def plotRawODPerStrain(p, concentrations, strains, medialist, strainColors=0, mediaColors=0):
    if mediaColors==0:
        mediaColors= randomColors(np.size(np.unique(cl.values[:,0])))
        mediaColors= dict(zip(np.unique(cl.values[:, 0]), mediaColors))
    for x in range(0,np.size(strains)):
        plt.figure()
        plt.ylim([0.2, 1.5])
        for y in range(0, np.size(medialist)):
            arr=p.d[medialist[y]][strains[x]]['OD']
            # #print((arr)
            plt.plot(p.t, arr, color=mediaColors[y]) 
        plt.title('Raw OD of '+strains[x]+' in all media')
        plt.ylabel('Raw OD')
        plt.xlabel('Time (Hrs)')

def plateReaderFluorescenceReport(p, concentrations, concvalues, plotvar=0, strainColors=0, mediaColors=0, plotByMaxOD=1, doseResponseTimes=[3]):
    #obtain a condition list of the experiment
    cl=conditionList(p, excludeNull=True).values
    strains=p.allstrains
    #assign colors to every media and every strain anyway
    if mediaColors==0:
        mediaColors= randomColors(np.size(np.unique(cl.values[:,0])))
        mediaColors= dict(zip(np.unique(cl.values[:, 0]), mediaColors))
    if strainColors==0:
        strainColors=randomColors(np.size(strains))  ###code also assumes that all concentrations are present for all strains.
        strainColors= dict(zip(strains,strainColors))
    markers= ['>', '^', 'o', 'D', '*', 's', '+', 'v', 'x',]
    ##print( strains
    medialist= p.d.keys()
    if mediaColors==0:
        mediaColors=randomColors(np.size(p.d.keys()))
    ##print( mediaColors
    if strainColors==0:
        strainColors=randomColors(np.size(p.d[concentrations[0]].keys()))  ###code also assumes that all concentrations are present for all strains.
    plotRawODPerStrain(p, concentrations, strainColors, mediaColors, strains, medialist)
    plotRawFLPerStrain(p, concentrations, strainColors, mediaColors, strains, medialist)
    #now we get rid of the WT strain and the null in order to compute fun things with fluorescence
    inds=plt.mlab.find(np.array([i != 'WT' for i in p.d[concentrations[0]].keys()]) - np.array([i != 'null' for i in p.d[concentrations[0]].keys()]) != True)
    p.d[concentrations[0]].keys()
    ##print( inds
    strainsFL= np.array(strains)[inds]
    strainColorsFL=np.array(strainColors)[inds]
    plotStatOverTimePerStrain(p, normalflperodName,concentrations, strainColorsFL, mediaColors, strainsFL, medialist, plotvar)
    #def plotStatOverTimePerStrain(p,stat,concentrations,strainColors,mediaColors,strains,medialist, plotvar=0)
    if plotByMaxOD==1:
        funcByMaxODAllMediaByStrain(p, maxExpressionLevel, strainsFL, concentrations, strainColorsFL, normalizeToExt=0, subtractBg=0)
        funcByMaxODAllMediaByStrain(p, maxGrowthRate, strainsFL, concentrations, strainColorsFL, normalizeToExt=0, subtractBg=0,)
        funcByMaxODAllMediaByStrain(p, maxExpressionLevelTime, strainsFL, concentrations, strainColorsFL, normalizeToExt=0, subtractBg=0)
        funcByMaxODAllMediaByStrain(p, integrateExpression, strainsFL, concentrations, strainColorsFL, normalizeToExt=0, subtractBg=0)
        #funcByMaxODAllMediaByStrain(p, expressionDerivative, strainsFL, concentrations, strainColorsFL, normalizeToExt=0, subtractBg=0)
        funcByMaxODAllMediaByStrain(p, lagTime, strainsFL, concentrations, strainColorsFL, normalizeToExt=0, subtractBg=0)
        funcByMaxODAllMediaByStrain(p, expLevelAtMaxGrowthRateTime, strainsFL, concentrations, strainColorsFL, normalizeToExt=0, subtractBg=0)
        plotByODAllMediaByStrain(p, normalflperodName, strainsFL, medialist, mediaColors, normalizeToSelf=0, normalizeToExt=0, subtractBg=0, supress=1)
    for j in doseResponseTimes:
        statByConcentrationManyStrains(p, normalflperodName, strainsFL, concvalues, concentrations, j, strainColorsFL)
    funcByConcentrationByStrain(p, maxGrowthRate, strainsFL, concvalues, concentrations, strainColorsFL)
    funcByConcentrationByStrain(p, lagTime, strainsFL, concvalues, concentrations, strainColorsFL)
    funcByConcentrationByStrain(p, maxExpressionLevelTime, strainsFL, concvalues, concentrations, strainColorsFL)
    funcByConcentrationByStrain(p, integrateExpression, strainsFL, concvalues, concentrations, strainColorsFL)
    funcByConcentrationByStrain(p, expLevelAtMaxGrowthRateTime, strainsFL, concvalues, concentrations, strainColorsFL)
    dic={"mediaColors":mediaColors, "strainColors":strainColors  }
    return dic

##### ROBUST REPORT FUNCTIONS: these functions will just plot everything that can be plotted, without assuming anything about the media, strain, kind of medium, or existing combinations.
##### the only mild assumption that could be made is that wt and null do not contain normalflperodName. 
#These functions generally include mediacolors and straincolors as optional inputs. these are dictionaries of the form d={'media1': 'hexcolor1', 'media2': 'hexcolor2'}, d={'strain1': 'hexcolor1', 'strain2': 'hexcolor2'}.
# if no colors are provided, a random set of colors will be generated. 
# the functions also output the specific color dictionary used for plotting, such that it can be input to other functions in the future (for example, if you liked a color combination a lot'. 


# def plotRawConditionRobust(p, media, strain):
# 
# plt.plot(p.t, arr, color= mediaColors[sdf['media'][x]], label=sdf['media'][x])
    
def plotStatOverTimePerStrainRobust(p,stat,mediaColors=0, plotvar=0, xlim=False, ylim=False):
    #mediaColors must be a dictionary where thee field names
    # are the media type and the values are colors. Starting from a media and a color vector, do dict(zip(media,color)).
    allarrs={}
    cl=conditionList(p, excludeNull=True)
    if mediaColors==0:
        mediaColors= randomColors(np.size(np.unique(cl.values[:,0])))
        mediaColors= dict(zip(np.unique(cl.values[:, 0]), mediaColors))
    strains= np.unique(cl.values[:, 1])
    
    for x in range(0,np.size(strains)):
        allarrs[strains[x]]=[]
        if containsFL(p, cl.values[x, 0], cl.values[x,1]):
            print(strains[x])
            plt.figure()
            plt.title(stat+' of '+strains[x])
            plt.ylabel(stat)
            plt.xlabel('Time (Hrs)')
            subconditions= cl[cl['strain']==strains[x]].values
            print('subconditions: ', subconditions)
            for y in range(0, np.size(subconditions[:,0])):
                arr=p.d[subconditions[y,0]][subconditions[y,1]][stat]
                arrvar=p.d[subconditions[y,0]][subconditions[y,1]][stat+'var']
                #if size(arr)>1:
                plt.plot(p.t, arr, color= mediaColors[subconditions[y,0]], label=subconditions[y,0], linewidth=2)
                allarrs[strains[x]].append(arr)
                if xlim !=False:
                    plt.xlim(xlim)
                if ylim != False:
                    plt.ylim(ylim)
                if plotvar!=0:
                    plt.plot(p.t, arr+np.sqrt(arrvar), color= mediaColors[y], label=medialist[y])
                    plt.plot(p.t, arr-np.sqrt(arrvar), color= mediaColors[y], label=medialist[y])
        else:
            continue
    return allarrs

def plotStatOverAlignedTimePerStrainRobust(p,stat,mediaColors=0, plotvar=0, alignBy='GrowthInflectionPoint', normalize=True, xlim=False, ylim=False):
    #mediaColors must be a dictionary where thee field names
    # are the media type and the values are colors. Starting from a media and a color vector, do dict(zip(media,color)).
    cl=conditionList(p, excludeNull=True)
    if mediaColors==0:
        mediaColors= randomColors(np.size(np.unique(cl.values[:,0])))
        mediaColors= dict(zip(np.unique(cl.values[:, 0]), mediaColors))
    strains= np.unique(cl.values[:, 1])
    for x in range(0,np.size(strains)):
        if stat==normalflperodName and (strains[x]=='null'):
            continue
        plt.figure()
        plt.title(stat+' of '+strains[x])
        plt.ylabel(stat)
        plt.xlabel('Time (Hrs)')
        subconditions= cl[cl['strain']==strains[x]].values
        print('subconditions: ', subconditions)
        for y in range(0, np.size(subconditions[:,0])):
            arr=p.d[subconditions[y,0]][subconditions[y,1]][stat]
            arrvar=p.d[subconditions[y,0]][subconditions[y,1]][stat+'var']
            #if size(arr)>1:
            plt.plot(alignByGrinfp(p, subconditions[y,0], subconditions[y,1]), (arr-mean(arr))/std(arr), color= mediaColors[subconditions[y,0]], label=subconditions[y,0], linewidth=2)
            if xlim !=False:
                plt.xlim(xlim)
            if ylim != False:
                plt.ylim(ylim)
            if plotvar!=0:
                plt.plot(p.t, arr+np.sqrt(arrvar), color= mediaColors[y], label=medialist[y])
                plt.plot(p.t, arr-np.sqrt(arrvar), color= mediaColors[y], label=medialist[y])
                
def plotByODAllMediaByStrainRobust(p,stat,mediaColors=0, plotvar=0,xlim= False, ylim=False):
    #mediaColors must be a dictionary where thee field names
    # are the media type and the values are colors. Starting from a media and a color vector, do dict(zip(media,color)).
    cl=conditionList(p)
    if mediaColors==0:
        mediaColors= randomColors(np.size(np.unique(cl.values[:,0])))
        mediaColors= dict(zip(np.unique(cl.values[:, 0]), mediaColors))
    strains= np.unique(cl.values[:, 1])
    
    for x in range(0,np.size(strains)):
        if stat==normalflperodName and (strains[x]==p.refstrain or strains[x]=='null'):
            continue
        plt.figure()
        plt.title(stat+' of '+strains[x])
        plt.ylabel(stat)
        plt.xlabel('OD595')
        subconditions= cl[cl['strain']==strains[x]].values
        for y in range(0, np.size(subconditions[:,0])):
            arr=p.d[subconditions[y,0]][subconditions[y,1]][stat]
            arrvar=p.d[subconditions[y,0]][subconditions[y,1]][stat+'var']
            #if size(arr)>1:
            plt.plot(p.d[subconditions[y,0]][subconditions[y,1]]['OD'].mean(1), arr, color= mediaColors[subconditions[y,0]], label=subconditions[y,0], linewidth=2)
            if plotvar!=0:
                plt.plot(p.d[subconditions[y,0]][subconditions[y,1]]['OD'].mean(1), arr+np.sqrt(arrvar), color= mediaColors[subconditions[y,0]], label=subconditions[y,0])
                plt.plot(p.d[subconditions[y,0]][subconditions[y,1]]['OD'].mean(1), arr-np.sqrt(arrvar), color= mediaColors[subconditions[y,0]], label=subconditions[y,0])
            if xlim !=False:
                plt.xlim(xlim)
            if ylim != False:
                plt.ylim(ylim)
            #plt.legend(loc=2)
            
            
def plotRawODPerStrainRobust(p, mediaColors=0, containsFL=True, xlim= False, ylim=False, prFormat=2):
    cl= conditionList(p, excludeNull=True)
    #print( cl)
    if mediaColors==0:
        mediaColors= randomColors(np.size(np.unique(cl['media'])))
        mediaColors= dict(zip(np.unique(cl['media']), mediaColors))
        #print( mediaColors)
    for  str in p.allstrains:
        plt.figure();
        sdf= cl[cl['strain']==str]
        sdf=sdf.reset_index()
        legendNames=[]
        patches=[]
        for x in range(0, np.size(sdf,0)):
            arr=p.d[sdf['media'][x]][sdf['strain'][x]]['OD']
            plt.plot(p.t, arr, color= mediaColors[sdf['media'][x]], label=sdf['media'][x])
            legendNames.append(sdf['media'][x])
            patches.append(pch.Patch(color=mediaColors[sdf['media'][x]], label=sdf['media'][x]))
            plt.title('Raw OD of '+sdf['strain'][x]+' in all media')
            plt.ylabel('Raw OD')
            plt.xlabel('Time (Hrs)')
        if xlim !=False:
            plt.xlim(xlim)
        if ylim != False:
            plt.ylim(ylim)
        plt.figlegend(patches, legendNames, 'upper right')
    return mediaColors
    
    
    
def plotRawStatPerStrainRobust(p, mediaColors=0, dtype='GFP', xlim=False, ylim=False, newfig=1):
    cl= conditionList(p, excludeNull=True)
    if mediaColors==0:
        mediaColors= randomColors(np.size(np.unique(cl['media'])))
        mediaColors= dict(zip(np.unique(cl['media']), mediaColors))
        #print( mediaColors)
    for  str in p.allstrains:
        plt.figure();
        sdf= cl[cl['strain']==str]
        sdf=sdf.reset_index()
        #print( sdf)
        legendNames=[];
        patches=[];
        for x in range(0, np.size(sdf,0)):
            legendNames.append(sdf['media'][x])
            patches.append(pch.Patch(color=mediaColors[sdf['media'][x]], label=sdf['media'][x]))
            arr=p.d[sdf['media'][x]][sdf['strain'][x]][dtype]
            if np.size(np.shape(arr))==1:
                plt.plot(p.t, arr, color= mediaColors[sdf['media'][x]], label=sdf['media'][x])
            else:
                plt.plot(p.t, arr, color= mediaColors[sdf['media'][x]], label=sdf['media'][x])
            if xlim !=False:
                plt.xlim(xlim)
            if ylim != False:
                plt.ylim(ylim)
            plt.title('Raw '+dtype+' of '+sdf['strain'][x]+' in all media')
            plt.ylabel('Raw '+dtype)
            plt.xlabel('Time (Hrs)')
        plt.figlegend(patches, legendNames, 'upper right')
    return mediaColors
    
def plotRawStatPerMediaRobust(p, strainColors=0, dtype='GFP', xlim=False, ylim=False, normalize=False, ignoreStrains=[]):
    cl= conditionList(p, excludeNull=True)
    if strainColors==0:
        strainColors= randomColors(np.size(np.unique(cl['strain'])))
        strainColors= dict(zip(np.unique(cl['strain']), strainColors))
        #print( strainColors)
    for  str in list(p.d.keys()):
        plt.figure();
        sdf= cl[cl['media']==str]
        sdf=sdf.reset_index()
        #print( sdf)
        legendNames=[];
        patches=[];
        for x in range(0, np.size(sdf,0)):
            if sdf['strain'][x] in ignoreStrains:
                continue
            arr=p.d[sdf['media'][x]][sdf['strain'][x]][dtype]
            #print(sdf['media'][x])
            legendNames.append(sdf['strain'][x])
            patches.append(pch.Patch(color=strainColors[sdf['strain'][x]], label=sdf['strain'][x]))
            if np.size(np.shape(arr))==1:
                if normalize==True:
                    arr=(arr-np.min(arr))/np.max(arr-np.min(arr))
                plt.plot(p.t, arr, color= strainColors[sdf['strain'][x]], label=sdf['media'][x])
            else:
                if normalize==True:
                    for j in range(1,np.size(arr,1)):
                        arr[:, j]=(arr[:,j]-np.min(arr[:,j]))/np.max(arr[:,j]-np.min(arr[:,j]))
                plt.plot(p.t, arr, color= strainColors[sdf['strain'][x]], label=sdf['media'][x])
            if xlim !=False:
                plt.xlim(xlim)
            if ylim != False:
                plt.ylim(ylim)
            plt.title('Raw '+dtype+' of '+ 'all strains in '+sdf['media'][x])
            plt.ylabel('Raw '+dtype)
            plt.xlabel('Time (Hrs)')
        plt.figlegend(patches, legendNames, 'upper right')
    return strainColors
    
def plateReaderFluorescenceReportRobust(p, FL='GFP', plotvar=0, strainColors=0, mediaColors=0, AutoFL=True, tlim=False, odlim=False, FLlim=False, concentrations=0):
    markers= ['>', '^', 'o', 'D', '*', 's', '+', 'v', 'x']
    cl=conditionList(p, excludeNull=True).values
    medialist=np.unique(cl[:,0]) ##these strains preserve are listed by condition which means that media may appear more than once. same for strains.
    strains=np.unique(cl[:,1])
    if mediaColors==0:
        mediaColors=randomColors(np.size(medialist))
        mediaColors= dict(zip(medialist,mediaColors))
    if strainColors==0:
        strainColors=randomColors(np.size(strains))  ###code also assumes that all concentrations are present for all strains.
        strainColors= dict(zip(strains,strainColors))
    experimentOverview(p, dtype= 'OD', colormap='cool', colorMapRange=odlim, timeRange=False)
    experimentOverview(p, dtype= FL, colormap='cool', addFL={FL:FLlim}, timeRange=False)
    plotRawStatPerStrainRobust(p, mediaColors, dtype='OD', xlim=tlim, ylim=odlim)
    plotRawStatPerStrainRobust(p, mediaColors, dtype=FL, xlim=tlim, ylim=FLlim)
    plotRawStatPerStrainRobust(p, mediaColors=mediaColors, dtype='FLperod', xlim=tlim, ylim=FLlim)
    plotRawStatPerMediaRobust(p, strainColors=strainColors, dtype='FLperod', xlim=tlim, ylim=FLlim)
    if(AutoFL==True):
        plotRawStatPerStrainRobust(p, mediaColors, dtype='AutoFL', xlim=tlim, ylim=FLlim)
    #plotStatOverTimePerStrainRobust(p,normalflperodName,mediaColors, xlim=tlim, ylim=FLlim)
    #plotStatOverTimePerStrainRobust(p,'gr',mediaColors)
    plotByODAllMediaByStrainRobust(p,normalflperodName,mediaColors, plotvar=plotvar, xlim=odlim, ylim=FLlim) ###as this is OD on the x axis
    dic={"mediaColors":mediaColors, "strainColors":strainColors  }
    ###concentrations is expected to be a dictionary where all media are associated with a numeric value. if a particular media is not of interest, then the value should be nan.
    #if concentrations != 0: ### if concentrations is not 0, then do the following
        ###retrieve all media whose concentration is not null
        ### order the concentrations from lowest to highest
        ###produce relevant stat by concentration plots
    return dic

######### MISCELLANEUOUS SMALL FUNCTIONS: these functions are needed throughout and are heavily used by all the above functions. They also have value and usefulness on their own. 


def hasKey(obj, k):
    if k in obj.keys():
        return True
    else:
        return False

def containsCondition(p, media, strain):
     if hasKey(p.d, media) and hasKey(p.d[media], strain):
        return 1
     else:
        return 0

def getAllStrains(p):
    strains=[]
    for straindict in p.d.values():
        strains.append(list(straindict.keys()))
    strains= list(np.unique(np.concatenate(strains)))
    return strains

def conditionList(p, excludeNull=False):
    finalmedia=[]
    finalstrains=[]
    for media in p.d.keys():
        for strain in p.d[media].keys():
            #print('media:', media)
            #print('strains:', strain)
            finalmedia.append(media)
            finalstrains.append(strain)
            cl=pd.DataFrame(list(zip(finalmedia, finalstrains)), columns= ['media', 'strain'])
    if excludeNull==True:
        cl=cl[cl['strain']!='null']
    return cl
        


def savePDF(fig, filename='tempfig'):
    import matplotlib.backends.backend_pdf
    ## saves only existing figures regardless of their number or order
    pdf = matplotlib.backends.backend_pdf.PdfPages(filename+'.pdf')
    pdf.savefig( fig )
    pdf.close()

def savePDFs(filename='tempfig'):
    import matplotlib.backends.backend_pdf
    for fig in plt.get_fignums(): ## saves only existing figures regardless of their number or order
        pdf = matplotlib.backends.backend_pdf.PdfPages(filename+str(fig)+'.pdf')
        pdf.savefig( fig )
        pdf.close()
        
def saveEPSs(filename='tempfig'):
    import matplotlib.backends.backend_pdf
    for fig in plt.get_fignums(): ## saves only existing figures regardless of their number or order
        plt.figure(fig)
        plt.savefig(filename+str(fig)+'.eps', format='eps', dpi=1000)

def savePNGs(filename='tempfig'):
    import matplotlib.backends.backend_pdf
    for fig in plt.get_fignums(): ## saves only existing figures regardless of their number or order
        plt.figure(fig)
        plt.savefig(filename+str(fig)+'.png', format='png', dpi=1000)
        
def saveJPGs(filename='tempfig'):
    import matplotlib.backends.backend_pdf
    for fig in plt.get_fignums(): ## saves only existing figures regardless of their number or order
        plt.figure(fig)
        plt.savefig(filename+str(fig)+'.jpg', format='jpg', dpi=1000)


def saveAllPDF(filename='output.pdf'):
    import matplotlib.backends.backend_pdf
    pdf = matplotlib.backends.backend_pdf.PdfPages(filename)
    for fig in plt.get_fignums(): ## saves only existing figures regardless of their number or order
        pdf.savefig( fig )
    pdf.close()


def randomColor():
    hex_digits = np.array(['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'])
    digit_array=(hex_digits[np.random.randint(0,15,6)])
    joined_digits = ''.join(digit_array)
    color = '#' + joined_digits
    return color
	
def randomColors(num):
	colorVector=[]
	for i in range(0,num):
		colorVector.append(randomColor())
	return colorVector
	
	
def sign2color(values, positive='red', negative='cyan'):

    colors= [];
    for j in range(0, np.size(values)):
        if np.sign(values[j])==1:
            colors.append(positive)
        else:
            colors.append(negative)
    
    return(colors)
    
def colorArray(arr, cmap='cool'):
    c=matplotlib.cm.get_cmap(cmap)
    return [c(x) for x in arr]
    
    
def plateCoordinates(plateloc):
    mp= {'A':0, 'B':1, 'C':2, 'D':3, 'E':4, 'F':5,'G':6, 'H':7}
    s=''
    if np.size(list(plateloc))>2:
        coordinates=[mp[list(plateloc)[0]], int(s.join((list(plateloc)[1], list(plateloc)[2])))-1] 
    else:
        coordinates=[mp[list(plateloc)[0]], int(list(plateloc)[1])-1]
    return coordinates
    

def scat3d(x,y,z, new=1, markerNum=1):
    markers= ['o', '^', 's', '+']
    fig=plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c='blue', marker='o')
    
def statTSFrames(p, media, strain, stat , finaltp, col='black', tpPerFrame=2, labelStrain=1):
    for j in range(tpPerFrame,finaltp,tpPerFrame):
        plt.figure();
        plt.scatter(p.t[0:j], p.d[media][strain][stat][0:j], color=col, marker='s', label=strain)
        plt.xlim([p.t[0],p.t[finaltp]])
        plt.ylim([0, max(p.d[media][strain][stat][0:j])])

def experimentOverview(p, dtype= 'OD', colormap='cool', colorMapRange=False, timeRange=False, addFL=False):
    defaultRange={'OD': [0,1.5], 'GFP':[0,8000], 'AutoFL': [0,1000], 'mCherry': [0, 5000]}
    if addFL!= False:
        defaultRange.update(addFL)
    if colorMapRange==False:
        colorMapRange= defaultRange[dtype]
    legends=[colorMapRange[1], '', '', '', colorMapRange[0]]
    f, axarr = plt.subplots(8, 12)
    axarr[0,6].set_title(dtype+' in all wells')
    refaxis = f.add_axes([0.07, .1, 0.05, 0.7])
    gradient = np.linspace(colorMapRange[0], colorMapRange[1], 256)
    gradient = np.flipud(np.vstack((gradient, gradient)).T)
    refaxis.imshow(gradient, aspect='auto', cmap=colormap)
    refaxis.get_xaxis().set_visible(False)
    refaxis.get_yaxis().set_ticks( np.linspace(256,0, 6) )
    refaxis.set_yticklabels(np.linspace(colorMapRange[0], colorMapRange[1], 6))
    #plt.figlegend(patches, legends, 'upper left')
    cl=conditionList(p)
    problematicWells=[]
    for x in range(0, np.size(cl,0)):
        media=cl.values[x,0]
        strain=cl.values[x,1]
        conditionCoordinates=[plateCoordinates(a) for a in p.d[media][strain]['plateloc']]
        print('plateloc: ', p.d[media][strain]['plateloc'])
        print('conditionCoordinates: ', conditionCoordinates)
        for cc in range(0,np.size(conditionCoordinates,0)):
            print('cc is ', cc )
            if p.d[media][strain]['plateloc'][cc] in p.ignoredwells:
                print('well ', p.d[media][strain]['plateloc'][cc], ' is in ignoredwells')
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].plot([1,2,3,4,5,6], [1,2,3,4,5,6], c='red')
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].plot([1,2,3,4,5,6], [6,5,4,3,2,1], c='red')
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].get_xaxis().set_visible(False)
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].get_yaxis().set_visible(False)
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].legend()
                continue
            try:
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].scatter(p.t, p.d[media][strain][dtype][:,cc], c=p.d[media][strain][dtype][:,cc], cmap=colormap, vmin=colorMapRange[0], vmax=colorMapRange[1], edgecolor='None')
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].get_xaxis().set_visible(False)
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].get_yaxis().set_visible(False)
                axarr[conditionCoordinates[cc][0], conditionCoordinates[cc][1]].legend()
            except IndexError:
                #print('time length: ', np.size(p.t))
                #print('vector length: ', np.size(p.d[media][strain][dtype],0))
                print('Index error:  well ', p.d[media][strain]['plateloc'][cc], ', media ', media, ', strain ', strain)

def plateLocMap(p):
    strains=[]
    media= list(p.d.keys())
    for m in media:
        strains.append(list(p.d[m].keys()))
    strains= np.unique(strains)
    finalDict=dict()
    for i in range(0, np.size(media)):
        for j in range(0,np.size(strains)):
            if containsCondition(p, media[i], strains[j]):
                for pl in p.d[media[i]][strains[j]]['plateloc']:
                    finalDict[pl]={'media': media[i], 'strain': strains[j]}
    return finalDict
        

#######This section is work in progress. Please ignore
# def initialFinalFL(p, strain, media, finalTime):
# def initialMaxFL(p, strain, media)
# def lagFinalFL(p, strain, media, finalTime)
# def lagMaxFL(p, strain, media)

# def extractConcentrations(medialist):
#     out=np.zeros(array(size(medialist)))
#     #print( out
#     for i in range(0, size(concentrations)): 
#         out[i]= float(medialist[i].split(" ")[1].split("%")[0])
#     return out
# 
# def replicateAnalysis(platereaders, strain, media)
#     cols= randomColors(3)
#     for i in range(0, size(platereaders)):
#         plotalign(platereaders{1}, strain, media, col=cols[i])
#         iFFL[i, :]= initialFinalFL(platereaders{i}, strain, media, finalTime)
#         iMFL[i,:]=initialMaxFL(platereaders{i}, strain, media)
#         lagFFL[i,:]=lagFinalFL(platereaders{i}, strain, media, finalTime)
#         lagMFL[i,:]=lagMaxFL(platereaders{i}, strain, media)
        

def statByGrowthStage(p, stat, threshold=0.01):
# for this to work there already has to be  a calculation of the derivative of the growth rate
    cl=conditionList(p, excludeNull=True)
    cl.reset_index()
    upMean=[]
    downMean=[]
    flatMean=[]
    upSd=[]
    downSd=[]
    flatSd=[]
    for x in range(0, np.size(cl, 0)):
        media= cl.iloc[x].media
        strain=cl.iloc[x].strain
        print('media'+media+'strain '+strain) 
        a=scipy.interpolate.interp1d(p.t, p.d[media][strain]['d/dtgr'])
        interpTime=np.linspace(0, p.t[-1],600)
        interpol= a(interpTime)
        grUpIndices= np.where(interpol>threshold)
        grDownIndices=np.where(interpol< -threshold)
        grFlatIndices=np.where(abs(interpol)<=threshold)
        b=scipy.interpolate.interp1d(p.t, p.d[media][strain][stat])
        interpolb=b(interpTime)
        upMean.append(np.mean(interpolb[grUpIndices]))
        downMean.append(np.mean(interpolb[grDownIndices]))
        flatMean.append(np.mean(interpolb[grFlatIndices]))
        upSd.append(np.std(interpolb[grUpIndices]))
        downSd.append(np.std(interpolb[grDownIndices]))
        flatSd.append(np.std(interpolb[grFlatIndices]))
    cl[stat+'risingGR']=np.array(upMean)
    cl[stat+'fallingGR']=np.array(downMean)
    cl[stat+'stationaryPhase']=np.array(flatMean)
    return(cl)

#determining the places where  there is a sign change
#zeroCrossingIndices= np.where(np.gradient(sign(a(interpTime)))) 
#in some tests, the gradient retrieves consecutive numbers, which doesn't make much sense
#therefore we remove any consecutive points as probably noise.
#zeroCrossingIndices= nonConsecutive(zeroCrossingIndices)


def nonConsecutive(lst, i=0): # list of potentially consecutive numbers. we toss the next consecutive number
    '''
    nonConsecutive screens a list of numbers, and if there are any consecutive ones, retrieves only
    the first one. Useful when looking sharp boundaries and consecutive numbers are a result of numerical glitches.
    ''' 
    inds= [x for x in range(0, np.size(lst))]
    if i==np.size(lst)-1 or i+1==np.size(lst)-1:
        return(lst)
    else: 
        print(' i is not penultimate\n')
        if lst[i+1]==lst[i]+1:  #if the next index is the next number#return the full array except for the next index (because it is a consecutive)
            print(' consecutive found\n')
            inds.remove(i+1)
            return(nonConsecutive([lst[x] for x in inds], i+1))
        else:
            print(' no consecutive found\n')
            return(nonConsecutive(lst, i+1))






	
	
	
	