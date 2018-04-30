expt=p
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
	rawdata[dt]=pd.DataFrame(data=actualdata, index=welllabels,columns=timelabels).replace('OVER', value=NaN) 
##now that we have got the raw data we can apply functions to the matrices directly
pickle.dump(rawdata, open('rawdatadict.pkl', 'wb'))



stat1='GFP60'
stat2='GFP80'
stat3='AutoFL'
notnans=np.isnan(rawdata[stat1])==False
reg=scipy.stats.linregress(rawdata[stat1].values[notnans].flatten(),rawdata[stat2].values[notnans].flatten())
ln=np.linspace(0, np.nanmax(rawdata[stat1]), 300)
ln2=np.linspace(0, np.nanmax(rawdata[stat1]), 300)*reg.slope+reg.intercept
plot=True
if plot==True:
	plt.scatter(rawdata[stat1],rawdata[stat2])
	plt.xlabel(stat1)
	plt.ylabel(stat2)
	plt.plot(ln, ln2, color='red')
	plt.legend([ 'y='+str(reg.slope)+'x'+stringSign(reg.intercept)+str(reg.intercept)])
#fixing the nans in stat2 by applying a regression from values in stat 1
substitution=rawdata[stat1][np.isnan(rawdata[stat2])]*reg.slope+reg.intercept     
rawdata[stat2][np.isnan(substitution)==False]=substitution[np.isnan(substitution)==False]



#now we standardise the data to be comparable across all gains and all machines. for this we need:
## a zero, which could well be the media
##a one, which will be a strain that is always measured in the same od for the same measurement. for example, every measurement is done for both medium and wt at OD 0.4. 
## hoping that the values are very consistent for that measurement across experiments, then the ratio of both could be a good basis for standardisation

nullindices=[d for d in range(84, 96)]  #find something more elborate in the future
wtindices=[d for d in range(0,11)]
#assuming the differences are minor (noisy) between one media measurement and another and allowing some time differences
#in the measurement, we subtract the minimum media value at each timepoint
rawdata['OD']-rawdata['OD'].iloc[nullindices,: ].min(0)




###
plt.figure();
plt.subplot(131)
plt.scatter(rawdata['OD'].iloc[wtindices,: ],rawdata['GFP80'].iloc[wtindices,: ])
plt.title('GFP80 of WT')
plt.xlabel('OD595')
plt.ylabel('GFP80')
plt.subplot(132)
plt.scatter(rawdata['OD'].iloc[wtindices,: ],rawdata['GFP60'].iloc[wtindices,: ])
plt.title('GFP60 of WT')
plt.xlabel('OD595')
plt.ylabel('GFP60')
plt.subplot(133)
plt.scatter(rawdata['GFP60'].iloc[wtindices,: ],rawdata['GFP80'].iloc[wtindices,: ])
plt.title('GFP60 vs GFP80 of WT')
plt.xlabel('GFP60')
plt.ylabel('GFP80')
plt.xlim([0,100]) 
plt.ylim([0,800])
reg=scipy.stats.linregress(rawdata['GFP60'].iloc[wtindices,: ].values.flatten(),rawdata['GFP80'].iloc[wtindices,: ].values.flatten())
ln=np.linspace(0, np.nanmax(rawdata['GFP60']), 300)
ln2=np.linspace(0, np.nanmax(rawdata['GFP60']), 300)*reg.slope+reg.intercept
plt.plot(ln, ln2, color='red');plt.legend([ 'y='+str(reg.slope)+'x'+stringSign(reg.intercept)+str(reg.intercept)])


##script to do internal data standardisation based on the WT fluorescence increase with respect to the OD
def internalStd(rawdata, stat, wtindices, fr=0.4, to=0.7, len=30):
	vals=np.zeros([np.size(wtindices),len])
	xvals=np.matlib.repmat(np.linspace(fr, to, len), np.size(vals, 0),1)     
	for d in range(0, size(wtindices)):
		func= interp1d(rawdata['OD'].iloc[wtindices[d],: ].values, rawdata[stat].iloc[wtindices[d],: ].values)
		vals[d, :]=func(linspace(fr, to, len))
	reg=scipy.stats.linregress(xvals.flatten(),vals.flatten())
	transformeddata= (rawdata[stat]-reg.intercept)/reg.slope
	return reg, transformeddata
#plt.figure();
#plt.plot(xvals.T, vals.T)
#plt.plot(xvals.T, xvals.T*reg1.slope+reg1.intercept, color='red', linewidth=2)

plt.figure()
c=1
tdata={}
for j in np.linspace(0.5, 1, 5):
	plt.subplot(1,5,c)
	reg1, tdata1=internalStd(rawdata, stat1, wtindices, fr=0.4, to=j)
	reg2, tdata2=internalStd(rawdata, stat2, wtindices, fr=0.4, to=j)
	reg3, tdata3=internalStd(rawdata, stat3, wtindices, fr=0.4, to=j)
	plt.plot((rawdata[stat1].iloc[wtindices,: ].values.T-reg1.intercept)/reg1.slope, color='blue')
	plt.plot((rawdata[stat2].iloc[wtindices,: ].values.T-reg2.intercept)/reg2.slope, color='red')
	plt.plot((rawdata[stat3].iloc[wtindices,: ].values.T-reg3.intercept)/reg3.slope, color='green')
	c+=1
plt.suptitle('Increasing datapoints increases match')


##free code for standardisation
for d in range(0, size(wtindices)):
	print('go')
	func= interp1d(rawdata['OD'].iloc[wtindices[d],: ].values, rawdata[stat2].iloc[wtindices[d],: ].values)
	vals[d, :]=func(linspace(fr, to, len))
reg2=scipy.stats.linregress(xvals.flatten(),vals.flatten())
plt.figure();
plt.plot(xvals.T, vals.T)
plt.plot(xvals.T, xvals.T*reg2.slope+reg2.intercept, color='red', linewidth=2)

plt.figure()	
plt.plot((rawdata[stat1].iloc[wtindices,: ].values.T-reg1.intercept)/reg1.slope, color='blue')
plt.plot((rawdata[stat2].iloc[wtindices,: ].values.T-reg2.intercept)/reg2.slope, color='red')


tdata={}
stat1='GFP60'
stat2='GFP80'
stat3='AutoFL'
tdata['OD']=rawdata['OD']
reg1, tdata[stat1]=internalStd(rawdata, stat1, wtindices, fr=0.4, to=j)
reg2, tdata[stat2]=internalStd(rawdata, stat2, wtindices, fr=0.4, to=j)
reg3, tdata[stat3]=internalStd(rawdata, stat3, wtindices, fr=0.4, to=j)



#### trying to create a new excel file with the new preprocessed data
df2=df #create a copy of the original data
for i in itime:
	d= i+2
	#rawdata[df.ix[i-2][0]]=pd.DataFrame(df.ix[d:d+matrixrows-1, 1:1+np.size(timepoints)], index=df.ix[d:d+matrixrows-1,0], columns=timepoints) #getting rid of incomplete timepoints
	dt=df.ix[i-2][0] ##the datatype
	df2.ix[d:d+matrixrows-1, 1:1+np.size(timepoints)]=tdata[dt].values #measurements
##now that we have got the raw data we can apply functions to the matrices directly
df2.to_excel(expt.name+'preprocessed.xlsx', header=True, index=False)



