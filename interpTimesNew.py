def interpTimesNew(self, replicateMatrix=False, dtype='OD', centeringVariable='time', upperLim=16, exceptionShift=0.01, ignoreExps=False, experiments='all'):    
	'''
	interpTimes(self, media, strain, dtype='FLperod', centeringVariable='time', upperLim=16, exceptionShift=0.01, ignoreExps=False, experiments='all')   
	interpolate all replicates across the same time scale.
	replicateMatrix is a matrix of the format of self.allReplicates.

	'''
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
	experimentList=replicateMatrix['experiment'].values
	adjustedTimes[dtype]=[]
	adjustedTimes[centeringVariable]=[]
	for j in range(0, size(replicateMatrix, 0)): ###retireving the limiting values for interpolation amongst all experiments.
		expt, media, strain=replicateMatrix.values[j, [0,1,2]]
		#we make a list to store all the replicates for dtype
		maxLengths.append(np.around(self.data[expt].d[media][strain][centeringVariable][-1],2))
		startPoints.append(np.around(self.data[expt].d[media][strain][centeringVariable][0],2))
		#startpoints=np.array(startPoints)
		#maxLengths=np.array(maxLengths)
	#print('startpoints: ', startPoints)
	#print('maxlengths: ', maxLengths)
	interpRange=[np.max(np.array(startPoints)), np.min(np.array(maxLengths))]
	#print('Interprange: ', interpRange)
	func= lambda x: fitsBounds(x, interpRange) ### this lambda creates function func which will evaluate whether x falls within the interpRange
	for j in range(0, size(replicateMatrix, 0)):
		expt, media, strain, plateloc=replicateMatrix.values[j, [0,1,2,3]]
		#finding the points that fit the range
		fitPoints=np.where([func(x) for x in self.data[expt].d[media][strain][centeringVariable]])
		#print(np.shape(self.data[expt].d[media][strain][dtype])[1])
		#Index in the data matrix that corresponds to this specific well
		platelocIndex=np.where(self.data[expt].d[media][strain]['plateloc']==plateloc)
		if not(dtype.endswith('mn')) and not(dtype.endswith('gr')) and not(dtype.endswith('perod')) and not(dtype.endswith('var')):
			adjustedTimes[dtype].append(self.data[expt].d[media][strain][dtype][fitPoints, platelocIndex])
		else:
			adjustedTimes[dtype].append(self.data[expt].d[media][strain][dtype][fitPoints])
		adjustedTimes[centeringVariable].append(np.around(self.data[expt].d[media][strain][centeringVariable][fitPoints],2))
	finalDict={};
	finalDict[dtype]=np.empty([np.size(adjustedTimes[centeringVariable][0]), len(adjustedTimes[dtype])], dtype=None)
	finalDict[centeringVariable]=np.around(adjustedTimes[centeringVariable][0],2) #arbitrarily taking the times of the first condition as a reference
	for j in len(adjustedTimes[dtype]): 
		try:
			fint=scint.interp1d(adjustedTimes[centeringVariable][j],adjustedTimes[dtype][j]) #interpolate times j and response j
			finalDict[dtype][:, j]=fint(finalDict[centeringVariable])
		except:
			finalDict[dtype][:, j]=np.empty([size(finalDict[centeringVariable],0),1])*np.nan
	return finalDict
		
		
		
		
		
		
		
		
		
		
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
