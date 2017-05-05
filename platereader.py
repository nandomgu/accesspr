#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import gaussianprocess as gp
from fitderiv import fitderiv
import genutils as gu
import pandas as pd
import datetime
try:
    import seaborn as sns
except ImportError:
    pass
plt.rcParams['figure.max_open_warning']= 0



class platereader:
    '''
    for analyzing platereader data, correcting for autofluorescence, and extracting growth rates

    A typical work flow is:

    import platereader as pr
    p= pr.platereader('GALdata.xls', 'GALcontents.xlsx')
    p.plot('labels')
    p.plot('OD',  plate= True)
    p.plot('OD')
    p.plot(conditionincludes= ['Gal', 'Glu'], strainexcludes= 'HXT7')
    p.correctauto()
    p.plot('c-GFPperod')
    p.savefigs()
    p.getstats('1% Gal', 'GAL2')
    p.getstats('1% Gal', 'GAL2', dtype= 'FLperod')

    See also http://swainlab.bio.ed.ac.uk/software/platereader/platereader.html

    Wells are removed using, for example,

    p.ignorewells(['B1', 'B2'])

    but the analysis routines should then be re-run.

    All wells can be reinstated with:

    p.ignorewells()

    General properties of the data are shown with:

    p.info()
    p.attributes()

    Pandas dataframes can be created to create new plots and to export to .xlsx or .csv file.
    For example,

    p.makedateframe(['OD', 'GFP'], strains= ['WT'])
    p.exportdataframe(ftype= 'csv')
    p.rmdataframe()

    Bounds on the hyperparameters for all Gaussian processes can be specified. For example,

    p.correctauto(bd= {1: (-6, -4)})

    changes the bounds on the second hyperparameter for the Gaussian process used to fit data for the reference strain to be [10^-6, 10^-4]. This hyperparameter describes the flexibility of the fitted curve.

    As a citation for this software, please use:

    Lichten CA, White R, Clark IBN, Swain PS. Unmixing of fluorescence spectra to resolve quantitative time-series measurements of gene expression in plate readers. BMC Biotechnol 14 (2014) 11

    The technique for estimating growth rates is given in:

    Swain PS, Stevenson K, Leary A, Montano-Gutierrez LF, Clark IB, Vogel J, Pilizota T. Inferring time derivatives including cell growth rates using Gaussian processes. Nat Commun 7 (2016) 13766

    '''

    #####
    def __init__(self, dname, aname= 'default', prtype= 'Tecan', wdir= '', dsheetname= 0, asheetname= 0,
                 ODfname= 'ODcorrection_Glucose_Haploid.txt', warn= False, info= True, standardgain= False):
        '''
        Requires a data file, fname, (.csv, .txt, .xls, or .xlsx). An annotation file, aname, giving the contents of each well is optional. The plate reader is assumed to be a Tecan machine (but alternatives can in principle be specified using the platereader argument). A working directory can be set using wdir, such as wdir= 'data/'.

        Arguments
        --
        dname: file name of data (if no extension is given, '.xlsx' is assumed)
        aname: file name for annotation (if 'default', assumed to be the root of dname + '_contents.xlsx)
        prtype: type of plate reader ('Tecan')
        dsheetname: specifies the sheet of the data Excel file, e.g. dsheetname= 'Sheet2'
        asheetname: specifies the sheet of the annotation Excel file
        wdir: name of directory where the data files are stored
        ODfname: file name for dilution data for corrected OD (default is 'ODcorrection.txt')
        warn: if False (default), warnings created by covariance matrices that are not positive semi-definite are stopped
        info: if True (default), display information on the plate after loading
        standardgain: if defined, fluorescence measurements are corrected to this gain
        '''
        self.version= '4.6'
        self.dsheetname= dsheetname
        self.asheetname= asheetname
        if '.' not in dname: dname += '.xlsx'
        # specify working directory
        self.wdir= wdir
        self.dname= dname
        self.ODfname= ODfname

        # general parameters
        self.gamma= 0.114   # ratio of 585 to 525 for eGFP
        self.nosamples= 100  # for estimating error through sampling
        self.consist= 2   # number of stds of corrected reference strain data that is considered measurement noise
        rows= 'ABCDEFGH'
        self.overflow= -999.99

        # correction has not been performed
        self.standardgain= standardgain
        self.ODcorrected= False
        self.processedref1= False
        self.processedref2= False
        self.autocorrected= {}
        self.mediacorrected= {}
        self.ignoredwells= []
        self.negativevalues= False

        if not warn:
            # warning generated occasionally when sampling from the Gaussian process likely because of numerical errors
            import warnings
            warnings.simplefilter("ignore", RuntimeWarning)

        print('Platereader', self.version, ': loading data')
        # define annotation
        if aname == 'default': aname= dname.split('.')[0] + '_contents.xlsx'
        if aname:
            # import annotation
            try:
                r= pd.ExcelFile(self.wdir + aname).parse(asheetname)
            except FileNotFoundError:
                raise(SystemExit("\nError: Can't find " + self.wdir + aname))
        else:
            # no annotation: use wells
            r= [[a + str(j) + ' in unspecified' for j in np.arange(1,13)] for a in rows]
            r= pd.DataFrame(r, index= [a for a in rows], columns= [j for j in np.arange(1,13)])


        # import data
        df, dlabels, gains, datatypes, t, idstartindex, dstartindex= self.importdata(prtype, dname, dsheetname)
        # check data type names are defined otherwise use defaults
        defaults= ['OD'] + ['FL'+ str(i) for i in range(len(datatypes)-1)]
        for j, dtn in enumerate(datatypes):
            if dtn == '' or dtn == 'nan' or (type(dtn) == float and np.isnan(dtn)):
                datatypes[j]= defaults[j]


        # create data structure
        nodata= len(t)
        S= {}
        platelabels= {}
        allconditions= []
        allstrains= []
        alldata= []
        for let in rows:
            for no in np.arange(1,13):
                plateloc= let + str(no)
                try:
                    alabel= r[no][let]
                    if ' in ' in alabel:
                        # analyse annotation
                        ells= alabel.split('in')
                        strain= ells[0].strip()
                        condition= 'in'.join(ells[1:]).strip()
                        alldata.append(strain + ' in ' + condition)
                        platelabels[plateloc]= strain + ' in ' + condition
                        # get data
                        id= np.nonzero(dlabels == plateloc)[0][idstartindex:]
                        if np.any(id):
                            dsub= np.empty((nodata, len(datatypes)))
                            for jd in range(len(id)):
                                dsub[:,jd]= np.transpose(df.ix[id[jd]][dstartindex:].dropna().values)
                                # replace overflow values with NaN
                                dsub[:,jd][dsub[:,jd] == self.overflow]= np.nan
                            # add to data structure
                            if condition not in S.keys():
                                # new condition
                                S[condition]= {}
                                # record all conditions in data set
                                if condition != 'media': allconditions.append(condition)
                            if strain in S[condition].keys():
                                # remember location on plate
                                S[condition][strain]['plateloc'].append(plateloc)
                                S[condition][strain]['originalplateloc'].append(plateloc)
                                # add data to structure
                                S[condition][strain]['data']= np.append(S[condition][strain]['data'],
                                                                        dsub, axis= 1)
                            else:
                                # new strain
                                S[condition][strain]= {'data' : dsub, 'plateloc' : [plateloc], 'time': t,
                                                       'originalplateloc' : [plateloc]}
                                # record all strains in data set
                                if strain not in allstrains and strain != 'null':
                                    allstrains.append(strain)

                        else:
                            print('Data missing for', plateloc)
                    else:
                        print('Ignoring', alabel, '(wrong notation)')
                except TypeError:
                    print('No label specified for well', plateloc)
        self.d= S
        self.t= t
        self.nodata= nodata
        self.datatypes= datatypes
        self.gains= gains
        for dn in datatypes:
            self.autocorrected[dn]= False
            self.mediacorrected[dn]= False
        self.nooutchannels= len(datatypes)
        self.extractdata()
        self.platelabels= platelabels
        self.allstrains= allstrains
        self.allconditions= allconditions
        self.alldata= list(np.unique(alldata))
        self.importtime= '{:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now())
        if info: self.info()


    #####
    def importdata(self, prtype, dname, dsheetname):
        '''
        Internal function: Creates and parses dataframe from plate reader Excel file.
        '''
        try:
            xl= pd.ExcelFile(self.wdir + dname)
        except FileNotFoundError:
            raise(SystemExit("\nError: Can't find " + self.wdir + dname))
        df= xl.parse(dsheetname)
        if prtype == 'Tecan':
            # for Tecan plate reader post 2015
            df= df.replace('OVER', self.overflow)
            dlabels= df[df.columns[0]].values
            # find plate
            snField= df.keys()[['Tecan i-control ' in j for j in df.keys()]][0]
            self.serialnumber= df[snField][0].split('Serial number: ')[1]
            if self.serialnumber == '1006006292':
                self.machine= 'Plate Reader 1'
            elif self.serialnumber == '1411011275':
                self.machine= 'Plate Reader 2'
            else:
                self.machine= 'unknown'
            # find date of experiment
            idate= np.nonzero(df[df.columns[0]] == 'Date:')[0][0]
            self.expdate= df.loc[idate, df.columns[1]]
            # time in hours
            itime= np.nonzero(dlabels == 'Time [s]')[0]
            t= df.ix[itime[0]][1:].dropna().values.astype('float')/3600.0
            # find gains
            gains= [df.ix[igain][4] for igain in np.nonzero(dlabels == 'Gain')[0]]
            # indices to extract data
            idstartindex= 0
            dstartindex= 1
            # data types
            datatypes= [str(df.ix[it-2][0]) for it in itime]
        elif prtype == 'Hidex':
            # Hidex plate reader
            dlabels= df[df.columns[0]].values
            # time in hours
            itime= np.nonzero(df[df.columns[1]].values == 'Time (s)')[0]
            t= df.ix[itime[0]][2:].dropna().values.astype('float')/3600.0
            gains= []
            # indices to extract data
            idstartindex= 1
            dstartindex= 2
            # use default datatypes
            datatypes= ['nan']*(np.diff(itime)[0]-1)
        else:
            print('Platereader not recognized')
            return
        return df, dlabels, gains, datatypes, t, idstartindex, dstartindex


    #####
    def extractdata(self):
        '''
        Internal function: Creates individual data fields for OD, GFP, etc., for each strain from raw data.
        '''
        S, noout, datatypes= self.d, self.nooutchannels, self.datatypes
        for c in self.getcons('all'):
            for s in self.getstrains('all', c):
                noreps= int(S[c][s]['data'].shape[1]/noout)
                od, ignoredwells= [], []
                for i in range(noout): exec('f' + str(i) + '= []')
                for p, i in zip(S[c][s]['plateloc'], range(0, noreps*noout, noout)):
                     # check for ignored wells
                    if p in self.ignoredwells:
                        print('Ignoring well', p)
                        ignoredwells.append(p)
                    else:
                        # extract data
                        od.append(S[c][s]['data'][:,i])
                        for j in range(noout-1):
                            if self.standardgain and self.standardgain != self.gains[j]:
                                # correct fluorescence for gain
                                cdata= gcm*S[c][s]['data'][:,i+j+1] + gcc
                                exec('f' + str(j) + ".append(cdata)")
                            else:
                                # keep original gain
                                exec('f' + str(j) + ".append(S[c][s]['data'][:,i+j+1])")
                # remove ignored wells
                for p in ignoredwells:
                    S[c][s]['plateloc'].pop(S[c][s]['plateloc'].index(p))
                S[c][s]['ignoredwells']= list(set(S[c][s]['originalplateloc']) - set(S[c][s]['plateloc']))
                # store results and means and variances
                S[c][s][datatypes[0]]= np.transpose(od)
                for i in range(noout-1):
                    exec("S[c][s]['" + datatypes[i+1] + "']= np.transpose(f" + str(i) + ")")
        self.updatemeans()


    def updatemeans(self):
        '''
        Internal function: Calculates means and variances of all datatypes
        '''
        S= self.d
        for c in self.getcons('all'):
            for s in self.getstrains('all', c):
                for dn in self.datatypes:
                    if S[c][s][dn].shape[1] > 1:
                        S[c][s][dn + 'mn']= np.mean(S[c][s][dn], 1)
                        S[c][s][dn + 'var']= np.var(S[c][s][dn], 1)
                    else:
                        S[c][s][dn + 'mn']= gu.makerow(S[c][s][dn])
                        S[c][s][dn + 'var']= np.zeros(len(S[c][s][dn]))


    #####
    def ignorewells(self, exclude= [], clearall= False):
        '''
        Allows wells to be ignored in any future processing. If called several times, the default behaviour is for any previously ignored wells not to be re-instated.

        Arguments
        --
        exclude: list of labels (platereader locations) of wells to be excluded
        clearall: if True, all previously ignored wells are re-instated.
        '''
        exclude= gu.makelist(exclude)
        S= self.d
        if clearall and hasattr(self, 'ignoredwells'):
            # forget any previously ignoredwells
            delattr(self, 'ignoredwells')
        else:
            # wells cannot be ignored twice
            exclude= list(set(exclude) - set(self.ignoredwells))
        # set plateloc to initial values if necessary
        for c in self.getcons('all'):
            for s in self.getstrains('all', c):
                S[c][s]['plateloc']= list(S[c][s]['originalplateloc'])
        # store ignoredwells
        if hasattr(self, 'ignoredwells'):
            self.ignoredwells += exclude
        else:
            self.ignoredwells= exclude
        # extract data
        self.extractdata()
        print('Warning: analysis routines should be re-run')


    #####
    def info(self):
        '''
        Displays conditions, strains, datatypes, and corrections made.

        Arguments
        --
        None
        '''
        print('\n'+self.dname.split('.')[0])
        print('---')
        print('Conditions:')
        for c in self.allconditions: print('\t', c)
        print('Strains:')
        for s in self.allstrains: print('\t', s)
        # check if null present
        for s in self.alldata:
            if 'null' == s.split()[0]:
                print('\t', 'null')
                break
        print('Data types:')
        for i, d in enumerate(self.datatypes):
            if d == 'OD':
                print('\t', d)
            else:
                print('\t {:10} (gain {:3d})'.format(d, self.gains[i-1]))
        if self.standardgain: print('Gains normalized to gain', self.standardgain)
        print('Ignored wells:')
        if self.ignoredwells:
            for d in self.ignoredwells: print('\t', d)
        else:
            print('\t', 'None')
        print('Corrected for nonlinearities in OD:', self.ODcorrected)
        for dn in self.datatypes:
            print(dn + ' corrected for media:', self.mediacorrected[dn])
            if dn != 'OD':
                print(dn + ' corrected for autofluorescence:', self.autocorrected[dn])


    def attributes(self):
        '''
        Displays the names of the attributes available.

        Arguments
        --
        None
        '''
        ignore= ['d', 'consist', 't', 'nosamples', 'gamma', 'ODfname', 'overflow', 'nooutchannels', 'nodata']
        for a in self.__dict__:
            if 'corrected' not in a and 'processed' not in a and a not in ignore: print(a)




    #####
    # Internal functions
    #####
    def getcons(self, conditions, nomedia= False, includes= False, excludes= False):
        '''
        Internal function: Creates a list of conditions
        '''
        if conditions == 'all' or includes or excludes:
            cons= list(self.d.keys())
            if nomedia and 'media' in cons:
                cons.pop(cons.index('media'))
            # find those conditions containing keywords given in 'includes'
            if includes:
                includes= gu.makelist(includes)
                newcons= []
                for condition in cons:
                    gotone= 0
                    for item in includes:
                        if item in condition: gotone += 1
                    if gotone == len(includes): newcons.append(condition)
                cons= newcons
            # remove any conditions containing keywords given in 'excludes'
            if excludes:
                excludes= gu.makelist(excludes)
                exconds= []
                for condition in cons:
                    for item in excludes:
                        if item in condition:
                            exconds.append(condition)
                            break
                for ex in exconds:
                    cons.pop(cons.index(ex))
        else:
            cons= gu.makelist(conditions)
        if cons:
            return sorted(cons)
        else:
            if includes:
                raise(SystemExit('No conditions have ' + ' and '.join(includes)))
            else:
                raise(SystemExit('No conditions found'))



    #####
    def getstrains(self, strains, condition, nonull= False, includes= False, excludes= False):
        '''
        Internal function: Creates a list of strains
        '''
        if strains == 'all' or includes or excludes:
            ss= list(self.d[condition].keys())
            if nonull and 'null' in ss:
                ss.pop(ss.index('null'))
            # find those strains containing keywords given in 'includes'
            if includes:
                includes= gu.makelist(includes)
                newss= []
                for strain in ss:
                    sc= 0
                    for item in includes:
                         if item in strain: sc += 1
                    if sc == len(includes): newss.append(strain)
                ss= newss
            # remove any strains containing keywords given in 'excludes'
            if excludes:
                excludes= gu.makelist(excludes)
                exsts= []
                for strain in ss:
                    for item in excludes:
                        if item in strain:
                            exsts.append(strain)
                            break
                for ex in exsts:
                    ss.pop(ss.index(ex))
        else:
            ss= gu.makelist(strains)
        if ss:
            return sorted(ss)
        else:
            if includes:
                print('Warning: No strains have ' + ' and '.join(includes) + ' for ' + condition)
            else:
                print('Warning: No strains found for ' + condition)
            return []


    #####
    def getconditionsandstrains(self, conditions= 'all', strains= 'all', conditionincludes= False,
                                strainincludes= False, conditionexcludes= False, strainexcludes= False,
                                nomedia= True, nonull= True):
        '''
        Returns a list of (condition, strain) tuples for iteration.

        Arguments
        --
        conditions: list of experimental conditions to include (default is 'all')
        strains: list of strains to include (default is 'all')
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        nomedia: ignores condition= media if True (default)
        nonull: ignores strain= null if True (default)
        '''
        return [(c,s) for c in self.getcons(conditions, nomedia= nomedia, includes= conditionincludes, excludes= conditionexcludes)
                for s in self.getstrains(strains, c, nonull= nonull, includes= strainincludes, excludes= strainexcludes)]


    #####
    def getkeys(self, conditions= 'all', strains= 'all', conditionincludes= False, strainincludes= False,
             conditionexcludes= False, strainexcludes= False, title= True):
        '''
        Returns the keys (the items) calculated for each strain.

        Arguments
        --
        conditions: list of experimental conditions to include (default is 'all')
        strains: list of strains to include (default is 'all')
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        title: display condition and strain
        '''
        S= self.d
        for c in self.getcons(conditions, True, conditionincludes, conditionexcludes):
            for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                if title: print(c, 'in', s, '\n---')
                for k in sorted(list(S[c][s].keys())): print('\t'+k)
                print('\n')


    #####
    # OD correction
    #####
    def correctOD(self, ODfname= False, figs= True, correctmedia= True, conditionincludes= False, strainincludes= False,
             conditionexcludes= False, strainexcludes= False):
        '''
        Corrects OD data for a non-linear relationship between OD and cell number.
        Requires a dilution data set in the file ODfname. The default dilution curve is 'ODcorrection.txt' (measured by C Lichten in 2010).

        Arguments
        --
        ODfname: file name for dilution data
        figs: if True, a plot of the fit to the dilution data is produced
        correctmedia: if True (default), correct OD measurements by the OD of the media
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        if not self.ODcorrected:
            # correct for media
            if correctmedia: self.correctmedia(datatypes= 'OD', conditionincludes= conditionincludes,
                                               conditionexcludes= conditionexcludes)
            # if no ODfname use the default
            if not ODfname: ODfname= self.wdir + self.ODfname
            S, noout= self.d, self.nooutchannels
            # fit dilution data
            if not hasattr(self, 'gc'): self.findODcorrection(ODfname, figs)
            # correct all wells containing strains
            for c in self.getcons('all', True, conditionincludes, conditionexcludes):
                for s in self.getstrains('all', c, False, strainincludes, strainexcludes):
                    for i in range(S[c][s]['OD'].shape[1]):
                        # perform correction
                        self.gc.predict(S[c][s]['OD'][:,i])
                        S[c][s]['OD'][:,i]= self.gc.f
            self.updatemeans()
            self.ODcorrected= True
        else:
            print('OD is already corrected')


    #####
    def findODcorrection(self, ODfname, figs):
        '''
        Internal function: Uses a Gaussian process to fit serial dilution data to correct for non-linearities in the relationship between OD and cell density. The data are expected in a file ODfname.
        '''
        print('Fitting dilution data for OD correction for non-linearities')
        try:
            od, dilfac= np.loadtxt(ODfname, unpack= True)
            print('Using', ODfname)
        except FileNotFoundError:
            raise(SystemExit("\nError: Can't find " + ODfname))
        # process data
        dilfac= dilfac[np.argsort(od)]
        od= np.sort(od)
        # run Gaussian process
        gc= gp.sqexplinGP({0: (-4, 2), 1: (-3, 1), 2: (-6, 1), 3: (-6, 1)}, od, dilfac)
        gc.findhyperparameters(noruns= 5, exitearly= True, quiet= True)
        gc.predict(od)
        if figs:
            plt.figure()
            gc.sketch('.')
            plt.xlim([0, 1.05*np.max(od)])
            plt.ylim([0, 1.05*np.max(od)])
            plt.gca().set_aspect('equal', adjustable='box')
            plt.draw()
            plt.xlabel('OD')
            plt.ylabel('relative cell density')
            plt.title('function for correcting OD')
            plt.show(block= False)
        self.gc= gc


    #####
    # Media correction
    #####
    def correctmedia(self, datatypes= 'OD', conditions= 'all', figs= True, noruns= 3, exitearly= False, bd= False,
                     results= False, conditionincludes= False, conditionexcludes= False, mean= False):
        '''
        Corrects OD or fluorescence for the OD or fluorescence of the media. Uses a Gaussian process to fit the time-series of measurements of all replicates of the media and subtracts this time series from the raw data.

        Arguments
        --
        conditions: list of experimental conditions to be corrected
        datatypes: data types to be corrected (default is 'OD')
        figs: if True, display fits of the media
        noruns: number of attempts used to fit the media
        exitearly: if True, stop at the first successful fit; if False, take the best fit from all successful fits
        bd: can be used to change the limits on the hyperparameters for the Gaussian process, which has a squared exponential kernel, used to fit the media
        results: if True, display best-fit parameters
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        mean: if True, the mean over time of the media values is used for the correction (rather than fitting with a Gaussian process)
        '''
        S= self.d
        cons= self.getcons(conditions, False, conditionincludes, conditionexcludes)
        datatypes= gu.makelist(datatypes)
        for dn in datatypes:
            if not self.mediacorrected[dn]:
                if 'media' in S.keys():
                    # one set of media measurements for all wells
                    print('Correcting', dn, 'for media')
                    mc= self.findmediacorrection(dn, 'media', figs, noruns, exitearly, bd, results,
                                                 'media correction for ' + dn)
                    for c in cons: self.performmediacorrection(dn, c, mc)
                else:
                    # media measurements for each condition
                    for c in cons:
                        print('Correcting', dn, 'for media for', c)
                        if 'null' in S[c]:
                            self.performmediacorrection(dn, c, figs, noruns, exitearly, bd, results, mean)
                        else:
                            print(' No well annotated "null" was found and media correction abandoned')
                self.updatemeans()
                self.mediacorrected[dn]= True
            else:
                print(dn + ' is already corrected for the media')

    #####
    def performmediacorrection(self, dtype, condition, figs, noruns, exitearly, bd, results, mean):
        '''
        Internal function: Uses a Gaussian process to fit the media over time and subtracts the best-fit media values from the data, as well as storing the uncorrected data.
        '''
        S= self.d
        t, data= S[condition]['null']['time'], S[condition]['null'][dtype]
        if mean:
            f= np.mean(data)
        else:
            # d is data and x is time
            d= data.flatten('F')
            x= np.tile(t, data.shape[1])
            b= {0: (-6,4), 1: (-5,-2), 2: (-10,0)}
            if bd: b= gu.mergedicts(original= b, update= bd)
            # fit with Gaussian process
            try:
                f, fvar= gu.smoothGP(x, d, xp= t, bd= b, results= results)
                if figs:
                    plt.figure()
                    plt.plot(x, d, 'ro', t, f, 'b-')
                    plt.xlabel('time (hours)')
                    plt.title('media correction for ' + dtype + ' in ' + condition)
                    plt.show(block= False)
            except gp.gaussianprocessException:
                raise(SystemExit('Fitting media failed'))
        # perform correction
        for s in self.getstrains('all', condition, True):
            sh= S[condition][s][dtype].shape
            S[condition][s][dtype] -= np.reshape(np.tile(f, sh[1]), (-1, sh[1]), order= 'F')
            if np.any(S[condition][s][dtype] < 0):
                print(' Warning: negative data for', s , 'in', condition)
                wstr= '\t' + dtype + ': ' + s + ' in ' + condition + '\n'
                if not self.negativevalues:
                    self.negativevalues= wstr
                else:
                    self.negativevalues += wstr


    #####
    # Fluorescence corrections
    #####
    def correctauto(self, f= ['GFP', 'AutoFL'], conditions= 'all', strains= 'all', refstrain= 'WT', figs= True,
                    correctOD= True, noruns= 2, bd= False, no1samples= 100, conditionincludes= False,
                    strainincludes= False, conditionexcludes= False, strainexcludes= False, results= False,
                    correctmedia= True, mediausemean= False, forceprocessref= False):
        '''
        Corrects fluorescence data for autofluorescence in comparison with a reference strain. Emissions can be measured at one or two wavelengths.

        Arguments
        --
        f: fluorescence measurements to be used, typically either ['mCherry'] or  ['GFP', 'AutoFL'] (default)
        conditions: list of experimental conditions to be corrected (default: 'all')
        strains: list of strains to be corrected (default: 'all')
        refstrain: the name of the reference strain ('WT' is the default)
        figs: if True (default), display fits used in applying the correction
        correctOD: if True (default), corrects OD measurements for non-linearities between OD and the number of cells
        noruns: number of fit attempts to be used
        bd: to change the limits on the hyperparameters for the Gaussian process, which has a squared exponential kernel and is used for fitting data from the reference strain
        no1samples: number of samples used to estimate errors when correcting fluorescence data measured at one wavelength
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        results: if True, display best-fit parameters
        correctmedia: if True (default), correct for media, which is only strictly necessary when using spectral unmixing (two fluoresence measurements) to correct autofluorescence
        mediausemean: if True (default: False), use the mean of the media for correcting (rather than fitting a Gaussian process)
        forceprocessref: if True (default: False), re-run the processing of the reference strain
        '''
        f= gu.makelist(f)
        S= self.d
        # correct OD
        if correctOD: self.correctOD(correctmedia= correctmedia, conditionincludes= conditionincludes,
                                     conditionexcludes= conditionexcludes)
        # correct autofluorescence
        print('Using', refstrain, 'as the reference')
        if len(f) == 2:
            # check have enough replicates
            go= True
            for c in self.getcons(conditions, True, conditionincludes, conditionexcludes):
                for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                    if S[c][s]['OD'].shape[1] < 2:
                        print('Not enough OD replicates to correct autofluorescence for', s, 'in', c)
                        go= False
                    if S[c][s][f[0]].shape[1] < 2:
                        print('Not enough fluorescence replicates to correct autofluorescence for', s, 'in', c)
                        go= False
            if go == False:
                print('Try specifying just one fluorescence measurement')
            else:
                self.correctauto2(f, conditions, strains, refstrain, figs, correctOD, noruns, bd,
                                  conditionincludes, strainincludes, conditionexcludes, strainexcludes,
                                  results, correctmedia, mediausemean, forceprocessref)
                if self.negativevalues:
                    print('Warning: correcting media has created negative data values.')
                    print('These values have been ignored for:')
                    print(self.negativevalues)
        elif len(f) == 1:
            self.correctauto1(f, conditions, strains, refstrain, figs, correctOD, noruns, bd, no1samples,
                              conditionincludes, strainincludes, conditionexcludes, strainexcludes,
                              results, correctmedia, mediausemean, forceprocessref)
        else:
            print('f must be a list of length 1 or 2')



    #####
    def correctauto1(self, f, conditions, strains, refstrain, figs, correctOD, noruns, bd, nosamples,
                     conditionincludes, strainincludes, conditionexcludes, strainexcludes, results,
                     correctmedia, mediausemean, forceprocessref):
        '''
        Internal function: Corrects for autofluorescence for experiments with measured emissions at one wavelength using the fluorescence of the wild-type interpolated to the OD of the tagged strain.
        '''
        S= self.d
        if forceprocessref: self.processedref1= False
        # correct for media
        if correctmedia:
            self.correctmedia(datatypes= f, conditions= conditions, conditionincludes= conditionincludes,
                              conditionexcludes= conditionexcludes, results= results, mean= mediausemean)
        # process reference strain
        if not self.processedref1:
            self.processref1(f, conditions, refstrain, figs, noruns, bd, conditionincludes, conditionexcludes,
                             results)
        else:
            print('Reference strain is already processed')

        # correct autofluorescence
        print('Correcting autofluorescence')
        # run through all conditions
        for c in self.getcons(conditions, True, conditionincludes, conditionexcludes):
            gfr= S[c][refstrain]['gpref']
            for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                if s != refstrain:
                    noreps= S[c][s][f[0]].shape[1]
                    fl, fl2, flperod, flperod2= (0,)*4
                    # run over all replicates
                    for i in range(noreps):
                        # interpolate WT fluorescence errors to OD values of strain
                        from scipy.interpolate import interp1d
                        interpf= interp1d(S[c][refstrain]['ODmn'], S[c][refstrain]['fmerr'])
                        try:
                            men= interpf(S[c][s]['OD'][:,i])
                        except ValueError:
                            men= np.median(S[c][refstrain]['fmerr'])*np.ones(len(S[c][refstrain]['fmerr']))
                        # predict WT fluorescence at OD values of strain
                        gfr.predict(S[c][s]['OD'][:,i], merrorsnew= men, addnoise= True)
                        # sample for estimating errors
                        fs= gfr.sample(nosamples)
                        if correctOD:
                            # sample ODs corrected for non-linearities
                            self.gc.predict(S[c][s]['OD'][:,i])
                            cods= self.gc.sample(nosamples)
                        else:
                            cods= gu.tilec(S[c][s]['OD'][:,i], nosamples)
                        # find fluorescence per cell
                        for j in range(nosamples):
                            d= S[c][s][f[0]][:,i] - fs[:,j]
                            dperod= d/cods[:,j]
                            fl += d
                            fl2 += d**2
                            flperod += dperod
                            flperod2 += dperod**2
                    bname= 'c-' + f[0]
                    S[c][s][bname]= fl/(noreps*nosamples)
                    S[c][s][bname + 'var']= fl2/(noreps*nosamples) - S[c][s][bname]**2
                    S[c][s][bname + 'perod']= flperod/(noreps*nosamples)
                    S[c][s][bname + 'perodvar']= flperod2/(noreps*nosamples) - S[c][s][bname + 'perod']**2
        self.autocorrected[f[0]]= True
        print('Created:\n' + '\tc-' + f[0] + 'perod (corrected fluorescence per cell)\n' + '\tc-' + f[0]
              + ' (corrected fluorescence)')



    #####
    def processref1(self, f, conditions, refstrain, figs, noruns, bd, conditionincludes, conditionexcludes, results):
        '''
        Processes reference strain for data with one fluorescence measurement. Uses a Gaussian process to fit the fluorescence as a function of OD.

        Arguments
        --
        f: fluorescence to be corrected, such as ['mCherry']
        conditions: list of experimental conditions to be corrected
        refstrain: the name of the reference strain, such as 'WT'
        figs: if True, display fits of the fluorescence
        noruns: number of fit attempts used
        bd: to change the limits on the hyperparameters for the Gaussian process, which has a squared exponential kernel and is used to fit the fluorescence
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        results: if True, display best-fit parameters
        '''
        S= self.d
        # bounds for Gaussian process for fitting reference strain
        b= {0: (3,10), 1: (-4,4), 2: (-4,2)}
        if bd: b= gu.mergedicts(original= b, update= bd)

        # run through all conditions
        for c in self.getcons(conditions, True, conditionincludes, conditionexcludes):
            print('Processing reference strain for', c)
            # fit reference strain's fluorescence as a function of OD
            x= S[c][refstrain]['OD'].flatten('F')
            y= S[c][refstrain][f[0]].flatten('F')
            me= gu.findsmoothvariance(S[c][refstrain][f[0]])
            ys= y[np.argsort(x)]
            mes= np.tile(me, S[c][refstrain][f[0]].shape[1])[np.argsort(x)]
            xs= np.sort(x)
            # fit with Gaussian process
            gfr= gp.sqexpGP(b, xs, ys, merrors= mes)
            try:
                gfr.findhyperparameters(noruns= noruns)
                if results: gfr.results()
                gfr.predict(xs)
            except gp.gaussianprocessException:
                raise(SystemExit('Fitting reference strain failed'))
            if figs:
                # plot fit
                plt.figure()
                gfr.sketch('o')
                plt.xlabel('OD')
                plt.ylabel(f[0])
                plt.title('fitting ' + refstrain + ' for ' + c)
                plt.show(block= False)
            S[c][refstrain]['gpref']= gfr
            S[c][refstrain]['fmerr']= me
        # store reference strain information
        for c in self.getcons('all', True, conditionincludes, conditionexcludes):
            for s in self.getstrains('all', c, True):
                S[c][s]['refstrain']= refstrain
        self.processedref1= True



    #####
    def correctauto2(self, f, conditions, strains, refstrain, figs, correctOD, noruns, bd, conditionincludes,
                     strainincludes, conditionexcludes, strainexcludes, results, correctmedia, mediausemean,
                     forceprocessref):
        '''
        Internal function: Corrects for autofluorescence using spectral unmixing for experiments with measured emissions at two wavelengths (following Lichten et al.)
        '''
        S, noout= self.d, self.nooutchannels
        if forceprocessref: self.processedref2= False
        # correct for media
        if correctmedia:
            self.correctmedia(conditions= conditions, datatypes= f, conditionincludes= conditionincludes,
                              conditionexcludes= conditionexcludes, results= results, mean= mediausemean)
        # process reference strain
        if not self.processedref2:
            self.processref2(f, conditions, refstrain, figs, noruns, bd, conditionincludes, conditionexcludes, results)
        else:
            print('Reference strain is already processed')

        # correct for autofluorescence
        print('Correcting autofluorescence')
        for c in self.getcons(conditions, True, conditionincludes, conditionexcludes):
            gr= S[c][refstrain]['gpref']
            for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                if s != refstrain:
                    nodata, noreps= S[c][s][f[0]].shape
                    # remove autofluorescence
                    fl= self.applyautoflcorrection(S[c][refstrain]['fratio'], S[c][s][f[0]], S[c][s][f[1]])
                    # estimate error
                    varf= np.var(S[c][s][f[0]], 1)
                    varcf= np.var(S[c][s][f[1]], 1)
                    if np.array_equal(S[c][refstrain]['time'], S[c][s]['time']):
                        gr.predict(S[c][refstrain]['time'], merrorsnew= S[c][refstrain]['qrerr'])
                    else:
                        raise(SystemExit('Error: the reference strain is measured at different time points from '
                                         + s + ' in ' + c))
                    rs= np.transpose(gr.sample(self.nosamples))
                    varfl= np.zeros(nodata)
                    # average variance over samples of ra
                    for ra in rs:
                        varfl += (varf + ra**2*varcf)/(self.gamma - ra)**2
                    varfl /= self.nosamples
                    bname= 'c-' + f[0]
                    S[c][s][bname]= np.mean(fl,1)
                    S[c][s][bname + 'var']= varfl
                    # fluorescence per cell
                    S[c][s][bname + 'perod']= S[c][s][bname]/S[c][s]['ODmn']
                    flperodvar= np.zeros(nodata)
                    for i in np.random.randint(noreps, size= self.nosamples):
                        if correctOD:
                            self.gc.predict(S[c][s]['OD'][:,i])
                            cc= self.gc.sample(1)
                            cc= cc.flatten('F')
                            flperodvar += S[c][s][bname + 'var']/cc**2
                        else:
                            flperodvar += S[c][s][bname + 'var']/S[c][s]['OD'][:,i]**2
                    flperodvar /= self.nosamples
                    S[c][s][bname + 'perodvar']= flperodvar
                    # calculate corrected levels that are above the expectation from reference strain
                    reflevel= self.consist*np.sqrt(S[c][refstrain][bname + 'var'])
                    keep= np.nonzero(S[c][s][bname] > reflevel)[0]
                    S[c][s]['s' + bname]= np.zeros(np.size(S[c][s][bname]))
                    S[c][s]['s' + bname][keep]= S[c][s][bname][keep]
                    S[c][s]['s' + bname + 'var']= np.zeros(np.size(S[c][s][bname + 'var']))
                    S[c][s]['s' + bname + 'var'][keep]= S[c][s][bname + 'var'][keep]
                    reflevel /= S[c][refstrain]['ODmn']
                    keep= np.nonzero(S[c][s][bname+ 'perod'] > reflevel)[0]
                    S[c][s]['s' + bname + 'perod']= np.zeros(np.size(S[c][s][bname + 'perod']))
                    S[c][s]['s' + bname + 'perod'][keep]= S[c][s][bname + 'perod'][keep]
                    S[c][s]['s' + bname + 'perodvar']= np.zeros(np.size(S[c][s][bname + 'perodvar']))
                    S[c][s]['s' + bname + 'perodvar'][keep]= S[c][s][bname + 'perodvar'][keep]
        self.autocorrected[f[0]]= True
        print('Created:\n' + '\tc-' + f[0] + 'perod (corrected fluorescence per cell)\n' + '\tc-' + f[0]
              + ' (corrected fluorescence)\n' + '\tsc-' + f[0]
              + 'perod (only statistically significant corrected fluorescence per cell)\n' + '\tsc-'
              + f[0] + ' (only statistically significant corrected fluorescence)')

    #####
    def processref2(self, f, conditions, refstrain, figs, noruns, bd, conditionincludes, conditionexcludes, results):
        '''
        Processes reference strain data for spectral unmixing (for experiments with two fluorescence measurements). Uses a Gaussian process to fit the ratio of emitted fluorescence measurements and checks that reference strain data is itself corrected to zero.

        Arguments
        --
        f: fluorescence to be corrected, such as ['GFP', 'AutoFL']
        conditions: list of experimental conditions to be corrected
        refstrain: the name of the reference strain, such as 'WT'
        figs: if True, display fits of the ratios and, as a check, the correction of the reference strain by itself
        noruns: number of attempts used to fit the ratio of fluorescences
        bd: can be used to change the limits on the hyperparameters for the Gaussian process, which has a squared exponential kernel and is used to fit the ratio of fluorescences
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        results: if True, display best-fit parameters
        '''
        S= self.d
        # bounds for Gaussian process for fitting reference strain
        b= {0: (-5,3), 1: (-4,2), 2: (-4, 4)}
        if bd: b= gu.mergedicts(original= b, update= bd)

        # run through all conditions
        for c in self.getcons(conditions, True, conditionincludes, conditionexcludes):
            print('Processing reference strain for', c)
            noreps= S[c][refstrain][f[0]].shape[1]
            f0, f1, t= S[c][refstrain][f[0]], S[c][refstrain][f[1]], S[c][refstrain]['time']
            # find negative values
            dels= np.unique(np.append(np.nonzero(f0 <= 0)[0], np.nonzero(f1 < 0)[0]))
            # find values with zero variance
            dels= np.unique(np.append(dels, np.nonzero(np.var(f1/f0, 1) == 0)[0]))
            # remove offending values
            f0r= np.delete(f0, dels, axis= 0)
            f1r= np.delete(f1, dels, axis= 0)
            tr= np.delete(t, dels)
            # fit ratio of fluorescence
            qr= (f1r/f0r).flatten('F')
            qrerr= np.var(f1r/f0r, 1)
            x= np.tile(tr, noreps)
            xs= np.sort(x)
            qrs= qr[np.argsort(x)]
            qrerrs= np.repeat(qrerr, noreps)[np.argsort(x)]
            # fit with Gaussian process
            gr= gp.nnGP(b, xs, qrs, merrors= qrerrs)
            try:
                # add back missing time points
                if np.any(dels): qrerr= np.interp(t, tr, qrerr)
                # optimize hyperparameters of Gaussian process
                gr.findhyperparameters(noruns)
                gr.predict(t, merrorsnew= qrerr)
            except gp.gaussianprocessException:
                raise(SystemExit('Fitting reference strain failed'))
            # store results
            S[c][refstrain]['qr']= qr
            S[c][refstrain]['gpref']= gr
            S[c][refstrain]['qrerr']= qrerr
            S[c][refstrain]['fratio']= gr.f
            if figs:
                # plot fit
                plt.figure()
                gr.sketch('o')
                plt.xlabel('time (hours)')
                plt.ylabel(f[1] + '/' + f[0])
                plt.title('fitting ' + refstrain + ' for ' + c)
                plt.show(block= False)
            # check autofluorescence correction for reference strain
            flref= self.applyautoflcorrection(S[c][refstrain]['fratio'], S[c][refstrain][f[0]], S[c][refstrain][f[1]])
            bname= 'c-' + f[0]
            S[c][refstrain][bname]= np.mean(flref, 1)
            S[c][refstrain][bname + 'perod']= np.mean(flref/S[c][refstrain]['ODmn'][:,None], 1)
            S[c][refstrain][bname + 'var']= np.var(flref, 1)
            S[c][refstrain][bname + 'perodvar']= np.var(flref/S[c][refstrain]['ODmn'][:,None], 1)
            if figs:
                # plot correction for reference strain
                plt.figure()
                plt.plot(S[c][refstrain]['time'], flref, '.')
                plt.plot(S[c][refstrain]['time'], self.consist*np.sqrt(S[c][refstrain][bname + 'var']), 'r:')
                plt.plot(S[c][refstrain]['time'], -self.consist*np.sqrt(S[c][refstrain][bname + 'var']), 'r:')
                plt.plot(S[c][refstrain]['time'], np.zeros(np.size(S[c][refstrain]['time'])), 'k')
                plt.ylabel('corrected ' + refstrain + ' fluorescence')
                plt.xlabel('time (hours)')
                plt.title(c + ': consistency check for reference strain')
                plt.show(block= False)
        # store reference strain information
        for c in self.getcons('all', True, conditionincludes, conditionexcludes):
            for s in self.getstrains('all', c, True):
                S[c][s]['refstrain']= refstrain
                if s != refstrain:
                    S[c][s][bname + 'refstd']= np.sqrt(S[c][refstrain][bname + 'var'])
                    S[c][s]['refODmn']= S[c][refstrain]['ODmn']
        self.processedref2= True


    #####
    def applyautoflcorrection(self, ra, fdata, cfdata):
        '''
        Internal function: Corrects for autofluorescence returning an array of replicates.
        '''
        noreps= fdata.shape[1]
        raa= np.reshape(np.tile(ra, noreps), (np.size(ra), noreps), order= 'F')
        return (raa*fdata - cfdata)/(raa - self.gamma*np.ones(np.shape(raa)))

    #####
    def reset(self):
        '''
        Reset so that all data processing can be re-run.

        Arguments
        --
        None
        '''
        print('Resetting...')
        self.processedref1= False
        self.processedref2= False
        self.ODcorrected= False
        for dn in self.datatypes:
            self.autocorrected[dn]= False
            self.mediacorrected[dn]= False
        self.extractdata()


    #####
    # Plotting routines
    #####
    def plot(self, dtype= 'OD', conditions= 'all', strains= 'all', plate= False, onefig= False,
             plotod= True, nonull= False, conditionincludes= False, strainincludes= False,
             conditionexcludes= False, strainexcludes= False, includeref= False):
        '''
        Plots data for specified strains in specified conditions.

        Arguments
        --
        dtype: data type, either
           'labels' (shows labels for each well in a single figure)
            'OD', 'GFP', etc.,
            'c-GFPperod' (corrected fluorescence per cell),
            'c-GFP' (corrected fluorescence),
            'sc-GFPperod' (only statistically significant corrected fluorescence per cell),
            'sc-GFP' (only statistically significant corrected fluorescence)
        conditions: list of experimental conditions to be plotted (default is 'all')
        strains: list of strains to be plotted (default is 'all')
        plate: if set to True, data for each well are plotted in one figure
        onefig: if True, data for all strains are plotted in the same figure
        plotod: if True (default), mean(OD) is plotted on a right-hand y-axis for processed data
        nonull: if True, do not plot 'null' strains
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        includeref: whether or not to include the reference strain in fluorescence plots
        '''
        S= self.d
        if dtype == 'labels': plate= True
        if 'c-' in dtype: nonull= True
        if plate: onefig= True
        if onefig:
            # one figure for all plots
            plt.figure()
            ax= plt.subplot(111)
            # count number of lines
            nlines= 0
            for c in self.getcons(conditions, False, conditionincludes, conditionexcludes):
                for s in self.getstrains(strains, c, nonull, strainincludes, strainexcludes):
                    nlines += S[c][s][self.datatypes[0]].shape[1]
            # define colour map
            colormap= plt.cm.nipy_spectral
            from cycler import cycler
            ax.set_prop_cycle(cycler('color', [colormap(i) for i in np.linspace(0, 1, nlines)]))
        # draw all plots
        for c in self.getcons(conditions, False, conditionincludes, conditionexcludes):
            for s in self.getstrains(strains, c, nonull, strainincludes, strainexcludes):
                if 'c-' in dtype and not includeref and 'refstrain' in list(S[c][s].keys()) and s == S[c][s]['refstrain']:
                    # ignore reference strain for fluorescence data
                    ignorestrain= True
                else:
                    ignorestrain= False
                if plate:
                    # plot as a whole plate
                    i= 0
                    for pl in S[c][s]['plateloc']:
                        rowl= 'ABCDEFGH'.index(pl[0])
                        coll= int(pl[1:])
                        sindex= coll + 12*rowl
                        plt.subplot(8, 12, sindex)
                        if dtype == 'labels':
                            # plot label
                            plt.axis([0, 10, 0, 10])
                            plt.text(5, 3, s + '\nin\n' + c, horizontalalignment='center',
                                     color= 'navy', fontsize= 'small')
                        else:
                            # plot data
                            plt.plot(S[c][s]['time'], S[c][s][dtype][:,i], '-')
                            i += 1
                        plt.tick_params(labelbottom= 'off', labelleft= 'off')
                        # label well locations
                        for j in range(12):
                            if sindex == j+1: plt.title(j+1)
                        for j, k in enumerate(np.arange(1, 96, 12)):
                            if sindex == k: plt.ylabel('ABCDEFGH'[j] + ' ', rotation= 0)
                else:
                    # plot data (not whole plate)
                    if 'c-' in dtype:
                        # plotting processed data
                        if not ignorestrain:
                            if not onefig: plt.figure()
                            ax= plt.subplot(1,1,1)
                            ax.plot(S[c][s]['time'], S[c][s][dtype], '.', label= s + ' in ' + c)
                            ax.plot(S[c][s]['time'], S[c][s][dtype]+np.sqrt(S[c][s][dtype+'var']), 'k:', alpha=0.4)
                            ax.plot(S[c][s]['time'], S[c][s][dtype]-np.sqrt(S[c][s][dtype+'var']), 'k:', alpha=0.4)
                            rname= 'c-' + dtype.split('-')[1].split('perod')[0]
                            if rname + 'refstd' in S[c][s]:
                                if 'perod' in dtype:
                                    ax.plot(S[c][s]['time'], self.consist*S[c][s][rname + 'refstd']/S[c][s]['refODmn'], 'r:')
                                else:
                                    ax.plot(S[c][s]['time'], self.consist*S[c][s][rname + 'refstd'], 'r:')
                            if plotod and not onefig:
                                # add OD to plot
                                ax2= ax.twinx()
                                ax2.plot(S[c][s]['time'], S[c][s]['ODmn'] , 'DarkOrange', linewidth=5, alpha=0.2)
                    else:
                        # plot raw data
                        if not onefig: plt.figure()
                        # remove ignored wells for the legend
                        pls= list(S[c][s]['plateloc'])
                        for igwell in self.ignoredwells:
                            if igwell in pls: pls.remove(igwell)
                        # plot each replicate
                        if S[c][s][dtype].ndim > 1:
                            for i in range(S[c][s][dtype].shape[1]):
                                if S[c][s][dtype].shape[1] == 1:
                                    label= s + ' in ' + c
                                else:
                                    label= pls[i] + ': ' + s + ' in ' + c
                                plt.plot(S[c][s]['time'], S[c][s][dtype][:,i], '.-', label= label)
                        else:
                            plt.plot(S[c][s]['time'], S[c][s][dtype], '.-', label= s + ' in ' + c)
                if not onefig and not ignorestrain:
                    # display and add labels for each plot
                    plt.title(dtype + ' of ' + s + ' in ' + c)
                    plt.xlabel('time (hours)')
                    if 'c-' not in dtype:
                        plt.ylabel(dtype)
                        plt.legend(pls, loc= 'lower right')
                    else:
                        ax.set_ylabel(dtype)
                        ax2.set_ylabel('mean(OD)')
                    plt.ylim(ymin= 0)
                    plt.show(block= False)
        if onefig:
            # display and add labels for single figure
            if plate:
                if dtype == 'labels':
                    plt.suptitle(self.dname.split('.')[0])
                else:
                    plt.suptitle(dtype)
            else:
                plt.title(dtype)
                plt.xlabel('time (hours)')
                # shrink current axis by 20%
                box= ax.get_position()
                ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
                # put legend to the right of the current axis
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.ylim(ymin= 0)
            plt.show(block= False)



    #####
    def savefigs(self, ftype= 'png', onefile= False):
        '''
        Saves all current figures, either each to a separate file (default) or all to one file.

        Arguments
        --
        ftype: type of file ('png' is the default)
        onefile: if True, all figures are saved to on PDF file
        '''
        if onefile:
            savename= self.wdir + self.dname.split('.')[0] + '.pdf'
            from matplotlib.backends.backend_pdf import PdfPages
            with PdfPages(savename) as pdf:
                for i in plt.get_fignums():
                    plt.figure(i)
                    pdf.savefig()
        else:
            for i in plt.get_fignums():
                plt.figure(i)
                savename= str(plt.getp(plt.gcf(), 'axes')[0].title).split("'")[1]
                savename= savename.replace(' ', '_')
                if savename == '': savename= 'Whole_plate_Figure_' + str(i)
                print('Saving', savename)
                plt.savefig(self.wdir + savename + '.' + ftype)


    #####
    def close(self):
        '''
        Close all figures.

        Arguments
        --
        None
        '''
        plt.close('all')


    #####
    # Statistical analysis
    #####
    def getstats(self, conditions= 'all', strains= 'all', dtype= 'OD', bd= False, cvfn= 'sqexp',
                 esterrs= False, noruns= 3, stats= True, plotodgr= False, conditionincludes= False,
                 strainincludes= False, conditionexcludes= False, strainexcludes= False):
        '''
        Calls fitderiv.py to estimate the first time-derivate.
        The data calculated by fitderiv is added to the overall data structure. For example, p.d['1% Gal']['GAL2']['flogOD'] gives the fit to the log OD curve if OD data is used.
        For OD data, 'flogOD' (the fitted log OD), 'fODvar' (the variance - error - in the fitted log OD), 'gr' (the inferred growth rate), 'grvar' (the variance in the inferred growth rate), and statistics, such as 'max growth rate' and 'lag time', are added for each strain.
        Equivalent variables are added for other types of data.

        Arguments
        --
        conditions: list of experimental conditions to be included
        strains: list of strains to be included
        dtype: type of data ('OD' - default, 'GFP', 'c-GFPperod', or 'c-GFP')
        bd: can be used to change the limits on the hyperparameters for the Gaussian process used for the fit. For example, p.odstats('1% Gal', 'GAL2', bd= {1: [-2,-2])}) fixes the flexibility to be 0.01
        cvfn: covariance function used for fit, either 'sqexp' (default) or 'nn'
        esterrs: if True, measurement errors are empirically estimated from the variance across replicates at each time point; if False, the size of the measurement error is fit from the data assuming that this size is the same at all time points
        noruns: number of attempts made for each fit (default is 3), each run is made with random initial estimates of the parameters
        stats: calculate statistics if True
        plotodgr: for OD data, plots growth rate versus log(OD) if True (default)
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        S, noout= self.d, self.nooutchannels
        for c in self.getcons(conditions, True, conditionincludes, conditionexcludes):
            for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                figtitle= s + ' in ' + c
                print('\nFitting', dtype, 'for', figtitle)
                try:
                    d= S[c][s][dtype]
                    if dtype == 'OD':
                        snames= ['max growth rate', 'time of max growth rate',
                                    'doubling time', 'max OD', 'lag time']
                        ylabels= ['log(OD)', 'growth rate']
                        logs= True
                    else:
                        esterrs= S[c][s][dtype + 'var']
                        snames= ['max derivative', 'time of max derivative',
                                 'inverse of max derivative',
                                 'max ' + dtype, 'lag time']
                        ylabels= [dtype, 'derivative of ' + dtype]
                        logs= False
                    f= fitderiv(S[c][s]['time'], d, figs= False, cvfn= cvfn, logs= logs,
                                bd= bd, esterrs= esterrs, statnames= snames, noruns= noruns,
                                exitearly= False, linalgmax= 5)
                    plt.figure()
                    plt.subplot(2,1,1)
                    f.plotfit('f', ylabel= ylabels[0], figtitle= figtitle + ' : growth rate')
                    plt.subplot(2,1,2)
                    f.plotfit('df', ylabel= ylabels[1])
                    plt.show(block= False)
                    if dtype == 'OD':
                        if plotodgr:
                            gu.plotxyerr(np.exp(f.f), f.df,
                                        (np.exp(f.f + np.sqrt(f.fvar)) - np.exp(f.f - np.sqrt(f.fvar)))/2.0,
                                        np.sqrt(f.dfvar), 'OD', 'growth rate', figtitle + ' : growth rate vs OD')
                        S[c][s]['flogOD']= f.f
                        S[c][s]['flogODvar']= f.fvar
                        S[c][s]['OD ratio']= np.exp(S[c][s]['flogOD'][-1] - S[c][s]['flogOD'][0])
                        S[c][s]['gr']= f.df
                        S[c][s]['grvar']= f.dfvar
                        S[c][s]['d/dtgr']= f.ddf
                        S[c][s]['d/dtgrvar']= f.ddfvar
                        # find local maximum derivative
                        from scipy.signal import argrelextrema
                        tpts= argrelextrema(f.df, np.greater)[0]
                        if np.any(tpts):
                            S[c][s]['local max growth rate']= np.max(f.df[tpts])
                            S[c][s]['time of local max growth rate']= S[c][s]['time'][tpts[np.argmax(f.df[tpts])]]
                        else:
                            S[c][s]['local max growth rate']= np.nan
                            S[c][s]['time of local max growth rate']= np.nan
                        # find area under gr vs OD
                        from scipy import integrate, interpolate
                        sod= np.exp(f.f)
                        integrand= lambda x: interpolate.interp1d(sod, f.df)(x)
                        S[c][s]['area under gr vs OD']= integrate.quad(integrand, np.min(sod), np.max(sod),
                                                                       limit= 100, full_output= 1)[0]
                    else:
                        S[c][s]['f' + dtype]= f.f
                        S[c][s]['f' + dtype + 'var']= f.fvar
                        S[c][s]['d/dt' + dtype]= f.df
                        S[c][s]['d/dt' + dtype + 'var']= f.dfvar
                        S[c][s]['d2/dt2' + dtype]= f.ddf
                        S[c][s]['d2/dt2' + dtype + 'var']= f.ddfvar
                    if stats:
                        for sname in f.ds.keys(): S[c][s][sname]= f.ds[sname]
                except KeyError:
                    print(dtype, 'is either not yet calculated or not recognized')
                    return
            print('\nProcessed strains now have these keys:')
            self.getkeys(c, s, title= False)

    #####
    def comparegrarea(self, ref, com, figs= True):
        '''
        Calculates the area between two growth rate versus OD curves, normalized by the length
        along the OD-axis where they overlap:

        e.g. p.comparegrarea(['1.9% Raffinose 0.0µg/ml cycloheximide', '77.WT'],
                             ['1.9% Raffinose 0.5µg/ml cycloheximide', '77.WT'])

        Arguments
        --
        ref: [condition, strain] array for the reference
        com: [condition, strain] array to compare to the reference
        '''
        # extract data
        od0= np.exp(self.d[ref[0]][ref[1]]['flogOD'])
        gr0= self.d[ref[0]][ref[1]]['gr']
        od1= np.exp(self.d[com[0]][com[1]]['flogOD'])
        gr1= self.d[com[0]][com[1]]['gr']
        # interpolate data
        from scipy import interpolate
        i0= interpolate.interp1d(od0, gr0)
        i1= interpolate.interp1d(od1, gr1)
        # find common OD range
        odinit= np.max([od0[0], od1[0]])
        odfin= np.min([od0[-1], od1[-1]])
        if figs:
            plt.figure()
            plt.plot(od0, gr0, '.-', od1, gr1, '.-')
            x= np.linspace(odinit, odfin, np.max([len(od0), len(od1)]))
            plt.fill_between(x, i0(x), i1(x), facecolor= 'red')
            plt.xlabel('OD')
            plt.ylabel('growth rate')
            plt.legend([ref[0] + ' in ' + ref[1], com[0] + ' in ' + com[1]],
                       loc= 'upper left', bbox_to_anchor= (0.5, 1.05))
            plt.title('area between ' + ref[1] + ' in ' + ref[0] + ' and ' + com[1] + ' in ' + com[0])
            plt.show(block= False)
        # perform integration
        from scipy import integrate
        igrand= lambda x: np.abs(i0(x) - i1(x))
        a= integrate.quad(igrand, odinit, odfin, limit= 100, full_output= 1)[0]
        # return normalized area between curves
        return a/(odfin-odinit)





    #####
    # Exporting
    #####
    def save(self, savename= False):
        '''
        Use pickle to save current instance.

        Arguments
        --
        savename: the name of the output file
        '''
        import pickle
        if not savename:
            savename= self.wdir + self.dname.split('.')[0] + '.pkl'
        pickle.dump(self, open(savename, 'wb'))


    #####
    def makedataframe(self, dnames, conditions= 'all', strains= 'all', dfname= 'df', nomedia= True, nonull= True,
                      conditionincludes= False, strainincludes= False, conditionexcludes= False, strainexcludes= False):
        '''
        Makes a data frame (from the Python panda package) and adds as .df to the current instance.
        For example, either for variables that are measured at each time point (a time column is automatically added)

        p.makedataframe(['gr', 'ODmn'])
        p.df.T

        or for other variables

        p.makedataframe(['max growth rate', 'lag time'])
        p.df

        Arguments
        --
        dnames : list of data fields to add (each must either refer to field with a single number or all fields should reference data that is a function of time)
        conditions : conditions of interest (default: 'all')
        strains : strains of interest (default: 'all')
        dfname : name to use for dataframe (default: 'df')
        nomedia : if True, wells labelled media are excluded
        nonull : if True, wells labelled null are excluded
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        S= self.d
        dnames= gu.makelist(dnames)
        cs= self.getcons(conditions, nomedia, conditionincludes, conditionexcludes)
        # test if data fields are variables that change with time
        timedata= True
        for dname in dnames:
            trialdata= S[cs[0]][self.getstrains(strains, cs[0], True, strainincludes,
                                                strainexcludes)[0]][dname]
            if type(trialdata) != np.ndarray or len(trialdata) != len(self.t):
                timedata= False
                break
        if timedata:
            # data fields are all a function of time
            data= {}
            data['Time (hrs)']= self.t
            for c in cs:
                for s in self.getstrains(strains, c,nonull, strainincludes, strainexcludes):
                    for dname in dnames:
                        try:
                            # multi-dimensional array
                            for i in range(S[c][s][dname].shape[1]):
                                data[dname +' for '+ s+' in '+ c+' replicate '+ str(i)]= S[c][s][dname][:,i]
                        except IndexError:
                            # 1-dimensional array
                            data[dname+' for '+ s+' in '+c]= S[c][s][dname]
            df= pd.DataFrame(data)
        else:
            # data fields are not a function of time
            snames, data= [], []
            for c in cs:
                for s in self.getstrains(strains, c, nonull, strainincludes, strainexcludes):
                    snames.append(s + ' in ' + c)
                    data.append([S[c][s][dname] for dname in dnames])
            df= pd.DataFrame(data, index= snames, columns= dnames)
            df.sort_index(inplace= True)
        if hasattr(self, dfname): print('Redefining', dfname)
        setattr(self, dfname, df)


    def rmdataframe(self, dfname= 'df'):
        '''
        Delete a dataframe.

        Arguments
        --
        dfname: dataframe to be deleted (default: 'df')
        '''
        if hasattr(self, dfname): delattr(self, dfname)


    def exportdataframe(self, dfname= 'df', ftype= 'xlsx', savename= False):
        '''
        Export a data frame to Excel.

        Arguments
        --
        dfname : name of dataframe to export
        ftype : 'xlsx' (default) or 'csv'
        savename : name of file to be saved (default is to add _df to the name of the data set)
        '''
        if hasattr(self, dfname):
            if not savename:
                savename= self.wdir + self.dname.split('.')[0] + '_' + dfname + '.' + ftype
            if ftype == 'xlsx':
                print('Writing Excel file:', savename)
                getattr(self, dfname).to_excel(savename, sheet_name= 'Sheet1', index= False)
            else:
                print('Writing CSV file:', savename)
                getattr(self, dfname).to_csv(savename, sep= ',', index= False)
        else:
            print('No dataframe', dfname, 'has been created')


##########################
# specialized platereader
##########################


class slpr(platereader):

    def __init__(self, dname, aname= 'default', prtype= 'Tecan', wdir= '', dsheetname= 0, asheetname= 0,
             ODfname= 'ODcorrection.txt', warn= False, info= True, standardgain= False):

        super().__init__(dname, aname, prtype, wdir, dsheetname, asheetname, ODfname,
                             warn, info, standardgain)
        self.__doc__= super().__doc__



#####

if __name__ == '__main__': print(platereader.__doc__)
