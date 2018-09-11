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
    p.plot('c-GFPperod', onefig= True)
    p.savefigs()
    p.getstats('OD', '1% Gal', 'GAL2')
    p.getstats('FLperod', '1% Gal', 'GAL2')

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

    p.makedateframe('timecol', dfname= 'df')
    p.exportdataframe(ftype= 'csv')
    p.rmdataframe(dfname= 'df')

    Bounds on the hyperparameters for all Gaussian processes can be specified. For example,

    p.correctauto(bd= {1: (-6, -4)})

    changes the bounds on the second hyperparameter for the Gaussian process used to fit data for the reference strain to be [10^-6, 10^-4]. This hyperparameter describes the flexibility of the fitted curve.

    As a citation for this software, please use:

    Lichten CA, White R, Clark IBN, Swain PS. Unmixing of fluorescence spectra to resolve quantitative time-series measurements of gene expression in plate readers. BMC Biotechnol 14 (2014) 11

    The technique for estimating growth rates is given in:

    Swain PS, Stevenson K, Leary A, Montano-Gutierrez LF, Clark IB, Vogel J, Pilizota T. Inferring time derivatives including cell growth rates using Gaussian processes. Nat Commun 7 (2016) 13766

    '''

    #####
    def __init__(self, dname, aname= 'default', platereadertype= 'Tecan', wdir= '', dsheetnumber= 0,
                 asheetnumber= 0, ODfname= 'ODcorrection_Glucose_Haploid.txt', warn= False, info= True,
                 standardgain= False, ignorestrains= []):
        '''
        Requires a data file, fname, (.csv, .txt, .xls, or .xlsx). An annotation file, aname, giving the contents of each well is optional. The plate reader is assumed to be a Tecan machine (but alternatives can in principle be specified using the platereader argument). A working directory can be set using wdir, such as wdir= 'data/'.

        Arguments
        --
        dname: file name of data (if no extension is given, '.xlsx' is assumed)
        aname: file name for annotation (if 'default', assumed to be the root of dname + '_contents.xlsx)
        platereadertype: type of plate reader ('Tecan')
        dsheetnumber: specifies the sheet of the data Excel file
        asheetnumber: specifies the sheet of the annotation Excel file
        wdir: name of directory where the data files are stored
        ODfname: file name for dilution data for corrected OD (default is 'ODcorrection.txt')
        warn: if False (default), warnings created by covariance matrices that are not positive semi-definite are stopped
        info: if True (default), display information on the plate after loading
        standardgain: if defined, fluorescence measurements are corrected to this gain
        ignorestrains: a list of strains in particular conditions to be ignored, such as ['GAL10 in 0.1% Gal']
        '''
        self.version= '4.84'
        self.dsheetnumber= dsheetnumber
        self.asheetnumber= asheetnumber
        if '.' not in dname: dname += '.xlsx'
        # specify working directory
        self.wdir= wdir
        self.name= dname.split('.')[0]
        self.ODfname= ODfname
        self.aname=aname
        # general parameters
        self.gamma= 0.114   # ratio of 585 to 525 for eGFP
        self.nosamples= 100  # for estimating error through sampling
        self.consist= 2   # number of stds of corrected reference strain data that is considered measurement noise
        self.overflow= -999.99

        # correction has not been performed
        self.standardgain= standardgain
        self.ODcorrected= False
        self.processedref1= {}
        self.processedref2= {}
        self.autocorrected= {}
        self.mediacorrected= {}
        self.ignoredwells= []
        self.ignored= []
        self.negativevalues= False
        self.mediaGP= {}

        if not warn:
            # warning generated occasionally when sampling from the Gaussian process likely because of numerical errors
            import warnings
            warnings.simplefilter("ignore", RuntimeWarning)

        print('Platereader', self.version, ': loading data')
        # define annotation
        if aname == 'default': aname= self.name + '_contents.xlsx'
        if aname:
            # import annotation
            try:
                r= pd.ExcelFile(self.wdir + aname).parse(asheetnumber)
            except FileNotFoundError:
                raise(SystemExit("\nError: Can't find " + self.wdir + aname))
        else:
            # no annotation: use wells
            r= [[a + str(j) + ' in unspecified' for j in np.arange(1,13)] for a in 'ABCDEFGH']
            r= pd.DataFrame(r, index= [a for a in 'ABCDEFGH'], columns= [j for j in np.arange(1,13)])


        # import data
        df, dlabels, gains, datatypes, t, idstartindex, dstartindex= self.importdata(platereadertype,
                                                                                     dname, dsheetnumber)
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
        for let in 'ABCDEFGH':
            for no in np.arange(1,13):
                plateloc= let + str(no)
                try:
                    alabel= r[no][let]
                    if ' in ' in alabel:
                        # analyse annotation
                        ells= alabel.split('in')
                        strain= ells[0].strip()
                        condition= 'in'.join(ells[1:]).strip()
                        sinc= strain + ' in ' + condition
                        if sinc not in ignorestrains:
                            alldata.append(sinc)
                            platelabels[plateloc]= sinc
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
                                print('Data missing for well', plateloc)
                        else:
                            print('Ignoring', sinc, 'at well', plateloc)
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
            self.mediaGP[dn]= {}
            for c in allconditions:
                self.mediacorrected[c + ' for ' + dn]= False
                if dn not in ['OD', 'AutoFL']:
                    self.autocorrected[c + ' for ' + dn]= False
                    self.processedref1[c + ' for ' + dn]= False
                    self.processedref2[c + ' for ' + dn]= False
        self.nooutchannels= len(datatypes)
        self.extractdata()
        self.platelabels= platelabels
        self.allstrains= allstrains
        self.allconditions= allconditions
        self.alldata= list(np.unique(alldata))
        self.importtime= '{:%Y-%b-%d %H:%M:%S}'.format(datetime.datetime.now())
        if info: self.info(full= False)


    #####
    def importdata(self, platereadertype, dname, dsheetnumber):
        '''
        Internal function: Creates and parses dataframe from plate reader Excel file.
        Returns:
            df: dataframe (arranged in rows by well label; data in form of time-series of OD then GFP, etc.)
            dlabels: used to find location of well in data (row number)
            gains: gains used
            datatypes: data types present
            t: time
            idstartindex: index in row for label of well
            dstartindex: index in row where data starts
        '''
        try:
            xl= pd.ExcelFile(self.wdir + dname)
        except FileNotFoundError:
            raise(SystemExit("\nError: Can't find " + self.wdir + dname))
        df= xl.parse(dsheetnumber)
        if platereadertype == 'Tecan':
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
        elif platereadertype == 'Sunrise':
            # Sunrise plate reader
            datatypes= ['OD']
            gains= ['']
            # rearrange data
            tdic= {}
            for well in [a + str(j) for j in np.arange(1,13) for a in 'ABCDEFGH']:
                tdic[well]= []
            for index, row in df.iterrows():
                for ir, rd in enumerate(row[1:]):
                    if row[0] in 'ABCDEFGH':
                        tdic[row[0] + str(ir+1)].append(rd)
            df= pd.DataFrame(tdic).transpose()
            dlabels= np.array(df.index.tolist())
            # time
            if not hasattr(self, 't'):
                # not specified by Magellan: assign default
                print('Warning: could not find any measurements of time')
                t= np.arange(0, df.shape[1])
            else:
                t= self.t
            # indices to extract data
            idstartindex= 0
            dstartindex= 0
        elif platereadertype == 'Hidex':
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
        for c in self.getconditions('all'):
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
        for c in self.getconditions('all'):
            for s in self.getstrains('all', c):
                for dn in self.datatypes:
                    if S[c][s][dn].size == 0:
                        # all wells for a condition and strain have been ignored
                        self.ignored.append(s + ' in ' + c)
                    elif S[c][s][dn].shape[1] > 1:
                        S[c][s][dn + ' mean']= np.mean(S[c][s][dn], 1)
                        S[c][s][dn + ' var']= np.var(S[c][s][dn], 1)
                    else:
                        S[c][s][dn + ' mean']= gu.makerow(S[c][s][dn])
                        S[c][s][dn + ' var']= np.zeros(len(S[c][s][dn]))
        if self.ignored:
            self.ignored= list(np.unique(self.ignored))
            for igns in self.ignored:
                print('Warning: all wells for', igns, 'have been ignored')
            raise(SystemExit('\nYou should reload the data specifying "ignorestrains" because all wells for a strain in a particular condition have been ignored.'))


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
            self.ignoredwells= []
        else:
            # wells cannot be ignored twice
            exclude= list(set(exclude) - set(self.ignoredwells))
        # set plateloc to initial values if necessary
        for c in self.getconditions('all'):
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
    def info(self, full= False):
        '''
        Displays conditions, strains, datatypes, and corrections made.

        Arguments
        --
        full: if True, display information on correcting for media and autofluorescence
        '''
        print('\n' + self.name)
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
        if full:
            print('---')
            print('Corrected for nonlinearities in OD:', self.ODcorrected)
            for dn in self.datatypes:
                print('---')
                for c in self.allconditions:
                    print(dn, 'in', c, 'corrected for media:', self.mediacorrected[c + ' for ' + dn])
                if dn != 'OD' and dn != 'AutoFL':
                    print('---')
                    for c in self.allconditions:
                        print(dn, 'corrected for autofluorescence in', c, ':', self.autocorrected[c + ' for ' + dn])


    def attributes(self):
        '''
        Displays the names of the attributes in the instance of platereader.

        Arguments
        --
        None
        '''
        ignore= ['d', 'consist', 't', 'nosamples', 'gamma', 'ODfname', 'overflow', 'nooutchannels', 'nodata', '__doc__']
        for a in self.__dict__:
            if 'corrected' not in a and 'processed' not in a and a not in ignore: print(a)


    #####
    def datavariables(self, individual= False, conditions= 'all', strains= 'all',
                      conditionincludes= False, strainincludes= False,
                      conditionexcludes= False, strainexcludes= False,
                      title= True, out= False):
        '''
        Returns the keys (the variables) calculated from the data for each strain.

        Arguments
        --
        individual: if True, separately list the variables for each strain
        conditions: list of experimental conditions to include (default is 'all')
        strains: list of strains to include (default is 'all')
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        title: display condition and strain
        out: if True, return the names of datavariables, timevariables, and othervariables
        '''
        S= self.d
        ignore= ['data', 'originalplateloc']
        if individual:
            for c in self.getconditions(conditions, True, conditionincludes, conditionexcludes):
                for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                    if title: print(s, 'in', c, '\n---')
                    for k in sorted(list(S[c][s].keys())):
                        if k not in ignore: print('\t'+k)
        else:
            datavariables= np.unique([var if var not in ignore else 'time'
                                      for c, s in self.getconditionsandstrains(conditions,
                                                                               strains,
                                                                               conditionincludes,
                                                                               strainincludes,
                                                                               conditionexcludes,
                                                                               strainexcludes,
                                                                               True, True)
                                      for var in S[c][s]])
            timevariables= np.unique([var if isinstance(S[c][s][var], np.ndarray) and
                                      S[c][s][var].shape[0] == S[c][s]['time'].shape[0]
                                      and var not in ignore
                                      else 'time'
                                      for c, s in self.getconditionsandstrains(conditions,
                                                                               strains,
                                                                               conditionincludes,
                                                                               strainincludes,
                                                                               conditionexcludes,
                                                                               strainexcludes,
                                                                               True, True)
                                      for var in S[c][s]])
            othervariables= np.setdiff1d(datavariables, timevariables)
            if out:
                return datavariables, timevariables, othervariables
            else:
                for dv in datavariables: print('\t'+dv)



    #####
    # Internal functions
    #####
    def getconditions(self, conditions= 'all', nomedia= False, conditionincludes= False,
                      conditionexcludes= False):
        '''
        Returns a list of conditions.

        Arguments
        --
        conditions: list of experimental conditions to include (default is 'all')
        nomedia: ignores condition= media if True (default)
        conditionincludes: selects only conditions with conditionincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        '''
        if conditions == 'all' or conditionincludes or conditionexcludes:
            cons= list(self.d.keys())
            if nomedia and 'media' in cons:
                cons.pop(cons.index('media'))
            # find those conditions containing keywords given in 'includes'
            if conditionincludes:
                conditionincludes= gu.makelist(conditionincludes)
                newcons= []
                for condition in cons:
                    gotone= 0
                    for item in conditionincludes:
                        if item in condition: gotone += 1
                    if gotone == len(conditionincludes): newcons.append(condition)
                cons= newcons
            # remove any conditions containing keywords given in 'excludes'
            if conditionexcludes:
                conditionexcludes= gu.makelist(conditionexcludes)
                exconds= []
                for condition in cons:
                    for item in conditionexcludes:
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
            if conditionincludes:
                raise(SystemExit('No conditions have ' + ' and '.join(conditionincludes)))
            else:
                raise(SystemExit('No conditions found'))



    #####
    def getstrains(self, strains, condition, nonull= False, strainincludes= False, strainexcludes= False):
        '''
        Returns list of strains

        Arguments
        --
        strains: list of strains to include (default is 'all')
        condition: the experimental condition for which strains are required
        nonull: ignores strain= null if True
        strainincludes: selects only strains with strainincludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        '''
        strainsincondition= list(self.d[condition].keys())
        if strains == 'all' or strainincludes or strainexcludes:
            ss= strainsincondition
            if nonull and 'null' in ss:
                ss.pop(ss.index('null'))
            # find those strains containing keywords given in 'includes'
            if strainincludes:
                strainincludes= gu.makelist(strainincludes)
                newss= []
                for strain in ss:
                    sc= 0
                    for item in strainincludes:
                         if item in strain: sc += 1
                    if sc == len(strainincludes): newss.append(strain)
                ss= newss
            # remove any strains containing keywords given in 'excludes'
            if strainexcludes:
                strainexcludes= gu.makelist(strainexcludes)
                exsts= []
                for strain in ss:
                    for item in strainexcludes:
                        if item in strain:
                            exsts.append(strain)
                            break
                for ex in exsts:
                    ss.pop(ss.index(ex))
        else:
            ss= gu.makelist(strains)
            for s in ss:
                if s not in strainsincondition:
                    ss.pop(ss.index(s))
        if ss:
            return sorted(ss)
        else:
            if strainincludes:
                print('Warning: No strains have ' + ' and '.join(strainincludes) + ' for ' + condition)
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
        return [(c,s) for c in self.getconditions(conditions, nomedia= nomedia,
                                                  conditionincludes= conditionincludes,
                                                  conditionexcludes= conditionexcludes)
                for s in self.getstrains(strains, c, nonull= nonull, strainincludes= strainincludes,
                                         strainexcludes= strainexcludes)]




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
            if correctmedia: self.correctmedia(datatypes= 'OD', figs= figs,
                                               conditionincludes= conditionincludes,
                                               conditionexcludes= conditionexcludes)
            # if no ODfname use the default
            if not ODfname: ODfname= self.wdir + self.ODfname
            S, noout= self.d, self.nooutchannels
            # fit dilution data
            if not hasattr(self, 'gc'): self.findODcorrection(ODfname, figs)
            # correct all wells containing strains
            for c in self.getconditions('all', True, conditionincludes, conditionexcludes):
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
            plt.grid(True)
            plt.xlabel('OD')
            plt.ylabel('relative cell density')
            plt.title('function for correcting OD')
            plt.show(block= False)
        self.gc= gc



    #####
    # Media correction
    #####
    def correctmedia(self, datatypes= 'all', conditions= 'all', figs= True, noruns= 3, exitearly= False, bd= False,
                     results= False, conditionincludes= False, conditionexcludes= False, mean= False,
                     commonmedia= False):
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
        commonmedia: condition containing a 'null' strain that is to be used to correct multiple other conditions
        '''
        S= self.d
        cons= self.getconditions(conditions, False, conditionincludes, conditionexcludes)
        if datatypes == 'all': datatypes= self.datatypes
        datatypes= gu.makelist(datatypes)
        for dn in datatypes:
            # correct for media
            for c in cons:
                if c != 'media':
                    if self.mediacorrected[c + ' for ' + dn]:
                        print(dn, 'in', c, 'is already corrected for the media')
                    elif commonmedia and 'null' in S[commonmedia]:
                        print('Correcting', dn, 'in', c, 'for media')
                        self.performmediacorrection(dn, c, figs, noruns, exitearly, bd, results, mean,
                                                    commonmedia)
                        self.mediacorrected[c + ' for ' + dn]= True
                    elif 'null' in S[c] or ('media' in cons and 'null' in S['media']) :
                        print('Correcting', dn, 'in', c, 'for media')
                        self.performmediacorrection(dn, c, figs, noruns, exitearly, bd, results, mean, c)
                        self.mediacorrected[c + ' for ' + dn]= True
                    else:
                        print('No well annotated "null" was found')
                        print('Correcting for media abandoned for', dn, 'in', c)
        self.updatemeans()
        if self.negativevalues:
            print('Warning: correcting media has created negative values for')
            print(self.negativevalues)

    #####
    def performmediacorrection(self, dtype, condition, figs, noruns, exitearly, bd, results, mean, commonmedia):
        '''
        Internal function: Uses a Gaussian process to fit the media over time and subtracts the best-fit media values from the data, as well as storing the uncorrected data.
        '''
        S= self.d
        if type(condition) is list:
            # `media' is specified for all conditions
            cons= condition[1]
            condition= condition[0]
        else:
            # data avaliable to correct per condition: 'null in c'
            cons= condition
        # find correction to use (either for 'media' and so all conditions or just for specific condition)
        if 'media' in self.getconditions():
            t, data= S['media']['null']['time'], S['media']['null'][dtype]
        else:
            t, data= S[commonmedia]['null']['time'], S[commonmedia]['null'][dtype]
        if mean:
            f= np.mean(data)
        elif commonmedia in self.mediaGP[dtype]:
            # Gaussian process has already been fit
            f, fvar= self.mediaGP[dtype][commonmedia].f, self.mediaGP[dtype][commonmedia].fvar
        else:
            # d is data and x is time
            d= data.flatten('F')
            x= np.tile(t, data.shape[1])
            b= {0: (-6,4), 1: (-5,-2), 2: (-10,0)}
            if bd: b= gu.mergedicts(original= b, update= bd)
            # fit with Gaussian process
            try:
                print('Fitting media for', dtype, 'in', commonmedia)
                f, fvar, self.mediaGP[dtype][commonmedia]= gu.smoothGP(x, d, xp= t, bd= b, results= results)
                if figs:
                    plt.figure()
                    plt.plot(x, d, 'ro', t, f, 'b-')
                    plt.xlabel('time (hours)')
                    plt.title('media correction for ' + dtype + ' in ' + condition)
                    plt.show(block= False)
            except gp.gaussianprocessException:
                raise(SystemExit('Fitting media failed'))
        # perform correction
        for c in gu.makelist(cons):
            for s in self.getstrains('all', c, nonull= True):
                sh= S[c][s][dtype].shape
                S[c][s][dtype] -= np.reshape(np.tile(f, sh[1]), (-1, sh[1]), order= 'F')
                if np.any(S[c][s][dtype] < 0):
                    wstr= '\t' + dtype + ': ' + s + ' in ' + c + '\n'
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
                    correctmedia= True, mediausemean= False, ignoreneg= False, minqrerr= 1.0e-6):
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
        ignoreneg: if True (default: False), proceed with correction despite negative fluorescence values
        minqrerr: minimum value allowed for the estimated error in the ratio of fluorescence (AutoFL/GFP) - too small values can cause instabilities in the fitting
        '''
        f= gu.makelist(f)
        S= self.d
        # check have enough replicates
        if len(f) == 2:
            go= True
            for c in self.getconditions(conditions, True, conditionincludes, conditionexcludes):
                for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                    if S[c][s]['OD'].shape[1] < 2:
                        print('Not enough OD replicates to correct autofluorescence for', s, 'in', c)
                        go= False
                    if S[c][s][f[0]].shape[1] < 2:
                        print('Not enough fluorescence replicates to correct autofluorescence for', s, 'in', c)
                        go= False
            if go == False:
                print('Try specifying just one fluorescence measurement')
                return
        # correct OD
        if correctOD:
            self.correctOD(correctmedia= correctmedia, conditionincludes= conditionincludes,
                           conditionexcludes= conditionexcludes)
        # correct for media
        if correctmedia:
            self.correctmedia(f, conditions, conditionincludes= conditionincludes,
                              conditionexcludes= conditionexcludes, results= results,
                              mean= mediausemean)
            if self.negativevalues:
                print('There are negative values for fluorescence measurements')
                print('Use\n\tignoreneg= True\nto force correction for autofluorescence')
                print('Alternatively, specify one fluorescence measurement, such as:')
                print("\tp.reset()\n\tp.correctauto('GFP', correctmedia= False)")
                if ignoreneg:
                    print('Going ahead...')
                else:
                    return
        # correct autofluorescence
        print('Using', refstrain, 'as the reference')
        if len(f) == 1:
            self.correctauto1(f, conditions, strains, refstrain, correctOD, figs, noruns, bd, no1samples,
                              conditionincludes, strainincludes, conditionexcludes, strainexcludes,
                              results)
        elif len(f) == 2:
            self.correctauto2(f, conditions, strains, refstrain, correctOD, figs, noruns, bd,
                              conditionincludes, strainincludes, conditionexcludes, strainexcludes,
                              results, minqrerr)
        else:
            print('f must be a list of length 1 or 2')



    #####
    def correctauto1(self, f, conditions, strains, refstrain, correctOD, figs, noruns, bd, nosamples,
                     conditionincludes, strainincludes, conditionexcludes, strainexcludes, results):
        '''
        Internal function: Corrects for autofluorescence for experiments with measured emissions at one wavelength using the fluorescence of the wild-type interpolated to the OD of the tagged strain.
        '''
        S= self.d
        # correct autofluorescence
        print('Correcting autofluorescence')
        # run through all conditions
        for c in self.getconditions(conditions, True, conditionincludes, conditionexcludes):
            # process reference strain
            if not self.processedref1[c + ' for ' + f[0]]:
                self.processref1(f, c, refstrain, figs, noruns, bd, conditionincludes, conditionexcludes,
                                 results)
            else:
                print('Reference strain is already processed for', f[0], 'in', c)
            gfr= S[c][refstrain]['gp1ref for ' + f[0]]
            for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                if s != refstrain:
                    noreps= S[c][s][f[0]].shape[1]
                    fl, fl2, flperod, flperod2= (0,)*4
                    # run over all replicates
                    for i in range(noreps):
                        # interpolate WT fluorescence errors to OD values of strain
                        from scipy.interpolate import interp1d
                        interpf= interp1d(S[c][refstrain]['OD mean'], S[c][refstrain]['fmerr for ' + f[0]])
                        try:
                            men= interpf(S[c][s]['OD'][:,i])
                        except ValueError:
                            men= np.median(S[c][refstrain]['fmerr for ' + f[0]])*np.ones(len(S[c][refstrain]['fmerr for ' + f[0]]))
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
                    S[c][s][bname + ' var']= fl2/(noreps*nosamples) - S[c][s][bname]**2
                    S[c][s][bname + 'perod']= flperod/(noreps*nosamples)
                    S[c][s][bname + 'perod var']= flperod2/(noreps*nosamples) - S[c][s][bname + 'perod']**2
            self.autocorrected[c + ' for ' + f[0]]= True
        print('Created:\n' + '\tc-' + f[0] + 'perod (corrected fluorescence per cell)\n' + '\tc-' + f[0]
              + ' (corrected fluorescence)')



    #####
    def processref1(self, f, conditions, refstrain, figs, noruns, bd, conditionincludes, conditionexcludes,
                    results):
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
        for c in self.getconditions(conditions, True, conditionincludes, conditionexcludes):
            print('Processing reference strain for', f[0], 'in', c)
            # fit reference strain's fluorescence as a function of OD
            try:
                x= S[c][refstrain]['OD'].flatten('F')
                y= S[c][refstrain][f[0]].flatten('F')
            except KeyError:
                print(refstrain, 'not found in', c)
                raise(SystemExit('Running correctauto failed'))
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
            # store reference strain information
            S[c][refstrain]['gp1ref' + ' for ' + f[0]]= gfr
            S[c][refstrain]['fmerr' + ' for ' + f[0]]= me
            self.processedref1[c + ' for ' + f[0]]= True
            for s in self.getstrains('all', c, True):
                S[c][s]['refstrain']= refstrain




    #####
    def correctauto2(self, f, conditions, strains, refstrain, correctOD, figs, noruns, bd, conditionincludes,
                     strainincludes, conditionexcludes, strainexcludes, results, minqrerr):
        '''
        Internal function: Corrects for autofluorescence using spectral unmixing for experiments with measured emissions at two wavelengths (following Lichten et al.)
        '''
        S, noout= self.d, self.nooutchannels
        # correct for autofluorescence
        print('Correcting autofluorescence')
        for c in self.getconditions(conditions, True, conditionincludes, conditionexcludes):
            # process reference strain
            if not self.processedref2[c + ' for ' + f[0]]:
                self.processref2(f, conditions, refstrain, figs, noruns, bd, conditionincludes,
                                 conditionexcludes, results, minqrerr)
            else:
                print('Reference strain is already processed for', f[0], 'in', c)
            gr= S[c][refstrain]['gp2ref' + ' for ' + f[0]]
            for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                if s != refstrain:
                    nodata, noreps= S[c][s][f[0]].shape
                    # remove autofluorescence
                    fl= self.applyautoflcorrection(S[c][refstrain]['fratio for ' + f[0]],
                                                   S[c][s][f[0]], S[c][s][f[1]])
                    # estimate error
                    varf= np.var(S[c][s][f[0]], 1)
                    varcf= np.var(S[c][s][f[1]], 1)
                    if np.array_equal(S[c][refstrain]['time'], S[c][s]['time']):
                        gr.predict(S[c][refstrain]['time'], merrorsnew= S[c][refstrain]['qrerr for ' + f[0]])
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
                    S[c][s][bname + ' var']= varfl
                    # fluorescence per cell
                    S[c][s][bname + 'perod']= S[c][s][bname]/S[c][s]['OD mean']
                    flperodvar= np.zeros(nodata)
                    for i in np.random.randint(noreps, size= self.nosamples):
                        if correctOD:
                            self.gc.predict(S[c][s]['OD'][:,i])
                            cc= self.gc.sample(1)
                            cc= cc.flatten('F')
                            flperodvar += S[c][s][bname + ' var']/cc**2
                        else:
                            flperodvar += S[c][s][bname + ' var']/S[c][s]['OD'][:,i]**2
                    flperodvar /= self.nosamples
                    S[c][s][bname + 'perod var']= flperodvar
                    # calculate corrected levels that are above the expectation from reference strain
                    reflevel= self.consist*np.sqrt(S[c][refstrain][bname + ' var'])
                    keep= np.nonzero(S[c][s][bname] > reflevel)[0]
                    S[c][s]['s' + bname]= np.zeros(np.size(S[c][s][bname]))
                    S[c][s]['s' + bname][keep]= S[c][s][bname][keep]
                    S[c][s]['s' + bname + ' var']= np.zeros(np.size(S[c][s][bname + ' var']))
                    S[c][s]['s' + bname + ' var'][keep]= S[c][s][bname + ' var'][keep]
                    reflevel /= S[c][refstrain]['OD mean']
                    keep= np.nonzero(S[c][s][bname+ 'perod'] > reflevel)[0]
                    S[c][s]['s' + bname + 'perod']= np.zeros(np.size(S[c][s][bname + 'perod']))
                    S[c][s]['s' + bname + 'perod'][keep]= S[c][s][bname + 'perod'][keep]
                    S[c][s]['s' + bname + 'perod var']= np.zeros(np.size(S[c][s][bname + 'perod var']))
                    S[c][s]['s' + bname + 'perod var'][keep]= S[c][s][bname + 'perod var'][keep]
            self.autocorrected[c + ' for ' + f[0]]= True
        print('Created:\n' + '\tc-' + f[0] + 'perod (corrected fluorescence per cell)\n' + '\tc-' + f[0]
              + ' (corrected fluorescence)\n' + '\tsc-' + f[0]
              + 'perod (only statistically significant corrected fluorescence per cell)\n' + '\tsc-'
              + f[0] + ' (only statistically significant corrected fluorescence)')

    #####
    def processref2(self, f, conditions, refstrain, figs, noruns, bd, conditionincludes, conditionexcludes,
                    results, minqrerr):
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
        minqrerr: if values for the estimated error in the fluorescence ratio fall below this value replace by this minimum value
        '''
        S= self.d
        # bounds for Gaussian process for fitting reference strain
        b= {0: (-5,3), 1: (-4,-1), 2: (-4, 4)}
        if bd: b= gu.mergedicts(original= b, update= bd)
        # run through all conditions
        for c in self.getconditions(conditions, True, conditionincludes, conditionexcludes):
            print('Processing reference strain for', f[0], 'in', c)
            try:
                noreps= S[c][refstrain][f[0]].shape[1]
            except KeyError:
                print(refstrain, 'not found in', c)
                raise(SystemExit('Running correctauto failed'))
            f0, f1, t= S[c][refstrain][f[0]], S[c][refstrain][f[1]], S[c][refstrain]['time']
            # find values with zero variance
            dels= np.unique(np.nonzero(np.var(f1/f0, 1) == 0)[0])
            # remove offending values
            f0r= np.delete(f0, dels, axis= 0)
            f1r= np.delete(f1, dels, axis= 0)
            tr= np.delete(t, dels)
            # find ratio of fluorescence
            qr= (f1r/f0r).flatten('F')
            # find error
            qrerr= np.var(f1r/f0r, 1)
            # check no errors too small
            if np.min(qrerr) < minqrerr:
                print('Warning: replacing small estimates for the error in the fluorescence ratio')
                qrerr[qrerr < minqrerr]= minqrerr
            # sort
            x= np.tile(tr, noreps)
            xs= np.sort(x)
            qrs= qr[np.argsort(x)]
            qrerrs= np.repeat(qrerr, noreps)[np.argsort(x)]
            # fit with Gaussian process
            gr= gp.sqexpGP(b, xs, qrs, merrors= qrerrs)
            try:
                # optimize hyperparameters of Gaussian process
                gr.findhyperparameters(noruns)
                # add back missing time points
                if np.any(dels): qrerr= np.interp(t, tr, qrerr)
                gr.predict(t, merrorsnew= qrerr)
            except gp.gaussianprocessException:
                raise(SystemExit('Fitting reference strain failed'))
            # store results
            S[c][refstrain]['(xs, qrs, qrerrs) for ' + f[0]]= (xs, qrs, qrerrs)
            S[c][refstrain]['qrerr for ' + f[0]]= qrerr
            S[c][refstrain]['gp2ref for ' + f[0]]= gr
            S[c][refstrain]['fratio for ' + f[0]]= gr.f
            if figs:
                # plot fit
                plt.figure()
                gr.sketch('o')
                plt.xlabel('time (hours)')
                plt.ylabel(f[1] + '/' + f[0])
                plt.title('fitting ' + refstrain + ' for ' + c)
                plt.show(block= False)
            # check autofluorescence correction for reference strain
            flref= self.applyautoflcorrection(S[c][refstrain]['fratio for ' + f[0]], S[c][refstrain][f[0]],
                                              S[c][refstrain][f[1]])
            bname= 'c-' + f[0]
            S[c][refstrain][bname]= np.mean(flref, 1)
            S[c][refstrain][bname + 'perod']= np.mean(flref/S[c][refstrain]['OD mean'][:,None], 1)
            S[c][refstrain][bname + ' var']= np.var(flref, 1)
            S[c][refstrain][bname + 'perod var']= np.var(flref/S[c][refstrain]['OD mean'][:,None], 1)
            if figs:
                # plot correction for reference strain
                plt.figure()
                plt.plot(S[c][refstrain]['time'], flref, '.')
                plt.plot(S[c][refstrain]['time'], self.consist*np.sqrt(S[c][refstrain][bname + ' var']), 'r:')
                plt.plot(S[c][refstrain]['time'], -self.consist*np.sqrt(S[c][refstrain][bname + ' var']), 'r:')
                plt.plot(S[c][refstrain]['time'], np.zeros(np.size(S[c][refstrain]['time'])), 'k')
                plt.ylabel('corrected ' + refstrain + ' fluorescence')
                plt.xlabel('time (hours)')
                plt.title(c + ': consistency check for reference strain')
                plt.show(block= False)
            # store reference strain information
            for s in self.getstrains('all', c, True):
                S[c][s]['refstrain']= refstrain
                if s != refstrain:
                    S[c][s][bname + 'refstd']= np.sqrt(S[c][refstrain][bname + ' var'])
                    S[c][s]['refODmn']= S[c][refstrain]['OD mean']
            self.processedref2[c + ' for ' + f[0]]= True


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
        ---
        None
        '''
        self.ODcorrected= False
        for dn in self.datatypes:
            for c in self.allconditions:
                self.mediacorrected[c + ' for ' + dn]= False
                if dn not in ['OD', 'AutoFL']:
                    self.processedref1[c + ' for ' + dn]= False
                    self.processedref2[c + ' for ' + dn]= False
                    self.autocorrected[c + ' for ' + dn]= False
        self.negativevalues= False
        self.ignoredwells= []
        self.extractdata()
        print('Reset: media and fluorescence corrections can now be re-run')


    #####
    # Plotting routines
    #####
    def plot(self, dtype= 'OD', conditions= 'all', strains= 'all', plate= False, onefig= False,
             plotod= True, nonull= False, conditionincludes= False, strainincludes= False,
             conditionexcludes= False, strainexcludes= False, includeref= False, errors= False):
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
        errors: if True, standard deviations over replicates are included for 'mn' dtypes
        '''
        S= self.d
        fldata= False
        if dtype == 'labels':
            plate= True
            oldparams= list(plt.rcParams["figure.figsize"])
            # make figure bigger
            plt.rcParams["figure.figsize"]= [17,17]
        elif 'c-' in dtype:
            fldata= True
            nonull= True
        if plate: onefig= True
        if onefig:
            # one figure for all plots
            plt.figure()
            ax= plt.subplot(111)
            # count number of lines
            nlines= 0
            for c in self.getconditions(conditions, False, conditionincludes, conditionexcludes):
                for s in self.getstrains(strains, c, nonull, strainincludes, strainexcludes):
                    if dtype in S[c][s]:
                        nlines += S[c][s][dtype].ndim
            # define colour map
            colormap= plt.cm.nipy_spectral
            from cycler import cycler
            colors= [colormap(i) for i in np.linspace(0, 1, nlines)]
            ax.set_prop_cycle(cycler('color', colors))
            icol= 0
        # draw all plots
        for c in self.getconditions(conditions, False, conditionincludes, conditionexcludes):
            for s in self.getstrains(strains, c, nonull, strainincludes, strainexcludes):
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
                            plt.text(5, 3, ',\n'.join(s.split(',')) + '\nin\n' + c,
                                     horizontalalignment='center',
                                     color= 'navy', fontsize= 'small')
                        else:
                            # plot data
                            plt.plot(S[c][s]['time'], S[c][s][dtype][:,i], '-')
                            i += 1
                        plt.tick_params(labelbottom= False, labelleft= False)
                        # label well locations
                        for j in range(12):
                            if sindex == j+1: plt.title(j+1)
                        for j, k in enumerate(np.arange(1, 96, 12)):
                            if sindex == k: plt.ylabel('ABCDEFGH'[j] + ' ', rotation= 0)
                else:
                    if (fldata and not includeref and 'refstrain' in list(S[c][s].keys())
                        and s == S[c][s]['refstrain']):
                        # ignore reference strain for fluorescence data
                        continue
                    else:
                        # plot data (not whole plate)
                        if not onefig:
                            plt.figure()
                            ax= plt.subplot(1,1,1)
                        # remove ignored wells for the legend
                        pls= [well for well in list(S[c][s]['plateloc']) if well not in self.ignoredwells]
                        if not errors and not fldata:
                            if S[c][s][dtype].ndim > 1:
                                # plot each replicate
                                for i in range(S[c][s][dtype].shape[1]):
                                    label= pls[i] + ': ' + s + ' in ' + c
                                    ax.plot(S[c][s]['time'], S[c][s][dtype][:,i], '.-', label= label)
                            else:
                                # plot single replicate
                                label= s + ' in ' + c
                                ax.plot(S[c][s]['time'], S[c][s][dtype], '.-', label= label)
                        else:
                            # plot with errors
                            var= S[c][s][dtype[:-5] + ' var'] if 'mean' in dtype else S[c][s][dtype + ' var']
                            if fldata:
                                # corrected fluorescence
                                ax.plot(S[c][s]['time'], S[c][s][dtype], '.', label= s + ' in ' + c)
                                ax.plot(S[c][s]['time'], S[c][s][dtype] + np.sqrt(var), 'k:', alpha=0.4)
                                ax.plot(S[c][s]['time'], S[c][s][dtype] - np.sqrt(var), 'k:', alpha=0.4)
                                # add control to fluorescence data
                                rname= 'c-' + dtype.split('-')[1].split('perod')[0]
                                if rname + 'refstd' in S[c][s]:
                                    if 'perod' in dtype:
                                        ax.plot(S[c][s]['time'], self.consist*S[c][s][rname + 'refstd']/S[c][s]['refODmn'], 'r:')
                                    else:
                                        ax.plot(S[c][s]['time'], self.consist*S[c][s][rname + 'refstd'], 'r:')
                                if plotod and not onefig:
                                    # add OD to plot
                                    ax2= ax.twinx()
                                    ax2.plot(S[c][s]['time'], S[c][s]['OD mean'] , 'DarkOrange', linewidth=5, alpha=0.2)
                            else:
                                # everything else
                                ax.errorbar(S[c][s]['time'], S[c][s][dtype], np.sqrt(var), fmt= '.',
                                            color= colors[icol], capsize= 0, label= s + ' in ' + c)
                                icol += 1
                    if not onefig:
                        # display and add labels for each plot
                        plt.title(dtype + ' of ' + s + ' in ' + c)
                        plt.xlabel('time (hours)')
                        if fldata:
                            ax.set_ylabel(dtype)
                            ax2.set_ylabel('mean(OD)')
                            ax.set_ylim(ymin= 0)
                            ax2.set_ylim(ymin= 0)
                        else:
                            plt.ylabel(dtype)
                            if len(plt.gca().get_legend_handles_labels()[1]) == len(pls):
                                plt.legend(pls, loc= 'lower right')
                            plt.ylim(ymin= 0)
                        plt.show(block= False)
        if onefig:
            # display and add labels for single figure
            if plate:
                if dtype == 'labels':
                    plt.suptitle(self.name, fontsize= 20)
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
        if dtype == 'labels': plt.rcParams["figure.figsize"]= oldparams



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
            gu.figs2pdf(self.wdir + self.name + '.pdf')
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
    def getstats(self, dtype= 'OD', conditions= 'all', strains= 'all', bd= False, cvfn= 'matern',
                 esterrs= False, noruns= 5, stats= True, plotodgr= False, conditionincludes= False,
                 strainincludes= False, conditionexcludes= False, strainexcludes= False,
                 keysmessage= True, figs= True, nosamples= 100):
        '''
        Calls fitderiv.py to estimate the first time-derivate.
        The data calculated by fitderiv is added to the overall data structure. For example, p.d['1% Gal']['GAL2']['flogOD'] gives the fit to the log OD curve if OD data is used.
        For OD data, 'flogOD' (the fitted log OD), 'fODvar' (the variance - error - in the fitted log OD), 'gr' (the inferred growth rate), 'grvar' (the variance in the inferred growth rate), and statistics, such as 'max growth rate' and 'lag time', are added for each strain.
        Equivalent variables are added for other types of data.

        Arguments
        --
        dtype: type of data ('OD' - default, 'GFP', 'c-GFPperod', or 'c-GFP')
        conditions: list of experimental conditions to be included
        strains: list of strains to be included
        bd: can be used to change the limits on the hyperparameters for the Gaussian process used for the fit. For example, p.odstats('1% Gal', 'GAL2', bd= {1: [-2,-2])}) fixes the flexibility to be 0.01
        cvfn: covariance function used for fit, either 'matern' (default) or 'sqexp' or 'nn' or, for example, 'sqexp : matern' to pick the covariance function with the highest maximum likelihood
        esterrs: if True, measurement errors are empirically estimated from the variance across replicates at each time point; if False, the size of the measurement error is fit from the data assuming that this size is the same at all time points
        noruns: number of attempts made for each fit (default is 5), each run is made with random initial estimates of the parameters
        stats: calculate statistics if True
        plotodgr: for OD data, plots growth rate versus log(OD) if True (default)
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        keysmessage: if True, displays keys for the strains that have been processed
        figs: if True (default), shows the fit and inferred derivative
        nosamples: number of samples to calculate errors in statistics for growth curve
        '''
        S, noout= self.d, self.nooutchannels
        linalgmax= 5
        for c in self.getconditions(conditions, True, conditionincludes, conditionexcludes):
            for s in self.getstrains(strains, c, True, strainincludes, strainexcludes):
                figtitle= s + ' in ' + c
                print('\nFitting', dtype, 'for', figtitle)
                try:
                    d= S[c][s][dtype]
                    if dtype == 'OD':
                        snames= ['max gr', 'time of max gr',
                                    'doubling time', 'max OD', 'lag time']
                        ylabels= ['log(OD)', 'gr']
                        logs= True
                    else:
                        esterrs= S[c][s][dtype + ' var']
                        snames= ['max derivative of ' + dtype,
                                 'time of max derivative of ' + dtype,
                                 'inverse of max derivative of ' + dtype,
                                 'max ' + dtype,
                                 'lag time of ' + dtype]
                        ylabels= [dtype, 'derivative of ' + dtype]
                        logs= False
                    if ' : ' in cvfn:
                        # multiple covariance functions
                        cvfns= cvfn.split(' : ')
                        fs= [fitderiv(S[c][s]['time'], d, figs= False, cvfn= cfn,
                                        logs= logs, bd= bd, esterrs= esterrs, statnames= snames,
                                        noruns= noruns, exitearly= False, linalgmax= linalgmax,
                                        nosamples= nosamples)
                             for cfn in cvfns]
                        # pick the covariance function with the highest maximum likelihood
                        ci= np.argmax([f.logmaxlike for f in fs])
                        f= fs[ci]
                        print(' Using', cvfns[ci])
                        if figs:
                            plt.figure()
                            for i, cfn in enumerate(cvfns):
                                ax= plt.subplot(2, len(cvfns), i+1)
                                fs[i].plotfit('f', ylabel= ylabels[0], figtitle= figtitle)
                                ax.set_title(cfn.upper()) if i == ci else ax.set_title(cfn)
                                plt.subplot(2, len(cvfns), i+1+len(cvfns))
                                fs[i].plotfit('df', ylabel= ylabels[1])
                            plt.suptitle(figtitle, fontweight='bold')
                            plt.show(block= False)
                    else:
                        # single covariance function
                        f= fitderiv(S[c][s]['time'], d, figs= False, cvfn= cvfn, logs= logs,
                                    bd= bd, esterrs= esterrs, statnames= snames, noruns= noruns,
                                    exitearly= False, linalgmax= linalgmax, nosamples= nosamples)
                        if figs:
                            plt.figure()
                            plt.subplot(2,1,1)
                            f.plotfit('f', ylabel= ylabels[0], figtitle= figtitle)
                            plt.subplot(2,1,2)
                            f.plotfit('df', ylabel= ylabels[1])
                            plt.show(block= False)
                    if dtype == 'OD':
                        if plotodgr:
                            gu.plotxyerr(np.exp(f.f), f.df,
                                        (np.exp(f.f + np.sqrt(f.fvar)) - np.exp(f.f - np.sqrt(f.fvar)))/2.0,
                                        np.sqrt(f.dfvar), 'OD', 'growth rate', figtitle + ' : growth rate vs OD')
                        S[c][s]['flogOD']= f.f
                        S[c][s]['flogOD var']= f.fvar
                        S[c][s]['gr']= f.df
                        S[c][s]['gr var']= f.dfvar
                        S[c][s]['d/dt gr']= f.ddf
                        S[c][s]['d/dt gr var']= f.ddfvar
                        # extra statistics for OD
                        fs, gs, hs= f.sample(nosamples)
                        # log2 OD ratio
                        da= np.log2(np.exp(fs[-1,:] - fs[0,:]))
                        S[c][s]['log2 OD ratio']= np.mean(da)
                        S[c][s]['log2 OD ratio var']= np.var(da)
                        # find local maximum derivative
                        from scipy.signal import argrelextrema
                        da, dt= [], []
                        for gsample in np.transpose(gs):
                            tpts= argrelextrema(gsample, np.greater)[0]
                            if np.any(tpts):
                                da.append(np.max(gsample[tpts]))
                                dt.append(f.t[tpts[np.argmax(gsample[tpts])]])
                        if np.any(da):
                            S[c][s]['local max gr']= np.mean(da)
                            S[c][s]['local max gr var']= np.var(da)
                            S[c][s]['time of local max gr']= np.mean(dt)
                            S[c][s]['time of local max gr var']= np.var(dt)
                        else:
                            S[c][s]['local max gr']= np.nan
                            S[c][s]['local max gr var']= np.nan
                            S[c][s]['time of local max gr']= np.nan
                            S[c][s]['time of local max gr var']= np.nan
                        # find area under gr vs OD
                        from scipy import integrate, interpolate
                        da, dna= [], []
                        for fsample, gsample in zip(np.transpose(fs), np.transpose(gs)):
                            sod= np.exp(fsample)
                            integrand= lambda x: interpolate.interp1d(sod, gsample)(x)
                            iresult= integrate.quad(integrand, np.min(sod), np.max(sod),
                                                     limit= 100, full_output= 1)[0]
                            da.append(iresult)
                            dna.append(iresult/(np.max(sod) - np.min(sod)))
                        S[c][s]['area under gr vs OD']= np.mean(da)
                        S[c][s]['area under gr vs OD var']= np.var(da)
                        S[c][s]['normalized area under gr vs OD']= np.mean(dna)
                        S[c][s]['normalized area under gr vs OD var']= np.var(dna)
                    else:
                        S[c][s]['f' + dtype]= f.f
                        S[c][s]['f' + dtype + ' var']= f.fvar
                        S[c][s]['d/dt' + dtype]= f.df
                        S[c][s]['d/dt' + dtype + ' var']= f.dfvar
                        S[c][s]['d2/dt2' + dtype]= f.ddf
                        S[c][s]['d2/dt2' + dtype + ' var']= f.ddfvar
                    S[c][s][dtype + ' logmaxlike']= f.logmaxlike
                    S[c][s][dtype + ' gp']= cvfn
                    S[c][s][dtype + '_GP']= f
                    if stats:
                        for sname in f.ds.keys(): S[c][s][sname]= f.ds[sname]
                except KeyError:
                    print(dtype, 'is either not yet calculated or not recognized')
                    return
            if keysmessage:
                print('\nProcessed strains now have these keys:')
                self.datavariables(individual= True, conditions= c, strains= s, title= False)

    #####
    def getfitnesspenalty(self, ref, com, figs= True, nosamples= 100):
        '''
        Calculates the area between two growth rate versus OD curves, normalized by the length
        along the OD-axis where they overlap:

        e.g. p.getfitnesspenalty(['1.9% Raffinose 0.0µg/ml cycloheximide', '77.WT'],
                             ['1.9% Raffinose 0.5µg/ml cycloheximide', '77.WT'])

        Arguments
        --
        ref: [condition, strain] array for the reference
        com: [condition, strain] array to compare to the reference
        figs: if True, an example of the area between the curves is shown
        nosamples: for the bootstraps to estimate error (default: 100)
        '''
        fps= np.zeros(nosamples)
        # get and sample from Gaussian processes
        try:
            f0s, g0s, h0s= self.d[ref[0]][ref[1]]['OD_GP'].sample(nosamples)
            f1s, g1s, h1s= self.d[com[0]][com[1]]['OD_GP'].sample(nosamples)
        except KeyError:
            print("Failed: getstats('OD') needs to be run for these strains")
            return
        # process samples
        for j, (f0sample, gr0, f1sample, gr1) in enumerate(zip(np.transpose(f0s), np.transpose(g0s),
                                                               np.transpose(f1s), np.transpose(g1s))):
            od0, od1= np.exp(f0sample), np.exp(f1sample)
            # remove any double values because of OD plateau'ing
            from scipy.signal import argrelextrema
            imax= argrelextrema(od0, np.greater)[0]
            if np.any(imax):
                od0= od0[:imax[0]]
                gr0= gr0[:imax[0]]
            imax= argrelextrema(od1, np.greater)[0]
            if np.any(imax):
                od1= od1[:imax[0]]
                gr1= gr1[:imax[0]]
            # interpolate data
            from scipy import interpolate
            i0= interpolate.interp1d(od0, gr0)
            i1= interpolate.interp1d(od1, gr1)
            # find common OD range
            odinit= np.max([od0[0], od1[0]])
            odfin= np.min([od0[-1], od1[-1]])
            # perform integration to find normalized area between curves
            from scipy import integrate
            igrand= lambda x: i0(x) - i1(x)
            fps[j]= integrate.quad(igrand, odinit, odfin, limit= 100, full_output= 1)[0]/(odfin - odinit)
            # an example figure
            if figs and j == 1:
                plt.figure()
                plt.plot(od0, gr0, 'k-', od1, gr1, 'b-')
                x= np.linspace(odinit, odfin, np.max([len(od0), len(od1)]))
                plt.fill_between(x, i0(x), i1(x), facecolor= 'red', alpha= 0.5)
                plt.xlabel('OD')
                plt.ylabel('growth rate')
                plt.legend([ref[1] + ' in ' + ref[0], com[1] + ' in ' + com[0]], loc= 'upper left', bbox_to_anchor= (0.5, 1.05))
                plt.title('fitness penalty for ' + com[1] + ' in ' + com[0] + ' relative to '+ ref[1] + ' in ' + ref[0])
                plt.show(block= False)
        return np.mean(fps), np.var(fps)





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
            savename= self.wdir + self.name + '.pkl'
        gu.putpkl(savename, self)


    def makedataframe(self, type, conditions= 'all', strains= 'all', dfname= 'df',
                       conditionincludes= False, strainincludes= False,
                       conditionexcludes= False, strainexcludes= False,
                       nomedia= True, nonull= True, sortby= False, ascending= True):
        '''
        Creates a dataframe (using Pandas) of three different types:
            'timerow'   : all time-dependent variables as rows
            'timecol'   : all time-dependent variables as columns
            'notime'    : only time-independent variables

        For example,

            p.makedataframe('timecol')
            p.df
            p.exportdataframe()

        Arguments
        --
        type : 'timerow', 'timecol', or 'notime'
        conditions : conditions of interest (default: 'all')
        strains : strains of interest (default: 'all')
        dfname : name to use for dataframe (default: 'df')
        conditionincludes: selects only conditions with conditionincludes in their name
        strainincludes: selects only strains with strainincludes in their name
        conditionexcludes: ignore conditions with conditionexcludes in their name
        strainexcludes: ignore strains with strainexcludes in their name
        nomedia : if True (default), wells labelled media are excluded
        nonull : if True (default), wells labelled null are excluded
        sortby: an optional variable (or list of variables) to sort the data by
        ascending: whether to sort in ascending or descencing order (default: True)
        '''
        variables, nottimevariables= self.datavariables(out= True, conditions= conditions,
                                                        strains= strains,
                                                        conditionincludes= conditionincludes,
                                                        strainincludes= strainincludes,
                                                        conditionexcludes= conditionexcludes,
                                                        strainexcludes= strainexcludes)[1:]
        both= self.getconditionsandstrains(conditions= conditions, strains= strains,
                                           conditionincludes= conditionincludes,
                                           strainincludes= strainincludes,
                                           conditionexcludes= conditionexcludes,
                                           strainexcludes= strainexcludes,
                                           nomedia= nomedia, nonull= nonull)
        # remove 'time' from the time-dependent variables
        if type == 'timecol' or type == 'notime':
            variables= np.delete(variables, np.nonzero(variables == 'time')[0][0])
        dfIndex= 0
        S= self.d
        if type == 'timecol':
            # each row corresponds to a time-
            columns= ['condition', 'strain', 'timepoint'] + list(variables)
            dfs= []
            for c, s in both:
                ic= 0
                dft= pd.DataFrame(columns= columns)
                dft[columns[0]]= [c for t in self.t]
                dft[columns[1]]= [s for t in self.t]
                dft[columns[2]]= self.t
                for var in variables:
                    dft[var]= np.array(S[c][s][var]).tolist() if var in S[c][s] else [np.nan for t in self.t]
                dfs.append(dft)
            df= pd.concat(dfs)
        elif type == 'timerow':
            # each time-point has its own column
            columns= ['condition', 'strain', 'variable'] + ['timepoint: ' + str(t)
                                                            for t in range(len(self.t))]
            allseries= []
            for c, s in both:
                allseries += [pd.Series([c, s, var] + [list(S[c][s][var]) if var in S[c][s]
                                                     else [np.nan for t in self.t]][0], index= columns)
                              for var in variables]
            df= pd.DataFrame(allseries, columns= columns)
        elif type == 'notime':
            # variables that are independent of time
            df= pd.DataFrame(columns= ['condition', 'strain'] + list(nottimevariables))
            for c, s in both:
                df.loc[dfIndex]= [c, s] + [S[c][s][var] if var in S[c][s] else [np.nan]
                                           for var in nottimevariables]
                dfIndex += 1
        else:
            print(type, 'is not recognized')
            return

        # sort data
        if sortby:
            df= df.sort_values(by= sortby, ascending= ascending).reset_index(drop= True)
        # store result
        setattr(self, dfname, df)
        print(dfname, 'created')



    def exportdataframe(self, dfname= 'df', ftype= 'xlsx', savename= False):
        '''
        Export a data frame to an Excel or CSV file.

        Arguments
        --
        dfname : name of dataframe to export
        ftype : 'xlsx' (default) or 'csv'
        savename : name of file to be saved (default is to add _df to the name of the data set)
        '''
        if hasattr(self, dfname):
            if not savename:
                savename= self.wdir + self.name + '_' + dfname + '.' + ftype
            if ftype == 'xlsx':
                print('Writing Excel file:', savename)
                getattr(self, dfname).to_excel(savename, sheet_name= 'Sheet1', index= False)
            else:
                print('Writing CSV file:', savename)
                getattr(self, dfname).to_csv(savename, sep= ',', index= False)
        else:
            print('No dataframe', dfname, 'has been created')


    def rmdataframe(self, dfname= 'df'):
        '''
        Delete a dataframe.

        Arguments
        --
        dfname: dataframe to be deleted (default: 'df')
        '''
        if hasattr(self, dfname): delattr(self, dfname)


##########################
# specialized platereader
##########################


class slpr(platereader):

    def __init__(self, dname, aname= 'default', platereadertype= 'Tecan', wdir= '', dsheetnumber= 0,
                 asheetnumber= 0, ODfname= 'ODcorrection_Glucose_Haploid.txt', warn= False, info= True,
                 standardgain= False, ignorestrains= []):

        if platereadertype == 'Sunrise':
            # extract time from sheet 0 of Excel file
            timedf= pd.ExcelFile(wdir + dname).parse(0)
            t= []
            for index, row in timedf.iterrows():
                if '/' in row[0]:
                    t.append(int(row[0].split('/')[-2].split('s')[0]))
            self.t= np.array(t)/3600
            # data is in sheet 1
            dsheetnumber= 1

        # call platereader
        super().__init__(dname, aname, platereadertype, wdir, dsheetnumber, asheetnumber, ODfname,
                             warn, info, standardgain, ignorestrains)
        self.__doc__= super().__doc__

        # for Luis
        if 'media' in self.d:
            S= self.d
            for nd in self.datatypes:
                if nd != 'OD':
                    S['media']['null']['c-' + dn]= np.nan(np.size(self.t))
                    S['media']['null']['c-' + dn + 'perod']= np.nan(np.size(self.t))




#####

if __name__ == '__main__': print(platereader.__doc__)
