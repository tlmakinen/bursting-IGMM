"""Python script following Desponds et al (2016) method for 
   drosophila embryo MS2 gene tagging autocorrelation signal 
   analysis. 
   Inputs: trace data (.csv format) and parameters:
   
   Default Parameters (from Desponds et al) 
   tPol = 6 seconds                     #  RNA polII loading time
   k_elong = 25 bp / second             #  elongation rate of polII along gene
   sizePol = tPol * k_elong = 150 bp    #  effective 'footprint' of polII along gene """

import numpy as np
import matplotlib.pyplot as plt
import gfp_signal
from telegraph import exponential
from scipy.optimize import curve_fit
from astropy.table import Table

# import our analysis functions

from loopFunction import ms2Loops

# helpful function for cutting the autocorrelation function into two parts 
# (only want righthand side decay shape)
def findMiddleIndex(input_list):
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return int(middle - .5)
    else:
        return (int(middle))

class tracePackageAutocorrelation:

    def __init__(self, tracetable, start_index, stop_index, tPol, k_elong):
        # compute autocorrelation with long time traces
        self.start_index = start_index
        self.stop_index = stop_index
        self.tPol  = tPol
        self.k_elong = k_elong
        self.sizePol = tPol * k_elong

        despondsfile = 'standalone/therightL.mat'

        table = tracetable
        autolist = []
        tracelist = []
        tracelist_corr = []
        max_intens = []


    # calibrate fluorescence signal: F(t) = I(t) / Io, Io = <Imax> / sum(loop)
    # get maximum fluorescence from each trace in given region
    # ***Note***: even if a cell doesn't give off fluorescence, we still 
    # count it when computing max intensity
        for j in table.colnames:
            max_intens.append(np.max(table[j][start_index:stop_index])) # get maximum fluorescence from each trace in given region
        
        Imax = np.mean(np.asarray(max_intens))   # mean of the maxima of all traces in set
        loopfn = ms2Loops(despondsfile, self.tPol, self.k_elong).loopFn
        I_o = Imax / np.sum(loopfn) *150
        
        for name in table.colnames:
            trace = table[name][start_index:stop_index]
            # calibrate fluorescence signal: F(t) = I(t) / Io        
            calibrated_trace = trace / I_o   # CALIBRATED trace
            # subtract out mean signal of each trace to prevent nucleus-nucleus variability
            corrected_trace = calibrated_trace - np.nanmean(calibrated_trace)   # CORRECTED trace
            # check for weird zero traces and fill their autocorrelations with nans
            colsum = np.sum(trace)    
            if (colsum == 0):                                # don't plot the zero signal traces
                auto_norm = np.ones(len(trace)) * np.nan     # set zero signal cells to nan to be ignored in autocorrelation

            else: 
                auto = np.correlate(corrected_trace, corrected_trace, 'full')
                auto = auto[np.argmax(auto):]                # take half of the autocorrelation function

                # finally, normalize according to Desponds et al     
                for r in range(len(corrected_trace)):
                    norm = []
                    sm2 = 0.
                    for k in range(1,len(corrected_trace)):
                        sm2 += (corrected_trace[k]**2)
                    norm.append(((len(corrected_trace)-float(r)) / len(corrected_trace)) * sm2)
                auto_norm = auto / (np.asarray(norm))
            autolist.append(auto_norm) #/ auto_norm[1])   # normalize by second data point--> VERIFY
            tracelist.append(calibrated_trace)
            tracelist_corr.append(corrected_trace)
            
        autoav = np.nanmean(np.asarray(autolist), axis=0)
        autostd = np.nanstd(np.asarray(autolist), axis=0)

        self.autoav = autoav[1:]                # return the function starting from the second, normalized point
        self.autostd = autostd[1:]              # array of standard deviations of each autocorrelated time point
        self.tracelist_corr = tracelist_corr    # F(t) - mean(F(t)) : calibrated and corrected traces
        self.tracelist = tracelist              # calibrated traces (mean not subtracted.. for Pon calculation)



        def getAverageFlors(self):
            # compute average fluorescence for each raw simulated trace
            avg_flors = []
            for i in self.tracelist:
                avg_flors.append(np.mean(i))   # stack in an array
            return np.asarray(avg_flors)



        # fit for pon in a SIMPLE way
        def fitPon(self):
            despondsfile = 'standalone/therightL.mat'
            self.tPol = tPol;                      # RNA polII loading time
            self.k_elong = k_elong                 # Elongation rate
            self.sizePol = tPol * k_elong          # Footprint, in basepairs, of polII
                
            avg_flors = getAverageFlors(self.tracelist)
            loops = ms2Loops.loopFn(despondsfile, self.tPol, self.k_elong, self.sizePol)     # pull in loop function    

            pon = np.mean(avg_flors) / np.sum(loops) * 150
            pon_var = np.std(avg_flors)**2

            return pon,pon_var
            # print('pon          = ', pon, '\npon variance = ', pon_var)
