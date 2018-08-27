"""Python script following Desponds et al (2016) method for 
   drosophila embryo MS2 gene tagging autocorrelation signal 
   analysis. 
   
   This analysis cleans data traces to account for cell-to-cell
   signal variability and generates an AVERAGE autocorrelation 
   function for the set of cell traces. This eliminates wonky 
   autocorrelation shapes in any given trace. The cleaner average
   can then be fit for characteristic "burst" time. 

   Inputs: trace data (list format) and parameters:

        loop_function   #  known agglomeration MS2 loop function along transcribed gene
        start_index     #  start index of steady-state signal for PACKAGE
        stop_index      #  stop index of steady-state signal for PACKAGE
        tPol            #  RNA polII loading time
        k_elong         #  elongation rate of polII along gene
           
   """

import numpy as np
from scipy.optimize import curve_fit
from astropy.table import Table

# import our analysis functions


# helpful function for cutting the autocorrelation function into two parts 
# (only want righthand side decay shape)
def findMiddleIndex(input_list):
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return int(middle - .5)
    else:
        return (int(middle))


class tracePackageAutocorrelation:

    def __init__(self, tracelist, loop_function, start_index, stop_index):
        # compute autocorrelation with long time traces
        self.tracelist = tracelist     # list of spot intensities over time
        self.start_index = start_index
        self.stop_index = stop_index
        self.loop_function = np.asarray(loop_function)
        

    
        # compute autocorrelation with long time traces
        autolist = []
        calibrated_tracelist = []
        corrected_tracelist = []
        max_intens = []

        # calibrate fluorescence signal: F(t) = I(t) / Io, Io = <Imax> / sum(loop)
        for j in self.tracelist:
            max_intens.append(np.max(j[start_index:stop_index])) # get maximum fluorescence from each trace in given region

        Imax = np.mean(np.asarray(max_intens))   # mean of the maxima of all traces in set
        I_o = Imax / np.sum(loop_function)       # calibration constant depends on loop function over gene
        #print('Io = ', I_o)

        for trace in tracelist:
                    
            # calibrate fluorescence signal: F(t) = I(t) / Io        
            calibrated_trace = trace / I_o   

            # subtract out mean signal of each trace to prevent nucleus-nucleus variability
            corrected_trace = calibrated_trace - np.nanmean(calibrated_trace)
            
            # check for weird zero traces and fill their autocorrelations with nans
            colsum = np.sum(trace)    

            if (colsum == 0):                                      # don't plot the zero signal traces
                auto_norm = np.ones(len(trace)) * np.nan         # set zero signal cells to nan to be ignored

            else: 
                auto = np.correlate(corrected_trace, corrected_trace, 'full')
                #auto = autocorrelateSignal(corrected_trace)
                auto = auto[np.argmax(auto):]     # take half of the autocorrelation function

                # normalize each autocorrelation function
                auto_norm = auto / (auto[1])

            autolist.append(auto_norm)
            calibrated_tracelist.append(calibrated_trace)
            corrected_tracelist.append(corrected_trace)
            
        # compute average autocorrelation, ignoring inactive cells
        autoav = np.nanmean(np.asarray(autolist), axis=0)  
        autoav = autoav / autoav[1]                        # normalize averaged autocorrelation

        # compute standard error according to Desponds
        autovar = np.nanvar(np.asarray(autolist), axis=0)
        #std_error = autovar / (len(autolist) * (len(calibrated_trace) - (np.arange(len(autoav)))))
        std_error = np.nanstd(np.asarray(autolist), axis=0)

        # weighted standard error of average connected autocorr function
        err = std_error / ((len(calibrated_trace) - np.arange(len(autoav))) / len(calibrated_trace))

        # compute average fluorescence for each raw simulated trace
        avg_flors = []
        for i in calibrated_tracelist:
            avg_flors.append(np.mean(i))                   # stack in an array

        # now fit pon and save as an aspect of the processing package
        loops = loop_function                                            # pull in interpolated loop function    
        pon = np.mean(avg_flors) / np.sum(loops)
        pon_std = (np.std(avg_flors) 
                                / np.sqrt(len(avg_flors))) / np.sum(loops)     # standard error on pon for uncertainty estimations
        
        pon_low = pon - pon_std
        pon_hi = pon + pon_std                             # upper and lower bounds for pon fit
                                             
        self.pon = pon
        self.loop_function = loop_function                 # for fitting reference
        self.avgflors = np.asarray(avg_flors)              # return array of average fluorescence
        self.autoav = autoav[1:]                           # ignore the unfitted first point
        self.auto_err = err[1:]   
        self.calibrated_tracelist = calibrated_tracelist   # intensity-calibrated traces
        self.corrected_tracelist  = corrected_tracelist    # mean-subtracted traces
        self.autolist = autolist

