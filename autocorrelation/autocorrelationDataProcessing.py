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

                # finally, normalize        
                #auto_norm = auto / auto[0]
                for r in range(len(auto)):
                    norm = []
                    sm2 = 0.
                    for k in range(1,len(corrected_trace)):
                        sm2 += (corrected_trace[k]**2)
                    norm.append(((len(corrected_trace)-float(r)) / len(corrected_trace)) * sm2)

                auto_norm = auto / (np.asarray(norm))

            autolist.append(auto_norm) #/ auto_norm[1])   # normalize by second data point
            calibrated_tracelist.append(calibrated_trace)
            corrected_tracelist.append(corrected_trace)
            
        # compute average autocorrelation, ignoring inactive cells
        autoav = np.nanmean(np.asarray(autolist), axis=0)  
        autoav = autoav / autoav[1]                        # normalize averaged autocorrelation
        # compute standard deviation, weighting error by number of points correlated
        autostd = np.nanstd(autoav, axis=0, 
                            ddof=np.arange(len(autoav)))  
        # compute average fluorescence for each raw simulated trace
        avg_flors = []
        for i in calibrated_tracelist:
            avg_flors.append(np.mean(i))            # stack in an array

        self.loop_function = loop_function                 # for fitting reference
        self.avgflors = np.asarray(avg_flors)              # return array of average fluorescence
        self.autoav = autoav[1:]                           # ignore the unfitted first point
        self.autostd = autostd[1:]   
        self.calibrated_tracelist = calibrated_tracelist   # intensity-calibrated traces
        self.corrected_tracelist  = corrected_tracelist    # mean-subtracted traces


