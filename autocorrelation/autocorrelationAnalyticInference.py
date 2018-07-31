""" 
Curve fitting (least squares error minimization) routine for MS2 trace data
autocorrelation function. Analytic autocorrelation function for a two-state 
system coded as in Desponds et al (2016).

Function works as follows:
    1) Takes in trace data package object from autocorrelationDataProcessing.py
    2) Computes average autocorrelation function for trace set
    3) Using the average fluorescence of the set, Pon is estimated according to
       Desponds et al.  <avg fluorescence> = Pon * Sum(loop agglomeration signal)
    4) Takes Pon and feeds it to analytic model
    5) Curve-fits analytic model to average autocorrelation function
"""

import numpy as np
import matplotlib.pyplot as plt
import gfp_signal
from telegraph import exponential
from scipy.optimize import curve_fit
from astropy.table import Table

# import analysis functions
from autocorrelationDataProcessing import tracePackageAutocorrelation

class fitAutocorrelationFunction():

    def __init__(self, tracePackageAutocorrelation):
        tracelist = tracePackageAutocorrelation.tracelist
        corrected_tracelist = tracePackageAutocorrelation.corrected_tracelist
        calibrated_tracelist = tracePackageAutocorrelation.calibrated_tracelist

        autoav = tracePackageAutocorrelation.autoav
        autostd = tracePackageAutocorrelation.autostd
        avgflors = tracePackageAutocorrelation.avgflors



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
