""" 
Curve fitting (least squares error minimization) routine for MS2 trace data
autocorrelation function. Analytic autocorrelation function for a two-state 
system coded as in Desponds et al (2016).

Function works as follows:
    1) Takes in corrected, autocorrelated trace data package object from autocorrelationDataProcessing.py
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
        self.tracelist = tracePackageAutocorrelation.tracelist
        self.corrected_tracelist = tracePackageAutocorrelation.corrected_tracelist
        self.calibrated_tracelist = tracePackageAutocorrelation.calibrated_tracelist

        self.autoav = tracePackageAutocorrelation.autoav
        self.autostd = tracePackageAutocorrelation.autostd
        self.avgflors = tracePackageAutocorrelation.avgflors
        self.loop_function = tracePackageAutocorrelation.loop_function


# fit for pon in a SIMPLE way
    def fitPon(self):
        loops = self.loop_function    # pull in loop function    
        pon = np.mean(self.avgflors) / np.sum(loops) * 150
        pon_std = np.std(self.avgflors) / len(self.avgflors)

        return pon,pon_std
        # print('pon          = ', pon, '\npon variance = ', pon_var)

    """"
    Autocorrelation function analysis for Two-State Model using Despond et al's (2016)
    parameters:
        tPol: polII loading time
        k_elong:  MS2 GFP loop elongation rate along gene
        tPol:     polII loading time, assumed to load at every observation time
                ***NOTE: this is subject to change since there may exist a blocking time,
                or characteristic "traffic jam" time during polII loading ***
        sizePol:  "Footprint" or size of polII in basepairs, along gene
        
        k_gfp:    gfp loop agglomeration rate such that L(t) = k_gfp*t, L(t) !> maxloops
        maxloops: max number of gfp loops allowed in chain
        
        FITTING FOR:
        ratesum: the inverse characteristic decay timescale for autocorrelation. Tells us how
                similar the signal is to itself and can reveal burst characteristics
                chartime = 1 / (k_on + k_off)
                
        Pon:      Probability that system is in "ON" state.
        k_on:     given as 1/mu for exponential distribution of ON wait times
        k_off:    1/mu for exponential distribution of OFF wait times

        ***NOTE: k_on and k_off decoupled by Pon = k_on / (k_on + k_off)
        """


    def autocorrAnalytic(self, t, chartime):  # take in parameters and t (signal data array index in seconds)
        # define all needed parameters #
        despondsfile = 'standalone/therightL.mat'
        stepsize = 6            # time between observations, seconds
        tPol=6;                 # polII loading time
        k_elong=25;             # Elongation rate
        sizePol = tPol * k_elong     # Footprint, in basepairs, of polII
        
        time = np.arange(len(t)) * stepsize
        
        tracelength = len(t)        # calculate from trace inputs
        
        # get the loop function
        loops = self.loop_function
        
        # do the Pon fitting (once I've written the function)
        # pon,pon_upper,pon_lower = pon_fit()
        
        # FOR NOW:
        #p_on = fitPon(tracelist)
        p_on,pon_std = self.fitPon()
        p_off = 1-p_on
        
        # compute normalization constant
        avg_flors = self.avgflors
        
        # write analytic autocorrelation function according to Desponds et al:    
        # chartime is defined as k_on + k_off
        
        delta = 1 - chartime    

        # write a for loop to do the double sums to compute connected correlation:
        
        c_arr = []
        
        for t in range(len(time)):
            sm0 = 0
            for i in range(len(loops)):
                for j in range(len(loops)):
                    sm0 += p_on*p_off*(loops[i] * loops[j] * np.exp((delta-1)*np.abs(t - j + i)))
        
            c_arr.append(sm0)
        #plt.plot(c_arr)
            
        connected_corr = np.asarray(c_arr)    # the two-state connected correlation function
        #return connected_corr / np.max(connected_corr)
    
        N = len(time)         # CONSTANT trace length
        # Add in the finite trace correction for the Ornstein-Uhlenbeck process    
        # perform the summations
        # initialize corrected lists
        corrected_full = []
        Co = connected_corr[0]    # initial condition of connected correlation function
            
        # Now we're going to correct EVERY data point in the connected autocorrelation function
        for r in range(len(connected_corr)):
        
            sm1 = 0    
            for k in range(1,N):
                sm1 += 2*(N - k)*connected_corr[k]
        
            sm2 = 0   
            for j in range(1,r):
                sm2 += 2*(r-j)*connected_corr[j]

            sm3 = 0
            for m in range(1,N):
                sm3 += connected_corr[m] * (np.min(np.asarray([m+r, N])) - np.max(np.asarray([r, m])))
            # full correction
            corrected_full.append((connected_corr[r] + 
                                        (1/N) *((1/N) - (2./(N-r))) * (N*Co + sm1)) + ((2/(N*(N-r))) * (r*Co + sm2 + sm3)))
    
        normed = np.asarray(corrected_full)
        
        return (normed / normed[1])


    def leastSquaresAutoFit(self):
        """
        When fitting the autocorrelation function within the standard deviation error range, we need to exclude the first point
        since the standard error here is zero (all autocorr functions normalized there). Otherwise, we get a divide by zero and 
        can't fit the fuction
        """
        t = np.arange(len(self.autoav))
        popt,pcov = curve_fit(self.autocorrAnalytic, t[0:], self.autoav[0:], sigma=self.autostd[0:])

        ratesum = popt[0]

        kon_fit = ratesum*p_on
        koff_fit = ratesum - kon_fit
        chrtime = 1 / ratesum

        print("k_on + k_off              = ", ratesum, 's^-1')
        print("k_on                      = ", kon_fit)
        print("k_off                     = ", koff_fit)
        print("t_polII_block             =  6 seconds")
        print("characteristic timescale  = ", chrtime, 'time units')
        print("covariance                = ", pcov[0][0])

        return kon_fit,koff_fit,chrtime,pcov[0][0]
