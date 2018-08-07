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
from loopFunction import loopInterpolate   # get interpolation function to calculate Pon

class fitAutocorrelationFunction():

    def __init__(self, tracePackageAutocorrelation, tPol, k_elong, stepsize):
        self.tracelist = tracePackageAutocorrelation.tracelist
        self.corrected_tracelist = tracePackageAutocorrelation.corrected_tracelist
        self.calibrated_tracelist = tracePackageAutocorrelation.calibrated_tracelist

        self.autoav = tracePackageAutocorrelation.autoav
        self.autostd = tracePackageAutocorrelation.autostd
        self.avgflors = tracePackageAutocorrelation.avgflors
        self.loop_function = tracePackageAutocorrelation.loop_function

        # Assumed parameters
        self.stepsize = stepsize
        self.tPol = tPol
        self.k_elong = k_elong
        # interpolated loop function
        self.interploops = loopInterpolate(self.loop_function, self.k_elong, self.tPol)
    '''
    Fitting for Pon, the probability that our promoter is ON during the trace time window.
    This mini function takes in a loop function (by basepair along the gene) and interpolates
    it to get the number of loops indexed by the discrete polII location on the gene. 
    For example, a polII with loading time tPol=6sec and elongation rate k_elong=25bp/sec, the polII has an
    effective footprint on the gene of 150 bp. Thus the gene index i is
    len(gene) / (size of polII)
    '''
    def fitPon(self):
        loops = self.loop_function                     # pull in interpolated loop function    
        sizePol = self.tPol * self.k_elong           # polII footprint on gene
        pon = np.mean(self.avgflors) / np.sum(loops)
        pon_std = np.std(self.avgflors / np.sum(loops)) 

        return pon
        # print('pon          = ', pon, '\npon variance = ', pon_var)


    def tracePackageBootstrap(self):
    
        # let's get an idea of the spread of our data by using the bootstrap method
        # create a routine that bootstrap fits a set of data arrays
        
        tracelist = np.asarray(self.calibrated_tracelist)          # convert to a numpy array to play with indices
        n_traces = len(tracelist)                # number of traces in our dataset
        n_trials = 100                      # how many times we wish to compute the bootstrap loop
        trace_indx = np.arange(n_traces)      # the index range of the list of traces
        auto_averages = []                    # list of averaged autocorrelation functions (should be n_trials long)


        # sample randomly, with replacement, a new set of n_traces 10,000 times
        for i in range(n_trials):

            # from our list of traces, sample radomly the trace_index of these datum in the list
            random_indx = np.random.choice(trace_indx, size=n_traces, replace=True)  
            # then use this array of indices to create our random distribution of traces
            sample_set = tracelist[random_indx]
            autolist = []
            # next compute the autocorrelation function from this set of traces        
            for sample in sample_set:
                
                corrected_sample = sample - np.nanmean(sample)            
                # check for weird zero traces and fill their autocorrelations with nans
                colsum = np.sum(sample)    
                if (colsum == 0):                                      # don't plot the zero signal traces
                    auto_norm = np.ones(len(sample)) * np.nan         # set zero signal cells to nan to be ignored

                else: 
                    auto = np.correlate(corrected_sample, corrected_sample, 'full')
                    #auto = autocorrelateSignal(corrected_trace)
                    auto = auto[np.argmax(auto):]     # take half of the autocorrelation function        
                    auto_norm = auto
                autolist.append(auto_norm)
                
            # compute average autocorrelation, ignoring inactive cells
            autoav = np.nanmean(np.asarray(autolist), axis=0)  
            autoav = autoav / autoav[1]                        # normalize averaged autocorrelation
            # append average to list
            auto_averages.append(autoav[1:])


            # compute weighted standard errors from bootstrap package

        N = len(auto_averages[1]) # number of tau time delays in autocorrelation function
        pts = np.arange(N)        # index of points        
        std_dev_arr = np.nanstd(auto_averages, axis=0, ddof=0)
        weightedstderr = std_dev_arr / ((N-pts)/(np.sqrt(N)*n_trials))
        return weightedstderr



        # Now let's take the 68% confidence intervals of the set of 10000 autocorrelation functions. This is the STANDARD ERROR
        # on our dataset    
        #upperlim,lowerlim = np.percentile(a=np.asarray(auto_averages), axis=0, q=[84,16])
        #med = np.nanmean(auto_averages, axis=0)
        #return med,(med - lowerlim)/weights,(upperlim-med)/weights    # returns distance between average and up,low limits


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


    def autocorrAnalytic(self, t, ratesum):  # take in parameters and t (signal data array index in seconds)
        # define all needed parameters #
        stepsize = self.stepsize            # time between observations, seconds
        time = np.arange(len(t)) * stepsize
        # get the loop function
        loops = self.interploops
        # compute Pon
        p_on = self.fitPon()
        p_off = 1-p_on
                
        # write analytic autocorrelation function according to Desponds et al (2016):    
        # ratesum is defined as k_on + k_off
        delta = 1 - ratesum    

        # write a for loop to do the double sums to compute connected correlation:
        c_arr = []
        for t in range(len(time)):
            sm0 = 0
            for i in range(len(loops)):
                for j in range(len(loops)):
                    sm0 += p_on*p_off*(loops[i] * loops[j] * np.exp((delta-1)*np.abs(t - j + i)))
            c_arr.append(sm0)
        connected_corr = np.asarray(c_arr)    # the two-state connected correlation function   
        N = len(time)         # CONSTANT trace length

        # Add in the finite trace correction for the Ornstein-Uhlenbeck process    
        # perform the summations
        corrected_full = []
        Co = connected_corr[0]                # initial condition of connected correlation function      
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
            
            corrected_full.append((connected_corr[r] + 
                                        (1/N) *((1/N) - (2./(N-r))) * (N*Co + sm1)) + ((2/(N*(N-r))) * (r*Co + sm2 + sm3)))
    
        normed = np.asarray(corrected_full)
        return (normed / normed[1])             # normalize at second data point


    def leastSquaresAutoFit(self, printvals, upperbound, lowerbound):
        """
        When fitting the autocorrelation function within the standard deviation error range, we need to exclude the first point
        since the standard error here is zero (all autocorr functions normalized there). Otherwise, we get a divide by zero and 
        can't fit the fuction
        """
        if upperbound != None:
            bounds = ([lowerbound], [upperbound])
        else:
            bounds = None
        t = np.arange(len(self.autoav))
        weightedstd = self.tracePackageBootstrap()
        print(bounds)
        popt,pcov = curve_fit(f=self.autocorrAnalytic, xdata=t[1:], 
                        ydata=self.autoav[1:], bounds=bounds)#, sigma=weightedstd[1:])  # fit everything but first (pinned) data point, using lowerlim distance as standard error

        ratesum = popt[0]
        p_on = self.fitPon()
        kon_fit = ratesum*p_on
        koff_fit = ratesum - kon_fit
        chrtime = 1 / ratesum

        if printvals == True:
            print("Pon                       = ", p_on)
            print("k_on + k_off              = ", ratesum, 's^-1')
            print("k_on                      = ", kon_fit)
            print("k_off                     = ", koff_fit)
            print("t_polII_block             =  6 seconds")
            print("characteristic timescale  = ", chrtime, 'time units')
            print("covariance                = ", pcov[0][0])

        return kon_fit,koff_fit,chrtime,p_on,popt,pcov,weightedstd
