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
from telegraph import exponential
from scipy.optimize import curve_fit
from astropy.table import Table

# import analysis functions
from autocorrelationDataProcessing import tracePackageAutocorrelation
from loopFunction import loopInterpolate   # get interpolation function to calculate Pon

class autocorrelationAnalytic:
    """
    Create a class for the autocorrelation analytic model so that we can vary pon and still
    implement a ratesum-dependent curve fitting later.
    """
    def __init__(self, tPol, k_elong, stepsize, pon, interploops):
        self.tPol = tPol
        self.k_elong = k_elong
        self.stepsize = stepsize
        self.pon = pon
        self.interploops = interploops   # loop function interpolated by discrete polII location on gene


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

    def autocorrAnalyticFunction(self, t, ratesum):  # take in parameters and t (signal data array index in seconds)
        # define all needed parameters #
        time = t                                # our time array is what we pass to the function (usually in seconds)
        # get the loop function
        loops = self.interploops
        # take in pon as aspect of function
        p_on = self.pon
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
                    sm0 += p_on*p_off*(loops[i] * loops[j] * np.exp((delta-1)*np.abs(time[t] - j + i)))
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





class fitAutocorrelationFunction():

    def __init__(self, tracePackageAutocorrelation, tPol, k_elong, stepsize):
        self.tracelist = tracePackageAutocorrelation.tracelist
        self.corrected_tracelist = tracePackageAutocorrelation.corrected_tracelist
        self.calibrated_tracelist = tracePackageAutocorrelation.calibrated_tracelist

        self.autoav = tracePackageAutocorrelation.autoav
        self.auto_err = tracePackageAutocorrelation.auto_err
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
        loops = self.loop_function                                            # pull in interpolated loop function    
        sizePol = self.tPol * self.k_elong                                    # polII footprint on gene
        pon = np.mean(self.avgflors) / np.sum(loops)
        pon_std = (np.std(self.avgflors) 
                                / np.sqrt(len(self.avgflors))) / np.sum(loops)     # standard error on pon for uncertainty estimations
        
        pon_low = pon - pon_std
        pon_hi = pon + pon_std
        return pon,pon_low,pon_hi,pon_std                                     # returns mean, upper, and lower bounds for pon fit
        

    def tracePackageBootstrap(self):
    
        # let's get an idea of the spread of our data by using the bootstrap method
        # create a routine that bootstrap fits a set of data arrays
        
        tracelist = np.asarray(self.calibrated_tracelist)          # convert to a numpy array to play with indices
        n_traces = len(tracelist)                # number of traces in our dataset
        n_trials = 1000                      # how many times we wish to compute the bootstrap loop
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
        weightedstderr = std_dev_arr / ((N-pts)*(np.sqrt(N)))
        return weightedstderr



        # Now let's take the 68% confidence intervals of the set of 10000 autocorrelation functions. This is the STANDARD ERROR
        # on our dataset    
        #upperlim,lowerlim = np.percentile(a=np.asarray(auto_averages), axis=0, q=[84,16])
        #med = np.nanmean(auto_averages, axis=0)
        #return med,(med - lowerlim)/weights,(upperlim-med)/weights    # returns distance between average and up,low limits

    def leastSquaresAutoFit(self, upperbound, lowerbound, printvals=True, fitpon=True, pon=0.3):
        """
        Least Squares fitting using scipy.optimize module. 
        Inputs:
        Upper and Lower bounds denote the range in which we perform the inference for the ratesum = kon + koff
        
        printvals: if True, we print the inferred values and generate a plot of the fitted function on top of
        the data, with standard error bars shown

        fitpon: if chosen to be False, don't fit pon from trace package

        pon: optional, for testing inference method

        To access the autocorrelation function, we call in the class defined above
        """
        # first, fit for pon and get error
        if fitpon == True:
            pon,ponlow,ponhi,ponstd = self.fitPon()
        
                
        # then create the analytic model class object
        autocorrelationAnalyticPack = autocorrelationAnalytic(self.tPol, self.k_elong, 
                                                    self.stepsize, pon, self.interploops)   # our package with mean pon

        if upperbound != None:
            bounds = ([lowerbound], [upperbound])
        else:
            bounds = None

        t = np.arange(len(self.autoav)) * self.stepsize    # fit function in units of seconds
        self.autocorrFunc = autocorrelationAnalyticPack.autocorrAnalyticFunction
        popt,pcov = curve_fit(f=self.autocorrFunc, xdata=t[1:], 
                        ydata=self.autoav[1:], bounds=bounds)#, sigma=self.auto_err[1:])  # fit everything but first (pinned) data point, using lowerlim distance as standard error

        ratesum = popt[0]
        kon_fit = ratesum*pon
        koff_fit = ratesum - kon_fit

        # get error range on the decoupling of kon, koff
        # lower bound:
        #kon_low = ratesum*ponlow
        #koff_low = ratesum - kon_low
        # upper bound:
        #kon_up = ratesum * ponhi
        #koff_up = ratesum - kon_up

        # compute errors in kon, koff
        #kon_up_err = np.abs(kon_up - kon_fit)
        #kon_low_err = np.abs(kon_fit - kon_low)

        #koff_up_err = np.abs(koff_fit - koff_up)
        #koff_low_err = np.abs(koff_fit - koff_low)

        chrtime = 1 / ratesum

        if printvals == True:
            print("Pon                       = ", pon)#, '+/-'), ponstd)
            print("k_on + k_off              = ", ratesum, 's^-1')
            print("k_on                      = ", kon_fit)#, '+', kon_up_err, '-', kon_low_err)
            print("k_off                     = ", koff_fit)#, '+/-', koff_up_err)
            print("t_polII_block             =  6 seconds")
            print("characteristic timescale  = ", chrtime, 'seconds')
            print("covariance                = ", pcov[0][0])

            # incorporate a plot of the fitted autocorrelation with data
            fig,ax = plt.subplots(1, 1, figsize=(10,5), sharex=True)
            ax.plot(t, self.autocorrFunc(t,ratesum), label='analytic model best fit')
            # show the mean of our original dataset in red
            ax.scatter(t[::2],self.autoav[::2], marker='.', 
                            color='r', label = 'Mean Simulated Data Autocorrelation')
            
            ax.errorbar(x=t, y=self.autoav, 
                        yerr=(self.auto_err*1, self.auto_err*1), 
                                ecolor='b', alpha=0.2, label = r'1-$\sigma$ weighted standard error from mean')  

            ax.plot(chrtime, 0, marker='+', zorder=10, linestyle='none',
                                            color='#de2d26', label='Characteristic Time Point')

            plt.legend(loc="best")
            ax.set_ylim(-1, 1.3)
            ax.set_xlabel(r'Wait time $\tau$ (seconds)', fontsize=15)
            ax.set_ylabel(r'autocorrelation M($\tau$)', fontsize=15)
            # add in a summary of our fitted parameters
            plt.text(x=50, y=-.8, s='Analytic fitted parameters \n\n$k_{on}$ = ' 
                                            + str(kon_fit) + '\n$k_{off}$ = ' + str(koff_fit))
            plt.show()


        return kon_fit,koff_fit,chrtime,pon,popt,pcov

