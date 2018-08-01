"""
Python script to stochastically simulate (SS) a package of N MS2 GFP loop traces of common length M. 
Method will utilize a two-state promoter model, assuming exponentially-distributed
kon and koff values 

   MS2 Casette system input parameters:
   
   Default Parameters (from Desponds et al) 
   tPol = 6 seconds                     #  RNA polII loading time
   k_elong = 25 bp / second             #  elongation rate of polII along gene
   sizePol = tPol * k_elong = 150 bp    #  effective 'footprint' of polII along gene
   max_loops = 24   # max number of gfp loops that attach to each polII site

For simulating N traces to create a package that we'd expect to get from actual data.
    Parameters:
    k_on:         mean of exponential "on" wait time distribution
    k_off:        mean of exponential "off" wait time distribution
    duration:     M length trace, in seconds
    stepsize:     observation time, e.g. number of seconds / emission point.



num = 100           # number of traces in our simulated dataset
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.table import Table

# import our analysis functions
import pol_signal
from telegraph import exponential


class tracePackageSimulation():

    def __init__(self, num_traces, k_on, k_off, duration, stepsize, loop_function, k_elong, tPol):
        
        self.num_traces = num_traces
        self.tPol = tPol                      # RNA polII loading time
        self.k_elong = k_elong                # polII speed along gene
        self.sizePol = tPol * k_elong         # polII "footprint" on gene
        self.loop_function = loop_function    # array of loops indexed by gene basepair coordinate
        self.kon = k_on                       # ON time distribution rate constant
        self.koff = k_off                     # OFF time distribution rate constant

        tracelist = []                        # empty list of traces 
        max_list = []
        tel_list = []                         # list of telegraph signals

        for i in range(self.num_traces):
            tel = exponential(k_on, k_off, duration)  # create a new signal every time      
            pol = pol_signal.pol_signal(telegraph=tel.signal, k_elong=self.k_elong, 
                                        loop_function=self.loop_function, tPol=tPol, stepsize=stepsize)
            trace = np.asarray(pol.signal) 
            tracelist.append(trace)
            tel_list.append(tel.signal)
            max_list.append(np.max(np.asarray(trace)))    
        
        self.tracelist = tracelist                  # array of generated traces
        self.max_list = max_list                    # array of trace maxima for p_on fitting later 
        self.tel_list = tel_list