"""
Python script to generate a package of N MS2 GFP loop traces of common length M. 
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
from loopFunction import ms2Loops
import gfp_signal
from telegraph import exponential




class tracePackageSimulation():

    def __init__(self, k_on, k_off, duration, stepsize, max_loops, k_elong, tPol):
        self.tPol = tPol
        self.k_elong = k_elong
        self.sizePol = tPol * k_elong

        self.kon = k_on
        self.koff = k_off
        # load in Desponds et al's loop function
        #despondsfile = "standalone/therightL.mat"
        #ms2loops = ms2Loops(dspondsfile, self.tPol, self.k_elong)

auto_traces = []    # empty list of autocorrelation arrays
tracelist = []      # empty list of traces 
corrected_traces = []     # corrected traces
max_list = []



for i in range(num):
    tel = exponential(k_on, k_off, duration, stepsize)  # create a new signal every time      
    gfp = gfp_signal.gfp_signal(telegraph=tel.signal, k_elong=k_elong, max_loops=max_loops, tPol=tPol, stepsize=stepsize)
    
    trace = np.asarray(gfp.signal) 
    
    tracelist.append(trace)
    max_list.append(np.max(np.asarray(trace)))