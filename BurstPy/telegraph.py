"""This function generates a telegraph signal that oscillates from an
'ON' and 'OFF' position with wait times between each state distributed
according to a given distribution. This signal can be used to numerically
simulate promoter site activation for gene transcription"""


import numpy as np
class exponential:
    def __init__(self, k_on, k_off, duration, stepsize):
        """ 
        # k_on is exponential distribution clustering paramter for "ON" states, k_off for "OFF"
        # duration is the total simulation time, in seconds
        # stepsize is the number of seconds per step (e.g. 3 seconds per observation) 
        Note: We draw our wait times in terms of stepsize in order to most closely replicate
        real laboratory data. 
        """
        self.k_on = k_on*stepsize          # mean of "ON" time exponential distribution mu = 1./k_on
        self.k_off = k_off*stepsize        # mean of "OFF" time exponential distribution
        self.stepsize = stepsize

        on_off = np.zeros(int(duration//self.stepsize))                                      # create array of on-off durations in terms of observations
        on_off[0::2] = np.random.exponential(scale=(1./self.k_on), size=len(on_off[0::2]))   # fill every other value with exponentially-distributed on durations
        on_off[1::2] = np.random.exponential(scale=(1./self.k_off), size=len(on_off[1::2]))  # then fill in every other place with "off" durations

        promoter_arr = np.zeros(int(duration//self.stepsize))                                # empty array of ON/OFF promoter switches
            
        state = 0       # start in the off state
        t0 = 0          # start at t0
        t1 = 0          # step size to be chosen in loop
        
                                            
        for i in range(len(promoter_arr)):
            t1 = int(on_off[i])                         # find the time until which the signal changes state
            promoter_arr[t0:(t0 + t1)] = state          # for the given duration of the ith state, fill the value with either "ON" or "OFF" values
            state = 1 - state                           # change the signal by flipping the switch to either 1 or 0.
            t0 = t0 + t1                                # move on in time
            
        self.signal = promoter_arr                      # returns telegraph signal of promoter indexed by observation (discrete steps).
        


