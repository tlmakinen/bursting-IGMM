"""This function generates a telegraph signal that oscillates from an
'ON' and 'OFF' position with wait times between each state distributed
according to a given distribution. This signal can be used to numerically
simulate promoter site activation for gene transcription"""


import numpy as np
class exponential:
    def __init__(self, k_on, k_off, duration):
    # k_on is exponential distribution clustering paramter for "ON" states, k_off for "OFF"
    # duration is the total simulation time, in seconds
    # stepsize is the number of seconds per step (e.g. 3 seconds per observation)
        self.k_on = k_on            # mean of "ON" time exponential distribution mu = 1./k_on
        self.k_off = k_off          # mean of "OFF" time exponential distribution
        
        on_off = np.zeros((duration))                                                   # create array of on-off durations in seconds
        on_off[0::2] = np.random.exponential(scale=(1./k_on), size=len(on_off[0::2]))   # fill every other value with exponentially-distributed on durations
        on_off[1::2] = np.random.exponential(scale=(1./k_off), size=len(on_off[1::2]))  # then fill in every other place with "off" durations

        promoter_arr = np.zeros((duration))
            
        state = 0      # start in the off state
        t0 = 0          # start at t0
        t1 = 0          # step size to be chosen in loop
        
                                            
        for i in range(len(promoter_arr)):
            t1 = int(on_off[i])                         # find the time until which the signal changes state
            promoter_arr[t0:(t0 + t1)] = state         # for the given duration of the ith state, fill the value with either "ON" or "OFF" values
            state = 1 - state                         # change the signal by flipping the switch to either 1 or 0.
            t0 = t0 + t1                                # move on in time
            
        self.signal = promoter_arr                      # returns telegraph signal of promoter indexed by stepsize (discrete steps).
        


