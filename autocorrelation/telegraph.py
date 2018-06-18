"""This function generates a telegraph signal that oscillates from an
'ON' and 'OFF' position with wait times between each state distributed
according to a given distribution. This signal can be used to numerically
simulate promoter site activation for gene transcription"""


import numpy as np
class exponential:
    def __init__(self, k_on, k_off, duration, stepsize):
    # k_on is exponential distribution clustering paramter for "ON" states, k_off for "OFF"
    # duration is the total simulation time, in seconds
    # stepsize is the number of seconds per step (e.g. 3 seconds per observation)

        on_off = np.zeros((duration//stepsize))                                             # create array of on-off durations
        on_off[0::2] = np.random.exponential(scale=(1./k_off), size=len(on_off[0::2]))      # fill every other value with exponentially-distributed off durations
        on_off[1::2] = np.random.exponential(scale=(1./k_on), size=len(on_off[1::2]))       # then fill in every other place with "ON" durations

        promoter_arr = np.zeros((duration//stepsize))
            
        signal = 0           # start in the off state
        t0 = 0               # start at t0
        t1 = 0               # step size to be chosen in loop
        
        step_arr = np.zeros(duration//stepsize)         # create our time step array (not really needed)
                                    
        for i in range(len(promoter_arr)):
            t1 = int(on_off[i] / stepsize)              # find the time until which the signal remains in its state
            promoter_arr[t0:(t0 + t1)] = signal         # for the given duration of the ith state, fill the value with either "ON" or "OFF" values
            signal = 1 - signal                         # change the signal by flipping the switch to either 1 or 0.
            t0 = t0 + t1                                # move on in time
            step_arr[i] = t0                            # 

        self.signal = promoter_arr                      # returns telegraph signal of promoter indexed by stepsize (discrete steps).
        self.step_arr = step_arr
        

