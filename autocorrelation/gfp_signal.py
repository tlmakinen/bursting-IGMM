"""Simulates GFP step-loop agglomeration according to a given promoter state binary telegraph signal """
import numpy as np

class gfp_signal:

    def __init__(self, telegraph, k_gfp, cts_per_gfp, k_polII, max_loops, stepsize):
        
        self.telegraph = telegraph
        self.k_gfp = k_gfp
        self.cts_per_gfp = cts_per_gfp
        self.k_polII = k_polII
        self.max_loops = max_loops
        self.stepsize = stepsize

        
        loops_per_step = k_gfp / stepsize          # get k_gfp from molecules / sec --> molecules / step
        
        polII_per_step = k_polII / stepsize    # start rate of polII in terms of discrete timesteps

        molecule_signal = np.zeros(len(telegraph))    # empty array of mRNA-GFP loops
        
        polII_arr = []   # empty array of polII molecules that will grow according to the pol_per_step rate.
                        # This is a list of strands of molecules
        
        # max_loops = 24   # max number of GFP loops per strand
        counter = 0      # keep track of how many time steps we've taken since we added the last polII molecule
        loopnum = 0      # keep track of the number of steps we've taken since we added the last loop

        
        for s in range(len(telegraph)):                 # iterate by time step
            
            if telegraph[s] == 1:                       # When the telegraph signal is ON, start a new pol II molecule according to the rate pol_per_step
                counter += 1
                if (counter/polII_per_step) >= 1:       # if we passed the start point of a new polII, add a new pol II molecule
                    polII_arr = np.append(polII_arr, 0)
                    counter = 0                         # restart counter
                            
            loopnum += 1                                # start the GFP loop counter

            if (loopnum/loops_per_step) >= 1:           # add a new GFP loop to each polII molecule in the loop array according to GFP loop rate               
                polII_arr = [j+1 for j in polII_arr] 
                loopnum = 0                             # reset the GFP loop counter for the discretized rate.
                
            #for i in range(len(polII_arr)):
                #if polII_arr[i] > max_loops:
                    #polII_arr[i] = 0                

            keep = np.asarray(polII_arr) <= max_loops   # if we've gone over the max number of GFP loops, set polII element to zero since 
                                                        # the chain of GFP attached to that polII will decay really quickly     
                                                                         
            molecule_signal[s] = sum(keep*polII_arr)    # sum up the number of GFP molecules from the polII list
            loopnum=0                                   # reset the loop counter

        self.counts = molecule_signal*cts_per_gfp       # returns the simulated GFP photon counts
        self.norm = molecule_signal                     # signal without photon counts
            