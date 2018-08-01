"""
Simulates GFP step-loop agglomeration according to a given promoter state binary telegraph signal 
Assumes that the MS2 system loop agglomeration function is known

Inputs:

 
   MS2 Casette system input parameters:
   telegraph                            #  binary array of ON/OFF times indexed in seconds representing the 
                                        #  promoter state

   Default Parameters (from Desponds et al) 

   tPol = 6 seconds                     #  RNA polII loading time
   k_elong = 25 bp / second             #  elongation rate of polII along gene
   sizePol = tPol * k_elong = 150 bp    #  effective 'footprint' of polII along gene
   loop_function                        #  stepwise function (array indexed by bp) describing the number of MS2 loops accumulated
                                        #  at given position (in base pairs) along gene


"""
import numpy as np

class pol_signal:

    def __init__(self, telegraph, k_elong, tPol, loop_function, stepsize):

              
        self.telegraph = np.array(telegraph)           # promoter state signal in SECONDS
        self.k_elong = k_elong               # polII elongation rate; default is 25 bp/s
        self.tPol = tPol                     # polII loading time in SECONDS
                                             # default PolII loading time is 6s for desponds,3 for Lagha
        self.loop_function = np.asarray(loop_function)   # known loop agglomeration pattern for MS2 system
        self.stepsize = stepsize             # observation time

        # now let's load in the stepwise MS2 loop agglomeration. Desponds et al
        # calculated this and created a nice little function of MS2 loops by basepair
        # along the gene. We're going to put this in terms of seconds using the polII 
        # rate, k_elong.
        ms2loops = self.loop_function
        gene = np.arange(len(ms2loops))               # gene index, by basepair
        molecule_signal = np.zeros(len(telegraph))    # empty array of signal at every given time point
        polII_arr = []   # empty array of polII molecules that will grow according to the pol_per_step rate.
        pol_times = []   # array of the times each polII molecule on gene has been on gene
        counter = 0      # keep track of how many time steps we've taken since we added the last polII molecule
        
        
        for s in range(len(telegraph)):                    # iterate every second            
            if telegraph[s] == 1:                          # When the telegraph signal is ON, start a new pol II molecule according to the rate pol_per_step
                counter += 1
                if (counter/tPol) >= 1:                    # if we passed the start point of a new polII, add a new pol II molecule
                    polII_arr = np.append(polII_arr, 0)
                    pol_times = np.append(pol_times, 0)    # add each polII's coordinate to array
                    counter = 0                            # restart counter
                            
                                               
            #polII_arr = [j+1 for j in polII_arr]          # add a new GFP loop to each polII molecule in the loop array according to GFP loop rate
            
            pol_times = [j+1 for j in pol_times]                     # move each polII along in time
            pol_positions = [int(j * k_elong) for j in pol_times]    # update each polII molecule's position, in bp using k_elong        
            
            for j in range(len(pol_positions)):
                if pol_positions[j] > (len(gene)-1):
                     polII_arr[j] = 0                       # set to 0 when polII reaches end of gene
                                                            # simulates rapid decay of MS2 chain
                else:
                    polII_arr[j] = ms2loops[pol_positions[j]]    # add a GFP loop according to loop function

            #keep = (np.asarray(polII_arr) <= max_loops)   # if we've gone over the max number of GFP loops, set polII element to zero since 
                                                           # the chain of GFP attached to that polII will decay really quickly     
                                                                         
            #molecule_signal[s] = np.sum(keep*polII_arr)    # sum up the number of GFP molecules from the polII list
            molecule_signal[s] = np.sum(polII_arr)
                        
        self.signal = molecule_signal                      # signal in terms of number of loops
        