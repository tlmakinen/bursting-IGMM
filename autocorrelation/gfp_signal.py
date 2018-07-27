"""Simulates GFP step-loop agglomeration according to a given promoter state binary telegraph signal """
import numpy as np
from loopFunction import ms2Loops  # get Desponds MS2 loop function

class gfp_signal:

    def __init__(self, telegraph, k_elong, tPol, max_loops, stepsize):

        despondsfile = 'standalone/therightL.mat'
        
        self.telegraph = telegraph
        self.k_elong = k_elong               # MS2 chain elongation rate; default is 25 bp/s
        self.tPol = tPol                     # default PolII loading time is 6s for desponds,3 for Lagha
        self.max_loops = max_loops           # 24 loops for MS2 system
        self.stepsize = stepsize             # observation time

        
        sizePol = k_elong * tPol            # "footprint" size of polII molecule on gene, in bp
        
        # now let's load in the stepwise MS2 loop agglomeration. Desponds et al
        # calculated this and created a nice little function of MS2 loops by basepair
        # along the gene. We're going to put this in terms of seconds using the polII 
        # rate, k_elong.
        ms2loops = ms2Loops(despondsfile, self.tPol, self.k_elong).loopsByBp
        gene = np.arange(len(ms2loops))        # GFP loop index, by basepair
        

        polII_per_step = tPol / stepsize              # loading rate of polII onto gene in terms of discrete timesteps
        molecule_signal = np.zeros(len(telegraph))    # empty array of mRNA-GFP loops
        
        polII_arr = []   # empty array of polII molecules that will grow according to the pol_per_step rate.
                        # This is a list of strands of molecules
        
        
        counter = 0      # keep track of how many time steps we've taken since we added the last polII molecule
        loopnum = 0      # keep track of the number of steps we've taken since we added the last loop
        pol_times = []
        
        for s in range(len(telegraph)):                   # iterate by time step
            
            position = int(s*k_elong)



            if telegraph[s] == 1:                         # When the telegraph signal is ON, start a new pol II molecule according to the rate pol_per_step
                counter += 1
                if (counter/tPol) >= 1:         # if we passed the start point of a new polII, add a new pol II molecule
                    polII_arr = np.append(polII_arr, 0)
                    pol_times = np.append(pol_times, 0)  # add each polII's coordinate to array
                    counter = 0                                  # restart counter
                            
                                               
            #polII_arr = [j+1 for j in polII_arr]          # add a new GFP loop to each polII molecule in the loop array according to GFP loop rate
            
            pol_times = [j+1 for j in pol_times]   # move each polII along in time
            pol_positions = [int(j * k_elong) for j in pol_times]    # each polII molecule's position, in bp

            for j in pol_positions:
                if j > np.max(gene):
                    j = 0                               # set position to 0 when the polII reaches the end of the gene

            for j in range(len(pol_positions)):
                if pol_positions[j] > len(ms2loops):
                     polII_arr[j] = 0                       # set to 0 when polII falls off gene

                else:
                    polII_arr[j] = ms2loops[pol_positions[j]]             # add a GFP loop according to Desponds
            #keep = (np.asarray(polII_arr) <= max_loops)   # if we've gone over the max number of GFP loops, set polII element to zero since 
                                                          # the chain of GFP attached to that polII will decay really quickly     
                                                                         
            #molecule_signal[s] = np.sum(keep*polII_arr)    # sum up the number of GFP molecules from the polII list
            molecule_signal[s] = np.sum(polII_arr)
            loopnum=0                                      # reset the loop counter
            
        self.signal = molecule_signal                      # signal in terms of number of loops
        