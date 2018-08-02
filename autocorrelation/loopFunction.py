""" create the array of MS2 loops, L, as a function of position, i
    ***NOTE*** in Desponds et al: "Li" ~should really be written L(i) 
    (adapted from Desponds et al Matlab code) """
import scipy.io as spio
import numpy as np

class DespondsMs2Loops():

    def __init__(self, despondsfile, tPol, k_elong):
        # load in Desponds et al's loop function
        sizePol = tPol * k_elong
        loopFn = spio.loadmat(despondsfile)
        
        ms = loopFn['ms']
        ms = ms[0]              # stupid matLab puts an array inside an array
        
        
        # Interpolate the loop function over polII positions on gene to aid in computation time
        # L(i) function in Desponds paper with i = polII position
        Li_fn = []
        for i in range(len(ms)//sizePol):    
            Li_fn.append(np.sum(ms[(sizePol*(i-1)+1) : (sizePol*i)]) / sizePol)
        
        #if (i < len(ms) // sizePol):
            #Li_fn[i] = np.sum(ms[(sizePol*(i)+1):-1]) / sizePol
            
        self.loopsByPolPosition = np.asarray(Li_fn)
        self.loopsByBp = ms

class SnailPromoterMs2Loops():
    '''
    Create the loop agglomeration pattern for Snail promoter with
    known polII speed (k_elong) and polII loading time (tPol)
    '''
    def __init__(self):
        genelength = 5600    # in basepairs. Snail MS2 loads GFP loops from 5' end of gene
        maxloops = 24        # max loops in MS2 chain
        loadlength = 60      # how many bp it takes to start chain
        looplength = 19      # how many bp generates one loop 
        loopchain = np.zeros(genelength - loadlength)   # look at just the loop part of gene

        counter = 0
        loopnum = 0
        waitbp = 20
        for bp in range(len(loopchain)):
            if (counter / looplength >= 1):
                loopnum += 1
                counter = 0
                
                if (loopnum % 2 == 0):   # on even loop numbers, distance between loops is 50
                    waitbp = 50
                else:
                    waitbp = 20

                if loopnum > maxloops:    # When we reach the maximum loops in the chain, don't add any more.
                    loopnum = maxloops

                loopchain[bp:(bp + waitbp)] = loopnum
            else:
                counter += 1
        self.loop_function = np.append(np.zeros(loadlength), loopchain)   # add in the start region at beginning of gene
        


    # Now interpolate loop function by polII position to cut down on computation time
def loopInterpolate(loop_function, k_elong, tPol):
    sizePol = k_elong * tPol
    loop_fn_interp = []
    counter=0
    for i in range(len(loop_function)):
        if counter / sizePol >= 1:
            loop_fn_interp.append(loop_function[i])
            counter=0
        else:
            counter += 1
    return np.asarray(loop_fn_interp)
        
