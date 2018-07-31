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
        
        
        # make the L(i) function (array of values)
        Li_fn = []
        for i in range(len(ms)//sizePol):    
            Li_fn.append(np.sum(ms[(sizePol*(i-1)+1) : (sizePol*i)]) / sizePol)
        
        #if (i < len(ms) // sizePol):
            #Li_fn[i] = np.sum(ms[(sizePol*(i)+1):-1]) / sizePol
            
        self.loopsByPolPosition = np.asarray(Li_fn)
        self.loopsByBp = ms

#class SnailPromoterMs2Loops():
    '''
    Create the loop agglomeration pattern for Snail promoter with
    known polII speed (k_elong) and polII loading time (tPol)
    '''
    #def __init__(self, tPol, k_elong):
         #self.tPol = tPol
        # self.k_elong = k_elong
