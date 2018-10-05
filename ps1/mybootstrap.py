import numpy as np
import scipy.stats as spstat


def BootstrapResample(x,NSamples):
    #Does resampling
    #x: sample  1 X N array
    #NSamples: number of bootstrap samples
    sampled = np.random.choice(x, size=NSamples*x.shape[0], replace=True)
    sampled = np.reshape(sampled, (x.shape[0],NSamples))
    
    return sampled


def BootstrapTTableCIMean(x,NSamples,confidence=0.9):
    #Computes CI for mean using t-table of bootstrap means
    #x: sample 1 X N
    #NSamples: number of bootstrap samples
    #confidence: confidence level, 0.9 by default
    
    m0=np.mean(x)
    s0=np.std(x)

    resampled = BootstrapResample(x,NSamples) # computes resample, then mean
    
    ts = (np.mean(resampled,axis=0)-m0)/np.std(resampled,axis=0)
    ts = np.sort(ts)
    ci= [ m0+ts[np.floor((1-confidence)/2*NSamples).astype('int32')]*s0,\
          m0+ts[np.floor( (confidence+(1-confidence)/2)*NSamples).astype('int32')]*s0]
    
    return ci
  
