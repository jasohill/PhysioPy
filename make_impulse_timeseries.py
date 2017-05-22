# -*- coding: utf-8 -*-
#   Author:  Dr. Jason E. Hill (post-doc fellow with the CNT at TTU)
#   Updated: 17 MAY 2017
from pandas import *
from pylab import *
from random import *
import numpy as np


def make_impulse_timeseries(dt,Nt,TI,TI_sigma,rng_seed = set()):
    """ Generates an impulse time-series characterized by 
    the sampling time dt over a total of Nt time points,
    and average interval time TI with variance varTI and optional random number generator seed rng_seed.
    NOTE: the resulting time-series will have a Gaussian spectrum.
    """

    if (rng_seed != set()):
        seed(rng_seed)

    TT = dt*Nt
    Nk = ceil(TT*TI)
    Nk = int(ceil(Nk*(1.1 + 2*TI_sigma))) # initial estimation of the number of intervals
    TIs = TI_sigma*standard_normal((Nk,)) + TI
    TIts = TIs.cumsum(axis=0)
    (indices,) = (TIts>TT).nonzero()
    indices = list(indices)
    lastIndex = indices[0]
    Nk = lastIndex-1
    TIts = TIts[0:Nk]
    TItis = zeros(Nk)
    TItis = np.round(TIts/dt)

    # generate impulse time-series
    ts_data = zeros(Nt)
    for i in range(0,Nk):
        ts_data[int(TItis[i])] = 1
    
    ts_time = arange(dt,(Nt+1)*dt,dt)
    ts = Series(ts_data,ts_time)
    # ts.TimeInfo.Units = 'seconds';
    # ts = setuniformtime(ts,'StartTime',dt);
    # ts = setuniformtime(ts,'Interval',dt);
    
    return ts