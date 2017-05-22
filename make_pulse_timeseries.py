# -*- coding: utf-8 -*- 
#    author:         Dr. Jason E. Hill (post-doc fellow with CNT at TTU)
#    updated:        17 MAY 2017
from pandas import *
from pylab import *
from random import *
import numpy as np
from make_impulse_timeseries  import make_impulse_timeseries 


def make_pulse_timeseries(dt,Nt,TI = 4.0,TI_sigma = 0.25,p = 4.0,A_sigma = 0.005, impulse_seed = set(),dz_seed = set()):
    """
    Generates a timeseries of pulsation events characterized by 
    interval time TI between pulses with standard deviation TI_sigma
    normal relaxed position is baseline z = 0. 
    Shape of pulses are half-cosines raised to the p-power
    Output is normalized to average amplitude.
    Defaults assume a typical respiratory pulse for healthy adult in supine position.
    """

    # make the impulse timeseries
    impulse_ts = make_impulse_timeseries(dt,Nt,TI,TI_sigma,rng_seed = impulse_seed)
    (indices,) = (impulse_ts == 1).nonzero()
    N = indices.size

    # calculate pulse intervals
    PIs = zeros(N)
    PIs[0:N-1] = indices[1:N] - indices[0:N-1]
    PIs[N-1] = PIs[N-2]
    PIs[1:N-1] = np.round(0.5*(PIs[0:N-2] + PIs[1:N-1]))
    PIs = 2*np.round(0.5*PIs)+1; #ensure odd durations 

    z_ts_Data = zeros(Nt)

    if (dz_seed != set()):
        seed(dz_seed)
        
    As = A_sigma*standard_normal((N,)) + 1.0*(PIs*dt)/TI
    if (indices[0] - 0.5*(PIs[0]+1)) < 0:
        start_i = 0.5*(PIs[0]+1)-indices[0]-1
        phi = array(range(start_i,PIs[0]))
        z_phi = As[0]*cos(pi*(phi - 0.5*(PIs[0]+1))/PIs[0])**p
        z_ts_Data[:len(phi)] = z_phi + z_ts_Data[:len(phi)]
    else:
        phi = array(range(0,int(PIs[0])))
        z_phi = As[0]*cos(pi*(phi - 0.5*(PIs[0]+1))/PIs[0])**p
        z_ts_Data[int((indices[0] - 0.5*(PIs[0]-1)) - 1):int(indices[0] + 0.5*(PIs[0]-1))] = z_phi + z_ts_Data[int((indices[0] - 0.5*(PIs[0]-1)) - 1):int(indices[0] + 0.5*(PIs[0]-1))]
 
    for n in range(1,N-1):
        phi = array(range(0,int(PIs[n])))
        z_phi = As[n]*cos(pi*(phi - 0.5*(PIs[n]+1))/PIs[n])**p
        z_ts_Data[int((indices[n] - 0.5*(PIs[n]-1)) - 1):int(indices[n] + 0.5*(PIs[n]-1))] = z_phi + z_ts_Data[int((indices[n] - 0.5*(PIs[n]-1)) - 1):int(indices[n] + 0.5*(PIs[n]-1))]

    if (indices[N-1] + 0.5*(PIs[N-1]+1)) > Nt:
        end_i = int(Nt-0.5*(PIs[N-1]+1)-indices[N-1])
        phi = array(range(0,end_i))
        z_phi = As[N-1]*cos(pi*(phi - 0.5*(PIs[N-1]+1))/PIs[N-1])**p
        z_ts_Data[int((indices[N-1] - 0.5*(PIs[N-1]+1)) - 1):end_i] = z_phi + z_ts_Data[int((indices[N-1] - 0.5*(PIs[N-1]-1)) - 1):end_i]
    else:
        phi = array(range(0,int(PIs[N-1])))
        z_phi = As[N-1]*cos(pi*(phi - 0.5*(PIs[N-1]+1))/PIs[N-1])**p
        z_ts_Data[int((indices[N-1] - 0.5*(PIs[N-1]-1)) - 1):int(indices[N-1] + 0.5*(PIs[N-1]-1))] = z_phi + z_ts_Data[int((indices[N-1] - 0.5*(PIs[N-1]-1))-1):int(indices[N-1] + 0.5*(PIs[N-1]-1))]
    
    ts_time = arange(dt,(Nt+1)*dt,dt)
    
    z_ts = Series(z_ts_Data, index=ts_time)
    #z_ts.TimeInfo.Units = 'seconds';
    #z_ts = setuniformtime(z_ts,'StartTime',dt);
    #z_ts = setuniformtime(z_ts,'Interval',dt);

    return z_ts
