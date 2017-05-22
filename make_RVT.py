# -*- coding: utf-8 -*-
#      author:         Dr. Jason E. Hill (post-doc fellow with the CNT at TTU)
#      updated:        16 MAY 2017
from pandas import *
from pylab import *
from make_pulse_timeseries  import make_pulse_timeseries

def make_RVT(dt,Nt,RTI = 4.0,RTI_sigma = 0.25,weight = 75.0,A_sigma = 0.005, impulse_seed = set(),dz_seed = set()):
    """ Generates a respiratory volume time-series (RVT) characterized by 
    average respiratory time RTI with variance var_RTI for a person of 
    a specific weight [kg] (optional), variance in amplitude A_sigma, with 
    random seeds for the impulse time series generation, impulse_seed, and
    for the chest displacement time series generation, dz_seed.
    Outputs are in units of mL, listed as: 
    FRC - Functional Residual Capacity, the lung volume left after a normal exhalation
    TIV - Typical Inspiratory Volume, the lung volume after normal inhalation.
    IRV - Inspiratory Reserve Volume (IRV), the lung volume after taking the
          largest inhalation possible
    """

    """
    NOTE: in regular respiration, (A. E. Lujan, J. M. Balter and R. K. Ten Haken  
    "A method for incorporating organ motion due to breathing into 3D
    dose calculations in the liver: sensitivity to variations in motion" 
    published in Med. Phys. 30 2643–9 [2003]) 
    a single respiratory cycle is represented with a functional form:
    z(phi) = A*(cos(pi*phi)^4)
    """

    # generate the respiratory motion time-series

    z_ts = make_pulse_timeseries(dt,Nt,RTI,RTI_sigma,p = 4.0,A_sigma = A_sigma,impulse_seed = impulse_seed,dz_seed = dz_seed)

    """ Average lung volume parameters as detailed in the thesis 
    "Simulation of an Artificial Respiratory System"
    by A. Delawari & R. Doelman (2010). [Weight = 100kg]
    Weight proportions from 
    https://en.wikipedia.org/wiki/Lung_volumes
    https://en.wikipedia.org/wiki/Lung_volumes#/media/File:Lungvolumes_Updated.png
    -----------------------------------------------------------------------
    Model lungs as an ellipsoid with dimensions a = 11, b = 11, c = 6 for
    a  fter regular exhale V = ERV + RV ~ 3.0 L.
    NAME                                    VOLUME[L]@100kg  [mL/kg] DeltaR
    Residual Volume (RV)                          1.5         15     -1.75
    Expirational Reserve Volume (ERV)             1.5         15
    Functional Residual Capacity (FRC) = ERV + RV 3.0         30      0.00      
    Tidal volume (typical breath) (TV or V_T)     0.7         7 
    Typical inhalation volume = V_T + ERV + RV    3.7         37     +0.58
    Inspirational Reserve Volume (IRV)            4.3         43
    Vital Capacity (VC)  = IRV + TV + ERV         6.5         65
    Total Lung Capacity (TLC) = IRV+TV+ ERV+ RV = 8.0         80     +3.3
    NOTE: average weight is 75 kg -> TLC = 6.0 L
    """  

    a = 11.0;    b = 11.0;    c = 6.0
    deltaR = 0.58*z_ts.values

    FRC = (4.0*pi*weight/300.0)*a*b*c
    TIV = (4.0*pi*weight/300.0)*(a+0.58)*(b+0.58)*(c+0.58)
    RVT_Data = list((4.0*pi*weight/300.0)*(a+deltaR)*(b+deltaR)*(c+deltaR))

    ts_time = arange(dt,(Nt+1)*dt,dt)
    
    RVT = Series(RVT_Data, index=ts_time)
#    RVT.TimeInfo.Units = 'seconds'
#    RVT = setuniformtime(RVT,'StartTime',dt)
#    RVT = setuniformtime(RVT,'Interval',dt)

    return RVT, FRC, TIV