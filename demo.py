# -*- coding: utf-8 -*-
"""
Created on Mon Jun 05 14:42:40 2017

@author: Jason E. Hill
"""
from make_RVT import make_RVT
from pylab import *

RVT, FRC, TIV = make_RVT(0.05,1000,4.0,0.24,80)

figure(1)
plot(RVT)
ylabel('respiratory volume [L]')
xlabel('time [s]')
title('Respiratory Volume Time-series (RVT)')
show()

print 'The functional residual capacity is ', FRC, 'cc.'
print 'The typical inhaled volume is ', TIV, 'cc.'