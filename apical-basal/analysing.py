# -*- coding: utf-8 -*-
"""
Created on Tue May 28 12:25:52 2019

@author: isphording
"""

import os
import pickle
import time
import matplotlib as mpl
from matplotlib import pyplot as plt
from main import gif_from_data
import parameters as p

mpl.rcParams.update({'figure.max_open_warning': 0})

# For every folder:
#   Make a video
#   Plot #Stemcells over time
#   Plot

end = 0
for x in os.listdir('data'):
    if '.pickle' in x:
        end = max(int(x[:-7]), end)
end -= end%20

p.init_from_file('parameters.pickle')

## Timelapse
gif_from_data('data')

time.sleep(0.1)
