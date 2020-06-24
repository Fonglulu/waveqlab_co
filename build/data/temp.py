# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 14:31:06 2020

@author: shilu
"""
import numpy as np
import matplotlib.pyplot as plt
def load_data(filename1, filename2):
    
    
    data1 = np.loadtxt(filename1, dtype = float)
    data2 = np.loadtxt(filename2, dtype = float)
    
    # time steps
    steps1 = np.shape(data1)[0]
    steps2 = np.shape(data2)[0]
    print(steps1)
 
    time1 = round(data1[-1,0],1)
    time2 = round(data2[-1,0],1)


    time1 = np.linspace(0,time1, steps1 )
    time2 = np.linspace(0,time2, steps2 )
    
    # Choose data in x_direction
    x_dir1 = data1[:,2]
    x_dir2 = data2[:,2]
    
    fig, ax = plt.subplots()
    
    ax.plot(time1, x_dir1 ,  time2, x_dir2)

    ax.set(xlabel='time (s)', ylabel='magnitude',
       title='someting')

    plt.show()
    