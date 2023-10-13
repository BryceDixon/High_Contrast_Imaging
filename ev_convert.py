"""
Code to convert stellar evolutionary models to contrast curves
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as cs

def ev_plot(age, mag, mass, save_folder = None, savefilename = None):
    """Function to plot the stellar evolutionary models

    Parameters
    ----------
    age : list of floats
        list of the ages for each point in the same order as the mag array
    mag : array of floats
        array of the magnitude data for all tracks in the form n columns = mass, m rows = age, (one track should be one whole column)
    mass : list of floats
        list of the masses for each track in the same order as the mag array
    """
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1,1,1)
    
    for i in range(len(mag[0,:])):
        new_age = []
        new_mag = []
        pos = np.isnan(mag[:,i])
        for j in range(len(pos)):
            if pos[j] == False:
                new_age.append(age[j])
                new_mag.append(mag[j,i])
        
        spline = cs(new_age, new_mag, bc_type = 'natural')
        x = np.geomspace(new_age[0], new_age[-1], 1000)
        y = spline(x)
        ax.plot(x, y, label = mass[i])
    
    ax.set_xscale('log')
    ax.invert_yaxis()
    ax.legend(title = 'Mass ($M_s$)')
    ax.set_xlabel('Age (Gyr)')
    ax.set_ylabel('Absolute Magnitude ($M_{ll}$)')
    ax.set_title('Stellar Evolutionary Models')
    
    if save_folder is not None:
        assert savefilename is not None, "You need to give both save_folder and savefilename to save the figure"
        plt.savefig(str(save_folder)+"/"+str(savefilename)+".png", bbox_inches='tight')

    if savefilename is not None and save_folder is None:
        print("Figure not saved, you need to provide savefilename and save_folder both") 


