"""
Code to convert stellar evolutionary models to contrast curves
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as cs
from scipy.interpolate import interp1d as interp

class contrast_convert:
    
    def __init__(self, age, mag, mass, save_folder = None, savefilename = None):
        """Class to convert a stellar evolutionary model to a contrast curve
        
        Parameters
        ----------
        age : list of floats
            list of the ages for each point in the same order as the mag array
        mag : array of floats
            array of the magnitude data for all tracks in the form n columns = mass, m rows = age, (one track should be one whole column)
        mass : list of floats
            list of the masses for each track in the same order as the mag array
        """
        self.age = age
        self.mag = mag
        self.mass = mass
        
    def ev_plot(self):
        """Function to plot the stellar evolutionary models
        """
    
        age = self.age
        mag = self.mag
        mass = self.mass
        
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
    
        if self.save_folder is not None:
            assert self.savefilename is not None, "You need to give both save_folder and savefilename to save the figure"
            plt.savefig(str(self.save_folder)+"/"+str(self.savefilename)+".png", bbox_inches='tight')

        if self.savefilename is not None and self.save_folder is None:
            print("Figure not saved, you need to provide savefilename and save_folder both")
    
    
    def ev_convert(self, age_track):
        """Function to convert the evolutionary model to a plot of mass against magnitude for a given age

        Parameters
        ----------
        age_track : float
            desired stellar age to plot the mass-magnitude relationship for in Gyrs
        """
        age = self.age
        mag = self.mag
        mass = self.mass
        
        mag_track = []
        mass_track = []
        
        assert age_track <= max(age), "Cannot produce a plot of an age larger than the maximum model age"
        assert age_track >= min(age), "Cannot produce a plot of an age smaller than the minimum model age"
        
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
            x = np.geomspace(min(new_age), max(new_age), 1000)
            y = spline(x)
            
            try:
                f = interp(x, y, kind='cubic')
                new_mag = f(age_track)
                mag_track.append(new_mag)
                mass_track.append(mass[i])
            except:
                pass
        
        spline = cs(mass_track, mag_track, bc_type = 'natural')
        x = np.geomspace(min(mass_track), max(mass_track), 1000)
        y = spline(x)
        
        ax.plot(x,y)
            
        ax.set_xscale('log')
        ax.invert_yaxis()
        ax.set_xlabel('Mass ($M_s$)')
        ax.set_ylabel('Absolute Magnitude ($M_{ll}$)')
        ax.set_title('Absolute Magnitude against Mass for a {} Gyrs Star'.format(age_track))
        
        self.x = x
        self.y = y
        self.age_track = age_track
        
        if self.save_folder is not None:
            assert self.savefilename is not None, "You need to give both save_folder and savefilename to save the figure"
            plt.savefig(str(self.save_folder)+"/"+str(self.savefilename)+".png", bbox_inches='tight')

        if self.savefilename is not None and self.save_folder is None:
            print("Figure not saved, you need to provide savefilename and save_folder both")


    def mag_to_mass(self, mag_data, sep_data):
        """Function to take the absolute magnitude and separation data from a contrast plot and convert the absolute magnitudes to a mass then plot mass against seperation

        Parameters
        ----------
        mag_data : array of floats
            array of absolute magnitudes from the contrast data
        sep_data : array of floats
            array of separations in AU from the contrast data
        """
        
        f = interp(self.y, self.x, kind = 'cubic')
        mass_data = f(mag_data)
        
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1,1,1)
        
        ax.plot(sep_data, mass_data)
        
        ax.set_yscale('log')
        ax.set_ylabel('Mass ($M_s$)')
        ax.set_xlabel('Separation (AU)')
        ax.set_title('Mass against Separation for a {} Gyrs Star'.format(self.age_track))
        
        if self.save_folder is not None:
            assert self.savefilename is not None, "You need to give both save_folder and savefilename to save the figure"
            plt.savefig(str(self.save_folder)+"/"+str(self.savefilename)+".png", bbox_inches='tight')

        if self.savefilename is not None and self.save_folder is None:
            print("Figure not saved, you need to provide savefilename and save_folder both")