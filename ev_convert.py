"""
Code to convert stellar evolutionary models to contrast curves

Contains:
    contrast_convert class
        ev_plot function
        ev_convert function
        mag_to_mass function
    convert function
    contrast_to_mag function
    arc_to_au function
    mag_contrast_plot function

Author: Bryce Dixon
Version: 16/10/2023
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as cs
from scipy.interpolate import interp1d as interp
import astropy.io.fits as fits

class contrast_convert:
    
    def __init__(self, age, mag, mass, save_folder = None):
        """Class to convert a stellar evolutionary model to a contrast curve
        
        Parameters
        ----------
        age : list of floats
            list of the ages for each point in the same order as the mag array
        mag : array of floats
            array of the magnitude data for all tracks in the form n columns = mass, m rows = age, (one track should be one whole column)
        mass : list of floats
            list of the masses for each track in the same order as the mag array
        save_folder : str, optional
            directory of the folder to save the plots to
        """
        self.age = age
        self.mag = mag
        self.mass = mass
        self.save_folder = save_folder
        
    def ev_plot(self, savefilename = None):
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
            assert savefilename is not None, "You need to give both save_folder and savefilename to save the figure"
            plt.savefig(str(self.save_folder)+"/"+str(savefilename)+".png", bbox_inches='tight')

        if savefilename is not None and self.save_folder is None:
            print("Figure not saved, you need to provide savefilename and save_folder both")
    
    
    def ev_convert(self, age_track, savefilename = None):
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
            assert savefilename is not None, "You need to give both save_folder and savefilename to save the figure"
            plt.savefig(str(self.save_folder)+"/"+str(savefilename)+".png", bbox_inches='tight')

        if savefilename is not None and self.save_folder is None:
            print("Figure not saved, you need to provide savefilename and save_folder both")


    def mag_to_mass(self, mag_data, sep_data, stellar_mag, distance, savefilename = None):
        """Function to take the absolute magnitude and separation data from a contrast plot and convert the absolute magnitudes to a mass then plot mass against seperation

        Parameters
        ----------
        mag_data : array of floats
            array of absolute magnitudes from the contrast data
        sep_data : array of floats
            array of separations in AU from the contrast data
        stellar_mag : float
            absolute magnitude of the star
        distance : float
            distance to the target in parsecs
        """
        
        f = interp(self.y, self.x, kind = 'cubic')
        #mass_data = f(mag_data)
        mass_data = []
        new_sep_data = []
        for i in range(len(mag_data)):
            try:
                mass_data.append(f(mag_data[i]))
                new_sep_data.append(sep_data[i])
            except:
                pass
        
        # convert solar masses to jupiter masses
        mass_data = np.array(mass_data)*1047
        
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1,1,1)
        
        ax.plot(new_sep_data, mass_data)
        
        ax.set_ylabel('Mass ($M_J$)')
        ax.set_xlabel('Separation (AU)')
        ax.set_title('Mass contrast curve for a {} Gyrs, {} Mag Star at a distance of {} pc'.format(self.age_track, stellar_mag, distance))
        
        if self.save_folder is not None:
            assert savefilename is not None, "You need to give both save_folder and savefilename to save the figure"
            plt.savefig(str(self.save_folder)+"/"+str(savefilename)+".png", bbox_inches='tight')

        if savefilename is not None and self.save_folder is None:
            print("Figure not saved, you need to provide savefilename and save_folder both")



def convert(age, mag, mass, age_track, mag_data, sep_data, stellar_mag, distance, save_folder = None):
    """Function to plot the given evolutionary model, then convert the absolute magnitudes of a given contrast curve to masses to plot a curve of mass against seperation for a given age

    Parameters
    ----------
    age : list of floats
        list of the ages for each point in the same order as the mag array
    mag : array of floats
        array of the magnitude data for all tracks in the form n columns = mass, m rows = age, (one track should be one whole column)
    mass : list of floats
        list of the masses for each track in the same order as the mag array
    save_folder : str, optional
        directory of the folder to save the plots to
    age_track : float
        desired stellar age to plot the mass-magnitude relationship for in Gyrs
    mag_data : array of floats
        array of absolute magnitudes from the contrast data
    sep_data : array of floats
        array of separations in AU from the contrast data
    stellar_mag : float
        absolute magnitude of the star
    distance : float
        distance to the target in parsecs
    save_folder : str, optional
        directory of the folder to save the plots to, defaults to None
    """
    
    cont = contrast_convert(age, mag, mass, save_folder)
    cont.ev_plot('evolutionary_plot')
    cont.ev_convert(age_track, 'mag_vs_mass_plot')
    cont.mag_to_mass(mag_data, sep_data, stellar_mag, distance, 'mass_vs_sep_plot')
    



def contrast_to_mag(stellar_mag, contrasts):
    """Function to convert an array of contrasts to an array of magnitudes for a given stellar magnitude

    Parameters
    ----------
    stellar_mag : float
        absolute magnitude of the star
    contrasts : array of floats
        array of contrasts to be converted to magnitudes
    zero_point_lum : float, optional
        zero point luminosity, defaults to 3.0128e+28
    
    Returns
    -------
    M_array : array of floats
        array of absolute magnitudes converted from the contrasts
    """
    
    contrasts = np.array(contrasts)
    zero_point_lum = 3.0128e+28
    
    L_star = zero_point_lum * 10 ** (-stellar_mag/2.5)
    L_array = L_star * contrasts
    M_array = -2.5 * np.log10(L_array/zero_point_lum)
    
    return M_array



def arc_to_au(arcseconds, distance):
    """Function to convert an array of separations in arcseconds to an array of separations in AU for a given distance

    Parameters
    ----------
    arcseconds : array of floats
        array of separations in arcseconds
    distance : float
        distance to the target in parsecs
    
    Returns 
    -------
    AU_array : array of floats
        array of separations in AU
    """
    
    arcseconds = np.array(arcseconds)
    
    AU_array = distance * arcseconds
    
    return AU_array


def mag_contrast_plot(fits_data, stellar_mag, distance, band, residual_type):
    """Function to convert a .fits contrast curve from sensitivity (contrast) against separation in arcseconds to absolute magnitude against separation in AU and plot it

    Parameters
    ----------
    fits_data : str
        directory for the .fits file 
    stellar_mag : float
        absolute magnitude of the star
    distance : float
        distance to the target in parsecs
    band : str
        name of the bandwidth filter
    residual_type : str
        type of residuals present in the data

    Returns
    -------
    mag_data : array of floats
        array of the absolute magnitudes
    sep_data : array of floats
        array of the separations in parsecs
    """
    data = fits.getdata(fits_data, ext=0)
    separation = data[0,:]
    contrasts = data[1,:]
    mag_data = contrast_to_mag(stellar_mag, contrasts)
    sep_data = arc_to_au(separation, distance)
    
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1,1,1)
    
    ax.plot(sep_data,mag_data)
    ax.invert_yaxis()
    ax.set_xlabel('Angular Separation (AU)')
    ax.set_ylabel('Absolute Magnitude')
    ax.set_title('{}-band, {} mag, {}, at {} parsec'.format(band, stellar_mag, residual_type, distance))
    
    return mag_data, sep_data