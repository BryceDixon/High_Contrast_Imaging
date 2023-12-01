import numpy as np
import ev_convert as ev
import matplotlib.pyplot as plt
from exo_dmc_custom import *
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import species
species.SpeciesInit()
import astropy.io.fits as fits
import os
from scipy.interpolate import CubicSpline as cs
from scipy.interpolate import interp1d as interp

class det_prob:
    
    def __init__(self, model, model_tag):
        
        self.database = species.Database()
        self.database.add_isochrones(model=model)
        self.read_iso = species.ReadIsochrone(tag=model_tag)
        
    def filter_list(self):
        print(self.read_iso.get_filters())
    
    def mass_cont(self, foldername, file_list, file_mag_list, distance, age, stellar_mag, filter_name, band, savefolder, star_name):
        
        x, y = generate_contrast(file_list, foldername, file_mag_list, stellar_mag)
        x = ev.arc_to_au(x, distance)
        masses = self.read_iso.contrast_to_mass(age=age,
                                   distance=distance,
                                   filter_name=filter_name,
                                   star_mag=stellar_mag,
                                   contrast=y,
                                   use_mag=False)
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1,1,1)
        cont = np.array([x,masses])
        ax.plot(x, masses)
        ax.set_ylabel('Mass ($M_J$)')
        ax.set_xlabel('Separation (AU)')
        if star_name is None:
            ax.set_title('{}-Band Mass contrast curve for a {} Myrs, {} Mag Star at a distance of {} pc'.format(band, age, stellar_mag, distance))
            if filter_name.find('/') != -1:
                savefilename = '{}_{}mag_{}myrs_{}pc'.format(band,stellar_mag,age,distance)
            else:
                savefilename = '{}_{}mag_{}myrs_{}pc'.format(filter_name,stellar_mag,age,distance)
        else:
            ax.set_title('{}-Band Mass contrast curve for {}, a {} Myrs, {} Mag Star at a distance of {} pc'.format(band, star_name, age, stellar_mag, distance))
            if filter_name.find('/') != -1:
                savefilename = '{}_{}_band'.format(star_name,band)
            else: 
                savefilename = '{}_{}_band'.format(star_name,filter_name)
        if savefolder is not None:
            plt.savefig(str(savefolder)+'/'+str(savefilename)+'_curve.png', bbox_inches='tight')
            mass_cont = os.path.join(savefolder, str(savefilename)+'_cont_data.txt')
            np.savetxt(mass_cont, cont)
        
        self.band = band
        self.age = age
        self.stellar_mag = stellar_mag
        self.distance = distance
        self.x = x
        self.masses = masses
        self.savefolder = savefolder
        self.savefilename = savefilename
        self.star_name = star_name
    
    def det_prob(self, xmin, xmax, logx):
        
        if self.star_name is None:
            ID='{}-Band Detection Probability for a {} Myrs, {} Mag Star at {} pc'.format(self.band, self.age, self.stellar_mag, self.distance)
        else:
            ID='{}-Band Detection Probability for {}, a {} Myrs, {} Mag Star at {} pc'.format(self.band, self.star_name, self.age, self.stellar_mag, self.distance)
        dist=([self.distance])
        map=exodmc(ID, dist)
        map.set_grid(x_min=xmin, x_max=xmax, logx=logx)
        xlim=self.x
        ylim=self.masses
        prob = map.DImode(xlim, ylim, lxunit = 'au', savefig = False)
        prob = np.array(prob)
        if self.savefolder is not None:
            plt.savefig(str(self.savefolder)+'/'+str(self.savefilename)+'_detprob.png', bbox_inches='tight')
            det_prob = os.path.join(self.savefolder, str(self.savefilename)+'_detprob.txt')
            np.savetxt(det_prob, prob[0])
        


def plot(model, model_tag, foldername, file_list, file_mag_list, distance, age, stellar_mag_list, filter_name_list, band_list, savefolder = None, xmin = 0.1, xmax = 100, logx = True, star_list = None):
    """Function to plot detection probability plots for a given population of stars 

    Parameters
    ----------
    model : str
        name of the stellar evolution model
    model_tag : str
        tag of the stellar model
    foldername : str
        directory of the folder containing the .fits contrast curve files
    file_list : list of str
        list of the file names of the .fits contrast curve files
    file_mag_list : list of floats
        list of the magnitudes corresponding to the contrast curve fits files
    distance : list of floats
        list of the distances of each star
    age : list of floats
        list of the ages of each star
    stellar_mag_list : list of floats
        list of the magnitudes of each star
    filter_name_list : list of str
        list of the filter names to be used for each star
    band_list : list of str
        list of the bands to be used for each star
    savefolder : str, optional
        folder directory to save the plots to, defaults to None
    xmin : float, optional
        minimum x value for the detection probability plots, defaults to 0.1
    xmax : float, optional
        maximum x value for the detection probability plots, defaults to 100
    logx : bool, optional
        if True, plot x on a log scale, if False, do not, defaults to True
    star_list : list of str, optional
        list of the names of the stars, defaults to None
    """
    assert len(stellar_mag_list) == len(filter_name_list) == len(band_list), "All lists must be the same length"
    
    _ = det_prob(model, model_tag)
    
    for N in range(len(stellar_mag_list)):
        
        if star_list is None:
            _.mass_cont(foldername, file_list, file_mag_list, distance[N], age[N], stellar_mag_list[N], filter_name_list[N], band_list[N], savefolder, None)
        else:
            _.mass_cont(foldername, file_list, file_mag_list, distance[N], age[N], stellar_mag_list[N], filter_name_list[N], band_list[N], savefolder, star_list[N])
        _.det_prob(xmin, xmax, logx)



        
def generate_contrast(file_list, foldername, mag_list, desired_mag):
        
    sep_list = []
    cont_list = []
    for filename in file_list:
        data = fits.getdata(str(foldername)+'/'+str(filename), ext=0)
        sep = data[0,:]
        cont = data[1,:]
        if len(sep_list) == 0:
            sep_list = np.concatenate((sep_list, sep))
            cont_list = np.concatenate((cont_list, cont))
        else:
            sep_list = np.column_stack((sep_list, sep))
            cont_list = np.column_stack((cont_list, cont))
    
    new_cont_list = []
    for l in range(len(sep_list[:,0])):
        spline = cs(mag_list, cont_list[l,:], bc_type = 'natural')
        x = np.linspace(min(mag_list)-0.5, max(mag_list)+0.5, 1000)
        y = spline(x)
        
        f = interp(x,y, kind = 'cubic')
        new_cont = f(desired_mag)
        new_cont_list.append(new_cont)
    
    contrasts = np.array(new_cont_list)
    seperations = np.array(sep_list[:,0])
    
    #fig = plt.figure(figsize=(10, 10))
    #ax = fig.add_subplot(1,1,1)
    #ax.plot(sep_list[:,0], new_cont_list, label = desired_mag)
    #ax.set_yscale('log')
    #ax.legend()
            
    return seperations, contrasts
            