import numpy as np
import ev_convert as ev
import matplotlib.pyplot as plt
from exo_dmc_custom import *
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import species
species.SpeciesInit()
import astropy.io.fits as fits

class det_prob:
    
    def __init__(self, model, model_tag):
        
        self.database = species.Database()
        self.database.add_isochrones(model=model)
        self.read_iso = species.ReadIsochrone(tag=model_tag)
        
    def filter_list(self):
        print(self.read_iso.get_filters())
    
    def mass_cont(self, foldername, filename, distance, age, stellar_mag, filter_name, band, savefolder):
        
        data = fits.getdata(str(foldername)+'/'+str(filename), ext=0)
        x = data[0,:]
        y = data[1,:]
        x = ev.arc_to_au(x, distance)
        masses = self.read_iso.contrast_to_mass(age=age,
                                   distance=distance,
                                   filter_name=filter_name,
                                   star_mag=stellar_mag,
                                   contrast=y,
                                   use_mag=False)
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1,1,1)
        ax.plot(x, masses)
        ax.set_ylabel('Mass ($M_J$)')
        ax.set_xlabel('Separation (AU)')
        ax.set_title('{}-Band Mass contrast curve for a {} Myrs, {} Mag Star at a distance of {} pc'.format(band, age, stellar_mag, distance))
        savefilename = '{}_{}mag_{}myrs_{}pc'.format(filter_name,stellar_mag,age,distance)
        plt.savefig(str(savefolder)+'/'+str(savefilename)+'_curve.png', bbox_inches='tight')
        
        self.band = band
        self.age = age
        self.stellar_mag = stellar_mag
        self.distance = distance
        self.x = x
        self.masses = masses
        self.savefolder = savefolder
        self.savefilename = savefilename
    
    def det_prob(self):
        
        ID='{}-Band Detection Probability a {} Myrs, {} Mag Star at {} pc'.format(self.band, self.age, self.stellar_mag, self.distance)
        dist=([self.distance])
        map=exodmc(ID, dist)
        map.set_grid(x_min=0.1, x_max=100, logx=True)
        xlim=self.x
        ylim=self.masses
        prob = map.DImode(xlim, ylim, lxunit = 'au', savefig = False)
        plt.savefig(str(self.savefolder)+'/'+str(self.savefilename)+'_detprob.png', bbox_inches='tight')


def plot(model, model_tag, foldername, file_list, distance, age, stellar_mag_list, filter_name_list, band_list, savefolder):
    
    assert len(file_list) == len(stellar_mag_list) == len(filter_name_list) == len(band_list), "All lists must be the same length"
    
    _ = det_prob(model, model_tag)
    
    for N,i in enumerate(file_list):
        
        _.mass_cont(foldername, i, distance, age, stellar_mag_list[N], filter_name_list[N], band_list[N], savefolder)
        _.det_prob()