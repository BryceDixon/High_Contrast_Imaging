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

class det_prob:
    
    def __init__(self, model, model_tag):
        
        self.database = species.Database()
        self.database.add_isochrones(model=model)
        self.read_iso = species.ReadIsochrone(tag=model_tag)
        
    def filter_list(self):
        print(self.read_iso.get_filters())
    
    def mass_cont(self, foldername, filename, distance, age, stellar_mag, filter_name, band, savefolder, star_name):
        
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
        cont = np.array([x,masses])
        ax.plot(x, masses)
        ax.set_ylabel('Mass ($M_J$)')
        ax.set_xlabel('Separation (AU)')
        if star_name is None:
            ax.set_title('{}-Band Mass contrast curve for a {} Myrs, {} Mag Star at a distance of {} pc'.format(band, age, stellar_mag, distance))
            savefilename = '{}_{}mag_{}myrs_{}pc'.format(filter_name,stellar_mag,age,distance)
        else:
            ax.set_title('{}-Band Mass contrast curve for {}, a {} Myrs, {} Mag Star at a distance of {} pc'.format(band, star_name, age, stellar_mag, distance))
            savefilename = '{}_{}_band'.format(star_name,filter_name)
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
        plt.savefig(str(self.savefolder)+'/'+str(self.savefilename)+'_detprob.png', bbox_inches='tight')
        det_prob = os.path.join(self.savefolder, str(self.savefilename)+'_detprob.txt')
        np.savetxt(det_prob, prob[0])
        


def plot(model, model_tag, foldername, file_list, distance, age, stellar_mag_list, filter_name_list, band_list, savefolder, xmin = 0.1, xmax = 100, logx = True, star_list = None):
    
    assert len(file_list) == len(stellar_mag_list) == len(filter_name_list) == len(band_list), "All lists must be the same length"
    
    _ = det_prob(model, model_tag)
    
    for N,i in enumerate(file_list):
        
        if star_list is None:
            _.mass_cont(foldername, i, distance[N], age[N], stellar_mag_list[N], filter_name_list[N], band_list[N], savefolder, None)
        else:
            _.mass_cont(foldername, i, distance[N], age[N], stellar_mag_list[N], filter_name_list[N], band_list[N], savefolder, star_list[N])
        _.det_prob(xmin, xmax, logx)