from astropy.io import fits
import matplotlib.pyplot as plt
import os

#os.chdir(os.path.normpath(os.path.expandvars('$HOME/INSTRUMENTS/METIS/heeps_metis/cc_preview/exeter')))
os.chdir(os.path.normpath(os.path.expandvars('$HOME/dox/METIS/METIS_HCI_HEEPS_simulations/FDR/exeter/fits')))

xlabel = 'Angular separation $[arcsec]$'
ylabel_adi = '5-$\sigma$ sensitivity (contrast)'
ylabel_raw = 'raw contrast'
r_CLC = 0.06

adi1 = fits.getdata('cc_adi_bckg1_mag7_L_CVC_exeter_all_effects.fits')
adi2 = fits.getdata('cc_adi_bckg1_mag7.5_L_CVC_exeter_all_effects.fits')
adi3 = fits.getdata('cc_adi_bckg1_mag8_L_CVC_exeter_all_effects.fits')
adi4 = fits.getdata('cc_adi_bckg1_mag7_M_CVC_exeter_all_effects.fits')
adi5 = fits.getdata('cc_adi_bckg1_mag7.5_M_CVC_exeter_all_effects.fits')
adi6 = fits.getdata('cc_adi_bckg1_mag8_M_CVC_exeter_all_effects.fits')
adi7 = fits.getdata('cc_adi_bckg1_mag7_N2_IMG_exeter_all_effects.fits')
adi8 = fits.getdata('cc_adi_bckg1_mag7.5_N2_IMG_exeter_all_effects.fits')
adi9 = fits.getdata('cc_adi_bckg1_mag8_N2_IMG_exeter_all_effects.fits')

figsize=(8,6)

plt.figure(figsize=figsize)
plt.plot(adi1[0], adi1[1], label='mag L=7')
plt.plot(adi2[0], adi2[1], label='mag L=7.5')
plt.plot(adi3[0], adi3[1], label='mag L=8')
plt.loglog()
plt.grid(True), plt.grid(which='minor', linestyle=':')
plt.legend()
plt.xlabel(xlabel)
plt.ylabel(ylabel_adi)
plt.title('L-band CVC with all effects')
plt.xlim(0.02, 0.75)
plt.ylim(1e-8,1e-2)
plt.xticks([0.02, 0.05, 0.1, 0.2, 0.5, 0.75])
plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter())
plt.savefig('exeter_L_CVC.png', transparent=True, dpi=300)

plt.figure(figsize=figsize)
plt.plot(adi4[0], adi4[1], label='mag M=7')
plt.plot(adi5[0], adi5[1], label='mag M=7.5')
plt.plot(adi6[0], adi6[1], label='mag M=8')
plt.loglog()
plt.grid(True), plt.grid(which='minor', linestyle=':')
plt.legend()
plt.xlabel(xlabel)
plt.ylabel(ylabel_adi)
plt.title('M-band CVC with all effects')
plt.xlim(0.02, 0.75)
plt.ylim(1e-8,1e-2)
plt.xticks([0.02, 0.05, 0.1, 0.2, 0.5, 0.75])
plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter())
plt.savefig('exeter_M_CVC.png', transparent=True, dpi=300)

plt.figure(figsize=figsize)
plt.plot(adi7[0], adi7[1], label='mag N=7')
plt.plot(adi8[0], adi8[1], label='mag N=7.5')
plt.plot(adi9[0], adi9[1], label='mag N=8')
plt.loglog()
plt.grid(True), plt.grid(which='minor', linestyle=':')
plt.legend()
plt.xlabel(xlabel)
plt.ylabel(ylabel_adi)
plt.title('N-band IMG (no vortex) with all effects')
plt.xlim(0.06, 1.2)
plt.ylim(1e-8,1e-2)
plt.xticks([0.06, 0.1, 0.2, 0.5, 1, 1.2])
plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter())
plt.savefig('exeter_N_IMG.png', transparent=True, dpi=300)
