import astropy.io.fits as fits
import glob
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile
import os
import numpy as np

basedir = '/n/fink2/ameisner/CI/v0001'

plot_num_dict = {'CIE': 4, 'CIN': 2, 'CIC': 5, 'CIS': 8, 'CIW': 6}

def rotate_1psf(stamp, extname):
    print('blat')

def plot_1exp(fname, outdir='/n/fink2/www/ameisner/ci_psfs'):
    hdul = fits.open(fname)

    plt.figure(figsize=(6, 6))

    for hdu in hdul:
        h = hdu.header
        print h['EXTNAME']
        stamp = hdu.data
        print stamp.shape

        if h['EXTNAME'] == 'CIE':
            stamp = np.rot90(stamp)
        if h['EXTNAME'] == 'CIW':
            stamp = np.rot90(stamp, k=3)
        if (h['EXTNAME'] == 'CIN') or (h['EXTNAME'] == 'CIC'):
            stamp = np.rot90(stamp, k=2)

        plt.subplot(3, 3, plot_num_dict[h['EXTNAME']])
        vmin = scoreatpercentile(np.ravel(stamp), 3.0)
        vmax = scoreatpercentile(np.ravel(stamp), 97.0)
        
        plt.imshow(stamp, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
        plt.title(h['EXTNAME'])
        plt.axis('off')

    print('~'*80)
    outname = os.path.basename(fname)
    outname = outname.replace('.fits', '.png')
    outname = os.path.join(outdir, outname)
    print(outname)
    plt.savefig(outname) # , bbox_inches='tight')
    plt.cla()

def _loop():
    flist = glob.glob(basedir + '/*/*/*psf*.fits')

    for f in flist:
        plot_1exp(f)

def _jump():
    fname = '/n/fink2/ameisner/CI/v0001/20190412/ci-00006319/ci-00006319_psf-a.fits'
    plot_1exp(fname, outdir='/n/fink2/www/ameisner')

# CIS orientation is correct without doing any rotations/flips !!!
