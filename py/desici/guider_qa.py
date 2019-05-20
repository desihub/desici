#!/usr/bin/env python

"""
Make guider QA plots from the NIGHT/EXPID/centroids-EXPID.json files
"""

import os, sys
import json
import argparse

import numpy as np
from astropy.table import Table
import fitsio
from pylab import * # sorry...
from matplotlib import cm

def main():
    parser = argparse.ArgumentParser(
        usage = "guider_qa [options]",
        description = "Runs QA on guider centroids",
        epilog = "At NERSC, you can provide either --infile or --expid and --night; at KPNO you must provide --infile")
    parser.add_argument("-i", "--infile", type=str,  help="input centroids file")
    parser.add_argument("-e", "--expid", type=int, help="exposure ID ")
    parser.add_argument("-n", "--night", type=int, help="YEARMMDD night of sunset")
    parser.add_argument("-o", "--outfile", help="output plot")

    args = parser.parse_args()

    if args.infile is None:
        if args.night is None or args.expid is None:
            print('Must provide --infile or --night plus --expid')
            return 1

        args.infile = os.path.join(
            os.getenv('DESI_ROOT'), 'spectro', 'data',
            str(args.night), '{:08d}'.format(args.expid),
            'centroids-{:08d}.json'.format(args.expid))
        if not os.path.exists(args.infile):
            print('Centroids file not found: {}'.format(args.infile))
            return 1

    #- Is outfile actually a directory?  if so, be clever
    if os.path.isdir(args.outfile):
        with open(args.infile) as fx:
            cx = json.load(fx)

        expid = cx['expid']
        if args.night is not None:
            night = args.night
        else:
            origdir = os.path.dirname(cx['filename'])
            night = os.path.basename(os.path.dirname(origdir))

        outpath = os.path.join(args.outfile, str(night))
        if not os.path.isdir(outpath):
            os.makedirs(outpath, exist_ok=True)

        args.outfile = os.path.join(outpath, 'offsets-{:08d}.png'.format(expid))

    data = load_centroids(args.infile)
    plot_offsets(data, args.outfile)

def load_centroids(filename):
    #- Load centroids
    cx = json.load(open(filename))

    #- ra, dec aren't in centroids file, so see if there is a guide file
    expid = cx['expid']
    dirname = os.path.dirname(filename)
    guidefile = os.path.join(dirname, 'guide-{:08d}.fits.fz'.format(expid))
    if os.path.exists(guidefile):
        hdr = fitsio.read_header(guidefile, 0)
        skyra, skydec = hdr['SKYRA'], hdr['SKYDEC']
    else:
        skyra, skydec = None, None

    #- platescale values taken from Echo22 design in
    #- desimodel/data/focalplane/platescale.txt
    center_platescale = 67.5 # um/arcsec
    radial_platescale = 76.3 # um/arcsec at r=400 mm
    az_platescale = 70.3 # um/arcsec at r=400m
    pixsize = 9 # um

    #- Convert json into Table; convert pixel offsets to arcsec
    rows = list()
    nframes = len(cx['frames'])
    for i, frame in enumerate(cx['frames'].values()):
        #- Use CIX and ROI=1 for the combined guider offset values
        ci = 'CIX'
        roi = 1

        #- NOTE: it appears that combined_y is aligned with the DEC axis
        #- *not* with the CIC y-axis which is -dec
        x_error = frame['combined_x']
        y_error = frame['combined_y']
        ra_err = pixsize * x_error / center_platescale
        dec_err = pixsize * y_error / center_platescale

        rows.append(dict(frame=i, ci=ci, roi=roi,
                         x_error=x_error, y_error=y_error,
                         ra_err=ra_err, dec_err=dec_err))

        #- Now loop over whatever cameras were included in this file
        for key, value in frame.items():
            if not key.startswith('CI'):
                continue

            ci = key[0:3]
            roi = int(key.split('_')[1])
            x_error = value['x_error']
            y_error = value['y_error']
            
            if ci == 'CIC':
                ra_err = pixsize * x_error / center_platescale
                dec_err = -pixsize * y_error / center_platescale
            elif ci == 'CIS':
                ra_err = -pixsize * x_error / az_platescale
                dec_err = pixsize * y_error / radial_platescale
            elif ci == 'CIN':
                ra_err = pixsize * x_error / az_platescale
                dec_err = -pixsize * y_error / radial_platescale
            elif ci == 'CIE':
                dec_err = -pixsize * x_error / az_platescale
                ra_err = -pixsize * y_error / radial_platescale
            elif ci == 'CIW':
                dec_err = pixsize * x_error / az_platescale
                ra_err = pixsize * y_error / radial_platescale

            else:
                raise ValueError(ci)
            
            rows.append(dict(
                frame=i, ci=ci, roi=roi,
                x_error=x_error, y_error=y_error,
                ra_err=ra_err, dec_err=dec_err)
                )

    data = Table(rows,
             names=['frame','ci','roi','x_error','y_error','ra_err','dec_err'],
             dtype=(int, str, int, float, float, float, float))

    data.meta['expid'] = cx['expid']
    data.meta['status'] = cx['status']
    data.meta['started_at'] = cx['started_at']
    data.meta['ended_at'] = cx['ended_at']
    data.meta['summary'] = cx['summary']
    data.meta['filename'] = cx['filename']
    if skyra is not None:
        data.meta['skyra'] = skyra
        data.meta['skydec'] = skydec

    #- Derive NIGHT from filename since it isn't in JSON
    night = int(os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(cx['filename'])))))
    data.meta['night'] = night

    return data

def plot_offsets(data, outfile):
    figure(figsize=(12,12))

    ci_names = ['CIE', 'CIN', 'CIC', 'CIS', 'CIW']

    subplot(2,2,1)
    for ci in ci_names:
        ii = (data['ci'] == ci)
        plot(data['ra_err'][ii], data['dec_err'][ii], '.', label=ci, alpha=0.5, ms=10)

    ii = data['ci'] == 'CIX'
    plot(data['ra_err'][ii], data['dec_err'][ii], 'kx', ms=10, label='guider corr', alpha=0.7)
# plot(combined_ra, combined_dec, 'kx', ms=10)

    grid()
    axhline(0, color='k', alpha=0.5)
    axvline(0, color='k', alpha=0.5)
    legend()
    xlabel('RA*cos(dec) error [arcsec]\nMean={:.2f}, RMS={:.2f}'.format(
        np.mean(data['ra_err']), np.std(data['ra_err'])))
    ylabel('dec error [arcsec]\nMean={:.2f}, RMS={:.2f}'.format(
        np.mean(data['dec_err']), np.std(data['dec_err'])))

    if 'skyra' in data.meta:
        title('Night {} Exp {} at RA,dec=({:.2f},{:.2f})'.format(
            data.meta['night'], data.meta['expid'],
            data.meta['skyra'], data.meta['skydec']))
    else:
        title('Night {} Exp {}; Unknown RA,DEC'.format(
            data.meta['night'], data.meta['expid']))

    xlim(-1, 1); ylim(-1, 1)

    for i, ci in enumerate(ci_names):
        subplot(25,2,2*i+2)
        ii = (data['ci'] == ci)
        hist(data['ra_err'][ii], 20, (-1, 1), alpha=0.8,
             color=cm.tab10.colors[i])
        ylabel(ci)
    xlabel('ra_err')


    for i, ci in enumerate(ci_names):
        subplot(25,2,12+2*i+2)
        ii = (data['ci'] == ci)
        hist(data['dec_err'][ii], 20, (-1, 1), alpha=0.8,
             color=cm.tab10.colors[i])
        ylabel(ci)
    xlabel('dec_err')

    subplot(4,1,3)
    for i, ci in enumerate(ci_names):
        for roi in set(data['roi'][data['ci'] == ci]):
            ii = (data['ci'] == ci) & (data['roi'] == roi)
            plot(data['frame'][ii], data['ra_err'][ii], '-',
                 alpha=0.7, color=cm.tab10.colors[i])

    ii = (data['ci'] == 'CIX')
    plot(data['frame'][ii], data['ra_err'][ii], '-', color='k')
    ylabel('ra_err [arcsec]')
    ylim(-1.1, 1.1)
    grid()

    subplot(4,1,4)
    for i, ci in enumerate(ci_names):
        for roi in set(data['roi'][data['ci'] == ci]):
            ii = (data['ci'] == ci) & (data['roi'] == roi)
            plot(data['frame'][ii], data['dec_err'][ii], '-',
                 alpha=0.7, color=cm.tab10.colors[i])

    ii = (data['ci'] == 'CIX')
    plot(data['frame'][ii], data['dec_err'][ii], '-', color='k')
    ylabel('dec_err [arcsec]')
    xlabel('frame number')
    ylim(-1.1, 1.1)
    grid()

    savefig(outfile)
    print('Wrote {}'.format(outfile))

if __name__ == '__main__':
    sys.exit(main())
