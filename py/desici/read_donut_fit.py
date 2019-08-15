import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import re

# donutana wants a dictionary of donuts. it looks like it can take a
# list of headers from the donutfit output files

def readdonutfittodict(filename):
    # the header for HDU0 has all the info we need
    hdulist = fits.open(filename)
    hdr = hdulist[0].header
    outdict = {}
    outdict['IFILE'] = filename
    for key in ['ZERN4','ZERN5','ZERN6','ZERN7','ZERN8','ZERN9','ZERN10','ZERN11','ZERN12','ZERN13','ZERN14','ZERN15','rzero','IX','IY','CHI2','DOF','FITSTAT','NELE','XPOS','YPOS','MOUNTAZ','MOUNTEL','EXTNAME','EXPID','ZERN4E','ZERN5E','ZERN6E','ZERN7E','ZERN8E','ZERN9E','ZERN10E']:
        outdict[key] = hdr[key]
    return outdict


# get a set of donut fit output files from <dirname> and return a list
# of  headers that have the fit info. the header for hdu0 has everything needed
def getdonutfithdrs(dirname, fileprefix, namematch):

    import os
    import re

    flist = []
    outdictlist = []
    if (os.path.isdir(dirname)) != True:
            return outdictlist
    dirflist = os.listdir(dirname)
    for fname in dirflist:
        if (re.search(fileprefix, fname) != None):
            if (re.search(namematch, fname) != None):
                flist.append(fname)
    for fname in flist:
        hdul = fits.open(dirname+'/'+fname)
        outdictlist.append(hdul[0].header)
    return outdictlist

def getdonutfitresults(dirname, fileprefix, namematch):

    import os
    import re
    
    dirflist = os.listdir(dirname)
    flist = []
    fitdictlist = []
    for fname in dirflist:
        if (re.search(fileprefix, fname) != None):
            if (re.search(namematch, fname) != None):
                flist.append(fname)
    for fname in flist:
        onedict = readdonutfittodict(dirname+'/'+fname)
        fitdictlist.append(onedict)
    # now re-group info by zernike, not by chip
    ndict = len(fitdictlist)
    outdict = {}
    #outkeys is minus the EXTNAME
    outkeys = ['ZERN4','ZERN5','ZERN6','ZERN7','ZERN8','ZERN9','ZERN10','ZERN11','ZERN12','ZERN13','ZERN14','ZERN15','rzero','IX','IY','CHI2','DOF','FITSTAT','NELE','XPOS','YPOS','MOUNTAZ','MOUNTEL','EXPID','ZERN4E','ZERN5E','ZERN6E','ZERN7E','ZERN8E','ZERN9E','ZERN10E']
    # make output arrays
    for ikey in outkeys:
        outdict[ikey] = np.zeros(ndict)
    outdict['EXTNAME'] = []
    outdict['DONUTCUT'] = np.zeros(ndict,dtype=np.bool)
    # now stuff them
    for i in range(ndict):
        onedict = fitdictlist[i]
        for ikey in outkeys:
            outdict[ikey][i] = onedict[ikey]
        outdict['EXTNAME'].append(onedict['EXTNAME'])
        #cutval = np.sqrt(onedict['ZERN5']**2 + onedict['ZERN6']**2 + onedict['ZERN7']**2 + onedict['ZERN8']**2)
        # there are also z4 cuts, could add those
        if ((onedict['ZERN4']<10) & (np.abs(onedict['ZERN5'])<1.) & (np.abs(onedict['ZERN6'])<2.) & (np.abs(onedict['ZERN7'])<1) & (np.abs(onedict['ZERN8'])<1) & (np.abs(onedict['ZERN9'])<1.) & (np.abs(onedict['ZERN10'])<1.)):
            cutval = 1
        else:
            cutval = 0
        outdict['DONUTCUT'][i] = ((cutval > 0.5) & (onedict['NELE'] > 100000.))
        #outdict['DONUTCUT'][i] = ((cutval < 3.0) & (onedict['NELE'] > 100000.))
    return outdict

# from a list of dictionaries that come from fitting donuts to a frame
# make a scatter plot of the desired quantity
def scatplotdonutfits(fitdict, keyname, titlestr=None, oroot=None, dointerp=False,idx=None):
    
    #darr = np.array([idata[keyname] for idata in fitdict])
    #xpos = np.asarray([idata['xpos'] for idata in fitdict])
    #ypos = np.array([idata['ypos'] for idata in fitdict])
    darr = fitdict[keyname]
    xpos = fitdict['XPOS']
    ypos = fitdict['YPOS']

    if idx == None:
        idx = np.arange(len(darr))

    darr = darr[idx]
    xpos = xpos[idx]
    ypos = ypos[idx]
    
    if titlestr == None:
        fulltitlestr = keyname
    else:
        fulltitlestr = titlestr + '_' + keyname
    
    icin, = np.where(ypos > 1.)
    icis, = np.where(ypos < -1.)
    icie, = np.where(xpos < -1.)
    iciw, = np.where(xpos > 1.)
    icic, = np.where((ypos > -1) & (ypos < 1.) & (xpos > -1) & (xpos < 1))

    minval = np.min(darr)
    maxval = np.max(darr)
    bfig,((bax1,bax2,bax3),(bax4,bax5,bax6),(bax7,bax8,bax9)) = plt.subplots(3,3)    
    #axlist is in order with the order in the default for plotchips
    baxlist = [bax2,bax4,bax5,bax6,bax8]
    # no axes drawn in grid positions we will not use
    bax1.axis('off')
    bax3.axis('off')
    bax7.axis('off')
    bax9.axis('off')
    plotchips = ['CIN','CIE','CIC','CIW','CIS']
    idxlist = [icin,icie,icic,iciw,icis]
    for i in range(len(plotchips)):
        if (len(idxlist[i]) < 1):
            continue
        biax = baxlist[i]
        ciname = plotchips[i]
        if i > 0:
            tname = ciname
            biax.set_title(tname)
        xmin = np.min(xpos[idxlist[i]])
        xmax = np.max(xpos[idxlist[i]])
        ymin = np.min(ypos[idxlist[i]])
        ymax = np.max(ypos[idxlist[i]])
        nx = len(xpos[idxlist[i]])
        ny = len(ypos[idxlist[i]])
        dx = (xmax-xmin)/(nx-1)
        dy = (ymax-ymin)/(ny-1)
        xi = np.arange(nx)*dx+xmin
        yi = np.arange(ny)*dy+ymin

        bxplt,byplt = np.meshgrid(xi,yi)
        bplt = griddata((xpos[idxlist[i]],ypos[idxlist[i]]),darr[idxlist[i]],(bxplt,byplt),method='nearest')
        if (dointerp == True):
            bfoo = biax.imshow(bplt,vmin=minval,vmax=maxval,extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto')
        bfoo = biax.scatter(xpos[idxlist[i]],ypos[idxlist[i]],c=darr[idxlist[i]],vmin=minval,vmax=maxval,s=20)
        if i == 0:
            bcfoo = bfig.colorbar(bfoo)
        #tight_layout is magic to fix things like axis label overlap
        btitlestr = fulltitlestr 
        plt.suptitle(btitlestr,fontsize=8)
        plt.tight_layout()
        if oroot != None:
            beforeoutname = oroot+'.png'
            plt.savefig(beforeoutname)

# read donut fit results to get ndonuts per CI camera used in the fits
# implement cut from donutana, see getdonutfitresults and check cut!
# NOTE have to do this here, DONUTCUT is per donut, this condenses that info
def ndonutscamera(dirname, fileprefix, namematch,docut=True):

    rematch = re.search('[0-9]+',fileprefix)
    expnum = rematch.group(0)
    expdirname = dirname + '/' + expnum
    fitdict = getdonutfitresults(expdirname, fileprefix, namematch)
    xpos = fitdict['XPOS']
    ypos = fitdict['YPOS']
    cutval = fitdict['DONUTCUT']

    if docut is True:
        icin, = np.where((ypos > 1.) & (cutval == True))
        icis, = np.where((ypos < -1.) & (cutval == True))
        icie, = np.where((xpos < -1.) & (cutval == True))
        iciw, = np.where((xpos > 1.) & (cutval == True))
        icic, = np.where((ypos > -1) & (ypos < 1.) & (xpos > -1) & (xpos < 1) & (cutval == True))
    else:
        print('dont know what to do!')
        
    outdict = {}
    outdict['ndonuts_cin'] = len(icin)
    outdict['ndonuts_cis'] = len(icis)
    outdict['ndonuts_cie'] = len(icie)
    outdict['ndonuts_ciw'] = len(iciw)
    outdict['ndonuts_cic'] = len(icic)
    
    return(outdict)
    
