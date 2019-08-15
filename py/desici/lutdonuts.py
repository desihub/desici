import numpy as np
import re
import glob
import matplotlib.pyplot as plt
#from scipy.interpolate import griddata
#from donutlib.donutana import donutana
from read_donut_fit import *
from getdonutanadata import *
import table as lut
import matplotlib.pyplot as plt
import scipy.optimize as optfit

# code for measuring the hexapod sensitivity matrix and for computing
# lookup tables

# input is a dictionary of donuntana results
# apply the dodx,dody,doxt,doyt results from donutana to the hexapod
# positions hexx,hexy,hexxt,hexyt
# this function also gets the LUT values to compare with donutana results
# LUT usage: lutobj.interpolate_delaunay(<az>,<alt>)
def applydohexapod(dandict, lutfile = '/Users/crockosi/desi/desitrac/HexapodLUT/trunk/TelescopeLUT_20190521.txt'):
    hexlut = lut.LUT_table(lutfile)
    # coordinates of hexapod postion when data were taken
    hexx = dandict['hexx']
    hexy = dandict['hexy']
    hexz = dandict['hexz']
    hexxt = dandict['hexxt']
    hexyt = dandict['hexyt']
    # telescope position when data were taken
    alt = dandict['MOUNTEL']
    az = dandict['MOUNTAZ']
    # corrections as computed by donutana
    dodx = dandict['dodx']
    dody = dandict['dody']
    dodz = dandict['dodz']
    doxt = dandict['doxt']
    doyt = dandict['doyt']

    # should add the corrections to the current coordinates if the signs are
    # correct. 
    newx = hexx+dodx
    newy = hexy+dody
    newz = hexz+dodz
    newxt = hexxt+doxt
    newyt = hexyt+doyt

    # lut table settings for this telecope position
    npoint = len(alt)
    lutx = np.zeros(npoint)
    luty = np.zeros(npoint)
    lutz = np.zeros(npoint)
    lutxt = np.zeros(npoint)
    lutyt = np.zeros(npoint)
    for i in range(npoint):
        lutvec = hexlut.interpolate_delaunay(az[i],alt[i])
        if lutvec[6] is not 'SUCCESS':
            print('lut interpolation failed')
            print(foo)
        lutx[i] = lutvec[0]
        luty[i] = lutvec[1]
        lutz[i] = lutvec[2]
        lutxt[i] = lutvec[3]
        lutyt[i] = lutvec[4]

    #now difference between lut and corrected alignment
    diffx = newx-lutx
    diffy = newy-luty
    diffz = newz-lutz
    diffxt = newxt-lutxt
    diffyt = newyt-lutyt

    # package output into dandict
    # new, corrected hexapod position
    dandict['newx'] = newx
    dandict['newy'] = newy
    dandict['newz'] = newz
    dandict['newxt'] = newxt
    dandict['newyt'] = newyt
    # lut table settings for this telecope position
    dandict['lutx'] = lutx
    dandict['luty'] = luty
    dandict['lutz'] = lutz
    dandict['lutxt'] = lutxt
    dandict['lutyt'] = lutyt
    # difference between lut and corrected alignment
    dandict['diffx'] = diffx
    dandict['diffy'] = diffy
    dandict['diffz'] = diffz
    dandict['diffxt'] = diffxt
    dandict['diffyt'] = diffyt
    
    return 

########################################################
# code for measuring the hexapod sensitivty matrix
#########################################################

# fit a slope. yvec should be a hexapod aligment value as solved by donutana
# error is the error on that, who knows.
def fitslope(xvec,yvec,yerrvec):

    def linfunc(depvals, slope, yint):
        out = depvals*slope + yint
        return out

    pvals,pcov = optfit.curve_fit(linfunc,xvec,yvec,sigma=yerrvec)
    fityvals = xvec*pvals[0]+pvals[1]
    plt.plot(xvec,fityvals,'b')

    return pvals

# nominal is a 6-element vector of nomial values used to separate out
# the sweeps along each axix in the alignment sweep data.
# dict is a dictionary of donutana results from a multi-axis alignment sweep
# j is 0,1,2,3,4 for x,y,z,xt,yt
def checkhexmatrix(sdict,nomvec,j):

    nomx = nomvec[0]
    nomy = nomvec[1]
    nomz = nomvec[2]
    nomxt = nomvec[3]
    nomyt = nomvec[4]
    # idx is all the data where the hexapod was not  at its nominal x pos'n
    # same for idy, etc.
    idx,=np.where((sdict['hexx'] > nomx+1) | (sdict['hexx'] < nomx-1))
    idy,=np.where((sdict['hexy'] > nomy+1) | (sdict['hexy'] < nomy-1))
    idz,=np.where((sdict['hexz'] > nomz+1) | (sdict['hexz'] < nomz-1))
    idtx,=np.where((sdict['hexxt'] > nomxt+1) | (sdict['hexxt'] < nomxt-1))
    idty,=np.where((sdict['hexyt'] > nomyt+1) | (sdict['hexyt'] < nomyt-1))
    ilist = [idx,idy,idz,idtx,idty]
    paramlist = [sdict['hexx'][idx]-nomx,sdict['hexy'][idy]-nomy,sdict['hexz'][idz]-nomz, sdict['hexxt'][idtx]-nomxt,sdict['hexyt'][idty]-nomyt]
    xlabellist = ['delta x, microns','deltay y, microns', 'delta z, microns','delta xtilt, arcsec','delta ytilt, arcsec']
    titlestrlist=['delta x','delta y','delta z','delta xtilt','delta ytilt']
    
    # the rest of this proc modeled on collectzernslopes
    ivec = ilist[j]
    paramvals = paramlist[j]
    xlabel = xlabellist[j]
    titlestr = titlestrlist[j]
    
#########
# Yes!  This seems to be the mapping of donutana thetax,thetay slopes and
# ns/ew. donutana's thetax really is the CIN-CIS slope, thetay is the
# CIE-CIW slope. 
    nsslopez5 = sdict['z5thetax'][ivec]
    nsslopez5err = sdict['z5thetaxErr'][ivec]
    nsslopez6 = sdict['z6thetax'][ivec]
    nsslopez6err = sdict['z6thetaxErr'][ivec]
    ewslopez5 = sdict['z5thetay'][ivec]
    ewslopez5err = sdict['z5thetayErr'][ivec]
    ewslopez6 = sdict['z6thetay'][ivec]
    ewslopez6err = sdict['z6thetayErr'][ivec]
    deltaz7 = sdict['z7meanDeltaBefore'][ivec]
    deltaz7err = sdict['z7deltaErr'][ivec]
    deltaz8 = sdict['z8meanDeltaBefore'][ivec]
    deltaz8err = sdict['z8deltaErr'][ivec]
    # nsslope and ewslope defined in Roodman 2014
    # would divide by 252 to get same units as Aaron's Fig 7
    # NS slope is slope in Y in field angle
    nsslope = 0.5*(ewslopez5-nsslopez6)
    nsslopeerr = np.sqrt((ewslopez5err**2+nsslopez6err**2)/2.)
    # EW slope is slope in X in field angle
    ewslope = 0.5*(ewslopez6+nsslopez5)
    ewslopeerr = np.sqrt((ewslopez6err**2+nsslopez5err**2)/2.)
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    ax1.set_ylim([-0.5,0.5])
    ax1.scatter(paramvals,nsslope)
    ax1.errorbar(paramvals,nsslope,yerr=nsslopeerr)
    ax1.set_title('NS Y slope, waves per deg field angle',fontsize=10)
    ax1.set_ylabel('waves/degree')
    ax1.set_xlabel(xlabel)
    ax2.set_ylim([-0.5,0.5])
    ax2.scatter(paramvals,ewslope)
    ax2.errorbar(paramvals,ewslope,yerr=ewslopeerr)
    ax2.set_title('EW X slope, waves per deg field angle',fontsize=10)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel('waves/degree')
    ax3.set_ylim([-0.2,0.2])
    ax3.scatter(paramvals,deltaz7)
    ax3.errorbar(paramvals,deltaz7,yerr=deltaz7err)
    ax3.set_title('Delta z7',fontsize=10)
    ax3.set_ylabel('waves of z7')
    ax3.set_xlabel(xlabel)
    ax4.set_ylim([-0.2,0.2])
    ax4.scatter(paramvals,deltaz8)
    ax4.errorbar(paramvals,deltaz8,yerr=deltaz8err)
    ax4.set_title('Delta z8',fontsize=10)
    ax4.set_ylabel('waves of z8')
    ax4.set_xlabel(xlabel)
    plt.suptitle(titlestr,fontsize=10,y=0.99)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    # output
    outdict={}
    outdict['paramvals'] = np.asarray(paramvals)
    outdict['yslope'] = nsslope
    outdict['yslopeerr'] = nsslopeerr
    outdict['xslope'] = ewslope
    outdict['xslopeerr'] = ewslopeerr
    outdict['deltaz7'] = deltaz7
    outdict['deltaz8'] = deltaz8
    outdict['deltaz7err'] = deltaz7err
    outdict['deltaz8err'] = deltaz8err
    return outdict

# xtiltvecs,ytiltvecs,xdecentervecs,ydecentervecs are the inputs to makehexapodmatrix.  hexmatfits is the fitdict output of makehexapodmatrix
def plotcheckhexmat(xtiltvecs,ytiltvecs,xdecentervecs,ydecentervecs,hexmatfits):

    xlabellist = ['delta x, microns','deltay y, microns', 'delta z, microns','delta xtilt, arcsec','delta ytilt, arcsec']
    titlestrlist=['delta x','delta y','delta z','delta xtilt','delta ytilt']

    # first x decenter
    xlabel = xlabellist[0]
    titlestr = titlestrlist[0]
    paramvals = xdecentervecs['paramvals']
    pmin = np.min(paramvals)
    pmax = np.max(paramvals)
    paramintvals = np.arange(100)*(pmax-pmin)/100.+pmin
    xslope = xdecentervecs['xslope']
    xslopeerr = xdecentervecs['xslopeerr']
    yslope = xdecentervecs['yslope']
    yslopeerr = xdecentervecs['yslopeerr']
    deltaz7 = xdecentervecs['deltaz7']
    deltaz7err = xdecentervecs['deltaz7err']
    deltaz8 = xdecentervecs['deltaz8']
    deltaz8err = xdecentervecs['deltaz8err']
    fitxslopeslope = hexmatfits['Apar'][0]
    fitxslopeint = hexmatfits['Apar'][1]
    fitxslopevals = fitxslopeslope*paramintvals+fitxslopeint
    fityslopeslope = hexmatfits['Epar'][0]
    fityslopeint = hexmatfits['Epar'][1]
    fityslopevals = fityslopeslope*paramintvals+fityslopeint
    fitdz7slope = hexmatfits['Jpar'][0]
    fitdz7int = hexmatfits['Jpar'][1]
    fitdz7slopevals = fitdz7slope*paramintvals+fitdz7int
    fitdz8slope = hexmatfits['Npar'][0]
    fitdz8int = hexmatfits['Npar'][1]
    fitdz8slopevals = fitdz8slope*paramintvals+fitdz8int
    oneaxishexplot(xlabel, titlestr, paramvals, paramintvals, xslope, xslopeerr, yslope, yslopeerr, deltaz7, deltaz7err, deltaz8, deltaz8err, fitxslopeslope, fitxslopeint, fitxslopevals, fityslopeslope, fityslopeint, fityslopevals,fitdz7slope, fitdz7int, fitdz7slopevals, fitdz8slope, fitdz8int, fitdz8slopevals)

    # ydecenter
    xlabel = xlabellist[1]
    titlestr = titlestrlist[1]
    paramvals = ydecentervecs['paramvals']
    pmin = np.min(paramvals)
    pmax = np.max(paramvals)
    paramintvals = np.arange(100)*(pmax-pmin)/100.+pmin
    xslope = ydecentervecs['xslope']
    xslopeerr = ydecentervecs['xslopeerr']
    yslope = ydecentervecs['yslope']
    yslopeerr = ydecentervecs['yslopeerr']
    deltaz7 = ydecentervecs['deltaz7']
    deltaz7err = ydecentervecs['deltaz7err']
    deltaz8 = ydecentervecs['deltaz8']
    deltaz8err = ydecentervecs['deltaz8err']
    fitxslopeslope = hexmatfits['Bpar'][0]
    fitxslopeint = hexmatfits['Bpar'][1]
    fitxslopevals = fitxslopeslope*paramintvals+fitxslopeint
    fityslopeslope = hexmatfits['Fpar'][0]
    fityslopeint = hexmatfits['Fpar'][1]
    fityslopevals = fityslopeslope*paramintvals+fityslopeint
    fitdz7slope = hexmatfits['Kpar'][0]
    fitdz7int = hexmatfits['Kpar'][1]
    fitdz7slopevals = fitdz7slope*paramintvals+fitdz7int
    fitdz8slope = hexmatfits['Ppar'][0]
    fitdz8int = hexmatfits['Ppar'][1]
    fitdz8slopevals = fitdz8slope*paramintvals+fitdz8int
    oneaxishexplot(xlabel, titlestr, paramvals, paramintvals, xslope, xslopeerr, yslope, yslopeerr, deltaz7, deltaz7err, deltaz8, deltaz8err, fitxslopeslope, fitxslopeint, fitxslopevals, fityslopeslope, fityslopeint, fityslopevals,fitdz7slope, fitdz7int, fitdz7slopevals, fitdz8slope, fitdz8int, fitdz8slopevals)

    # xtilt
    xlabel = xlabellist[3]
    titlestr = titlestrlist[3]
    paramvals = xtiltvecs['paramvals']
    pmin = np.min(paramvals)
    pmax = np.max(paramvals)
    paramintvals = np.arange(100)*(pmax-pmin)/100.+pmin
    xslope = xtiltvecs['xslope']
    xslopeerr = xtiltvecs['xslopeerr']
    yslope = xtiltvecs['yslope']
    yslopeerr = xtiltvecs['yslopeerr']
    deltaz7 = xtiltvecs['deltaz7']
    deltaz7err = xtiltvecs['deltaz7err']
    deltaz8 = xtiltvecs['deltaz8']
    deltaz8err = xtiltvecs['deltaz8err']
    fitxslopeslope = hexmatfits['Cpar'][0]
    fitxslopeint = hexmatfits['Cpar'][1]
    fitxslopevals = fitxslopeslope*paramintvals+fitxslopeint
    fityslopeslope = hexmatfits['Gpar'][0]
    fityslopeint = hexmatfits['Gpar'][1]
    fityslopevals = fityslopeslope*paramintvals+fityslopeint
    fitdz7slope = hexmatfits['Lpar'][0]
    fitdz7int = hexmatfits['Lpar'][1]
    fitdz7slopevals = fitdz7slope*paramintvals+fitdz7int
    fitdz8slope = hexmatfits['Qpar'][0]
    fitdz8int = hexmatfits['Qpar'][1]
    fitdz8slopevals = fitdz8slope*paramintvals+fitdz8int
    oneaxishexplot(xlabel, titlestr, paramvals, paramintvals, xslope, xslopeerr, yslope, yslopeerr, deltaz7, deltaz7err, deltaz8, deltaz8err, fitxslopeslope, fitxslopeint, fitxslopevals, fityslopeslope, fityslopeint, fityslopevals,fitdz7slope, fitdz7int, fitdz7slopevals, fitdz8slope, fitdz8int, fitdz8slopevals)

    # ytilt
    xlabel = xlabellist[4]
    titlestr = titlestrlist[4]
    paramvals = ytiltvecs['paramvals']
    pmin = np.min(paramvals)
    pmax = np.max(paramvals)
    paramintvals = np.arange(100)*(pmax-pmin)/100.+pmin
    xslope = ytiltvecs['xslope']
    xslopeerr = ytiltvecs['xslopeerr']
    yslope = ytiltvecs['yslope']
    yslopeerr = ytiltvecs['yslopeerr']
    deltaz7 = ytiltvecs['deltaz7']
    deltaz7err = ytiltvecs['deltaz7err']
    deltaz8 = ytiltvecs['deltaz8']
    deltaz8err = ytiltvecs['deltaz8err']
    fitxslopeslope = hexmatfits['Dpar'][0]
    fitxslopeint = hexmatfits['Dpar'][1]
    fitxslopevals = fitxslopeslope*paramintvals+fitxslopeint
    fityslopeslope = hexmatfits['Hpar'][0]
    fityslopeint = hexmatfits['Hpar'][1]
    fityslopevals = fityslopeslope*paramintvals+fityslopeint
    fitdz7slope = hexmatfits['Mpar'][0]
    fitdz7int = hexmatfits['Mpar'][1]
    fitdz7slopevals = fitdz7slope*paramintvals+fitdz7int
    fitdz8slope = hexmatfits['Rpar'][0]
    fitdz8int = hexmatfits['Rpar'][1]
    fitdz8slopevals = fitdz8slope*paramintvals+fitdz8int
    oneaxishexplot(xlabel, titlestr, paramvals, paramintvals, xslope, xslopeerr, yslope, yslopeerr, deltaz7, deltaz7err, deltaz8, deltaz8err, fitxslopeslope, fitxslopeint, fitxslopevals, fityslopeslope, fityslopeint, fityslopevals,fitdz7slope, fitdz7int, fitdz7slopevals, fitdz8slope, fitdz8int, fitdz8slopevals)

# helper fuction for plotcheckhexmat above
def oneaxishexplot(xlabel, titlestr, paramvals, paramintvals, xslope, xslopeerr, yslope, yslopeerr, deltaz7, deltaz7err, deltaz8, deltaz8err, fitxslopeslope, fitxslopeint, fitxslopevals, fityslopeslope, fityslopeint, fityslopevals,fitdz7slope, fitdz7int, fitdz7slopevals, fitdz8slope, fitdz8int, fitdz8slopevals):
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    ax1.set_ylim([-0.5,0.5])
    ax1.scatter(paramvals,xslope)
    ax1.errorbar(paramvals,xslope,yerr=xslopeerr,fmt='none')
    ax1.plot(paramintvals,fitxslopevals,'r-')
    ax1.set_title('X slope, waves per deg field angle',fontsize=10)
    ax1.set_ylabel('waves/degree')
    ax1.set_xlabel(xlabel)
    ax2.set_ylim([-0.5,0.5])
    ax2.scatter(paramvals,yslope)
    ax2.errorbar(paramvals,yslope,yerr=yslopeerr,fmt='none')
    ax2.plot(paramintvals,fityslopevals,'r-')
    ax2.set_xlabel(xlabel)
    ax2.set_title('Y slope, waves per deg field angle',fontsize=10)
    ax2.set_ylabel('waves/degree')
    ax3.set_ylim([-0.2,0.2])
    ax3.scatter(paramvals,deltaz7)
    ax3.errorbar(paramvals,deltaz7,yerr=deltaz7err,fmt='none')
    ax3.plot(paramintvals,fitdz7slopevals,'r-')
    ax3.set_title('Delta z7',fontsize=10)
    ax3.set_ylabel('waves of z7')
    ax3.set_xlabel(xlabel)
    ax4.set_ylim([-0.2,0.2])
    ax4.scatter(paramvals,deltaz8)
    ax4.errorbar(paramvals,deltaz8,yerr=deltaz8err,fmt='none')
    ax4.plot(paramintvals,fitdz8slopevals,'r-')
    ax4.set_title('Delta z8',fontsize=10)
    ax4.set_ylabel('waves of z8')
    ax4.set_xlabel(xlabel)
    plt.suptitle(titlestr,fontsize=10,y=0.99)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)


##############
# code for making lookup tables
##############

def getdist(az1,el1,az2,el2):
        """
        Angular distance beween two points on sphere.
        Units for both input and output are degrees.
        AJR using formulae used for angular clustering
        """
        
        sa1=np.sin(az1*np.pi/180.0)
        sa2=np.sin(az2*np.pi/180.0)
        ca1=np.cos(az1*np.pi/180.0)
        ca2=np.cos(az2*np.pi/180.0)
        se1=np.sin(el1*np.pi/180.0)
        se2=np.sin(el2*np.pi/180.0)
        ce1=np.cos(el1*np.pi/180.0)
        ce2=np.cos(el2*np.pi/180.0)
        ac = ce1*ce2*(ca1*ca2 + sa1*sa2) + se1*se2
        dist=np.arccos(ac)*180.0/np.pi
        return dist

#avradius in degrees
# doesnt' deal with z for now, the temp differences mean it is hard to
# bundle all the data
def avpoints(gridaz,gridel,az,el,hexx,hexy,hexz,hexxt,hexyt,avradius,doplots=True):

    npoints = len(gridaz)
    ndata = len(hexx)

    outhexx = []
    outhexy = []
    outhexz = []
    outhexxt = []
    outhexyt = []

    goodgridaz = []
    goodgridel = []
    
    for i in range(npoints):
        igaz = gridaz[i]
        igel = gridel[i]

        # goofy special cases to use all the data
        # get data at el 60 az 180
        #if ((igaz > 150.) & (igaz < 190) & (igel > 55) & (igel < 65)):
        #    print(igaz)
        #    print(igel)
        #    igaz = 180.
        # data at 225,32
        #if ((igaz > 200.) & (igaz < 230) & (igel > 25) & (igel < 35)):
        #    print(igaz)
        #    print(igel)
        #    igaz = 225.
        
        dist = getdist(igaz,igel,az,el)
        itake, = np.where(dist < avradius)

        if len(itake) < 2:
            # print('less than two points: '+'{:.1f}'.format(igaz) + ' El ' + '{:.1f}'.format(igel))
            continue

        goodgridaz.append(igaz)
        goodgridel.append(igel)

        xav = np.median(hexx[itake])
        outhexx.append(xav)
        yav = np.median(hexy[itake])
        outhexy.append(yav)
        zav = np.median(hexz[itake])
        outhexz.append(zav)
        xtav = np.median(hexxt[itake])
        outhexxt.append(xtav)
        ytav = np.median(hexyt[itake])
        outhexyt.append(ytav)

        if doplots == True:
    #        fig,((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3)
            fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    #        ax6.axis('off')
            titlestr = 'Az ' + '{:.1f}'.format(igaz) + ' El ' + '{:.1f}'.format(igel)
            plt.suptitle(titlestr)
            ax1.hist(hexx[itake],10)
            ax1.axvline(x=xav)
            ax1.set_xlabel('hexapod x, microns')
            ax2.hist(hexy[itake],10)
            ax2.axvline(x=yav)
            ax2.set_xlabel('hexapod y, microns')
    #        ax3.hist(hexz[itake],10)
    #        ax3.axvline(x=zav)
    #        ax3.set_xlabel('hexapod z, microns')
            ax3.hist(hexxt[itake],10)
            ax3.axvline(xtav)
            ax3.set_xlabel('hexapod x tilt, arcsec')
            ax4.hist(hexyt[itake],10)
            ax4.axvline(ytav)
            ax4.set_xlabel('hexapod y tilt, arcsec')
            plt.tight_layout()
            
    goodgridaz= np.asarray(goodgridaz)
    goodgridel = np.asarray(goodgridel)
    outhexx=np.asarray(outhexx)
    outhexy=np.asarray(outhexy)
    outhexz=np.asarray(outhexz)
    outhexxt=np.asarray(outhexxt)
    outhexyt=np.asarray(outhexyt)
    
    outdict = {}
    outdict['allaz'] = gridaz
    outdict['allel'] = gridel
    outdict['hexx'] = outhexx
    outdict['hexy'] = outhexy
    outdict['hexz'] = outhexz
    outdict['hexxt'] = outhexxt
    outdict['hexyt'] = outhexyt
    outdict['goodaz'] = goodgridaz
    outdict['goodel'] = goodgridel

    return outdict


def oldmakespheregrid(delta):

    # the 18 and 32 degree values are so I can allow a 10 deg match
    # radius to LUT points at 20 and 30 degrees and get the best close
    # match. Probably hand edit the LUT back to 20 and 30 deg
    elvals = np.asarray([18., 32., 45., 60., 75.])
    
    nel = len(elvals)
    elvec = np.zeros(0)
    azvec = np.zeros(0)
    print(elvec)
    print(azvec)
    for i in range(nel):
        eldeg = elvals[i]
        elrad = eldeg*np.pi/180.
        idaz = delta/np.cos(elrad)
        naz = np.int(np.trunc(np.cos(elrad)*360./delta))
        iazvec = np.arange(naz)*idaz
        ielvec = np.ones(naz)*eldeg
        azvec = np.append(azvec,iazvec)
        elvec = np.append(elvec,ielvec)
    # add point at 90 degrees by hand
    azvec = np.append(azvec,0.0)
    elvec = np.append(elvec,90.)
    return azvec,elvec

# filename is ocs-format lut file
# zvec is [x,y,z,xt,yt] at zenith
# subtract zenithval - LUT_absval to get deltas
# zenithval-delta = LUT vals
def readlut(filename,zenithvec):

    data = np.genfromtxt(filename)
    nlines = len(data)
    azvec = np.zeros(nlines)
    elvec = np.zeros(nlines)
    hexx = np.zeros(nlines)
    hexy = np.zeros(nlines)
    hexz = np.zeros(nlines)
    hexxt = np.zeros(nlines)
    hexyt = np.zeros(nlines)
    for i in range(nlines):
        ivec = data[i]
        azvec[i] = ivec[0]
        elvec[i] = ivec[1]
        hexx[i] = ivec[2]
        hexy[i] = ivec[3]
        hexz[i] = ivec[4]
        hexxt[i] = ivec[5]
        hexyt[i] = ivec[6]

    hexdx = zenithvec[0]-hexx
    hexdy = zenithvec[1]-hexy
    hexdz = zenithvec[2]-hexz
    hexdxt = zenithvec[3]-hexxt
    hexdyt = zenithvec[4]-hexyt
    
    outdict = {}
    outdict['az'] = azvec
    outdict['el'] = elvec
    outdict['hexdx'] = hexdx
    outdict['hexdy'] = hexdy
    outdict['hexdz'] = hexdz
    outdict['hexdxt'] = hexdxt
    outdict['hexdyt'] = hexdyt
    return outdict

def oldwritelut(lutdict, filename):

    az = lutdict['goodaz']
    el = lutdict['goodel']
    hexx = lutdict['hexx']
    hexy = lutdict['hexy']
    hexz = lutdict['hexz']
    hexxt = lutdict['hexxt']
    hexyt = lutdict['hexyt']
    nel = len(az)
    
    with open(filename, 'w') as fout:
        for i in range(nel):
            outstr = '{:.2f}'.format(az[i]) + ' {:.2f}'.format(el[i]) + ' {:.2f}'.format(hexx[i]) + ' {:.2f}'.format(hexy[i])  + ' {:.2f}'.format(hexz[i]) + ' {:.2f}'.format(hexxt[i]) + ' {:.2f}'.format(hexyt[i]) +'\n'
            fout.write(outstr)
            

# get correction dz to focus from temperature LUT
def gettempcorrfocus(donutdict):

    tempvec = donutdict['TRUSTEMP']
    dztvec = 110.*(7.-tempvec)
    return dztvec

# lookup table dz values come from measuring slope of difference between
# best focus at zenith and best focus vs. elevation, taking out the different
# absolute focus values at zenith due to temperature. 
# the slope is defined so dz = dfocus/delevation 
# delevation is measured 90-el degrees and delta-focus is measured
# foc_zenith - foc_el. So if the two data points used to compute the
# slope are focus = -83.2 at elevation 20 and focus = +15 at elevation 90,
# then slope = (15 - -83.2)/(90-20) = 1.4. Focus gets more negative at
# lower elevations. dz gets larger toward lower elevations, and subtract
# dz to get the hexapod value for that elevation.
# newlut  is a lookup-table that needs dz values. it is a dictionary
# and has an element that gives the elevation value at each point in the
# LUT, currently key "lutel" and stuffs into key "lutdz" the dz that at
# gives the difference between the hexapod value at zenith and each
#  elevation in goodel 
def makedzlut(dzslope,newlut):

    elvals = newlut['lutel']
    deltaelvals = 90.-elvals
    dzvals = dzslope*deltaelvals
    newlut['lutdz'] = dzvals

# fit amplitude of sine function
# xvec is azimuth, radians
def fitsineamp(xvec,yvec,yerrvec,doplot=False):

    def sinefunc(depvals, amplitude, phase):
        out = amplitude*np.sin(depvals+phase)
        return out

    pvals,pcov = optfit.curve_fit(sinefunc,xvec,yvec,sigma=yerrvec)

    if doplot == True:
        fityvals = pvals[0]*np.sin(xvec+pvals[1])
        plt.plot(xvec*180./np.pi,fityvals,'ro')
        
    return pvals

# lutdataaz,el are the vectors of az,el for each data point for the lut
# deltavec is delta relative to zenith for dx,dy,xtilt,ytilt (dz is computed
# using the slope in focus vs. elevation). error is the error on each
# delta measurement
def fitlutdelta(lutdataaz, lutdatael, deltavec,deltaerrvec, lutazvec,lutelvec):
    lutdataazrad = lutdataaz*np.pi/180.
    lutazvecrad = lutazvec*np.pi/180.
    # slices in elevation
    i20,=np.where(lutdatael < 25)
    i30,=np.where((lutdatael > 25) & (lutdatael < 35))
    i45,=np.where((lutdatael > 40) & (lutdatael < 50))
    i60,=np.where((lutdatael > 55) & (lutdatael < 65))
    i75,=np.where((lutdatael > 70) & (lutdatael < 80))
    # fit each slice
    lut20 = fitsineamp(lutdataazrad[i20],deltavec[i20],deltaerrvec[i20])
    lut30 = fitsineamp(lutdataazrad[i30],deltavec[i30],deltaerrvec[i30])
    lut45 = fitsineamp(lutdataazrad[i45],deltavec[i45],deltaerrvec[i45])
    lut60 = fitsineamp(lutdataazrad[i60],deltavec[i60],deltaerrvec[i60])
    lut75 = fitsineamp(lutdataazrad[i75],deltavec[i75],deltaerrvec[i75])
    # output
    # the best fit on the regular LUT az grid
    # first split in slies in elevation
    ilut20 = np.where(lutelvec < 22.)
    ilut30 = np.where((lutelvec > 25) & (lutelvec < 35))
    ilut45 = np.where((lutelvec > 40) & (lutelvec < 50))
    ilut60 = np.where((lutelvec > 55) & (lutelvec < 65))
    ilut75 = np.where((lutelvec > 70) & (lutelvec < 80))
    # make output array
    lutdeltavals =  np.zeros(len(lutazvec))
    # stuff
    lutdeltavals[ilut20] = lut20[0]*np.sin(lutazvecrad[ilut20]+lut20[1])
    lutdeltavals[ilut30] = lut30[0]*np.sin(lutazvecrad[ilut30]+lut30[1])
    lutdeltavals[ilut45] = lut45[0]*np.sin(lutazvecrad[ilut45]+lut45[1])
    lutdeltavals[ilut60] = lut60[0]*np.sin(lutazvecrad[ilut60]+lut60[1])
    lutdeltavals[ilut75] = lut75[0]*np.sin(lutazvecrad[ilut75]+lut75[1])
    
    return lutdeltavals

def makespheregrid():

    elvals = np.asarray([20., 30., 45., 60., 75.])
    azvals = np.asarray([0,45,90,135,180,225,270,315])
    naz = len(azvals)
    nel = len(elvals)
    elvec = np.zeros(0)
    azvec = np.zeros(0)
    for i in range(nel):
        eldeg = elvals[i]
        ielvec = np.ones(naz)*eldeg
        iazvec = np.zeros(naz)+azvals
        azvec = np.append(azvec,iazvec)
        elvec = np.append(elvec,ielvec)
    # add point at 90 degrees by hand
    azvec = np.append(azvec,0.0)
    elvec = np.append(elvec,90.)
    return azvec,elvec

# remember, these are deltas of the hexapod position from zenith
# subtract these from the zenith value to get the absolote value
# of the hexapod position
def fitlut(lutdataaz, lutdatael, deltax, deltay, deltaxt, deltayt, dxerr, dyerr, dxterr, dyterr):

    lutaz,lutel = makespheregrid()
    fitdx = fitlutdelta(lutdataaz,lutdatael,deltax,dxerr,lutaz,lutel)
    fitdy = fitlutdelta(lutdataaz,lutdatael,deltay,dyerr,lutaz,lutel)
    fitdxt = fitlutdelta(lutdataaz,lutdatael,deltaxt,dxterr,lutaz,lutel)
    fitdyt = fitlutdelta(lutdataaz,lutdatael,deltayt,dyterr,lutaz,lutel)
    
    lutdict = {}
    lutdict['lutaz'] = lutaz
    lutdict['lutel'] = lutel
    lutdict['lutdx'] = fitdx
    lutdict['lutdy'] = fitdy
    lutdict['lutdxtilt'] = fitdxt
    lutdict['lutdytilt'] = fitdyt
    return lutdict


def writelut(lutdict, filename):

    az = lutdict['lutaz']
    el = lutdict['lutel']
    hexx = lutdict['lutdx']
    hexy = lutdict['lutdy']
    hexz = lutdict['lutdz']
    hexxt = lutdict['lutdxtilt']
    hexyt = lutdict['lutdytilt']
    nel = len(az)
    
    with open(filename, 'w') as fout:
        for i in range(nel):
            outstr = '{:.2f}'.format(az[i]) + ' {:.2f}'.format(el[i]) + ' {:.2f}'.format(hexx[i]) + ' {:.2f}'.format(hexy[i])  + ' {:.2f}'.format(hexz[i]) + ' {:.2f}'.format(hexxt[i]) + ' {:.2f}'.format(hexyt[i]) +'\n'
            fout.write(outstr)
            
