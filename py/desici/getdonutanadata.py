# this doesn't write the mesh objects

import numpy as np
import re
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from donutlib.donutana import donutana
from read_donut_fit import *

# use fileroot as outfileprefix
# specific to desi data directory and file structure
# firstnum and lastnum are sequence numbers, as ints, no leading zeros
# filedir includes the date part of the path
# of sumseq = True, collect all the donuts in the sequence range and then
# run donutana
def rundonutana(rootdir,firstnum,lastnum,outfileprefix,fileprefix='ci-',paramDict=None,dosecond=False,sumseq=False):
    
# change log10(nele)>3.5 to log10(nele)>5

    if paramDict == None:
        paramDict  = {"deltaZSetPoint":0.0,                  
                      #"donutCutString":"nele>0. and numpy.log10(nele)>5 and numpy.sqrt(zern5*zern5+zern6*zern6+zern7*zern7+zern8*zern8)<3.0 and abs(zern4)>3 and abs(zern4)<15",
                      #"donutCutString":"nele>0. and numpy.log10(nele)>5 and zern4<10. and abs(zern5)<1 and abs(zern6)<2 and abs(zern7)<1 and abs(zern8)<1 and abs(zern9)<1 and abs(zern10)<1",
                      "donutCutString":"nele>0. and numpy.log10(nele)>5 and zern4<10. and abs(zern5)<1 and abs(zern6)<2 and abs(zern7)<1 and abs(zern8)<1",
                      "sensorSet":"CI",    # for DESI, CI or GFA
                      "unVignettedOnly":False, # not used DESI
                      "doTrefoil":False,
                      "doSpherical":False,
                      "doQuadrefoil":False,
                      "doRzero":False,
                      "histFlag":False,
                      "debugFlag":True,
                      #"z4PointsFile":"Donut/python/Donut/donutlib/z4Mesh.dat",
                      "z4PointsFile":"/Users/crockosi/desi/desitrac/Donut/python/Donut/donutlib/z4Mesh.dat",
                      "z5PointsFile":"/Users/crockosi/desi/desitrac/Donut/python/Donut/donutlib/z5Mesh.dat",
                      "z6PointsFile":"/Users/crockosi/desi/desitrac/Donut/python/Donut/donutlib/z6Mesh.dat",
                      "z7PointsFile":"/Users/crockosi/desi/desitrac/Donut/python/Donut/donutlib/z7Mesh.dat",
                      "z8PointsFile":"/Users/crockosi/desi/desitrac/Donut/python/Donut/donutlib/z8Mesh.dat",
                      "z9PointsFile":"/Users/crockosi/desi/desitrac/Donut/python/Donut/donutlib/z9Mesh.dat",
                      "z10PointsFile":"/Users/crockosi/desi/desitrac/Donut/python/Donut/donutlib/z10Mesh.dat",
                      "z11PointsFile":"/Users/crockosi/desi/desitrac/Donut/python/Donut/donutlib/z11Mesh.dat"} 

    nseq = lastnum-firstnum+1
    idx = list(range(nseq))
    seqlist = [i + firstnum for i in idx]
#    for seqnum in range(firstnum,lastnum+1):
    if sumseq == False:
        nloop = nseq
        nsum = 1
    else:
        nloop = 1
        nsum = nseq
        
    for i in range(nloop):
        for j in range(nsum):
            # if we are not summing, j is only zero and we step in i
            # if we are summing, i is only zero and we step in j
            idx = j+i 
            seqnum = seqlist[idx]
            filedir  = rootdir  + '/' + '{:08d}'.format(seqnum) + '/'
            fileroot = fileprefix+'{:08d}'.format(seqnum)
            outfileroot = outfileprefix + fileroot
            if dosecond != False:
                jdatasec = getdonutfithdrs(filedir, fileroot, 'second')
            jdatafirst = getdonutfithdrs(filedir, fileroot, 'first')
            if j == 0:
                datafirst  = jdatafirst
                if dosecond != False:
                    datasec = jdatasec
            else:
                datafirst.extend(jdatafirst)
                if dosecond != False:
                    datasec.extend(jdatasec)
            print(len(datafirst))

        # now run donutana
        if dosecond != False:
            dan = donutana(**paramDict)
            if len(datasec) > 0:
                try:
                    outsec = dan.analyzeDonuts(datasec)
                    writeresultsdict(outsec,outfileroot+'_second')
                except Exception as e:
                    print(str(e))
        dan = donutana(**paramDict)
        if len(datafirst) > 0:
            try:
                outfirst = dan.analyzeDonuts(datafirst)
                writeresultsdict(outfirst,outfileroot+'_first')
                print(outfileroot)
                # useful stuff direct from donut fits headers for each exposure, not 
                # outherwise output by donutana
                dkeys = list(datafirst[0].keys())
                if 'HEXPOS' in dkeys:
                    hinfofile = outfileroot + '_headerinfo.txt'
                    with open(hinfofile,'w') as hout:
                        outstr = 'HEXPOS ' + datafirst[0]['HEXPOS']+'\n'
                        hout.write(outstr)
                        if 'FOCUS' in dkeys:
                            outstr = 'FOCUS ' + datafirst[0]['FOCUS']+'\n'
                            hout.write(outstr)
                        if 'MOUNTAZ' in dkeys:
                            outstr = 'MOUNTAZ ' + "{:f}".format(datafirst[0]['MOUNTAZ'])+'\n'
                            hout.write(outstr)
                        if 'MOUNTEL' in dkeys:
                            outstr = 'MOUNTEL ' + "{:f}".format(datafirst[0]['MOUNTEL'])+'\n'
                            hout.write(outstr)
                        if 'TRUSTEMP' in dkeys:
                            outstr = 'TRUSTEMP ' +"{:f}".format(datafirst[0]['TRUSTEMP'])+'\n'
                            hout.write(outstr)
            except Exception as e:
                print(str(e))
                
def writeresultsdict(resultsdict, outfileroot, outfileprefix=None):

    mainoutfilename = outfileroot + '_donut_summary.txt'
    with open(mainoutfilename, 'w') as fout:
        donutsum = resultsdict['donut_summary']
        for key in donutsum.keys():
            outstr = key + ' ' + "{:f}".format(donutsum[key])+'\n'
            fout.write(outstr)
    # now loop through and write out the results for each zernike
    # write the numbers in one file, and the arrays one per file 
    for zdictname in ['z4ResultDict','z5ResultDict','z6ResultDict','z7ResultDict','z8ResultDict']:
        zdict = resultsdict[zdictname]
        for zdictkey in ['thetax','thetaxErr','thetay','thetayErr','delta','deltaErr','meanDeltaBefore','rmsDeltaBefore','meanDeltaAfter','rmsDeltaAfter']:
            zoutfilename = outfileroot + '_' + zdictname + '.txt'
            with open(zoutfilename, 'w') as fout:
                outstr = zdictkey + ' ' + "{:f}".format(zdict[zdictkey])+'\n'
                fout.write(outstr)
        for zarrayname in ['deltaArrayX','deltaArrayY','deltaArrayBefore','deltaArrayAfter']:
            zarray = zdict[zarrayname]
            arroutname = outfileroot + '_' + zdictname + '_' + zarrayname + '.txt'
            hdrstr = outfileroot + ' ' + zdictname + ' ' + zarrayname
            np.savetxt(arroutname,zarray,header=hdrstr)

   
def readdeltaarrayxy(fileroot):
    getxfilename = fileroot + '_z4ResultDict_deltaArrayX' + '.txt'
    xarry = np.genfromtxt(getxfilename)
    getyfilename = fileroot + '_z4ResultDict_deltaArrayY' + '.txt'
    yarry = np.genfromtxt(getyfilename)
    return xarry,yarry

# zcoeff is one of z4, z5, etc., as a string
def readzarr(fileroot, zcoeff,getxy=True):
    getfilename = fileroot + '_' + zcoeff + 'ResultDict_deltaArrayBefore' + '.txt'
    before = np.genfromtxt(getfilename)
    getfilename = fileroot + '_' + zcoeff + 'ResultDict_deltaArrayAfter' + '.txt'
    after = np.genfromtxt(getfilename)
    if getxy:
        getxfilename = fileroot + '_' + zcoeff + 'ResultDict_deltaArrayX' + '.txt'
        getyfilename = fileroot + '_' + zcoeff + 'ResultDict_deltaArrayY' + '.txt'
        xarr = np.genfromtxt(getxfilename)
        yarr = np.genfromtxt(getyfilename)
        return before,after,xarr,yarr
    else:
        return before,after

# plot the "before" and "after" arrays for one zernike coefficient
# zernike coeff is 'z5' etc., as string
# oroot is outfilenameroot
def plotdarrays(fileroot, zcoeff, titlestr=None, oroot = None):

    if titlestr is None:
        titlestr = fileroot

    before,after,xarr,yarr = readzarr(fileroot, zcoeff,getxy=True)

    fulltitlestr = titlestr + ' ' + zcoeff
    # eeeew, labelsize for the ticks is global
    plt.rc('xtick',labelsize=6)
    plt.rc('ytick',labelsize=6)
    
    icin, = np.where(yarr > 1.)
    icis, = np.where(yarr < -1.)
    icie, = np.where(xarr < -1.)
    iciw, = np.where(xarr > 1.)
    icic, = np.where((yarr > -1) & (yarr < 1.) & (xarr > -1) & (xarr < 1))
    print('CIC '+str(len(icic)))
    print('CIN '+str(len(icin)))
    print('CIS '+str(len(icis)))
    print('CIW '+str(len(iciw)))
    print('CIE '+str(len(icie)))
    
    bminval = np.min(before)
    bmaxval = np.max(before)
    aminval = np.min(after)
    amaxval = np.max(after)

    # grid of subplots that reprents CI camera placement in focal plane
    # before plots
    bfig,((bax1,bax2,bax3),(bax4,bax5,bax6),(bax7,bax8,bax9)) = plt.subplots(3,3)    
    #axlist is in order with the order in the default for plotchips
    baxlist = [bax2,bax4,bax5,bax6,bax8]
    # no axes drawn in grid positions we will not use
    bax1.axis('off')
    bax3.axis('off')
    bax7.axis('off')
    bax9.axis('off')
    # after plots
    afig,((aax1,aax2,aax3),(aax4,aax5,aax6),(aax7,aax8,aax9)) = plt.subplots(3,3)    
    # axlist is in order with the order in the default for plotchips
    aaxlist = [aax2,aax4,aax5,aax6,aax8]
    # no axes drawn in grid positions we will not use
    aax1.axis('off')
    aax3.axis('off')
    aax7.axis('off')
    aax9.axis('off')

    plotchips = ['CIN','CIE','CIC','CIW','CIS']
    idxlist = [icin,icie,icic,iciw,icis]
    for i in range(len(plotchips)):
        biax = baxlist[i]
        aiax = aaxlist[i]
        ciname = plotchips[i]
        if i > 0:
            tname = ciname
            #if (flipEWnames is True):
            #    if re.search('CIE',ciname) != None:
            #        tname = 'CIW'
            #    if re.search('CIW',ciname) != None:
            #        tname = 'CIE'
            aiax.set_title(tname)
            biax.set_title(tname)
        if len(idxlist[i]) > 1: 
            xmin = np.min(xarr[idxlist[i]])
            xmax = np.max(xarr[idxlist[i]])
            ymin = np.min(yarr[idxlist[i]])
            ymax = np.max(yarr[idxlist[i]])
            nx = len(xarr[idxlist[i]])
            ny = len(yarr[idxlist[i]])
            dx = (xmax-xmin)/(nx-1)
            dy = (ymax-ymin)/(ny-1)
            xi = np.arange(nx)*dx+xmin
            yi = np.arange(ny)*dy+ymin
            plt.figure(afig.number)       
            axplt,ayplt = np.meshgrid(xi,yi)
            aplt = griddata((xarr[idxlist[i]],yarr[idxlist[i]]),after[idxlist[i]],(axplt,ayplt),method='nearest')
            afoo = aiax.imshow(aplt, vmin=aminval,vmax=amaxval,extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto')
            afoo = aiax.scatter(xarr[idxlist[i]], yarr[idxlist[i]], c=after[idxlist[i]],vmin=aminval,vmax=amaxval,s=20)
            if i == 0:
                acfoo = afig.colorbar(afoo)
            plt.figure(bfig.number)             
            bxplt,byplt = np.meshgrid(xi,yi)
            bplt = griddata((xarr[idxlist[i]],yarr[idxlist[i]]),before[idxlist[i]],(bxplt,byplt),method='nearest')
            bfoo = biax.imshow(bplt, vmin=bminval,vmax=bmaxval,extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto')
            bfoo = biax.scatter(xarr[idxlist[i]],yarr[idxlist[i]],c=before[idxlist[i]],vmin=bminval,vmax=bmaxval,s=20)
            if i == 0:
                bcfoo = bfig.colorbar(bfoo)
    #tight_layout is magic to fix things like axis label overlap
    plt.figure(afig.number)
    atitlestr = fulltitlestr + ' ' + ' delta after subtracting fit'
    plt.suptitle(atitlestr,fontsize=8)
    plt.tight_layout()
    if oroot != None:
        afteroutname = oroot+'_'+zcoeff+'_DeltaAfter.png'
        plt.savefig(afteroutname)
    plt.figure(bfig.number)
    btitlestr = fulltitlestr + ' ' + 'delta before fit'
    plt.suptitle(btitlestr,fontsize=8)
    plt.tight_layout()
    if oroot != None:
        beforeoutname = oroot+'_'+zcoeff+'_DeltaBeforeFit.png'
        plt.savefig(beforeoutname)
    
# read a single donut_summary file
# filename is complete path
def readdonutsum(filename):
    sumdict = {}
    with open(filename,'r') as dfile:
        for iline in dfile:
            ilist = iline.split()
            sumdict[ilist[0]] = ilist[1]
    return(sumdict)


# read in many donut summary files as output by donutana, gather the
# output from each exposure info into arrays of the quantities
# if it is necessary to choose between first or second donut fits
# include that in matchstr
def gathermanydonutsum(fdir, matchstr = '*first_donut_sum*'):

    flist = glob.glob(fdir+'/'+matchstr)
    nfiles = len(flist)
    outdict = {}
    for i in range(nfiles):
        isumdata = readdonutsum(flist[i])
        # set up output dict of arrays
        if i == 0:
            keylist = list(isumdata.keys())
            for j in range(len(keylist)):
                outdict[keylist[j]] = np.zeros(nfiles)
            outdict['expnum'] = np.zeros(nfiles)
            outdict['expname'] = []
        # put data from this file into each array
        # get sequence num so we know what data goes with what exposure
        fname = flist[i]
        expmatch = re.search('ci-[0-9]+',fname)
        iexp = expmatch.group(0)
        outdict['expname'].append(iexp)
        inum = re.search('[0-9]+',iexp)
        inum = np.int(inum.group(0))
        outdict['expnum'][i] = inum
        for j in range(len(keylist)):
            jkey = keylist[j]
            jarr = outdict[jkey]
            jarr[i] = isumdata[jkey]
        iname = flist[i]
        if re.search('first',flist[i]) != None:
            hname = re.sub('first_donut_summary.txt','headerinfo.txt',iname)
        else:
            hname = re.sub('second_donut_summary.txt','headerinfo.txt',iname)
        #print(hname)
        with open(hname,'r') as hfile:
            for iline in hfile:
                ilist = iline.split()
                hkey = ilist[0]
                if hkey != 'FOCUS':
                    if hkey == 'HEXPOS':
                        hexvallist = ilist[1].split(',')
                        hexkeylist = ['hexx','hexy','hexz','hexxt','hexyt'] 
                        for k in range(len(hexkeylist)):
                            hexkey = hexkeylist[k]
                            if i == 0:
                                outdict[hexkey] = np.zeros(nfiles)
                            harr = outdict[hexkey]
                            harr[i] = hexvallist[k]
                    else:
                        if (i == 0):
                            outdict[hkey] = np.zeros(nfiles)
                        harr = outdict[hkey]
                        harr[i] = ilist[1]
        
    return outdict

# for a set of exposure ID strings from gathermanydonutsum,
# getndonutscam reads the individual donut fit info for each exposure
# and gets the number of donuts used per CI camera.
# stampdir is the directory where the donut stamp fits are.
# summarydirdir is where the donut_summary files are from donutana.
# this is really slow for a lot of files. 
def getndonutscam(infodict, stampdirname, namematch='first'):

    expnamelist = infodict['expname']
    nexp = len(expnamelist)
    infodict['ndonuts_cin'] = np.zeros(nexp)
    infodict['ndonuts_cis'] = np.zeros(nexp)
    infodict['ndonuts_cie'] = np.zeros(nexp)
    infodict['ndonuts_ciw'] = np.zeros(nexp)
    infodict['ndonuts_cic'] = np.zeros(nexp)
    
    for i in range(nexp):
        iexp = expnamelist[i]
        ndonutdict = ndonutscamera(stampdirname, iexp,namematch, docut=True)
        for ikey in ndonutdict.keys():
            infodict[ikey][i] = ndonutdict[ikey]

    return(infodict)

def showdonutfitresults(infodict, stampdirname, namematch='first'):

    expnamelist = infodict['expname']
    nexp = len(expnamelist)
    for i in range(nexp):
        iexp = expnamelist[i]
        rematch = re.search('[0-9]+',iexp)
        expnum = rematch.group(0)
        expdirname = stampdirname + '/' + expnum
        fitdict = getdonutfitresults(expdirname, iexp, namematch)
        xpos = fitdict['XPOS']
        ypos = fitdict['YPOS']
        cutval = fitdict['DONUTCUT']
        icin, = np.where((ypos > 1.) & (cutval == True))
        icis, = np.where((ypos < -1.) & (cutval == True))
        icie, = np.where((xpos < -1.) & (cutval == True))
        iciw, = np.where((xpos > 1.) & (cutval == True))
        icic, = np.where((ypos > -1) & (ypos < 1.) & (xpos > -1) & (xpos < 1) & (cutval == True))
        rcin, = np.where(ypos > 1.)
        rcis, = np.where(ypos < -1.)
        rcie, = np.where(xpos < -1.)
        rciw, = np.where(xpos > 1.)
        rcic, = np.where((ypos > -1) & (ypos < 1.) & (xpos > -1) & (xpos < 1))
        print(iexp + ' N fit: ' + str(len(rcic)) + ' ' + str(len(rcin)) + ' ' + str(len(rcis)) + ' ' + str(len(rcie)) + ' ' + str(len(rciw)) + ' total: ' + str(len(rcic)+len(rcin)+len(rcis)+len(rcie)+len(rciw)) + ' after cut: ' + str(len(icic)) + ' ' + str(len(icin)) + ' ' + str(len(icis)) + ' ' + str(len(icie)) + ' ' + str(len(iciw))+ ' total: ' + str(len(icic)+len(icin)+len(icis)+len(icie)+len(iciw)))

def showdonutfitzern(infodict, stampdirname, zernname, namematch='first'):

    expnamelist = infodict['expname']
    nexp = len(expnamelist)
    for i in range(nexp):
        iexp = expnamelist[i]
        rematch = re.search('[0-9]+',iexp)
        expnum = rematch.group(0)
        expdirname = stampdirname + '/' + expnum
        fitdict = getdonutfitresults(expdirname, iexp, namematch)
        print(iexp)
        print(fitdict[zernname])
        
# for a dictionary of arrays of info from many donut summaries as read in by
# gathermanydonutsum, use getdonutscam to get the number of donuts per CI
# camera and cut on <ncut> per camera and on the CUTVAL string computed in
# <getdonutfitresults> which are meant to be the same as the cutstring passed
# to getdonutana.
def cutdonutsum(infodict, ncut, stampdir, qacuts=True, namematch='first'):
    newinfodict = getndonutscam(infodict, stampdir, namematch)
    itake,=np.where((newinfodict['ndonuts_cin'] > ncut) & (newinfodict['ndonuts_cis'] > ncut) & (newinfodict['ndonuts_cie'] > ncut) & (newinfodict['ndonuts_ciw'] > ncut))
    outdict = {}
    for ikey in newinfodict.keys():
        if ikey is not 'expname':
            outdict[ikey] = newinfodict[ikey][itake]
    outexplist = []
    for i in range(len(itake)):
        outexplist.append(newinfodict['expname'][itake[i]])
    outdict['expname'] = outexplist
    return outdict

# pretty-print the contents of the dictionaries of many donutana analysis
# of exposures as read in by gathermanydonutsum
def printmanydonutsum(idict, idx=None):

    if idx is None:
        idx = np.arange(len(idict['ndonuts_used']))
    
    print('exposure name')
    print(idict['expname'])
    print('ndonuts_used')
    print(idict['ndonuts_used'][idx])
    print(' ')

    print('dodx')
    print(idict['dodx'][idx])
    print('mean dx '+"{:f}".format(np.mean(idict['dodx'][idx]))+' std '+"{:f}".format(np.std(idict['dodx'][idx])))
    print('dodxerr ')
    print(+idict['dodxerr'][idx])
    print('mean err '+"{:f}".format(np.mean(idict['dodxerr'][idx])))
    print(' ')

    print('dody')
    print(idict['dody'][idx])
    print('mean dy '+"{:f}".format(np.mean(idict['dody'][idx]))+' std '+"{:f}".format(np.std(idict['dody'][idx])))
    print('dodyerr ')
    print(+idict['dodyerr'][idx])
    print('mean err '+"{:f}".format(np.mean(idict['dodyerr'][idx])))
    print(' ')

    print('dodz')
    print(idict['dodz'][idx])
    print('mean dodz '+"{:f}".format(np.mean(idict['dodz'][idx]))+' std '+"{:f}".format(np.std(idict['dodz'][idx])))
    print('dodzerr ')
    print(+idict['dodzerr'][idx])
    print('mean err '+"{:f}".format(np.mean(idict['dodzerr'][idx])))
    print(' ')


    print('doxt')
    print(idict['doxt'][idx])
    print('mean xt '+"{:f}".format(np.mean(idict['doxt'][idx]))+' std '+"{:f}".format(np.std(idict['doxt'][idx])))
    print('dotxerr ')
    print(+idict['doxterr'][idx])
    print('mean err '+"{:f}".format(np.mean(idict['doxterr'][idx])))
    print(' ')

    print('doyt')
    print(idict['doyt'][idx])
    print('mean yt '+"{:f}".format(np.mean(idict['doyt'][idx]))+' std '+"{:f}".format(np.std(idict['doyt'][idx])))
    print('doyterr ')
    print(+idict['doyterr'][idx])
    print('mean err '+"{:f}".format(np.mean(idict['doyterr'][idx])))
    print(' ')    

    print('z6meanDeltaAfter')
    print(idict['z6meanDeltaAfter'][idx])
    print('z6thetax')
    print(idict['z6thetax'][idx])
    print('z6thetay')
    print(idict['z6thetay'][idx])
    print('z5meanDeltaAfter')
    print(idict['z5meanDeltaAfter'][idx])
    print('z7meanDeltaAfter')
    print(idict['z7meanDeltaAfter'][idx])
    print('z8meanDeltaAfter')
    print(idict['z8meanDeltaAfter'][idx])


# plot donutana results against various error, rms and residuals returned
# by donutana
def plotverr(idict,idx=None):

    if idx is None:
        idx = np.arange(len(idict['ndonuts_used']))

    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z4meanDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z4meanDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z4meanDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z4meanDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z4meanDeltaAfter')
    fig.tight_layout()
    
    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z5meanDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z5meanDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z5meanDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z5meanDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z5meanDeltaAfter')
    fig.tight_layout()

    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z6meanDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z6meanDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z6meanDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z6meanDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z6meanDeltaAfter')
    fig.tight_layout()

    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z7meanDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z7meanDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z7meanDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z7meanDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z7meanDeltaAfter')
    fig.tight_layout()

    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z8meanDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z8meanDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z8meanDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z8meanDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z8meanDeltaAfter')
    fig.tight_layout()

    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z4rmsDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z4rmsDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z4rmsDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z4rmsDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z4rmsDeltaAfter')
    fig.tight_layout()
    
    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z5rmsDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z5rmsDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z5rmsDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z5rmsDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z5rmsDeltaAfter')
    fig.tight_layout()
    
    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z6rmsDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z6rmsDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z6rmsDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z6rmsDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z6rmsDeltaAfter')
    fig.tight_layout()

    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z7rmsDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z7rmsDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z7rmsDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z7rmsDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z7rmsDeltaAfter')
    fig.tight_layout()

    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['z8rmsDeltaAfter'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['z8rmsDeltaAfter'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['z8rmsDeltaAfter'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['z8rmsDeltaAfter'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'z8rmsDeltaAfter')
    fig.tight_layout()

    fig,ax = plt.subplots(2,2)
    ax[0,0].plot(idict['ndonuts_used'][idx],idict['dodx'][idx],'bo')
    ax[0,1].plot(idict['ndonuts_used'][idx],idict['dody'][idx],'bo')
    ax[1,0].plot(idict['ndonuts_used'][idx],idict['doxt'][idx],'bo')
    ax[1,1].plot(idict['ndonuts_used'][idx],idict['doyt'][idx],'bo')
    ax[0,0].set_ylabel('dodx')
    ax[0,1].set_ylabel('dody')
    ax[1,0].set_ylabel('doxt')
    ax[1,1].set_ylabel('doyt')
    fig.text(0.5,0.01,'ndonuts_used')
    fig.tight_layout()
