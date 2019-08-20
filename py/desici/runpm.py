import numpy as np
import re
from astropy.io import fits
import subprocess
import os
import sys
import pdb

# this copied from working version in rungfaproctests on nersc 20190818
def rerungfaproc(seqlist=None):    
    if seqlist is None:
        seqlist = os.listdir(os.environ['PMDIR'])
    os.chdir(os.environ['PMDIR'])
    nseq = len(seqlist)
    for i in range(nseq):
        iseq = seqlist[i]
        if re.search('^[0-9]',iseq) == None:
            continue
        os.chdir(iseq)
        cmdstr1 = 'confignext ' + iseq
        oneout = subprocess.run(cmdstr1,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if oneout.returncode > 0:
            print('ERROR comfignext '+iseq + ' '+cmdstr1)
            #print(oneout.stdout.decode('utf-8'))
            os.chdir('..')
            continue
        cmdstr2 = 'gfaproc ' + iseq
        twoout = subprocess.run(cmdstr2,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if twoout.returncode > 0:
            print('ERROR gfaproc '+iseq)
            #print(twoout.stdout.decode('utf-8'))
        os.chdir('..')

# also taken from gfaproctests.ipynb on nersc 20190818
# call as:
#successlist,fewmatchlist,badpointinglist,nomatchlist,missinghdulist,missinghduerrlist,badimagelist,badnfsproclist,otherlist = sortgfaprocout()
def sortgfaprocout(seqlist=None):
    if seqlist is None:
        seqlist = os.listdir(os.environ['PMDIR'])
    # set up lists for the seqids that correspond to each type of PM error, or success
    # list of successful PM sequences
    successlist = []
    # missing HDU in the data
    missinghdulist = []
    # save error message
    missinghduerrlist = []
    # no gfa.fits file but nsfProc was successful
    badimagelist = []
    # nsfproc not successful
    badnfsproclist = []
    # pointing correction too big
    badpointinglist = []
    # too few matches
    fewmatchlist = []
    # no matches running imgStarMatch
    nomatchlist = []
    # other
    otherlist = []
    os.chdir(os.environ['PMDIR'])
    nseq = len(seqlist)
    for i in range(nseq):
        iseq = seqlist[i]
        if re.search('(^([0-9]+)$)',iseq) == None:
            continue
        os.chdir(iseq)
        # get all gfaproc.out files
        alllist = os.listdir()
        gfaoutlist = []
        gfaimagefile = ''
        for i in range(len(alllist)):
            if re.search('gfaproc',alllist[i]) !=None:
                gfaoutlist.append(alllist[i])
            if re.search('(^gfa-)([0-9]+)(\.)([0-9])(\.)(fits$)',alllist[i]) != None:
                gfaimagefile = alllist[i]
        # now find the most recent gfaproc.out file
        configlist = []
        for i in range(len(gfaoutlist)):
            splitlist = gfaoutlist[i].split(sep='.')
            configlist.append(int(splitlist[1]))
        configarr = np.asarray(configlist)
        idxmax = np.argmax(configarr)
        # we will need the actual PM seq number later
        maxnum = configarr[idxmax]
        # get most recent gfaproc.out file name
        gfaprocfile = gfaoutlist[idxmax]
        # now add this seqid to the list corresponding to what the error was
        # first check of gfaproc ended with Glorious Success
        cmdstr = 'grep Glorious ' + gfaprocfile
        ret = subprocess.run(cmdstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if ret.returncode == 0:
            successlist.append(iseq)
            os.chdir('..')
            continue
        # now check for nfsproc failing
        nfsprocname = 'nfsproc-'+iseq+'.'+'{:d}'.format(maxnum-1)+'.out'
        cmdstr = 'grep Glorious ' + nfsprocname
        ret = subprocess.run(cmdstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if ret.returncode == 1:
            badnfsproclist.append(iseq)
            os.chdir('..')
            continue
        # check if the gfa.fits file from image acquisition is missing
        if re.search('fits',gfaimagefile) == None:
            badimagelist.append(iseq)
            os.chdir('..')
            continue
        # check for no missing image in the gfa.image file. error is "can't read HDU(<name>): no such element 
        #in array" but grep for HDU b/c the camera nameas are in successful gfaproc.out files
        cmdstr = 'grep HDU '+ gfaprocfile
        ret = subprocess.run(cmdstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if ret.returncode == 0:
            missinghdulist.append(iseq)
            missinghduerrlist.append(ret.stdout)
            os.chdir('..')
            continue
        # look for pointing correction too big
        cmdstr = 'grep \'Pointing corrections too big\' '+ gfaprocfile
        ret = subprocess.run(cmdstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if ret.returncode == 0:
            badpointinglist.append(iseq)
            os.chdir('..')
            continue
        # now look for too few matches to gaia stars
        cmdstr = 'grep \'too few matches\' ' + gfaprocfile
        ret = subprocess.run(cmdstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if ret.returncode == 0:
            fewmatchlist.append(iseq)
            os.chdir('..')
            continue
        # look for no matches running imgStarMatch
        cmdstr = 'grep \'No matches running imgStarMatch\' ' + gfaprocfile
        ret = subprocess.run(cmdstr,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if ret.returncode == 0:
            nomatchlist.append(iseq)
            os.chdir('..')
            continue
        # if we got this far, we have an error but not one we specifically search for
        otherlist.append(iseq)
        os.chdir('..')
    os.chdir('..')
    return successlist,fewmatchlist,badpointinglist,nomatchlist,missinghdulist,missinghduerrlist,badimagelist,badnfsproclist,otherlist
                

