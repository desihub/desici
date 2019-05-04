import numpy as np
import json
import re
import copy
import random

# test scripts:
# import json
# json.load(open('scriptfile.json'))

# helper function to write pretty output
def writescript(scriptdictlist, sfname, scriptdir='/data/exposure_scripts/'):
    sfnameout = scriptdir + '/' + sfname + '.json'
    with open(sfnameout,"w") as fx:
        json.dump(scriptdictlist,fx,indent=2)
    return sfnameout

# does not set pointing or ADC
basicexpfoc = {"sequence":"CI", "flavor":"science","exptime":0,"program":"SetMe","focus":0.0}

# set focus and take and exposure
def expfoc(programname, focval, sfname, exptimeval, scriptdir='/data/exposure_scripts/'):
    myscript = copy.deepcopy(basicexpfoc)
    myscript['program'] = programname
    myscript['exptime'] = exptimeval
    myscript['focus'] = focval
    explist = [myscript]
    sfnameout = scriptdir + '/' + sfname + '.json'
    with open(sfnameout,"w") as fx:
        json.dump(explist,fx,indent=2)
    return sfnameout

# take an exposure
def expscript(programname, exptimeval, sfname, scriptdir='/data/exposure_scripts/'):

    myscript={"sequence":"CI", "flavor":"science","exptime":0,"program":"SetMe"}
    myscript['program'] = programname
    myscript['exptime'] = exptimeval
    expscript = [myscript]
    sfnameout = scriptdir + '/' + sfname + '.json'
    with open(sfnameout,"w") as fx:
        json.dump(explist,fx,indent=2)
    return sfnameout

def slewall(reqra, reqdec, exptime, program, focarray, sfname, useadc=False, scriptdir='/data/exposure_scripts/'):
    slewdict = {'sequence':'CI','flavor':'science','exptime':exptime,'program':program,'reqra':reqra,'reqdec':reqdec, 'focus':focarray, 'correct_for_adc':useadc}
    explist = [slewdict]
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename

# exptime for the fvc should always be 1.0 unless there is a good reason
def fvcadcscript(adc1poslist,adc2poslist, sfname, fvcexptimeval=1.0, scriptdir='/data/exposure_scripts/'):
    explist = []
    basefvc = {'sequence':'FVC','flavor':'science','exptime':fvcexptimeval,'fiducials':'on','leave_fiducials':'off','program':'fvc adc','adc':[0,0]}
    npos = len(adc1poslist)
    for i in range(npos):
        iadc = copy.deepcopy(basefvc)
        iadc['adc'] = [adc1poslist[i],adc2poslist[i]]
        explist.append(iadc)
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename

def cfocussweep(dfoc, dfocstart, nfoc, exptime, sfname, scriptdir='/data/exposure_scripts/'):
    
    explist = []
    template = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'Alignment','relfocus':'True','focus':0.0}
    dfocarr = np.arange(nfoc)*dfoc+dfocstart
    for i in range(nfoc):
        idict = copy.deepcopy(template)
        idict['focus'] = dfocarr[i]
        explist.append(idict)
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename


# assumes you have already set telescope position
# can set the focus value in the defaultfocvec vector
# and then put None in element 2 of defaultfocvec
# does not restore focus at the end, but it does restore tilts
# NOTE focs args are [dx,dy,dz,tx,ty,rot]
# put None for axes you don't want to move
# defaultfoc is the set of default values for all hexapod axes
def tiltsweep(defaultfocvec, exptime, sfname, focoffset=None, tiltvals = [-150,-100,-50,0,50,100], scriptdir='/data/exposure_scripts/'):
    
    explist = []
    itx = 3
    ity = 4

    if focoffset != None:
        fdict = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'Alignment setup','relfocus':True,'focus':focoffset}
        explist.append(fdict)

    template = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'Alignment','focus':defaultfocvec}
    # first tilt x
    for i in range(len(tiltvals)):
        idict = copy.deepcopy(template)
        ifoc = idict['focus']
        ifoc[itx] = tiltvals[i]
        idict['focus'] = ifoc
        explist.append(idict)
    # now tilt in y
    for i in range(len(tiltvals)):
        idict = copy.deepcopy(template)
        ifoc = idict['focus']
        ifoc[ity] = tiltvals[i]
        idict['focus'] = ifoc
        explist.append(idict)
    idict = copy.deepcopy(template)
    explist.append(idict)
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename

# if you don't want to move in RA,Dec, give None
def donutpair(relfocoff, exptime, sfname, reqra=None, reqdec=None, scriptdir='/data/exposure_scripts/'):

    explist = []
    if (reqra != None):
        if (reqdec != None):
            slewdict = {'sequence':'CI','flavor':'science','exptime':exptime,'program':'donut setup','reqra':reqra,'reqdec':reqdec}
            explist.append(slewdict)
        else:
            print('Specify both reqra and reqdec!!')
            exit
    
    ddictplus = {'sequence':'CI','flavor':'science','exptime':exptime,'program':'Donut','relfocus':True,'focus':relfocoff}
    ddictminus = {'sequence':'CI','flavor':'science','exptime':exptime,'program':'Donut','relfocus':True,'focus':-2*relfocoff}
    cleanup = {'sequence':'CI','flavor':'science','exptime':exptime,'program':'Donut','relfocus':True,'focus':relfocoff}

    explist.append(ddictplus)
    explist.append(ddictminus)
    explist.append(cleanup)
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename

# take a set of data for hexpod sensitivity matrix
# assumes you have slewed the telescope
# it will come out of focus if focoffset is set to a relative value
# set element 2 of defaultfocvec to None to then leave focus at focoffset
# for the entire script. Note it does not restore focus at the end
# but it does restore the tilts and decenters
# it goes to ABSOLUTE values of tilt and decenter
def hexcal(defaultfocvec, exptime, sfname, focoffset=None, scriptdir='/data/exposure_scripts/'):
    
    explist = []
    decenvals = [-4000,-2000,0,2000]
    tiltvals = [-150,-75,0,75,150]
    
    if focoffset != None:
        fdict = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'Alignment setup','relfocus':True,'focus':focoffset}
        explist.append(fdict)
    
    tplt = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'Alignment','focus':defaultfocvec}
    # dx
    for i in range(len(decenvals)):
        idict = copy.deepcopy(tplt)
        idict['focus'][0] = decenvals[i]
        explist.append(idict)
    # dy
    for i in range(len(decenvals)):
        idict  = copy.deepcopy(tplt)
        idict['focus'][1] = decenvals[i]
        explist.append(idict)
    # xtilt
    for i in range(len(tiltvals)):
        idict  = copy.deepcopy(tplt)
        idict['focus'][3] = tiltvals[i]
        explist.append(idict)
     # ytilt
    for i in range(len(tiltvals)):
        idict  = copy.deepcopy(tplt)
        idict['focus'][4] = tiltvals[i]
        explist.append(idict)
    # put things back, but don't reset focus
    explist.append(tplt)
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename


# 40 points in current LUT
# lutfile is full path to LUT file, standard Klaus format
# lineidx can be a list of lines to subsample in the LUT file
# [1-] gets all the lines, horrible hack sorry. 
def hexlutdonuts(sfname, npoints=10,lineidx=[-1],exptime=120, dfoc=660, dsize=600, ndither=5, lutfile='/software/products/HexapodLUT-trunk/TelescopeLUT.txt', scriptdir='/data/exposure_scripts/'):
    
    tablevals = readlut(lutfile,lineidx=lineidx)
    nlines = len(tablevals)
    allidx = list(range(nlines))
    random.shuffle(allidx)
    allidx.insert(0,0)

    explist = []
    expdict = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'donuts','focus':None}
    expditherdict = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'donuts','deltara':dsize,'deltadec':dsize}
    changefoc = {'sequence':'Action','action':'slew','relfocus':True,'focus':dfoc}
    replacefoc = {'sequence':'Action','action':'slew','relfocus':True,'focus':-1*dfoc}
        
    explist.append(changefoc)   

    for i in range(npoints):
        
        oneidx = allidx[i]

        itable = tablevals[oneidx]
        ireqaz = itable[0]
        ireqel = itable[1]
        ifocvec = itable[2:8]
  
        # slew
        sdict = {'sequence':'Action','action':'slew','track':True,'reqaz':ireqaz,'reqel':ireqel}
        explist.append(sdict)
        iexpdict = copy.deepcopy(expdict)
        iexpdict['focus'] = ifocvec
        explist.append(iexpdict)
        for j in range(ndither-1):
            explist.append(expditherdict)

    explist.append(replacefoc)        
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename

# read a lookup table file in Official Format and make the list of lists
# to pass to hexluttest.
def readlut(filename, lineidx=[-1]):
    data = np.genfromtxt(filename)
    nlines = len(data)
    outlist = []
    # elements we want in each line
    # alt az X,Y,Xtilt,Ytilt
    validxlist = [0,1,2,3,5,6]
    # specific lines we want
    if lineidx[0] < 0:
        lineidx = list(range(nlines))
    for i in range(len(lineidx)):
        iline = data[lineidx[i]]
        ilist = []
        for idx in validxlist:
            ilist.append(iline[idx])
        # no spin
        ilist.append(None)
        # no focus
        ilist.insert(4,None)
        ilist[2] = ilist[2]
        ilist[3] = ilist[3]
        outlist.append(ilist)
    return outlist
        
def genspintest(sfname, scriptdir='/data/exposure_scripts/'):

    fieldlist = [[90.,30.],[90.,45.],[90.,60.],[90.,75.],[90.,90.],[270.,75.],[270.,60.],[270.,45.],[270.,30.],[0.0,30.],[0.,45.],[0.,60.],[0.,75.],[0.,90.],[180.,75.],[180.,60.],[180.,45.],[180.,30.]]
    fieldlist2 = [[135.,30.],[135.,45.],[135.,60.],[135.,75.],[135.,90.],[315.,75.],[315.,60.],[315.,45.],[315.,30.],[45.,30.],[45.,45.],[45.,60.],[45.,75.],[45.,90.],[225.,75.],[225.,60.],[225.,45.],[225.,30.]]
    fieldlist.extend(fieldlist2)

    expdict =  {'sequence':'CI', 'flavor':'science', 'exptime':fieldexptime, 'program':'LUT spin check'}
    donutdict = {'sequence':'CI', 'flavor':'science', 'exptime':donutexptime, 'program':'LUT spin check donut','relfocus':True,'focus':donutfocoffset}
    undonutdict = {'sequence':'Action','action':'slew','relfocus':True,'focus':-1.0*donutfocoffset}
    explist = []
    for coord in fieldlist:
        az = coord[0]
        el = coord[1]
        islew =  {'sequence':'Action','action':'slew','track':True,'reqaz':az,'reqel':el}
        explist.append(islew)
        explist.append(expdict)
        #explist.append(donutdict)
        #explist.append(undonutdict)
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename


# disize is dither size in arcsec
# you don't need to dither by too much, the CI chips are 4x6 arcmin
# moves by dfoc and then puts it back
def donutsatonepos(exptime, sfname, dfoc=660., ndither=5, dsize=-600., scriptdir='/data/exposure_scripts/'):
    
    explist = []
    expdict = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'donuts'}
    expditherdict = {'sequence':'CI', 'flavor':'science', 'exptime':exptime, 'program':'donuts','deltara':dsize,'deltadec':dsize}
    changefoc = {'sequence':'Action','action':'slew','relfocus':True,'focus':dfoc}
    replacefoc = {'sequence':'Action','action':'slew','relfocus':True,'focus':-1*dfoc}
    
    #didx = np.arange(ndither)
    #evenidx = np.where(didx%2 < 1)
    #oddidx = np.where(didx%2 != 0)
    #dravals = np.ones(ndither)*dsize
    #dravals[evenidx] = dravals[evenidx]*-1.0
    #ddecvals = np.ones(ndither)*dsize
    #ddecvals[oddidx] = ddecvals[oddidx]*-1.0

    explist.append(changefoc)
    explist.append(expdict)
    for i in range(ndither-1):
        explist.append(expditherdict)
    explist.append(replacefoc)
    scriptfilename = writescript(explist, sfname, scriptdir=scriptdir)
    return scriptfilename
        
    
