
# ZBFitter8.py
# Fits ZB distortion model to each of 25 ADC states. 
# using added noise to eliminate the fit singulairty.  Works. 
# M.LamptonUCB SSL  15 July 2019

import csv
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

PROGNAME = 'ZBFitter8.py'

LUT = [0,  1,  2,  5,  6,   9,   20,  27, 28, 29, 30]   # 11 polynomials
# parm 0,  1,  2,  3,  4,   5,   6,   7,  8,  9,  10    # parm numbers
# ZB:  S2, S3, S4, S7, S8,  S11, S22, T4, T7, T8, T11    # Zhao-Burge labels

NPARMS = len(LUT)


infilename = 'RT170_Nrings_8_Rmax_28.txt'  # 5425 records = 25 x 217

def getFileArray(filename):
    nums2D = list()    
    with open(filename) as f:
        data = f.readlines()
    for row in data:
        numrow = list()
        words = row.split()
        for word in words:
            try:
                x = float(word)
            except:
                x = -0.0
            numrow.append(x)
        nums2D.append(numrow)
    #---Make nums2D rectangular, lest asarray() will return a 2D list---
    nrows = len(nums2D)
    ncols = 0
    for i in range(nrows):
        ncols = max(ncols, len(nums2D[i]))
    print(' ncols = ' + str(ncols))
    for i in range(nrows):
        while len(nums2D[i]) < ncols:
            nums2D[i].append(-0.0)
    myArray = np.asarray(nums2D)
    return myArray


#-------ZERNIKE NORMALIZED FUNCTIONS & DERIVATIVES--------
#-------------Using {n,m} definitions---------------------

def factorial(n):
    if n > 1:
       return int(n*factorial(n-1))
    else:
       return 1
       
def convertNolltoBW(noll):
    # converts a Noll Zernike index to the B&W {n,m,t} triplet
    n = int(-0.5+np.sqrt(2*noll-1.75))
    m = 0
    base = int(0.5*n**2 +0.5*n + 1)
    diff = noll - base
    if n%2==0:
        m = 2*int(0.5*diff + 0.7)
    else:
        m = 2*int(0.5*diff + 1.2) - 1
    if noll%2>0:  
        m = -m
    return np.array([n, m])

def convertWyanttoBW(wyant):
    # converts a Wyant Zernike index to the B&W {n,m,t}   
    halfsum = int(np.sqrt(wyant))
    idif = int(wyant - halfsum**2)
    halfdif = int (idif/2)
    n = halfsum + halfdif
    m = halfsum - halfdif
    if idif%2 >0:
        m = -m
    return np.array([n, m]) 
    
def getZernFuncXY(nm, xnorm, ynorm):   # BIG NINE:  #1
    # Here, xnorm and ynorm must lie within the unit circle
    rnorm = np.sqrt(xnorm*xnorm + ynorm*ynorm)
    angle = np.arctan2(ynorm,xnorm)
    return getZernRadial(nm,rnorm) * getZernAngular(nm,angle)

    
def getZernRadial(nm, rnorm):    # BIG NINE: #4
    n = nm[0]             # B&W
    m = np.abs(nm[1])     # B&W
    halfsum = (n+m)/2
    idif = n-m
    halfdif = int(idif/2)
    # n = halfsum + halfdif   # or, halfsum = (n+m)/2
    # m = halfsum - halfdif   # or, halfdif = (n-m)/2
    # loop through the polynomial
    result = 0.
    for i in range(0, halfdif+1):
        expon = int(n-2*i)
        sign = 1 if i%2 == 0 else -1
        numer = sign * factorial(n-i)
        denom = factorial(i) * factorial(halfsum-i) * factorial(halfdif-i)
        coef = numer / denom
        term = coef*math.pow(rnorm, expon)
        result = result + term
    return result  
    
def getZernAngular(nm, theta): 
    m = nm[1]    # B&W
    if m==0:
        return 1.
    if m>0:
        return math.cos(m*theta)
    m = np.abs(m)         # note this abs() function
    return math.sin(m*theta)  

def zernFormulaText(nm):           # BIG NINE: #8
    #---generates a text representation of the specified Zernike function--
    n = nm[0]  # B&W
    m = nm[1]  # B&W
    # print 'New zernFormulaText() is using n, m: ', n, m
    needsine = True if m<0 else False
    m = np.abs(m)
    halfsum = (n+m)/2
    idif = n-m
    halfdif = int(idif/2)
    
    #--first do the radial part----
    nterms = 0
    s = ''
    #--evaluate the radial polynomial-----
    for i in range(0, halfdif+1):
        nterms = nterms + 1
        # print "Starting with n,  m, i, nterms = ", n, m, i, nterms
        expon = int(n-2*i)                # start with highest exponent
        # print "  expon = ", expon
        sign = 1 if i%2 == 0 else -1      # alternating signs in Zernike series
        strsign = '+' if i%2==0 else '-'
        numer = sign * factorial(n-i)
        denom = factorial(i) * factorial(halfsum-i) * factorial(halfdif-i)
        coef = numer / denom
        scoef = str(coef)
        if coef==1 and expon>0:          # suppress showing coef=1
            scoef = ''
        if coef > 0 and nterms > 1:
            scoef = '+' + scoef
        s = s + scoef
        if expon > 0:
            s = s + 'r'
        if expon > 1:
            s = s + '^' + str(expon)
    if nterms>1 and m!=0:
        s = '('+ s + ')'
    #--then do the azimuthal part, if any--------
    if m==0:
        return s
    strm = ''
    if m>1:    
        strm = str(m)
    if needsine:
        s = s + '*sin(' + strm + 't)'
    else:
        s = s + '*cos(' + strm + 't)'
    return s

#--END ZERNIKES with {n,m} B&W indexing------------




#-----Zhao-Burge functions built from Zernikes---------------------

rh = np.sqrt(0.5)
rt = np.sqrt(2.0)


NPARMS = len(LUT)
  
def getZ(noll, x, y):
    return getZernFuncXY(convertNolltoBW(noll), x, y)

def getZhaoBurgeTerm(whichparm, x, y):
    # Given a Lampton index "which" 0....19 and an object point x, y,
    # fetches the needed the Zernikes via their Noll numbers..
    # Returns the modeled image point deviation x, y.
    # Case numbers shown are from Zhao & Burge Tables 1 and 2. 
    which = LUT[whichparm]
    if which==0:   # case "S2"  r^0, keep; X translate
        return getZ(1,x,y), 0.0
        
    if which==1:   # case "S3"  r^0, keep; Y translate
        return 0.0, getZ(1,x,y)
        
    if which==2:   # case "S4"   r^1, keep: magnification
        return rh*getZ(2,x,y), rh*getZ(3,x,y)
        
    if which==3:   # case "S5"   r^1   -1ppm Mangled, 3ppm PartlyMangled, 4ppm ADC45
        return rh*getZ(3,x,y), rh*getZ(2,x,y)
        
    if which==4:   # case "S6"    r^1  24ppm Mangled,  6ppm PartlyMangled,-5ppm ADC45
        return rh*getZ(2,x,y), -rh*getZ(3,x,y)
        
    if which==5:   # case "S7"   r^2 125ppm Mangled  -908ppm PartlyMangled, ADC=423
        return 0.5*getZ(5,x,y), rh*getZ(4,x,y)-0.5*getZ(6,x,y)
        
    if which==6:   # case "S8"    r^2, 1223 Mangled, 653 PartlyMangled 10ppm ADC
        return rh*getZ(4,x,y)+0.5*getZ(6,x,y), 0.5*getZ(5,x,y)
        
    if which==7:   # case "S9"   r^3, 77ppm Mangled, 383 PartlyMangled, 116ppm ADC
        return rh*getZ(5,x,y), rh*getZ(6,x,y)
        
    if which==8:   # case "S10"   r^2, -6ppm Mangled  -1ppm PartlyMangled, zero ADC
        return rh*getZ(6,x,y), -rh*getZ(5,x,y)
        
    if which==9:   # case "S11"   r^3, huge
        return rh*getZ(8,x,y), rh*getZ(7,x,y)
        
    if which==10:  # case "S12"   r^3, 2ppm Mangled,zero PartlyMangled, zero ADC
        return 0.5*getZ(8,x,y)+0.5*getZ(10,x,y), -0.5*getZ(7,x,y)+0.5*getZ(9,x,y)
        
    if which==11:  # case "S13"    r^3, zero Mangled  zero PartlyMangled, 1ppm ADC
        return 0.5*getZ(7,x,y)+0.5*getZ(9,x,y), 0.5*getZ(8,x,y)-0.5*getZ(10,x,y)
        
    if which==12:  # case "S14"    r^3, 1ppm Mangled  zero PartlyMangled zero ppm ADC
        return rh*getZ(10,x,y), -rh*getZ(9,x,y)
        
    if which==13:  # case "S15"    r^3, zero Mangled   zeroPartlyMangled  zero ppm ADC
        return rh*getZ(9,x,y), rh*getZ(10,x,y)
        
    if which==14:  # Case "S16"  r^4  38ppm Mangled; 9ppm PartlyMangled, zeroADCero ADC
        return rh*getZ(11,x,y)+0.5*getZ(12,x,y), 0.5*getZ(3,x,y)
        
    if which==15:  # Case "S17":   r^4   15ppmMangled, -7ppm PartlyMangled,   8ppm ADC
        return 0.5*getZ(3,x,y), rh*getZ(11,x,y)-0.5*getZ(12,x,y)
        
    if which==16:  # Case "S18"  r^4   -6ppm Mangled, zero PartlyMangled, zero ADC
        return 0.5*getZ(12,x,y)+0.5*getZ(14,x,y), 0.5*getZ(15,x,y)-0.5*getZ(13,x,y)
        
    if which==17:  # Case "S19"   r^4   2ppm Mangled  -1ppm PartlyMangled  1ppm ADC;
        return 0.5*getZ(13,x,y)+0.5*getZ(15,x,y), 0.5*getZ(12,x,y)-0.5*getZ(14,x,y)
        
    if which==18:  # Case "S20"  r^4   zero Mangled, zero PartlyMangled, zero ADC
        return rh*getZ(14,x,y), -rh*getZ(15,x,y)
        
    if which==19:  # Case "S21"  r^4   zero Mangled  1ppm PartlyMangled, zero ADC
        return rh*getZ(15,x,y), rh*getZ(14,x,y)
        
    if which==20:  # Case "S22"  r^5  171 ppm Mangled  172ppm PartlyMangled, 172ppm ADC
        return rh*getZ(16,x,y), rh*getZ(17,x,y)
        
    if which==21:  # Case "S23"  r^5  zero Mangled  1ppm PartlyMangled,  zero ADC
        return 0.5*getZ(17,x,y)+0.5*getZ(19,x,y), 0.5*getZ(16,x,y)-0.5*getZ(18,x,y)
        
    if which==22:  # Case "S24"   r^5  -1ppm Mangled  zeroPartlyMangled,zero ADC
        return 0.5*getZ(16,x,y)+0.5*getZ(18,x,y), -0.5*getZ(17,x,y)+0.5*getZ(19,x,y)
        
    if which==23:  # Case "S25"  r^5   zero Mangled  zero PartlyMangled zero ADC
        return 0.5*getZ(19,x,y)+0.5*getZ(21,x,y), 0.5*getZ(18,x,y)-0.5*getZ(20,x,y)
        
    if which==24:  # Case "S26"   r^5  zero ppm  zero PartlyMangled  zero ADC
        return 0.5*getZ(18,x,y)+0.5*getZ(20,x,y), -0.5*getZ(19,x,y)+0.5*getZ(21,x,y)

    if which==25:  # Case "S27"   r^5   zero ppm  zero, PartlyMangled; zero ADC
        return rh*getZ(21,x,y), rh*getZ(20,x,y)

    if which==26:  # Case "S28"  r^5  zero ppm  zero PartlyMangled, zero ADC
        return rh*getZ(20,x,y), -rh*getZ(21,x,y)

    if which==27:  #  case "T4"   r^1,  huge, -6ppm ADC.  Roll.
        return rh*getZ(3,x,y), -rh*getZ(2,x,y)
        
    if which==28:  # case "T7"    r^2,  -265ppm Mangled, 131 PartlyMangled, -4ppm ADC
        return rh*getZ(4,x,y)-0.5*getZ(6,x,y), -0.5*getZ(5,x,y)
        
    if which==29:  # case "T8"   r^2, -105ppm Mangled.  -544ppm PartlyMangled, -163 ADC
        return 0.5*getZ(5,x,y), -rh*getZ(4,x,y)+0.5*getZ(6,x,y)
        
    if which==30:  # case "T11"   r^3, 358ppm Mangled, 358 PartlyMangled  zero ADC
        return rh*getZ(7,x,y), -rh*getZ(8,x,y)
        
    if which==31:  # case "T12"   r^3,  -1ppm Mangled, zero PartlyMangled, zero ADC
        return -0.5*getZ(7,x,y)+0.5*getZ(9,x,y), -0.5*getZ(8,x,y)-0.5*getZ(10,x,y)
        
    if which==32:  # case "T13"  r^3, +1ppm Mangled,  zero PartlyMangled, zero ADC
        return 0.5*getZ(8,x,y)-0.5*getZ(10,x,y), -0.5*getZ(7,x,y)-0.5*getZ(9,x,y)
        
    print("ZhaoBurgeTerm() is exitting because which = ", which)
    quit()


"""
def doRemapping(myParms):
    # feeds each sky point {u,v} into getZhaoSum()
    # will be called by optimizer to try out a mix of parameters.
    for i in range(0, len(SkyGrid)):
        u = SkyGrid[i, 0]
        v = SkyGrid[i, 1]
        x, y = getZhaoSum(myParms, u, v)
        # print 'result:  i, x, y ={:6d}{:12.6f}{:12.6f}'.format(i, x, y) 
        X[i] = x
        Y[i] = y
"""     

def zbFunc(uv, *parms): # uv = 1D array input concatenate(u[],v[])
    # Model: any number of Zhao & Burge parameters
    # note: curve_fit() calls its given func(args, *parms)
    # print("anyFunc() hasreceived parms = ")
    # print(["{0:0.6f}".format(xx) for xx in parms])
    pq = np.zeros(len(uv))
    half = len(uv)//2
    for i in range(0, half):  # EACH OBJECT POINT NOT AN ARRAY
        u = uv[i]
        v = uv[i+half]
        x, y = getZhaoSum(u, v, *parms)  # always args then parms
        pq[i] = x
        pq[i+half] = y
    return pq   

def getZhaoSum(u, v, *coefs):
    # This models the distortion of ONE SKY POINT onto ONE FOCAL POINT.
    # Sums over all terms.
    # coefs is the array of 20 coefficients
    # u and v are Cartesian sky coordinates mapped into unit radius circle.
    x = u   # initially undeviated
    y = v   # initially undeviated
    for index in range(0, len(coefs)):
        dx, dy = getZhaoBurgeTerm(index, u, v)
        x += dx * coefs[index]
        y += dy * coefs[index]
    return x, y

def radialFunc(uv, *parms):
    # 6 parms model: dx, dy, a1, a3, a5, a7
    # untested!
    pq = np.zeros(len(uv))
    half = len(uv)//2
    for i in range(0, half):
        u = uv[i]
        v = uv[i+half]
        r = np.sqrt(u*u+v*v)
        if r==0.:
            u = 1E-12
        t = np.arctan2(u,v)
        rpoly = (1.+p[2])*r + p[3]*r**3 + p[4]*r**5 + p[5]*r**7
        x = p[0] + np.cos(t)*rpoly
        y = p[1] + np.sin(t)*rpoly
        pq[i] = x
        pq[i+half] = y
    return pq
     
#-------------------Main program-------------------------    


alldata = getFileArray(infilename)
print(alldata.shape)  # 5208 x 9; 24 fields of 217 spots each
print('  ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms')
Nadcs     = 25
Neach     = 217
NRINGS    = 8
UVMAX     = 28     # milliradians
XYMAX     = 407    # millimeters
NOISE     = 1.0E-6

# first get the forward (sky to FP) fits to the 25 ADC tasks

outfilename = 'ZBF8fwd_Nrings_' + str(NRINGS) + '_Rmax_' + str(UVMAX) + '.txt'

with open(outfilename, 'w') as outfile:
    for adc in range(Nadcs):
        adc12 = np.copy(alldata[adc*Neach:(adc+1)*Neach, 0:2])    # rows,cols
        uv    = np.copy(alldata[adc*Neach:(adc+1)*Neach, 2:4])    # rows,cols
        xy    = np.copy(alldata[adc*Neach:(adc+1)*Neach, 5:7])    # rows,cols
        uv *= 1000./UVMAX  # for unit circle scaling
        xy *= 1./XYMAX     # for unit circle scaling
        xy += np.random.normal(0.0, NOISE, xy.shape)
        print()

        pinitial = np.zeros(NPARMS)
        args     = np.concatenate([uv[:,0], uv[:,1]])
        goals    = np.concatenate([xy[:,0], xy[:,1]])
        popt, covar = curve_fit(zbFunc, args, goals, pinitial)
        #  note: curve_fit() calls its given func(args, *parms)
        print('\npopt: \n', popt)
        perror = np.sqrt(covar.diagonal())
        print('\nerrors: \n', perror)
        combo = np.concatenate((popt, perror))  # double parentheses
        result = ''
        for item in combo:
            result += '{:12.6f}'.format(item)
        outfile.writelines(result + '\n') # yes, plural even for one line

# then get the reverse (FP to sky) fits to the 25 ADC tasks

outfilename = 'ZBF8rev_Nrings_' + str(NRINGS) + '_Rmax_' + str(UVMAX) + '.txt'

with open(outfilename, 'w') as outfile:
    for adc in range(Nadcs):
        adc12 = np.copy(alldata[adc*Neach:(adc+1)*Neach, 0:2])    # rows,cols
        uv    = np.copy(alldata[adc*Neach:(adc+1)*Neach, 2:4])    # rows,cols
        xy    = np.copy(alldata[adc*Neach:(adc+1)*Neach, 5:7])    # rows,cols
        uv *= 1000./UVMAX  # for unit circle scaling
        xy *= 1./XYMAX     # for unit circle scaling
        xy += np.random.normal(0.0, NOISE, xy.shape)
        print()

        pinitial = np.zeros(NPARMS)
        goals    = np.concatenate([uv[:,0], uv[:,1]])
        args     = np.concatenate([xy[:,0], xy[:,1]])
        popt, covar = curve_fit(zbFunc, args, goals, pinitial)
        #  note: curve_fit() calls its given func(args, *parms)
        print('\npopt: \n', popt)
        perror = np.sqrt(covar.diagonal())
        print('\nerrors: \n', perror)
        combo = np.concatenate((popt, perror))  # double parentheses
        result = ''
        for item in combo:
            result += '{:12.6f}'.format(item)
        outfile.writelines(result + '\n') # yes, plural even for one line

