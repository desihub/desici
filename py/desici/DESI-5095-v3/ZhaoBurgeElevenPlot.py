

# ZhaoBurgeElevenPlot.py
# 10 Sept 2018 with Area normalization, T8Y typo corrected
# shows the tadpole maps of Lampton polynomial fields numbered 0...10
# 
# 
# Uses an implicit subplot enumeration, as suggested in the StackOverflow answer
# https://stackoverflow.com/questions/13365617/large-number-of-subplots-with-matplotlib
# This nicely overcomes the discordance of subplot numbering 1,2 3.. and array indexing. 

import math
import numpy as np
import matplotlib.pyplot as plt


 
#-----------ZERNIKE NORMALIZED EDGE=1 PACKAGE-----------
#-----------Using BORN-WOLF {n,m} definitions-----------

def factorial(n):
    if n > 1:
       return int(n*factorial(n-1))
    else:
       return 1
       
def convertNolltoBW(noll):
    # converts a Noll Zernike index to the B&W {n,m} doublet
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
    # converts a Wyant Zernike index to the B&W {n,m}   
    halfsum = int(np.sqrt(wyant))
    idif = int(wyant - halfsum**2)
    halfdif = int (idif/2)
    n = halfsum + halfdif
    m = halfsum - halfdif
    if idif%2 >0:
        m = -m
    return np.array([n, m]) 
    
def getZernFuncXY(nm, xnorm, ynorm):   # nm is the BornWolf index; BIG NINE:  #1
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
        term = coef*np.power(rnorm, expon)
        result = result + term
    return result  
    
def getZernAngular(nm, theta): 
    m = nm[1]    # B&W
    if m==0:
        return 1.
    if m>0:
        return np.cos(m*theta)
    m = np.abs(m)         # note this abs() function
    return np.sin(m*theta)  

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

#----------End Zernike package with {n,m} B&W indexing---------------------


#-----Zhao-Burge functions built from Zernikes Package---------------------

rh = np.sqrt(0.5)
rt = np.sqrt(2.0)

NCOEFS = 33  # available item coefs


# LUT= [0,1,2, 4,5,6, 7,8,9, 14,15,16, 20,27,28, 29,30]    # 17 polynomials
# LUT= [0,1,2, 4,5,6, 7,9, 14,15, 20,27,28, 29,30]         # 15 polynomials
# LUT = [0,1,2, 5,6,7, 9, 14,15, 20, 27, 28,29,30]         # 14 polynomials
#   parm 0,1,2, 3,4,5, 6,  7, 8,  9, 10, 11,12,13          # parm numbers

# LUT = [0, 1, 2,  5, 6, 7,   9,  15,  20,  27, 28, 29, 30] # 13 polynomials out of 33
# parm 0, 1, 2,  3, 4, 5,   6,   7,   8,   9, 10, 11, 12    # parm numbers
#  ZB: S2,S3,S4, S7,S8,S9, S11, S17, S22, T4, T7, T8, T11   # Zhao-Burge labels

# LUT = [0, 1,  2,   5, 6,    9,  15,  20,  27, 28, 29, 30]   # 12 polynomials
# parm   0, 1,  2,   3, 4,    5,   6,   7,   8,  9, 10, 11    # parm numbers
#  ZB:  S2, S3, S4, S7,S8,   S11, S17, S22, T4, T7, T8, T11   # Zhao-Burge labels; drop S9


# LUT = [0, 1, 2,   5, 6,    9,  20,  27, 28, 29,  30]   # 11 polynomials
# parm 0, 1, 2,   3, 4,    5,   6,   7,  8,  9,  10    # parm numbers
# ZB:  S2,S3,S4,  S7,S8,   S11, S22, T4, T7, T8, T11    # Zhao-Burge labels; drop S17

# LUT = [2,  5,  6,   9,  20,  28, 29,  30]    # 8 polynomials
# parm 0,  1,  2,   3,   4,   5,  6,   7     # parm numbers
# ZB: S4, S7, S8, S11, S22,  T7, T8,  T11    # Zhao-Burge labels; drop S17

LUT = [0,  1,  2,  5,  6,   9,   20,  27, 28, 29, 30]   # 11 polynomials
# parm 0,  1,  2,  3,  4,   5,   6,   7,  8,  9,  10    # parm numbers
# ZB:  S2, S3, S4, S7, S8,  S11, S22, T4, T7, T8, T11    # Zhao-Burge labels;


NPARMS = len(LUT)


#--------Zernike-Noll Renormalization for AreaIntegral = Pi----------------

squares = np.array([0, 1, 4, 4, 3, 6, 6, 8, 8, 8, 8,  5, 10, 10, 10, 10, 12, 12, 12, 12, 12, 12, 7, 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 16, 16, 16, 16, 9])
# for Noll index =  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,12, 13, 14, 15, 16, 17, 18, 19, 20, 21,22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37 

def normalizeArea(noll):
    # converts an edge=1 Zernike to an area=Pi normalized Zernike
    return np.sqrt(squares[noll])
  
def getZ(noll, x, y):
    #  Note: area integral = Pi is this normalization
    return normalizeArea(noll) * getZernFuncXY(convertNolltoBW(noll), x, y)

def getZhaoBurgeTerm(whichparm, x, y):
    # Cartesian input x,y; cartesian output x, y, and label.
    # Uses getZ() and delivers area normalized ZB terms
    which = LUT[whichparm]
    if which==0:   # case "S2"  r^0, keep; X translate; 
        result =  getZ(1,x,y), 0.0,            'S2 = X translate'
    elif which==1:   # case "S3"  r^0, keep; Y translate; 
        result =  0.0, getZ(1,x,y),             'S3 = Y translate'
    elif which==2:   # case "S4"  magnify
        result =  rh*getZ(2,x,y), rh*getZ(3,x,y), 'S4 = Magnify'
    elif which==3:   # case "S5"
        result =  rh*getZ(3,x,y), rh*getZ(2,x,y), 'S5'  # n/a
    elif which==4:   # case "S6"
        result =  rh*getZ(2,x,y), -rh*getZ(3,x,y), 'S6'  # n/a
    elif which==5:   # case "S7" 
        result =  0.5*getZ(5,x,y), rh*getZ(4,x,y)-0.5*getZ(6,x,y), 'S7 = UpDownUp'
    elif which==6:   # case "S8"
        result =  rh*getZ(4,x,y)+0.5*getZ(6,x,y), 0.5*getZ(5,x,y), 'S8 = RightLeftRight'
    elif which==7:   # case "S9"
        result =  rh*getZ(5,x,y), rh*getZ(6,x,y), 'S9'      # parm=5/13
    elif which==8:   # case "S10"
        result =  rh*getZ(6,x,y), -rh*getZ(5,x,y), 'S10'    # n/a but is SISTER TO S9 !!!
    elif which==9:   # case "S11" 
        result =  rh*getZ(8,x,y), rh*getZ(7,x,y), 'S11 = InOut'
    elif which==10:  # case "S12"
        result =  0.5*getZ(8,x,y)+0.5*getZ(10,x,y), -0.5*getZ(7,x,y)+0.5*getZ(9,x,y), 'S12' # n/a
    elif which==11:  # case "S13" 
        result =  0.5*getZ(7,x,y)+0.5*getZ(9,x,y), 0.5*getZ(8,x,y)-0.5*getZ(10,x,y), 'S13'  # n/a
    elif which==12:  # case "S14"
        result =  rh*getZ(10,x,y), -rh*getZ(9,x,y),  'S14'   # n/a
    elif which==13:  # case "S15"
        result =  rh*getZ(9,x,y), rh*getZ(10,x,y),  'S15'   # n/a
    elif which==14:  # Case "S16"
        result =  rh*getZ(11,x,y)+0.5*getZ(12,x,y), 0.5*getZ(13,x,y), 'S16'  # n/a; TYPO! 13 not 3
    elif which==15:  # Case "S17"
        result =  0.5*getZ(13,x,y), rh*getZ(11,x,y)-0.5*getZ(12,x,y), 'S17'  # parm=7/13  TYPO! 13 not 3
    elif which==16:  # Case "S18"
        result =  0.5*getZ(12,x,y)+0.5*getZ(14,x,y), 0.5*getZ(15,x,y)-0.5*getZ(13,x,y), 'S18' # n/a
    elif which==17:  # Case "S19"
        result =  0.5*getZ(13,x,y)+0.5*getZ(15,x,y), 0.5*getZ(12,x,y)-0.5*getZ(14,x,y), 'S19'  # n/a
    elif which==18:  # Case "S20"
        result =  rh*getZ(14,x,y), -rh*getZ(15,x,y), 'S20'  # n/a
    elif which==19:  # Case "S21"
        result =  rh*getZ(15,x,y), rh*getZ(14,x,y), 'S21', # n/a
    elif which==20:  # Case "S22"
        result =  rh*getZ(16,x,y), rh*getZ(17,x,y), 'S22 = OutInOut'
    elif which==21:  # Case "S23"
        result =  0.5*getZ(17,x,y)+0.5*getZ(19,x,y), 0.5*getZ(16,x,y)-0.5*getZ(18,x,y), 'S23'  # n/a
    elif which==22:  # Case "S24"
        result =  0.5*getZ(16,x,y)+0.5*getZ(18,x,y), -0.5*getZ(17,x,y)+0.5*getZ(19,x,y), 'S24' # n/a
    elif which==23:  # Case "S25" 
        result =  0.5*getZ(19,x,y)+0.5*getZ(21,x,y), 0.5*getZ(18,x,y)-0.5*getZ(20,x,y), 'S25'  # n/a
    elif which==24:  # Case "S26"
        result =  0.5*getZ(18,x,y)+0.5*getZ(20,x,y), -0.5*getZ(19,x,y)+0.5*getZ(21,x,y), 'S26'  # n/a
    elif which==25:  # Case "S27"
        result =  rh*getZ(21,x,y), rh*getZ(20,x,y), 'S27'                 # n/a
    elif which==26:  # Case "S28"
        result =  rh*getZ(20,x,y), -rh*getZ(21,x,y), 'S28'                # n/a
    elif which==27:  #  case "T4"
        result =  rh*getZ(3,x,y), -rh*getZ(2,x,y), 'T4 = Roll'
    elif which==28:  # case "T7" 
        result =  rh*getZ(4,x,y)-0.5*getZ(6,x,y), -0.5*getZ(5,x,y), 'T7 = LeftrightCurl'
    elif which==29:  # case "T8"
        result =  0.5*getZ(5,x,y), -rh*getZ(4,x,y)-0.5*getZ(6,x,y), 'T8 = UpDownCurl'
    elif which==30:  # case "T11" 
        result =  rh*getZ(7,x,y), -rh*getZ(8,x,y), 'T11 = Roll-AntiRoll'
    elif which==31:  # case "T12"
        result =  -0.5*getZ(7,x,y)+0.5*getZ(9,x,y), -0.5*getZ(8,x,y)-0.5*getZ(10,x,y), 'T12' # n/a
    elif which==32:  # case "T13"
        result =  0.5*getZ(8,x,y)-0.5*getZ(10,x,y), -0.5*getZ(7,x,y)-0.5*getZ(9,x,y), 'T13' # n/a
    else: 
        print("ZhaoBurgeTerm() is exitting because unsatisfied which = ", which) 
        quit()
    return result   

#---------end of Zhao-Burge evaluator with Sanity correction------------











#---ray generators-------------

def buildCircles(radius, ncircles):
    # circular grid within a circle
    pairs = list()
    pairs.append([0.,0.])
    count = 1
    for icirc in range(1, ncircles+1):
        daz = 60./icirc
        off = 0.0 if icirc%2==0 else daz/2.0
        r = icirc * radius/ncircles
        for jaz in range(0, 6*icirc):
            a = np.deg2rad(off + jaz*daz)
            x = r * np.cos(a)
            y = r * np.sin(a)
            pairs.append([x,y])
            count += 1
    print("buildCircles() has created npairs = ", count)
    return np.array(pairs)
    
def buildSquares(radius, npoints):
    # square grid, npoints x npoints
    pairs = list()
    step = 2.0*radius/(npoints-1)   
    count = 0
    for row in range(0, npoints):
        y = step*row - radius
        for col in range(0, npoints):
            x = step*col - radius
            r = np.sqrt(x**2 + y**2)
            if r <= radius:
                pairs.append([x,y])       
                count += 1
    print("buildSquares() has created npairs = ", count)
    return np.array(pairs)


  

     
#---graphics tadpole support------------

"""
def doOneTadpoleMap():
    #  Note: arrays U,V are precomputed SkyGrid positions: arguments to getZhaoBurgeTerm().
    #  Note: arrays X,Y are ready to receive deviated rays.  
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))  # left,bottom,width,height
    coef = 0.15
    which= 9     # THIS SPECIFIES WHICH ONE MAP IS WANTED
    for i in range(0, len(SkyGrid)):
        u = SkyGrid[i, 0]
        v = SkyGrid[i, 1]
        x, y, label = getTerm(u, v, which)
        X[i] = coef * x
        Y[i] = coef * y
    magnif = 1.0
    clearParms()
    h = U + magnif*X
    v = V + magnif*Y
    ax.plot(U, V, 'ko', markersize=1)       # dots
    ax.plot([U,h],[V,v],'r', linewidth=0.8) # lines
    ax.axis('equal')                        # unneeded with xlim, ylim?
    ax.set_xlim(-1.1, +1.1)
    ax.set_ylim(-1.1, +1.1)

    title = str(which)+' '+label
    color = 'k'
    ax.set_title(title, fontsize=8)
    ax.spines['bottom'].set_color(color)
    ax.spines['top'].set_color(color) 
    ax.spines['right'].set_color(color)
    ax.spines['left'].set_color(color)
    # ax.set_xticks([])   # no ticks anywhere
    # ax.set_yticks([])

    fig.savefig('OneTadpoleMap.eps', format='eps', dpi=600)
    plt.show()
"""

#---------------0----1----2----3----4----5----6----7----8----9---10---
spinecolors = ['k', 'k', 'b', 'b', 'b', 'b', 'b', 'r', 'r', 'r', 'r']



def doEightPortrait():
    #--3x3 plot grid showing 8 tadpole distortion maps
    #  Note: arrays U,V are precomputed SkyGrid positions: arguments to getZhaoBurgeTerm().
    #  Note: arrays X,Y are ready to receive deviated rays.  
    fig, axarr = plt.subplots(3, 3, figsize=(5,5))  # 3 rows, 3 plots each row
    coef = 0.15                 # keeps the tadpoles reasonable size
    for iplot in range(8):     # 8 maps, 0...7
        col = iplot % 3         # col = 0, 1, 2.
        row = int(iplot/3)      # row = 0, 1, 2.
        for i in range(0, len(SkyGrid)):
            u = SkyGrid[i, 0]
            v = SkyGrid[i, 1]
            x, y, label = getZhaoBurgeTerm(iplot, u, v)    # Now area normalized     
            X[i] = coef * x
            Y[i] = coef * y
            
        magnif = 1               # or could control tadpole lengths here. 
        h = U + magnif*X
        v = V + magnif*Y
        axarr[row, col].plot(U, V, 'ko', markersize=1)       # dots
        axarr[row, col].plot([U,h],[V,v],'r', linewidth=1.0) # lines
        # axarr[row, col].axis('equal')
        axarr[row, col].set_xlim(-1.2, +1.2)
        axarr[row, col].set_ylim(-1.2, +1.2)
        # axarr[row, col].axis('equal')
        title = str(iplot) +'    '+ label
        color = spinecolors[iplot]
        axarr[row, col].set_title(title, fontsize=8)
        axarr[row, col].spines['bottom'].set_color(color)
        axarr[row, col].spines['top'].set_color(color) 
        axarr[row, col].spines['right'].set_color(color)
        axarr[row, col].spines['left'].set_color(color)
        axarr[row, col].set_xticks([])   # no ticks anywhere
        axarr[row, col].set_yticks([])

    # blank plot box 8
    axarr[2,2].axis('off')
    
    fig.text(0.68, 0.24, 'Black: Laplacian')
    fig.text(0.68, 0.20, 'Blue: zero curl')
    fig.text(0.68, 0.16, 'Red: zero div')

    fig.savefig('EightPortrait.eps', format='eps', dpi=1200)
    plt.show()



def doElevenPortrait():
    #--4x3 plot grid showing 8 tadpole distortion maps
    #  Note: arrays U,V are precomputed SkyGrid positions: arguments to getZhaoBurgeTerm().
    #  Note: arrays X,Y are ready to receive deviated rays.  
    fig, axarr = plt.subplots(nrows=4, ncols=3, figsize=(6,8)) 
    coef = 0.15                 # keeps the tadpoles reasonable size
    iplot = 0
    for irow, ax_row in enumerate(axarr):
        print('Starting row = ', irow)
        for jcol, ax in enumerate(ax_row):
            if iplot<11: 
                print('Starting col = ', jcol)
                for i in range(len(SkyGrid)):  # 20 circles = 469 skygrid points
                    u = SkyGrid[i, 0]
                    v = SkyGrid[i, 1]
                    x, y, label = getZhaoBurgeTerm(iplot, u, v)    # Now area normalized     
                    X[i] = coef * x
                    Y[i] = coef * y
                magnif = 1               # or could control tadpole lengths here. 
                h = U + magnif*X
                v = V + magnif*Y
                ax.plot(U, V, 'ko', markersize=1)       # dots
                ax.plot([U,h],[V,v],'r', linewidth=1.0) # lines
                # axarr[row, col].axis('equal')
                ax.set_xlim(-1.2, +1.2)
                ax.set_ylim(-1.2, +1.2)
                # axarr[row, col].axis('equal')
                title = str(iplot) +'    '+ label
                color = spinecolors[iplot]
                ax.set_title(title, fontsize=8)
                ax.spines['bottom'].set_color(color)
                ax.spines['top'].set_color(color) 
                ax.spines['right'].set_color(color)
                ax.spines['left'].set_color(color)
                ax.set_xticks([])   # no ticks anywhere
                ax.set_yticks([])
                print('drawing box whose iplot = ', iplot)
                iplot += 1

    # blank plot box 12
    axarr[3,2].axis('off')
    
    fig.text(0.68, 0.24, 'Black: Laplacian')
    fig.text(0.68, 0.20, 'Blue: zero curl')
    fig.text(0.68, 0.16, 'Red: zero div')

    fig.savefig('ElevenPortrait.eps', format='eps', dpi=1200)
    plt.show()

#-------------------Main program-------------------------    
#-------set up sky grid for tadpoles-----------

SkyGrid = buildCircles(1.00, 8)
U = SkyGrid[:,0]
V = SkyGrid[:,1]
X = np.zeros_like(U)
Y = np.zeros_like(V)

          
#--Finally: do something!-----------

doElevenPortrait()


