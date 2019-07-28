#  RT171.py  MLL UCB SSL 25 July 2019
#
#  Similar to RT165a.py but now including atmospheric refraction = correction cases
#  Want to see >>combined<< field distortion.
#
#
#-----------------------ZA------ADC1------ADC2-----
# tasks = np.array([[   0.0,     0.0,     0.0],   
#                   [  45.0,   36.74,  -36.74], 
#                   [  60.0,    90.0,  -90.00]])                  
#
#
#  Just like RT164.py but here I have separated the ray direction finder matrix Mtel
#  from the ray in/put tasks, so as to simplify the work per ray.  This matrix need
#  be calculated & stored just once per CenterAZEL setting. 
#
#  Exploring field distortion as a function of ADC prism settings.
#  No SAG tables;  no Zernikes, no MultiSpot.
#
#  Includes a routine to pre-refract rays according to elevation angle within field
#
#  Signs: increasing ADC moves the image in +Y direction, up to 5mm.  See MAXADC.png
#  But:   increasing atmosphere moves image in -X direction, a few mm. 
#  Let's try to rotate the atmosphere to the -Y direction. 
#  The trouble must lie in rayUVW2azel() and its inverse azel2rayUVW()
#
#  Reminder: after altering OTILT OPITCH or OROLL, be sure to call setEulers().
#
#  Adopting double retro format  DESI-7.OPT  DESI-7.RAY  DESI-7.MED
#  To eliminate individual ray noise: NO SPIDERS, maybe NO OBSTRUCTIONS.
#  ADC roll is implemented using four coordinate breaks
#    Row 12 = CBin ahead of ADC1: roll angle +A1; for max A1=+90, roll=+90
#    Row 16 = CBin between ADCs:  roll angle -A1; so then roll -90
#    Row 17 = CBout between ADCs: roll angle -A2; for max A2=-90, roll=+90
#    Row 21 = CBout after ADC2:   roll angle +A2; so then roll=-90
#  With max ADC, these give pure Y deviations BlueY=+5.1mm, RedY=+4.9mm
#  See graphic RT163-ADC-CBs.pdf  29 June 2019
#
#  Positive roll CBin: makes ray antiroll in downstream lab system
#  Positive roll CBout: makes ray roll in downstream lab system. 
#
#  ADC1back positive pitch => thickest on +X edge
#  Positive CBin roll ahead of ADC1 moves +X rays to -Y, +Y rays to +X
#  Max ADC1 roll = +90deg puts +Y rays through thickest part of prism.
#  
#  ADC2front positive pitch => thickest on -X edge
#  Positive CBout ahead of ADC2 moves +X rays to +Y and +Y rays to -X
#  Again this puts +Y rays through thickest part of ADC2.
#
#  So: maxADC: ADC1 CBin rolls = +90, -90 and ADC2 CBout rolls = +90, -90
#  These combine to place super thick prism into +Y lab rays
#  CONFIRMED by the dispersed image BlueYfinal=+5.1mm, RedYfinal=+4.9mm
#  CONFIRMED by adding a full aperture plano prism, tilt=+0.030deg BK7
#  which makes its thick side -Y like low altitude and compensates spectrum.
#  See figure: Junk-7-FullAperturePrism.png
#  On the sky, the thick +Y part must lie towards the zenith
#  On the sky, the +Z axis lies towards the PM as seen by the corrector.
#  So, for right handed XxY=+Z, +X must point Eastward. 
#
#  Signs: incoming blue starlight appears elevated owing to atmospheric dispersion.
#  Elevated rays entering the telescope appear below average on the focal surface.
#  The ADC must elevate blue rays on the focal surface to compensate atmosphere. 
#
# OUTLINE:
#   Library imports:    line 90...
#   Macro definitions:   line 100...
#   List definitions:    lines 260...
#   Table parsing tools:  lines 320..
#   Math helpers:         lines 570...  
#   Surface generators:  lines 670...
#   Slope generators:   lines 730..
#   Interceptors:      lines 870...
#   Validator:           lines 1000...
#   Redirectors:            lines 1050...
#   Coordinate changers:      lines 1150...
#   Text output formatters:   lines 1250...
#   Input file unpackers:    lines 1350...
#   AutoRay support:           lines 1870...
#   Ray tracing:                  lines 2000...
#   Ray pattern generator:        lines 2200...
#   Output graphics:               lines 2230.....
#   Sky coordinates & refractors:      lines 2330...
#   Main program with INPUT FILE NAMES:  lines 2500...
#  
# M.Lampton UCB SSL 2018, 2019
#
#

#-------------------IMPORT ZONE-----------------------

import math
import numpy as np
import matplotlib.pyplot as plt       # for showing PSF output
import csv                            # for interpreting CSV files; 'rU' = universal EOL recognizer
import matplotlib.colors as mcolors   # for colorizing pixels with specified gamma

PROGNAME = 'RT171.py'

#---SET UP MACROS FOR INDEXING THE ARRAYS-------------
#---SET UP MACROS FOR INDEXING THE ARRAYS-------------
#---SET UP MACROS FOR INDEXING THE ARRAYS-------------


#----group 0 "OBASIC"-----------
OABSENT     = -1
OINDEX      = 0   # user-supplied refraction index
OACTIONTYPE = 1   # user-supplied surface  action type: lens, mirror, etc
OODIAM      = 2   # user-supplied limiting Outer Diameter
OIDIAM      = 3   # user-supplied inner diameter
OCURVE      = 4   # user-supplied curvature
OASPH       = 5   # user-supplied asphericity

#----group 1 "OXYZ"---------
OX          = 6   # user-supplied X coordinate of surface vertex, lab frame
OY          = 7   # user-supplied Y coordinate of surface vertex, lab frame
OZ          = 8   # user-supplied Z coordinate of surface vertex, lab frame

#----group 2 "OTPR" be sure to call setEulers()---------------
OTILT     = 9   # user-supplied tilt angle about lab X axis, degrees
OPITCH    = 10  # user-supplied pitch angle about tilted Y axis, degrees
OROLL     = 11  # user-supplied roll angle about pitched and tilted Z axis, degrees

#---group 3 "OPOLY"------------
OA1       = 12  # user-supplied polynomial coefficients
OA2       = 13 
OA3       = 14
OA4       = 15 
OA5       = 16
OA6       = 17
OA7       = 18
OA8       = 19
OA9       = 20
OA10      = 21 
OA11      = 22
OA12      = 23
OA13      = 24
OA14      = 25

OSAGFILE  = 26  # not an index into Oarray, but into OsagFileNames[]
#---Sag group indexes into Oarray---------
OSAGMULT  = 27  # user-supplied sag table height multiplier; can be zero.
OSAGSTEP  = 28  # user-supplied sag table grid spacing, mm
OSAGNSIDE = 29  # machine supplied from sagtable[0].shape
# Note: sagstep is not OODIAM/(NSIDE-1), it is user-supplied independent of OODIAM

#----Zernike group ------
# Numbering is R.N.Wilson 2000; and J.C.Wyant 
# but these are not well standardized.
# Born & Wolf (n,m) are radial and angular indices.
OZ0 = 30   # (0,0) piston
OZ1 = 31   # (1,1) x tilt
OZ2 = 32   # (1,1) y tilt
OZ3 = 33   # (2,0) defocus
OZ4 = 34   # (2,2) astig 3rd order 0deg
OZ5 = 35   # (2,2) astig 3rd order 45deg
OZ6 = 36   # (3,1) coma 3rd order, x
OZ7 = 37   # (3,1) coma 3rd order, y
OZ8 = 38   # (4,0) spherical 3rd order
OZ9 = 39   # (3,3) trefoil 5th 0deg
OZ10 = 40  # (3,3) trefoil 5th 30deg
OZ11 = 41  # (4,2) astig 5th 0deg
OZ12 = 42  # (4,2) astig 5th 45deg
OZ13 = 43  # (5,1) coma 5th 0deg
OZ14 = 44  # (5,1) coma 5th 90deg
OZ15 = 45  # (6,0) spherical 5th order
OZ16 = 46  # (4,4) tetrafoil 7th 0deg
OZ17 = 47  # (4,4) tetrafois 7th 22.5deg
OZ18 = 48  # (5,3) trefoil 7th 0deg
OZ19 = 49  # (5,3) trefoil 7th 30deg
OZ20 = 50  # (6,2) astig 7th 0deg
OZ21 = 51  # (6,2) astig 7th 45deg
OZ22 = 52  # (7,1) coma 7th 0deg
OZ23 = 53  # (7,1) coma 7th 90deg
OZ24 = 54  # (8,0) spherical 7th
OZ25 = 55  # (5,5) pentafoil 9th 0deg
OZ26 = 56  # (5,5) pentafoil 9th 18deg
OZ27 = 57  # (6,4) tetrafoil 9th 0deg
OZ28 = 58  # (6,4) tetrafoil 9th 22.5deg
OZ29 = 59  # (7,3) trefoil 9th 0deg
OZ30 = 60  # (7,3) trefoil 9th 30deg
OZ31 = 61  # (8,2) astig 9th 0deg
OZ32 = 62  # (8,2) astig 9th 45deg
OZ33 = 63  # (9,1) coma 9th 0deg
OZ34 = 64  # (9,1) coma 9th 90deg
OZ35 = 65  # (10,0) spherical 9th order

#---spider macro definitions------------
ONSPIDER   = 66  # number of equally spaced spider legs
OWSPIDER   = 67  # width of each spider leg
OMAXINPUT  = 68

#----group 4 "EULER" internal access only---------
OE11      = 71  # internal Euler matrix element
OE12      = 72  # internal Euler matrix element
OE13      = 73  # internal Euler matrix element
OE21      = 74  # internal Euler matrix element
OE22      = 75  # internal Euler matrix element
OE23      = 76  # internal Euler matrix element
OE31      = 77  # internal Euler matrix element
OE32      = 78  # internal Euler matrix element
OE33      = 79  # internal Euler matrix element
OFINAL    = 80

#---optics action possible code values---------

OLENSACTION    = 0
OMIRRORACTION  = 1
ORETROACTION   = 2
OIRISACTION    = 3
OCBINACTION    = 4
OCBOUTACTION   = 5
actionLookup   = ['OLENSACTION', 'OMIRRORACTION', 'ORETROACTION', 'OIRISACTION', 'OCBINACTION', 'OCBOUTACTION']



#---RAY FIELD HEADER IDENTIFIERS----------------
#--these apply at every surface: 0=Raystarts; F=final surface.
#--parser must supply the surface number.-------------

RX      = 0  #uppercase: lab frame
RY      = 1
RZ      = 2
RU      = 3
RV      = 4
RW      = 5
RWAVE   = 6
RXG     = 7   # goal field
RYG     = 8   # goal field

RFINALINPUT = 8

Rx      = 7  # lowercase = vertex frame
Ry      = 8
Rz      = 9
Ru      = 10
Rv      = 11
Rw      = 12
RFINAL  = 13

#--------ray failure classifiers--------------------------    

OK       = 0   # no failure
BACK     = 1   # intercept failure; only roots are d<=0, i.e. backwards
MISS     = 2   # intercept failure
ODIAM    = 3   # outer diameter validate failure
IDIAM    = 4   # inner diameter validate failure
SPI      = 5   # spider leg hit killed ray
TIR      = 6   # refraction failure
EDGE     = 7   # interpolation beyond safe boundary
UNKNOWN  = 8   # no valid action 

failures = ['OK', 'BACK', 'MISS', 'ODIAM', 'IDIAM', 'SPI', 'TIR', 'EDGE', 'UNKNOWN']







#------LIST AND INDEX DEFINITIONS USED BY THE MAIN PROGRAM-----------------
#------LIST AND INDEX DEFINITIONS USED BY THE MAIN PROGRAM-----------------
#------LIST AND INDEX DEFINITIONS USED BY THE MAIN PROGRAM-----------------

#----Helpers that just read these need not declare them global---
#----Helpers that create or modify these MUST declare them global--------
#----Here, the .CSV files have no ruler line; parsing is by field-----

#----Optics params and lists------------------------

Nsurfs          = 0                 # how many optical surfaces, modified by guide number.
Onfields        = 0                 # how many data fields in the table

Ostrings        = []                # list of strings starting with title string
Otable          = []                # list of lists of fields: 2D, starting with row=2
Oheaders        = []                # list of strings, one per field
OhasZern        = []                # one-based list of booleans
OhasPoly        = []                # one-based list of booleans
OsagFileNames   = []                # one-based list of sagfile names;   parsed from .OPT table
OhasSag         = []                # one-based list of booleans; computed here
OsagArrays      = []                # one-based list of sag arrays;  computed here.
OglassNames     = []                # one-based list of refractive indices and glass Names. 
OneedsMedia     = False             # if True, nonnumerical glass names will require a validated .MED table

Oarray          = np.empty([0,0])   # Working, possibly modified, optical specification [allSurfaces, allParms]

#----Ray arrays--------------------------

Nrays         = 0                    # how many ray starts, modified by guide number
Rnfields      = 0                    # how many data fields in the table
Rstrings      = []                   # list of strings, 1D, starting with title string
Rtable        = []                   # list of lists of fields: 2D, starting with row=2
Rheaders      = []                   # list of strings
RItoF         = []                   # lookup index-to-field list ONLY FOR RAY INPUT FIELDS
RFtoI         = []                   # lookup field-to-index list.
Rarray        = np.zeros([0,0,0])    # 3D numpy array[nlines,nsurfs,nparms]
Raystarts     = np.zeros([0,0])      # 2D numpy array[nrays,nparms]; parms=X,Y,Z,U,V,W,wavelength
RwaveNames    = []                   # list of wavelengths from "@wave" field

#----Media lists----------------------------

Mnglasses   = 0                      # how many glasses, modified by guide number
Mnfields    = 0                      # how many data fields in the table
Mstrings    = []                     # list of strings starting with title string
Mtable      = []                     # list of lists of fields: 2D, starting with row=2
Mheaders    = []                     # list of strings = wavelength Names
MwaveNames  = []                     # list of wavelength Names from headers > 0
MglassNames = []                     # list of glass Names from column zero
Marray      = np.empty([0,0])        # 2D numpy array

#----SAG TABLE SUPPORT---------------------

SsagSurfs   = []
Sarray      = np.zeros([0,0])        # 2D array of sag values










#------TABLE PARSING TOOLS------------------
#------TABLE PARSING TOOLS------------------
#------TABLE PARSING TOOLS------------------

def cuberoot(x):
    # allows negative arguments
    return math.copysign(math.pow(abs(x), 1./3.), x)

def makeTen(s):
    while len(s)<10:
        s += ' '
    return s

def getFloatValue(s):
    try:
        x = float(s)
    except ValueError:
        x = -0.0
    return x


def fixup(u,v,w):
    arg = 1 - u*u - v*v
    if arg < 0.:
        print('fixup() has u,v too large. Quitting.', u, v)
        quit()
    if w>= 0:
        return np.sqrt(arg)
    else:
        return -np.sqrt(arg)

def suckInt(s):
    # extracts a positive integer from a messy string
    slen = len(s)
    if slen<1:
        return 0
    numstring = ""
    i = 0
    while i<slen and not s[i].isdigit():
        i += 1
    while i<slen and s[i].isdigit():
        numstring += s[i]
        i += 1
    if len(numstring) ==0:
        return 0
    try:
        result = int(numstring)
        return result
    except ValueError:
        return 0
    
    
def getActionType(snippet):
    # this number can be put right into Oarray[OACTIONTYPE field]
    length = len(snippet)
    if length < 1:
        return OLENSACTION
    c0 = snippet[0].upper()
    c2 = ' '
    if length > 2:
        c2 = snippet[2].upper()
    if c0 == 'M':
        return OMIRRORACTION
    if c0 == 'R':
        return ORETROACTION
    if c0 == 'I':     # includes spider
        return OIRISACTION
    if c0 == 'C':    #coordinate break
        if c2 == 'I':
            return OCBINACTION
        if c2 == 'O':
            return OCBOUTACTION
    return OLENSACTION


def getOpticsAttribute(header):
    # From a column header, return its math array column index.
    # Surface numbers are row numbers, they do not come from headers.
    # Includes Zernike coefficients, SagMaps.
    guidenum = suckInt(header)
    header += "      "
    c0 = header[0]
    c1 = header[1]
    c2 = header[2]
    c3 = header[3]
    c0up = c0.upper()
    c1up = c1.upper()
    c2up = c2.upper()
    c3up = c3.upper()
    
    if c0=='D':
        return OODIAM
    elif c0=='d':
        return OIDIAM
    elif (c0up=='T' and c1up=='Y') or (c0up=='L' and c1up=='E') or (c0up=='M' and c1up=='I'):
        # rint 'getOpticsIndex() header and OACTIONTYPE:', header, OACTIONTYPE
        return OACTIONTYPE
    elif c0up=='I' or c0up=='G':
        return OINDEX
    elif c0up=='X':
        return OX
    elif c0up=='Y':
        return OY
    elif c0up=='Z' and not c1up=='E':
        return OZ
    elif c0up=='P':
        return OPITCH
    elif c0up=='T' and c1up=='I':
        return OTILT
    elif c0up=='R':
        return OROLL
    elif c0up=='C':
        return OCURVE
    elif c0up=='A' and c1up=='S':
        return OASPH
    elif c0up=='A' and guidenum>=0 and guidenum<15:
        sum = OA1 - 1 + guidenum
        # print 'getOpticsIndex() finds polynomial field whose index = ', sum
        return sum
    elif c0up=='Z' and c1up=='E' and guidenum<36:
        sum = OZ0 + guidenum
        # print 'getOpticsIndex() finds Zernike field whose index = ', sum
        return sum
    elif c0up=='S' and c3up=='F':   # sag header
        return OSAGFILE
    elif c0up=='S' and c3up=='M':
        return OSAGMULT
    elif c0up=='S' and c3up=='S':
        return OSAGSTEP
    elif c0up=='N' and c1up=='S':
        return ONSPIDER
    elif c0up=='W' and c1up=='S':
        return OWSPIDER
    else:
        return -1

    
def getRayStartAttribute(header):
    # Returns the input field parameter, or -1 if not an input field.
    # print 'getRayStart(header) has header = ', header
    if len(header) < 1:
        return -1
    header = header.upper()
    c0 = header[0]
    if c0=='@':
        return RWAVE
    c1 = ' '
    if len(header) > 1:
        c1 = header[1]
    # For AutoRay, X0... are ray inputs, XG... are ray goals   
    # but sometimes I want to accept any Xxxxx or Yxxxx as a goal. 
    # if len(header) < 2:
    #    return -1
    # if header[1] != '0':  
    #     return -1
    
    if c0=='X':
        if c1 == 'G':
            # print '......XG detected; returning RXG = ', RXG
            return RXG
        if c1 == '0':
            return RX
        return -1
    if c0=='Y':
        if c1 == 'G':
            return RYG
        if c1 == '0':
            return RY
        return -1
    if c0=='Z':
        if c1 == '0':
            return RZ
        return -1
    if c0=='U':
        if c1 == '0':
            return RU
        return -1
    if c0=='V':
        if c1 == '0':
            return RV
        return -1
    if c0=='W':
        if c1 == '0':
            return RW
        return -1
    return -1
    

def findGlassRow(glassname):
    # Search MglassNames to find a given glassname.
    # Return -1 if not found. 
    for kglass in range(1, len(MglassNames)):
        if glassname == MglassNames[kglass]:
            # print 'findGlassRow has glassname, kglass = ', glassname, kglass
            return kglass
    return -1
    
    
def findWaveColumn(wavename):
    # Search MwaveNames trying to locate a given wavename.
    # Skip column zero: it is the glass name header.
    # Return -1 if not found.
    for col in range(1, len(MwaveNames)):  # ignore column zero.
        if wavename == MwaveNames[col]:  
            # print 'findWaveColumn has wavename, column = ',wavename, col
            return col
    return -1
    

def findRefraction(iray, jsurf):
    glassname = OglassNames[jsurf]      # numbering 1...Nsurfs
    if len(glassname) < 1:
        return 1.0                      # assumes blank = air = 1.0    
    try:
        result = float(glassname)
        return result
    except ValueError:                  # need a Media table     
        if Mnglasses < 1:
            print('Media table is needed for glassname = ', glassname)
            quit()
        wavename = RwaveNames[iray]         # numbering 1...Nrays
        result = 1.0
        # print('Starting findRefraction() lookup for glassname = ', glassname)
        mediarow = findGlassRow(glassname)
        if mediarow < 0:
            print('findRefraction() is quitting since no GlassName = ', glassname)
            quit()
        mediacol = findWaveColumn(wavename)
        if mediacol < 1:
            print('findRefraction() is quitting since no WaveName = ', wavename)
            quit()
        result = float(Marray[mediarow][mediacol])
        # print('findRefraction() for mediarow, mediacol gets result = ', mediarow, mediacol, result)
    return result
    













#--------MATH HELPERS----------------
#--------MATH HELPERS----------------
#--------MATH HELPERS----------------

def deg(radians):
    return math.degrees(radians)
    
def rad(degrees):
    return math.radians(degrees)

def isMinusZero(x):
    return x==0. and np.signbit(x)==True


def getBothRoots(A, B, C):
    # Solves for the real roots of a quadratic function.
    # Returns 0, 1, or two roots as a tuple: "-0.0" means failed.
    # Method is Press et al 'Numerical Recipes' 2nd edition p.183
    if A==0.0 and B!= 0.0:
        return -C/B, -0.0
    if B==0.0 and A!=0.0:
        if C/A>0:
            return np.sqrt(C/A), -np.sqrt(C/A)
        else:
            return -0.0, -0.0
    if C==0 and A!= 0.0:
        return -B/A, -0.0
    D = B*B-4*A*C
    if D <= 0.0:
        return -0.0, -0.0
    Q = -0.5*(B + np.sign(B)*np.sqrt(D))
    return Q/A, C/Q

def setEulers():  # call this after any change in OTILT, OPITCH, or OROLL
    for j in range(1, Nsurfs+1):      
        ct = np.cos(np.radians(Oarray[j, OTILT]))
        st = np.sin(np.radians(Oarray[j, OTILT]))
        cp = np.cos(np.radians(Oarray[j, OPITCH])) 
        sp = np.sin(np.radians(Oarray[j, OPITCH])) 
        cr = np.cos(np.radians(Oarray[j, OROLL])) 
        sr = np.sin(np.radians(Oarray[j, OROLL]))  
        Oarray[j,OE11] = cr*cp;               # X <- x; M11
        Oarray[j,OE12] = -sr*cp;              # X <- y; M12
        Oarray[j,OE13] = sp;                  # X <- z; M13
        Oarray[j,OE21] = cr*sp*st + sr*ct;    # Y <- x; M21
        Oarray[j,OE22] = cr*ct - sr*sp*st;    # Y <- y; M22
        Oarray[j,OE23] = -cp*st;              # Y <- z; M23
        Oarray[j,OE31] = -cr*sp*ct + sr*st;   # Z <- x; M31
        Oarray[j,OE32] = sr*sp*ct + cr*st;    # Z <- y; M32
        Oarray[j,OE33] = cp*ct;               # Z <- z; M33   
 
def dotproduct(abc, xyz):
    # returns the dot product of two triplets
    return abc[0]*xyz[0] + abc[1]*xyz[1] + abc[2]*xyz[2]
    
    
def crossproduct(abc, xyz):
    # returns the cross product of two triplets
    product = np.zeros(3)
    product[0] = abc[1]*xyz[2] - abc[2]*xyz[1]
    product[1] = abc[2]*xyz[0] - abc[0]*xyz[2]
    product[2] = abc[0]*xyz[1] - abc[1]*xyz[0]
    return product    
    
def normalize(norm):
    # modifies given host array.
    len = np.sqrt(norm[0]**2 + norm[1]**2 + norm[2]**2)
    if len==0:
        print("cannot normalize a zero vector")
        return
    norm[0] /= len
    norm[1] /= len
    norm[2] /= len

def testUVW(iray, jsurf):
    err = Rarray[iray, jsurf, RU]**2  \
        + Rarray[iray, jsurf, RV]**2  \
        + Rarray[iray, jsurf, RW]**2  \
        - 1.0
    if math.fabs(err) > 1E-14:
        print('UVW normalization error at iray, surf = ', iray, jsurf, err)

def isNegZero(x):
    return x==0. and np.signbit(x)==True

    













#-----SURFACE GENERATORS--------------------
#-----SURFACE GENERATORS--------------------
#-----SURFACE GENERATORS--------------------


def getZtotal(iray, jsurf, d):
    # "d" is a positive trial distance along current ray
    z = getZconic(iray, jsurf, d)
    if OhasPoly[jsurf]:
        z += getZpoly(iray, jsurf, d)
    if OhasZern[jsurf]:
        z += getZzern(iray, jsurf, d)
    if OhasSag[jsurf]:
        z += getZsag(iray, jsurf, d)
    return z   

    
def getZconic(iray, jsurf, d):
    #  coordinates here are local "vertex frame" values.
    # "d" is a positive trial distance along current ray being tested here.
    x = Rarray[iray, jsurf, Rx] + d * Rarray[iray, jsurf, Ru]
    y = Rarray[iray, jsurf, Ry] + d * Rarray[iray, jsurf, Rv]
    r2 = x*x + y*y
    s = Oarray[jsurf, OASPH] + 1.0
    c = Oarray[jsurf, OCURVE]
    numer = c*r2
    arg = 1 - s*c*c*r2
    if arg < 0:
        print('negative argument found by getZconic(); flange case failure code -0.0')
        return -0.0, MISS  # failure code
    denom = 1 + np.sqrt(arg)
    zconic = numer/denom
    return zconic, OK


def getZpoly(iray, jsurf, d):
    #  coordinates here are local vertex frame values.
    x = Rarray[iray, jsurf, Rx] + d * Rarray[iray, jsurf, Ru]
    y = Rarray[iray, jsurf, Ry] + d * Rarray[iray, jsurf, Rv]
    r = np.sqrt(x*x + y*y)
    product = 1.0
    sum = 0.0
    for attrib in range(OA1, OA14+1):   # OA1=11 ... OA14=24
        product *= r
        sum += Oarray[jsurf, attrib] * product
    return sum, OK
    
    









#----SURFACE SLOPES AND NORMALS---------    
#----SURFACE SLOPES AND NORMALS---------    
#----SURFACE SLOPES AND NORMALS---------    
    
def getNormal(iray, jsurf):
    # There are two surface normals. Should not matter. I always use the one with Nz>0.
    gx, gy = gradTotal(iray, jsurf)
    normal = np.array([-gx, -gy, 1.0])
    normalize(normal)
    return normal

def gradTotal(iray, jsurf):
    gx, gy = gradConic(iray, jsurf)
    if OhasPoly[jsurf]:
        px, py = gradPoly(iray, jsurf)
        gx += px
        gy += py
    if OhasZern[jsurf]:
        zx, zy = gradZern(iray, jsurf)
        gx += zx
        gy += zy
    if OhasSag[jsurf]:
        sx, sy = gradSag(iray, jsurf)
        gx += sx
        gy += sy
    # print 'gradTotal() is returning gx, gy = ', gx, gy
    return gx, gy
    
def gradConic(iray, jsurf):
    s = Oarray[jsurf, OASPH] + 1
    c = Oarray[jsurf, OCURVE]
    if c==0:               # plano case
        return 0.0, 0.0
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    r2 = x*x + y*y
    arg = 1.0 - s*c**2*r2
    if arg <= 0:           # flange case
        # print '   gradConic() is returning the flange case.'
        return -0.0, -0.0
    coef = c/np.sqrt(arg)  # conic case
    gx = x*coef
    return x*coef, y*coef
    
def gradPoly(iray, jsurf):
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    r = np.sqrt(x*x + y*y)
    if r==0:
        return 0.0, 0.0
    product = 1.0
    dzdr = 0.0
    for index in range(OA1, OA14+1):   # OA1=11 ... OA14=24
        coef = 1 + index - OA1
        dzdr += coef * Oarray[jsurf, index] * product
        product *= r
    return (x/r)*dzdr, (y/r)*dzdr 
    
def gradZern(iray, jsurf):    
    Radius = 0.5*Oarray[jsurf, OODIAM]
    if Radius <= 0:
        print("Zernikes require a specified Diameter.  Returning zero.")
        return 0.0
    xnorm = Rarray[iray, jsurf, Rx]/Radius
    ynorm = Rarray[iray, jsurf, Ry]/Radius
    rnorm = math.sqrt(xnorm*xnorm + ynorm*ynorm)  
    if rnorm > 1.0:
        print("Zernike normalized ray lies outside unit radius circle. Returning zero.")
        return 0.0
    # r=0 has special handling within getZernDerivXY().
    gx = 0.
    gy = 0.
    for attrib in range (OZ1, OZ35+1):
        # Note: OZZ0 has no gradient; skip it. 
        if Oarray[jsurf, attrib] != 0.0:
            index = attrib - OZ0
            tx, ty = getZernDerivXY(index, xnorm, ynorm)
            gx += tx * Oarray[jsurf, attrib]/Radius
            gy += ty * Oarray[jsurf, attrib]/Radius
    return gx, gy
    
def gradSag(iray, jsurf):
    coef = Oarray[jsurf, OSAGMULT]
    nside = int(Oarray[jsurf, OSAGNSIDE])
    step = Oarray[jsurf, OSAGSTEP]
    diam = step*(nside-1)    

    radius = 0.5*diam
    if coef==0.:
        return 0., 0.
    coef = Oarray[jsurf, OSAGMULT]
    # print 'gradSag() sees coef = ', coef
    if radius <= 0:
        print("Sags require a specified Diameter.  Returning zero.")
        return 0.0, 0.0
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    # print 'gradSag() sees optical x, y = ', x, y
    r2 = x**2 + y**2
    if r2 > SAFETY * radius**2:
        print('gradSag() is out of safety range. Returning zero.')
        return 0.0, 0.0
    ibin, xfrac = getBinAndFrac(x, radius, step)
    jbin, yfrac = getBinAndFrac(y, radius, step)
    if ibin<1 or ibin>nside-3 or jbin<1 or jbin>nside-3:
       print("gradSag() finds ray ", iray, " outside bounds.  Returning zero.")
       return 0.0, 0.0
       
    #--create the 4x4 interpolation grid-----
    # print 'gradSag() is using jbin-1, jbin+3, ibin-1, ibin+3: ', jbin-1, jbin+3, ibin-1, ibin+3
    # print 'gradSag() is using xfrac, yfrac = ', xfrac, yfrac
    
    Sarray = OsagArrays[jsurf]
    # print 'sag array shape = ', Sarray.shape
    
    pinterp = Sarray[jbin-1:jbin+3, ibin-1:ibin+3]   # jy=downward, ix=rightward
    # print 'gradSag() has generated pinterp = \n', pinterp
    gx, gy = cubicDerivs2D(pinterp, xfrac, yfrac)
    # print 'gradSag() has received UNIT gx, gy from cubicDerivs2D():  ', gx, gy
    gx *= coef/step
    gy *= coef/step
    # print 'gradSag() is returning gx, gy =', gx, gy
    return gx, gy

















#----------INTERCEPTOR TOOLS------------------------
#----------INTERCEPTOR TOOLS------------------------
#----------INTERCEPTOR TOOLS------------------------

def intercept(iray, jsurf):
    # Works entirely in local frame Rx,Ry,Rz,Ru,Rv,Rw.
    # These numbers come entirely out of labToVx() routine: lines 1330-1366.
    # zlocal = Rarray[iray, jsurf, Rz]   # ??? Yikes! Rz is sometimes negative!  FIX THIS BLUNDER PLEASE
    # wlocal = Rarray[iray, jsurf, Rw]   # always negative, as planned.
    # print '   intercept start iray, jsurf, zLocal, wLocal={:4d}{:4d}{:16.6f}{:16.6f}'.format(iray, jsurf, zlocal, wlocal) 
    
    if OhasPoly[jsurf] or OhasZern[jsurf] or OhasSag[jsurf]:
        return higherIntercept(iray, jsurf)
    return conicIntercept(iray, jsurf)


def func(iray, jsurf, d):  
    # Function whose root is to be found using Newton's method.
    zc, code = getZconic(iray, jsurf, d)
    if code!=OK:
        return -0.0, code
    # print 'zc = ', zc 
    zp, code = getZpoly(iray, jsurf, d)
    if code!=OK:
        return -0.0, code
    # print 'zp = ', zp
    z0 = Rarray[iray, jsurf, Rz]
    w = Rarray[iray, jsurf, Rw]
    sum = zc + zp - (z0 + w*d)
    # print 'func() gets total sum = ', sum
    return sum, OK
    
    
def deriv(iray, jsurf, d):
    # Estimator of the derivative of func() for Newton's method.
    DELTA = 0.00001   # should be made adaptive
    fplus, status = func(iray, jsurf, d+DELTA)
    fminus, status = func(iray, jsurf, d-DELTA)
    return (fplus - fminus)/(2.*DELTA)


def higherIntercept(iray, jsurf):
    # First do a conic intercept to get close to the correct root. 
    # print '\nHigherIntercept starting its conic intercept...'
    status  = conicIntercept(iray, jsurf)
    if status !=OK:
       return status
    # Set up for using the Newton method.
    # print '\nHigherIntercept() is setting up for Newton rootfinder'
    d = 0
    niters = 0
    while True:   # Newton rootfinder
        niters += 1
        f, status = func(iray, jsurf, d)
        if status!=OK:
            return status
        slope = deriv(iray, jsurf, d)
        d -= f/slope
        # print '    Newton rootfinder: niters, d, f: ', niters, d, f
        if abs(f) < 1E-12 or niters > 8:
            break;
    Rarray[iray, jsurf, Rx] = Rarray[iray, jsurf, Rx] + d*Rarray[iray, jsurf, Ru]
    Rarray[iray, jsurf, Ry] = Rarray[iray, jsurf, Ry] + d*Rarray[iray, jsurf, Rv]
    Rarray[iray, jsurf, Rz] = Rarray[iray, jsurf, Rz] + d*Rarray[iray, jsurf, Rw]
    return status
    
    
def conicIntercept(iray, jsurf):   
    s = Oarray[jsurf, OASPH] + 1.0
    c = Oarray[jsurf, OCURVE]
    
    # Note: labtovx() will have already set the vertex-frame ray starts.
    
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    z = Rarray[iray, jsurf, Rz]
    u = Rarray[iray, jsurf, Ru]
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw]
    sos = u*u + v*v + w*w
    err = sos - 1.0
    if abs(err) > 1E-12:
        print('Yikes, faulty direction cosines received by conicIntercept at jsurf =', jsurf)
    # print '   conic intercept input  x,y,z,u,v,w: {:12.4f}{:12.4f}{:12.4f}{:12.6f}{:12.6f}{:12.6f}'.format(x,y,z,u,v,w)
    
    A = c*(u*u + v*v + s*w*w)
    B = 2*c*(x*u + y*v + s*z*w) - 2*w
    C = c*(x*x + y*y + s*z*z) - 2*z
    r1, r2 = getBothRoots(A, B, C)
    rLessPositive = min(r1, r2)
    rMorePositive  = max(r1, r2)
    if rMorePositive<=0. and rLessPositive<0.:
        return BACK
        
    d=-0.0
    sheetLessPositive = s*c*(z + w * rLessPositive)
    sheetMorePositive  = s*c*(z + w * rMorePositive)
    
    # Always try the shorter path first...
    if sheetLessPositive<1.0 and rLessPositive>0.0:  # if OK? this is our winner.
        d=rLessPositive
    elif sheetMorePositive<1.0 and rMorePositive>0.0:  # else? maybe this is our winner. 
        d=rMorePositive
    if d<=0.0:                             # Neither? well then... failed.
        # print 'intercept failed: d1, d2 = ', rLessPositive, rMorePositive
        return MISS
        
    # Now we have a good ray segment length "d"  -- SO PROPAGATE.
    # print 'interceptor has d = ', d
    Rarray[iray, jsurf, Rx] = Rarray[iray, jsurf, Rx] + d*u
    Rarray[iray, jsurf, Ry] = Rarray[iray, jsurf, Ry] + d*v
    Rarray[iray, jsurf, Rz] = Rarray[iray, jsurf, Rz] + d*w  
    return OK
    

    
    
    
    
    
    
    
    
    
    

    
#---------VALIDATOR with spider----------------------
#---------VALIDATOR with spider----------------------
#---------VALIDATOR with spider----------------------

def validate(iray, jsurf):
    # Is intercept within {Diameter, diameter} pupil annulus?
    # print 'Starting validation for jsurf = ', jsurf
    x = Rarray[iray, jsurf, Rx]       # Rx is index for local frame ray x
    y = Rarray[iray, jsurf, Ry]       # Ry is index for local frame ray y
    r = np.sqrt(x*x + y*y)            # ray's local radius off axis
    Diam = Oarray[jsurf, OODIAM]      # outer Diameter; -0.0 if absent
    if Diam > 0.0:
        if r > 0.5*Diam:
            # print("clobbered by OODIAM.")
            return ODIAM
    if r < 0.5*Oarray[jsurf, OIDIAM]:  # inner diameter; -0.0 if absent
        # print("clobbered by IDIAM")
        return IDIAM

    nlegs = int(Oarray[jsurf, ONSPIDER])
    halfwidth = 0.5*Oarray[jsurf, OWSPIDER]
    if nlegs < 1 or halfwidth <= 0.0:   # no spider declared
        return OK
        
    # Test for rectangular spider leg blockage. Don't rotate the spider: rotate the ray. 
    rolldeg = Oarray[jsurf, OROLL]                      # spider roll pattern, degrees CCW from +X
    raydeg = deg(np.arctan2(y,x))                       # degrees CCW from +X;  0/0 is OK here
    # print('Spider at iray, jsurf: nlegs, halfwidth, rolldeg, raydeg = {:6d}{:3d}{:3d}{:8.2f}{:8.2f}{:8.2f}'.format(iray,jsurf, nlegs, halfwidth, rolldeg, raydeg))
    for ileg in range(0, nlegs):
        diff = ileg*360.0/nlegs + rolldeg - raydeg      # degrees between ray and leg
        xtemp = math.fabs(r * math.sin(rad(diff)))
        ytemp = r * math.cos(rad(diff))
        # print('...spider ileg, rolldeg, raydeg, diffdeg, xtemp, ytemp = {:3d}{:8.3f}{:8.3f}{:8.3f}{:8.1f}{:8.1f}'.format(ileg, rolldeg, raydeg, diff, xtemp, ytemp))
        if xtemp < halfwidth and ytemp > 0.:
            # print('....Finding that xtemp<halfwidth and ytemp>0, so declaring ray failure of type SPI')
            return SPI
    return OK



















#-------REDIRECTORS-------------------
#-------REDIRECTORS-------------------
#-------REDIRECTORS-------------------


def redirect(iray, jsurf):
    # This is the switchyard to the detail redirectors..
    if jsurf == Nsurfs:
        return OK                 # no redirection at final surface.
        
    u = Rarray[iray, jsurf, Ru]   # input directions to be modified
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw]
    action = int(Oarray[jsurf, OACTIONTYPE])
    # print '================redirect starting with u,v,w,action = {:12.6f}{:12.6f}{:12.6f}{:4d}'.format(u,v,w,action)+"  "+actionLookup[action]
    if action == OIRISACTION:     # cannot fail.
        return OK

    if action == OMIRRORACTION:   # cannot fail.
        normal = getNormal(iray, jsurf)
        dotp = u*normal[0] + v*normal[1] + w*normal[2]
        u = Rarray[iray, jsurf, Ru] = Rarray[iray, jsurf, Ru] - 2*dotp*normal[0]
        v = Rarray[iray, jsurf, Rv] = Rarray[iray, jsurf, Rv] - 2*dotp*normal[1]
        w = Rarray[iray, jsurf, Rw] = Rarray[iray, jsurf, Rw] - 2*dotp*normal[2]
        return OK

    if action == OLENSACTION:   # can fail via TIR.    
        numer = findRefraction(iray, jsurf)
        if numer==0.0:
            numer = 1.0
        denom = findRefraction(iray, jsurf+1)
        if denom==0.0:
            denom = 1.0
        # print 'numer, denom = ', numer, denom
        mu = numer/denom
        normal = getNormal(iray, jsurf)
        kinput = np.array([u, v, w])
        kparallel = crossproduct(normal, crossproduct(kinput, normal))
        
        kparallel = mu * kparallel                         # vector equation 
        kparallelSQ = dotproduct(kparallel, kparallel)     # scalar equation 
        kperpSQ = 1.- kparallelSQ                          # Pythagoras
        if kperpSQ <= 0.0:
            return TIR
        kperpmagnitude = np.sqrt(kperpSQ)                  # scalar equation 
        perpsign = np.sign(dotproduct(normal, kinput))
        kperp = perpsign*kperpmagnitude * normal           
        kout = kparallel + kperp                           # vector equation 
        Rarray[iray, jsurf, Ru] = kout[0]
        Rarray[iray, jsurf, Rv] = kout[1]
        Rarray[iray, jsurf, Rw] = kout[2]
        return OK
        
    if action == ORETROACTION:    # cannot fail
        Rarray[iray, jsurf, Ru] *= -1.0
        Rarray[iray, jsurf, Rv] *= -1.0
        Rarray[iray, jsurf, Rw] *= -1.0
        return OK
        
    if action == OCBINACTION or action == OCBOUTACTION:
         return OK
         
    return UNKNOWN
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#---------COORDINATE CHANGERS-------------------
#---------COORDINATE CHANGERS-------------------
#---------COORDINATE CHANGERS-------------------

def labtovx(iray, jsurf):
    global Rarray
    # Assumes that raystarts have been loaded into Rarray[iray, 0].
    # Coordinate frame changer, moves from jsurf-1 LAB to current jsurf LOCAL.
    # Matrix OE converts local to lab coordinates; must transpose it here.
    # M.Lampton STELLAR SOFTWARE (C) 1989, 2003, 2017

    Xprev = Rarray[iray, jsurf-1, RX]
    Yprev = Rarray[iray, jsurf-1, RY]
    Zprev = Rarray[iray, jsurf-1, RZ]
    Uprev = Rarray[iray, jsurf-1, RU]
    Vprev = Rarray[iray, jsurf-1, RV]
    Wprev = Rarray[iray, jsurf-1, RW] # if forward: Rw is positive.  Reverse: Rw is negative.

    xLocal = Xprev - Oarray[jsurf, OX]
    yLocal = Yprev - Oarray[jsurf, OY]
    zLocal = Zprev - Oarray[jsurf, OZ]  # if forward: Rz is negative.  Reverse: Rz is positive.  
    
    x = Rarray[iray,jsurf,Rx] = xLocal*Oarray[jsurf,OE11] + yLocal*Oarray[jsurf,OE21] + zLocal*Oarray[jsurf,OE31]
    y = Rarray[iray,jsurf,Ry] = xLocal*Oarray[jsurf,OE12] + yLocal*Oarray[jsurf,OE22] + zLocal*Oarray[jsurf,OE32]
    z = Rarray[iray,jsurf,Rz] = xLocal*Oarray[jsurf,OE13] + yLocal*Oarray[jsurf,OE23] + zLocal*Oarray[jsurf,OE33]

    u = Rarray[iray,jsurf,Ru] = Uprev*Oarray[jsurf,OE11] + Vprev*Oarray[jsurf,OE21] + Wprev*Oarray[jsurf,OE31]
    v = Rarray[iray,jsurf,Rv] = Uprev*Oarray[jsurf,OE12] + Vprev*Oarray[jsurf,OE22] + Wprev*Oarray[jsurf,OE32]
    w = Rarray[iray,jsurf,Rw] = Uprev*Oarray[jsurf,OE13] + Vprev*Oarray[jsurf,OE23] + Wprev*Oarray[jsurf,OE33]

    return OK

def vxtovx(iray, jsurf):
    # used only by CBout coordinate break. No math; it just copies locals.
    Rarray[iray, jsurf, Rx] = Rarray[iray, jsurf-1, Rx]
    Rarray[iray, jsurf, Ry] = Rarray[iray, jsurf-1, Ry]
    Rarray[iray, jsurf, Rz] = Rarray[iray, jsurf-1, Rz]    
    Rarray[iray, jsurf, Ru] = Rarray[iray, jsurf-1, Ru]
    Rarray[iray, jsurf, Rv] = Rarray[iray, jsurf-1, Rv]
    Rarray[iray, jsurf, Rw] = Rarray[iray, jsurf-1, Rw]    

def vxtolab(iray, jsurf):
    # Coordinate frame changer at a single surface.
    # Here the Euler matrix is used directly, local to lab conversion. 
    # M.Lampton STELLAR SOFTWARE (C) 1989, 2003, 2017

    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    z = Rarray[iray, jsurf, Rz]    
    u = Rarray[iray, jsurf, Ru]
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw]
    
    Rarray[iray, jsurf, RU] = u*Oarray[jsurf,OE11] + v*Oarray[jsurf,OE12] + w*Oarray[jsurf,OE13]
    Rarray[iray, jsurf, RV] = u*Oarray[jsurf,OE21] + v*Oarray[jsurf,OE22] + w*Oarray[jsurf,OE23]
    Rarray[iray, jsurf, RW] = u*Oarray[jsurf,OE31] + v*Oarray[jsurf,OE32] + w*Oarray[jsurf,OE33]

    Rarray[iray, jsurf, RX] = x*Oarray[jsurf,OE11] + y*Oarray[jsurf,OE12] + z*Oarray[jsurf,OE13]
    Rarray[iray, jsurf, RY] = x*Oarray[jsurf,OE21] + y*Oarray[jsurf,OE22] + z*Oarray[jsurf,OE23]
    Rarray[iray, jsurf, RZ] = x*Oarray[jsurf,OE31] + y*Oarray[jsurf,OE32] + z*Oarray[jsurf,OE33]
    
    Rarray[iray, jsurf, RX] = Rarray[iray, jsurf, RX] + Oarray[jsurf, OX]
    Rarray[iray, jsurf, RY] = Rarray[iray, jsurf, RY] + Oarray[jsurf, OY]
    Rarray[iray, jsurf, RZ] = Rarray[iray, jsurf, RZ] + Oarray[jsurf, OZ]
    return OK



































#---------NUMERICAL TEXT DISPLAY TOOLS----------------
#---------NUMERICAL TEXT DISPLAY TOOLS----------------
#---------NUMERICAL TEXT DISPLAY TOOLS----------------


def showEuler(jsurf):
    eulerlist = [[Oarray[jsurf, OE11], Oarray[jsurf, OE12], Oarray[jsurf, OE13]],
                 [Oarray[jsurf, OE21], Oarray[jsurf, OE22], Oarray[jsurf, OE23]],
                 [Oarray[jsurf, OE31], Oarray[jsurf, OE32], Oarray[jsurf, OE33]]]
    euler = np.array(eulerlist)
    print(euler)
    
def displayInput(iray):
    X = Raystarts[iray, RX]
    Y = Raystarts[iray, RY]
    Z = Raystarts[iray, RZ]    
    U = Raystarts[iray, RU]
    V = Raystarts[iray, RV]
    W = Raystarts[iray, RW] 
    print('+++ Input:   iray,   XYZUVW: {:3d}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}'.format(iray,X,Y,Z,U,V,W))

def displayLocal(iray, jsurf):
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    z = Rarray[iray, jsurf, Rz]    
    u = Rarray[iray, jsurf, Ru]
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw] 
    print('*** Output: howfar,  xyzuvw:{:3d}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}'.format(jsurf,x,y,z,u,v,w))

def displayLabs(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    Z = Rarray[iray, jsurf, RZ]    
    U = Rarray[iray, jsurf, RU]
    V = Rarray[iray, jsurf, RV]
    W = Rarray[iray, jsurf, RW] 
    print('*** Output: howfar,  XYZUVW:{:3d}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}'.format(jsurf,X,Y,Z,U,V,W))

def displayLongOutput(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    Z = Rarray[iray, jsurf, RZ]
    print('*** Output: howfar, X,Y,Z:{:4d}{:24.12f}{:24.12f}{:24.12f}'.format(jsurf,X,Y,Z))

def displayXYUV(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    U = Rarray[iray, jsurf, RU]
    V = Rarray[iray, jsurf, RV]
    print(' X,Y,U,V{:22.14f}{:22.14f}{:22.14f}{:22.14f}'.format(X,Y,U,V))

def displayHXYZ(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    Z = Rarray[iray, jsurf, RZ]    
    print('*** Output: howfar,  XYZ:{:3d}{:18.12f}{:18.12f}{:20.12f}'.format(jsurf,X,Y,Z))

def displayXYZ(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    Z = Rarray[iray, jsurf, RZ]
    print(' X,Y,Z: {:24.16f}{:24.16f}{:24.16f}'.format(X,Y,Z))

def displayXY(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    print(' X, Y: {:24.16f}{:24.16f}'.format(X,Y))

def doMonsterListing():       # List user's ray table results
    prepTableRayGroup(False)  # True=wantListing
    runTableRayGroup(True)    # True=wantListing
    ngood = len(xpoints)
    print("\nNgoodRays = {:6d}".format(ngood) + '  ' + xRMSstring + '  ' + yRMSstring)
    print('Successive ray listing...')
    for iray in range(1, Nrays+1):
        X0 = Rarray[iray, 0, RX]
        U0 = Rarray[iray, 0, RU]
        X1 = Rarray[iray, 1, RX]
        U1 = Rarray[iray, 1, RU]
        X2 = Rarray[iray, 2, RX]
        U2 = Rarray[iray, 2, RU]
        Xf = Rarray[iray, Nsurfs, RX]
        Uf = Rarray[iray, Nsurfs, RU] 
        print('iray, X0, U0, X1, U1, X2, U2, Xf, Uf = {:4d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:21.15f}{:12.6f}'.format(iray, X0, U0, X1, U1, X2, U2, Xf, Uf))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#-------FILE READERS--------------
#-------FILE READERS--------------
#-------FILE READERS--------------

TITLEROW = 0   # top row of .CSV table
HEADERROW = 1  # next row of .CSV table

def isEmpty(anyList):
    if len(anyList) < 1:
        return True
    width = 0
    for i in range(0, len(anyList)):
        width = max(width, len(anyList[i]))
    return True if width==0 else False
    
    
    
def getColorFromNumber(wavename):
    # From an RwaveName string, get a ray color.
    # Python lacks switch() or case() ability.  Use "if" chain.
    wavel = getFloatValue(wavename)
    if wavel < 0.5:
        return 'b'
    if wavel < 0.7:
        return 'g'
    return 'r'



def getColorFromName(wavename):
    # from an RwaveName string, get a ray color.
    # Seven wavelengths:  i=365,  g=435,  F=486,   d=587,   C=656,   s=852,     t=1014
    # Seven colors:       b=blue, c=cyan, g=green, y=yellow, r=red,  m=magenta, k=black
    # Python lacks switch() or case() ability.  Use "if" chain.
    wavename = wavename.strip()    # remove all leading and trailing blanks
    if wavename=='i':  # 365 nm
        return 'b'     # blue 
    if wavename=='g':  # 435 nm
        return 'c'     # cyan
    if wavename=='F':  # 486 nm
        return 'g'     # green
    if wavename=='d':  # 587 nm
        return 'y'     # yellow
    if wavename=='C':  # 656 nm
        return 'r'     # red
    if wavename=='s':  # 852 nm
        return 'k'     # black 
    return 'm'         # magenta for "t" and all other inputs
    
    
    
def unpackCSV(fname):
    data = list()       # initially empty.
    print('\nunpackCSV() Trying: ', fname)
    try:
        data = list(csv.reader(open(fname)))  # 2D list of snippets
    except IOError:
        print("Could not open that file. Quitting this file.")
        return data     # empty.
    if len(data) < 3:
        print("Fewer than three CSV records are found. Quitting this file.")
        return data
    for irow in range(0, len(data)):    
        for jcol in range(0, len(data[irow])):
            data[irow][jcol] = data[irow][jcol].strip()  # unnecessary from Excel
    # print "Initial nrows = ", len(data)

    # delete all empty rows: crucial if empty rows are mixed into the .csv file
    initialLength = len(data)
    for irow in range(initialLength-1, 0, -1):  # don't eliminate title row even if empty
        if isEmpty(data[irow]):
            del data[irow]
    # print "After removing empties, nrows = ", len(data)
    if len(data) < 1:
        print("Nothing left to return. Quitting this file.")
        return data

    # get number of fields from longest row that holds any actual data
    nfields = 0
    for irow in range(0, len(data)):
       for jcol in range(0, len(data[irow])):
           if len(data[irow][jcol]) > 0:
               nfields = max(nfields, jcol+1)
    # print "Nfields = ", nfields
    if nfields < 1:
        print("Nfields is < 1.  Quitting this file.")
        return data

    # make all rows have nfields by appending empty fields where needed.
    for irow in range(len(data)):
        data[irow] = data[irow][:nfields]  # truncate beyond nfields
        while len(data[irow]) < nfields:   # append empty fields
            data[irow].append("")

    # print "Retaining the original title, nrows = ", len(data)
    # Cleanup is now complete.
    # data[0] is the original title row
    # data[1] is the header row
    # data[2].... are the actual data
    # failure is indicated by len(data)=0.
    return data
    
    
    
    
def crosshair(ix, iy, rad):  
    # units here are plot index integers: ix=right, iy=down
    plt.plot([ix-rad, ix+rad], [iy, iy], 'k-')
    plt.plot([ix, ix], [iy-rad, iy+rad], 'k-')
    








#==========FILE UNPACKERS===================
#==========FILE UNPACKERS===================
#==========FILE UNPACKERS===================


def getOpticsCSV(optname):
    # puts user CSV data into a global list "Odata" 
    global Onfields, Nsurfs, Oarray, OglassNames, OneedsMedia, OhasAnySags, Oheaders, Ojmirror, Ojfocal
    if len(optname) < 1:
        print("Optics table was not found, but is mandatory.  Quitting.")
        quit()
        
    Odata = unpackCSV(optname)
    # print 'Showing OpticsCSV data...'
    # for row in range(0, len(data)):
    #     print data[row]
    
    if len(Odata) < 1:
        print(optname, " has returned no data. Quitting.")
        quit()

    Onfields = len(Odata[HEADERROW])
    Nsurfs = len(Odata) - 2
    guideNumber = suckInt(Odata[TITLEROW][0])
    # print "Optics guideNumber = ", guideNumber
    if guideNumber>0:
        Nsurfs = min(Nsurfs, guideNumber)
    Ojfocal = Nsurfs
    print("  getOpticsCSV() is setting Ojfocal = Nsurfs = ", Ojfocal)
    Oheaders = Odata[HEADERROW]
    if Onfields<1 or Nsurfs < 1:
        print(optname, ' has no data available.  Quittimg.')
        quit()
    Oarray = np.zeros([Nsurfs+1, OFINAL])        # rebuild host Oarray
    Oarray.fill(-0.0)
    
    #---set up complete empty lists------------
    OneedsMedia = False
    del OglassNames[:]
    OglassNames.append("base=1")
    for k in range(1, Nsurfs+1):
        OglassNames.append("")

    #------recall definitions-------------

    # OsagFileNames   = []                # one-based list of sagfile names;    parsed from .OPT table
    # OsagMultipliers = []                # one-based list of sag multipliers;    parsed from .OPT table
    # OsagNside       = []                # one-based list of Nsides;           computed from shape of sagtable
    # OhasSag         = []                # one-based list of booleans;            computed here
    # OsagStepSize    = []                # one-based list of sag grid steps, mm;  computed here
    # OsagArrays      = []                # one-based list of sag arrays;          computed here.
    # remember: Nside is evaluated from parser; StepSize=OODIAM/(Nside-1).

    del OsagFileNames[:]           # empty the list
    OsagFileNames.append("base=1") # start the list
    for j in range(1, Nsurfs+1):   # build the list
        OsagFileNames.append("")   # will be replaced during parsing, below
    
    del OsagArrays[:]              # empty the list
    OsagArrays.append("base=1")    # start the list
    for j in range(1, Nsurfs+1):   # build the list
        OsagArrays.append(np.empty(shape=(0,0)))  # will be generated from filename
    
    del OhasSag[:]                 # empty the list
    OhasSag.append("base=1")       # start the list
    for j in range(1, Nsurfs+1):   # build the list
        OhasSag.append(False)      # will be generated from filename

    #-----set literal and numerical OpticsData() field by field--------
    for ifield in range(0, Onfields):
        header = Odata[HEADERROW][ifield]
        attrib = getOpticsAttribute(header)
        
        if attrib == OINDEX:                      # literal not numerical
            for jsurf in range(1, Nsurfs+1):      # jsurf=1, 2, ...Nsurfs
                snippet = Odata[jsurf+1][ifield]
                OglassNames[jsurf] = snippet 
                if len(snippet) > 0:
                    try:
                        x = float(snippet)        # if not numeric,
                    except:                       # will need a .MED lookup table.
                        OneedsMedia = True
                        
        elif attrib == OACTIONTYPE:
            for jsurf in range(1, Nsurfs+1):
                snippet = Odata[jsurf+1][ifield]
                iaction = getActionType(snippet)  # returns an action code number
                Oarray[jsurf, OACTIONTYPE] = iaction
                if iaction == OMIRRORACTION:
                    print('  getOpticsCSV() is setting Ojmirror = ', jsurf)
                    Ojmirror = jsurf
                
        elif attrib == OSAGFILE:                  # literal not numerical
            for jsurf in range(1, Nsurfs+1):      # jsurf=1, 2, ..Nsurfs
                snippet = Odata[jsurf+1][ifield]
                OsagFileNames[jsurf] = snippet    # can be null or empty
                if len(snippet) > 0:
                    print('  getOpticsCSV() finds sag file name = ', snippet)

        elif attrib>=0 and attrib<OMAXINPUT:      # numerical data fields including SAGMULT etc
            # if attrib in (26,27,28,29):           # sag attribute range
            #     print "++++++++Parsing attribute number ", attrib
            for jsurf in range(1, Nsurfs+1):      # jsurf=1, 2, ...Nsurfs
                snippet = Odata[jsurf+1][ifield]
                x = -0.0
                if len(snippet) > 0:
                    try:
                        x = float(snippet)
                    except ValueError:
                        x = -0.0
                Oarray[jsurf, attrib] = x 

    #---For each surface: any need polynomials? or Zernikes?---
    del OhasPoly[:]                           # empty the host list
    del OhasZern[:]
    OhasPoly.append("base=1")                 # 1-based to match jsurf
    OhasZern.append("base=1")                 # 1-based to match jsurf
    for jsurf in range(1, Nsurfs+1):          # jsurf = 1, 2, ...Nsurfs
        numPoly = 0                           # poly search
        for index in range(OA1, OA14+1):
            if Oarray[jsurf, index] != 0.0: 
                numPoly += 1   
        if numPoly>0:
            OhasPoly.append(True)             # 1-based list like jsurf
        else:
            OhasPoly.append(False)
        numZern = 0                           # Zern search
        for index in range(OZ1, OZ35+1):
            if Oarray[jsurf, index] != 0.0:
                numZern += 1
        if numZern >0:
            OhasZern.append(True)             # 1-based list like jsurf
        else:
            OhasZern.append(False)
            
    #---capture sag data if present----------
    # print '  getOpticsCSV() is detailing sag list...'
    for jsurf in range(1, Nsurfs+1):
        fname = OsagFileNames[jsurf]
        if len(fname) > 2:
            print("  At surface = "+str(jsurf)+ "  Found sagfile name = ", fname)
            nskip = 7   # ought to be generalized!
            saglist = np.loadtxt(fname, skiprows=nskip)
            print("  Parsing "  + fname + ' saglist.shape = ', saglist.shape)
            sagtable = np.array(saglist)
            print("  sagtable.shape = ", sagtable.shape)
            nside, mside = sagtable.shape
            Oarray[jsurf,OSAGNSIDE] = nside
            OsagArrays[jsurf] = sagtable
            OhasSag[jsurf] = True
            step = Oarray[jsurf,OSAGSTEP]
            diam = step * (nside-1)
            print('  Sag data for ', jsurf, ' is nside, step, diam : ', nside, step, diam)
            print('  --------DISMISS GRAPHIC TO PROCEED--------------')
            #---show the sag graphic-----
            fig = plt.figure(figsize = (6,5))         # inches square
            ax = fig.add_axes((0.08, 0.1, 0.8, 0.8))  # left,bottom,width,height
            plt.imshow(sagtable, cmap='jet')
            ax.set_xlim(0, nside-1)
            # ax.set_ylim(nside-1, 0)  # row number increases downward
            ax.set_ylim(0, nside-1)    # row number increases upward
            plt.colorbar(fraction=0.046, pad=0.04)
            plt.title(fname)
            if nside==820:    # ditch at edge of C2F
                crosshair(90, 615, 20)
            if nside==261:    # bird in C3F
                crosshair(115, 58, 10)
            fig.savefig(fname+'.eps', format='eps', dpi=600)
            plt.show()
            
    #----evaluate all the Euler matrices----
    setEulers()
    



def smartUVW(uvw):
    # given a triplet "uvw" this routine imposes |uvw|=1
    # But it is smart, maintaining u and v if w is unspecified.
    # fixes up host uvw triplet
    sum = uvw[0]**2 + uvw[1]**2
    if sum == 0.:
        if uvw[2] < 0.:
            uvw[2] = -1.
        else:
            uvw[2] = +1. 
        return uvw
    # from here onward, sum > 0.
    if isMinusZero(uvw[2]) or uvw[2] == 1.:   # want w>0
        if sum > 1:
            root = np.sqrt(sum)
            uvw[0] /= root
            uvw[1] /= root
            uvw[2] = 0.
            sum = 1.
        uvw[2] = np.sqrt(1.-sum)
        return uvw
    if uvw[2] == 0.:   # want w=0
        root = np.sqrt(sum)
        uvw[0] /= root
        uvw[1] /= root
        return uvw
    if uvw[2] == -1.:  # want w<0
        if sum > 1.:
            root = np.sqrt(sum)
            uvw[0] /= root
            uvw[1] /= root
            uvw[2] = 0.
            sum = 1.
        uvw[2] = -np.sqrt(1.-sum)
        return uvw
    # general case: full 3D renormalization
    sum += uvw[2]**2
    root = np.sqrt(sum)
    uvw[0] /= root
    uvw[1] /= root
    uvw[2] /= root
    # return uvw # unnecessary; use arg instead.
        
        
        
        
def getRaysCSV(rayname):
    global Rnfields, Nrays, Rarray, Raystarts, Rheaders, RwaveNames
    if len(rayname) < 1:
        print("  Ray table was not found, but is mandatory.  Quitting.")
        quit()
        
    data = unpackCSV(rayname)
    if len(data) < 1:
        print(rayname, " has returned no data. Quitting.")
        quit()
        
    Rnfields = len(data[HEADERROW])
    Rheaders = data[HEADERROW]
    Nrays = len(data) - 2 
    
    guideNumber = suckInt(data[TITLEROW][0])
    print("  Rays guideNumber = ", guideNumber)
    if guideNumber>0:
        Nrays = min(Nrays, guideNumber)

    if Rnfields<1 or Nrays < 1:
        print(rayname, ' has no data available.  Quitting.')
        quit()
        
    #---set Raystarts[iray,attrib] field by field-------
    # print 'Creating and Setting Raystarts[iray,attrib]'
    Raystarts = np.zeros([Nrays+1, RFINALINPUT+1])   # base=1, base=0
    Raystarts.fill(-0.0)
    del RwaveNames[:]
    del RFtoI[:]
    del RItoF[:]
    RwaveNames.append("base=1")  # numbering 1...Nrays inclusive
    for iray in range(1, Nrays+1):
        RwaveNames.append("")
    
    for ifield in range(0, Rnfields):
        header = data[HEADERROW][ifield]
        attrib = getRayStartAttribute(header)
        RFtoI.append(attrib)
        if attrib == RWAVE:
            for iray in range(1, Nrays+1):
                snippet = data[iray+1][ifield]
                RwaveNames[iray] = snippet  # replace the empty name with snippet
                
        if attrib>=0 and attrib<=RFINALINPUT:      #X0,Y0,Z0,U0,V0,W0,@wavel,XG,YG
            # found a ray start column in ray table.
            for iray in range(1, Nrays+1):
                snippet = data[iray+1][ifield]
                Raystarts[iray, attrib] = getFloatValue(snippet)
                
    #---evaluate RItoF[]-----------
    del RItoF[:]
    for attrib in range(RX, RFINALINPUT+1):
        try:
            field = RFtoI.index(attrib)
        except ValueError:
            field = -1
        RItoF.append(field)
        
    fixupRW()
    
    print('...getRaysCSV() fixup complete.')

    #---display what we have read and fixed up---------------
    for iray in range(1, Nrays+1):
        x = Raystarts[iray,RX]
        y = Raystarts[iray,RY]
        z = Raystarts[iray,RZ]
        u = Raystarts[iray,RU]
        v = Raystarts[iray,RV]
        w = Raystarts[iray,RW]
        # print('...getRaysCSV() Raystarts iray,x,y,z,u,v,w: {:4d}{:9.3f}{:9.3f}{:12.3f}{:10.6f}{:10.6f}{:10.6f}'.format(iray,x,y,z,u,v,w))


def fixupRW():
    global Raystarts
    for iray in range(1, Nrays+1): 
        uvw = np.zeros(3)
        uvw[0] = Raystarts[iray,RU]
        uvw[1] = Raystarts[iray,RV]
        uvw[2] = Raystarts[iray,RW]
        smartUVW(uvw) 
        Raystarts[iray,RU] = uvw[0]
        Raystarts[iray,RV] = uvw[1]
        Raystarts[iray,RW] = uvw[2]    
    
    
    

def getMediaCSV(medianame):
    global Mnfields, Mnglasses, Marray, MglassNames, MwaveNames
    if len(medianame) < 1:
        return
    data = unpackCSV(medianame)
    if len(data) < 1:
        print(medianame, " has returned no data. Quitting.")
        quit()
    Mnfields = len(data[HEADERROW])
    Mnglasses = len(data) -2
    guideNumber = suckInt(data[TITLEROW][0])
    if guideNumber > 0:
        Mnglasses = min(Mnglasses, guideNumber)
    if Mnfields<1 or Mnglasses<1:
        print(medianame, " has no data available.  Quitting.")
        quit()
    print('  Media guideNumber = ' + str(guideNumber) + ' Mnglasses = ' + str(Mnglasses)) 
    del MglassNames[:]  # empty the host list
    MglassNames.append("base=1")
    del MwaveNames[:]   # empty the host list
    Marray = np.zeros([Mnglasses+1, Mnfields])
    Marray.fill(-0.0)

    #----set Glass data() field by field---------
    for ifield in range(0, Mnfields): 
        header = data[HEADERROW][ifield]
        MwaveNames.append(header)   # do not use MwaveNames[0] it is not a wavename
        for kglass in range(1, Mnglasses+1):
            snippet = data[kglass+1][ifield]
            x = -0.0
            if len(snippet) > 0:
                try:
                    x = float(snippet)
                except ValueError:
                    x = -0.0
            Marray[kglass, ifield] = x
            
    #---gather glass Names from column zero----
    for kglass in range(1, Mnglasses+1):
        snippet = data[kglass+1][0]
        MglassNames.append(snippet)

        
def getMedia(myMedFileName):
    global Mnglasses, MwaveNames, MglassNames, OglassNames, RwaveNames    
    getMediaCSV(myMedFileName)
    if Mnglasses <1:
        print('  Media table is needed but no glasses are found. Quitting.')
        quit()
    # print("  MwaveNames: ", MwaveNames)
    # print("  MglassNames: ", MglassNames)
    for iglass in range(1, Nsurfs+1):   # skip base=0
        glassname = OglassNames[iglass]
        if len(glassname)>0 and not glassname in MglassNames:
            print('  Failed to find .OPT glass name ', glassname, ' in .MED table.  Quitting.')
            quit()
    for i in range(1, Nrays+1):   #skip base=0
        wavename = RwaveNames[i]
        if not wavename in MwaveNames:
            print('  Failed to find .RAY wave name ', wavename, ' in .MED fable.  Quitting.')
            quit()
                    
        
        
def getFraunhofer(s):         
    # Converts a Fraunhofer wavelength into a numerical wavelength in microns
    if s == 'i': return 0.365
    if s == 'g': return 0.436
    if s == 'F': return 0.486
    if s == 'd': return 0.588
    if s == 'C': return 0.656
    if s == 's': return 0.852
    if s == 't': return 1.014
    return -0.0
        
  
        
        
        
        
        
        
        
        
        
        
        
        
        
        


#----AutoRay Support-----------------------
#----AutoRay Support-----------------------
#----AutoRay Support-----------------------

def zeroOneRayGroup():
    print('Starting zeroOneRayGroup() chief rays...')
    for iray in range(1, Nrays+1): 
        x = Rarray[iray, 0, RX] = 0.
        y = Rarray[iray, 0, RY] = 0.
        z = Rarray[iray, 0, RZ] = 0.
        u = Rarray[iray, 0, RU] = 0.
        v = Rarray[iray, 0, RV] = 0.
        if Raystarts[iray,RW] >= 0.:
            Rarray[iray,0,RW] = 1.
        else:
            Rarray[iray,0,RW] = -1. 
        Rarray[iray,0,RWAVE]  = Raystarts[iray, RWAVE]           
        print("zeroOneRay() has completed iray = ", iray)

def rayfunc(iray, utrial, vtrial):
    # For a given chief ray "iray", tries out a {U0trial,V0trial} pair
    # If successful, it returns {Xf, Yf}, else {-0., -0.}
    # prepOneRay(iray, utrial, vtrial)   # installs desired utrial,vtrial
    howfar = runOneRay(iray)
    if howfar == Nsurfs:
        return Rarray[iray,Nsurfs,RX], Rarray[iray,Nsurfs,RY]
    else:
        return -0., -0.
    
def rayderivNumer(iray, u, v):
    DELTA = 1E-6
    xp, yp = rayfunc(iray, u+DELTA, v+DELTA)
    x, y = rayfunc(iray, u, v)
    sx = (xp-x)/DELTA
    sy = (yp-y)/DELTA
    # print "rayderiv() = {:24.12f}{:24.12f}".format(sx, sy)
    return sx, sy
   
def rayderiv(iray, u, v):
    return 15000.0, 15000.0
    
def doAutoRay():      # list chief ray starts {U0,V0}'s that achieve {Xgoal,Ygoal}'s
    MAXITERS = 20
    TOLFRAC = 1E-14
    zeroOneRayGroup()     # start all rays with {x0,y0,z0,u0,v0,w0} = {0.,0.,0.,0.,0.,-1.}
    runOneRayGroup(False) # verify that all the rays get safely to final surface.
    ngood = len(xpoints)
    if ngood < Nrays:
        print('doAutoRay(): Some rays have failed. Quitting. Ngood = ', ngood)
        quit()
    print('\ndoAutoRay() has built once and tested once its group of zero ray starts:')
    for iray in range(1, Nrays+1):   # base=1
        x0 = Rarray[iray, 0, RX]
        y0 = Rarray[iray, 0, RY]
        z0 = Rarray[iray, 0, RZ]
        u0 = Rarray[iray, 0, RU]
        v0 = Rarray[iray, 0, RV]
        w0 = Rarray[iray, 0, RW]
        print('{:4d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}'.format(iray,x0,y0,z0,u0,v0,w0))
        
    for iray in range(1, Nrays+1):   # base=1
        xgoal = Raystarts[iray, RXG]
        ygoal = Raystarts[iray, RYG]
        print('\ndoAutoRay() now starting iterations; iray, xgoal, ygoal = ', iray, xgoal, ygoal)
        # print 'Expect  intercept() to deliver POSITIVE zLocal and NEGATIVE wLocal......'
        # print 'This is because we are forward ray tracing from PM Vertex to each lens surface.'
        if isNegZero(xgoal) or isNegZero(ygoal):
            print('Missing ray start Xgoal or Ygoal data; quitting.')
            quit()
        niters = 0
        u = 0.   # initial guess
        v = 0.   # initial guess
        while True:
            niters += 1
            x, y = rayfunc(iray, u, v )
            if isNegZero(x):
                print('rayfunc() has failed to run iray, niter :', iray, niters)
            xerr = x - xgoal   # precision is limited here; using xfrac.
            yerr = y - ygoal   # precision is limited here; using yfrac.
            xfrac = xerr/(np.abs(x) + np.abs(xgoal) + 1)
            yfrac = yerr/(np.abs(y) + np.abs(ygoal) + 1)
            # print 'AutoRay(): xfrac, yfrac =          {:24.18f}{:24.18f}'.format(xfrac, yfrac)
            if np.abs(xfrac)<TOLFRAC and np.abs(yfrac)<TOLFRAC:
                print('AutoRay() tolerance OK: Exiting Newton loop, niters = ', niters)
                break
                
            slopex, slopey = rayderiv(iray, u, v)
            # print 'AutoRay(): slopex, slopey = ', slopex, slopey
            du = -xerr/slopex
            dv = -yerr/slopey
            # print 'AutoRay(): du, dv = ', du, dv
            u += du
            v += dv
            # print 'AutoRay(): niters, xerr, yerr: {:4d}{:24.18f}{:24.18f}'.format(niters, xerr, yerr)
            if niters > MAXITERS:
                print('Exceeded MAXITERS at iray = ', iray)
                break
        print('Done with iters; iray, niters, U, V =  {:4d}{:4d}{:18.12f}{:18.12f}'.format(iray, niters, u, v))


































#----RAY TRACE METHODS-----------------
#----RAY TRACE METHODS-----------------
#----RAY TRACE METHODS-----------------

def showLocals(iray, jsurf):
    # this diagnostic reveals what is going on in the local ray table.
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    z = Rarray[iray, jsurf, Rz]
    u = Rarray[iray, jsurf, Ru]
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw]
    print("    x,y,z,u,v,w = {:12.4f}{:12.4f}{:12.4f}{:12.6f}{:12.6f}{:12.6f}".format(x,y,z,u,v,w))

def modifyADC(adc1, adc2):  
    global Oarray
    # Modifies DESI-8.OPT Oarray table to specified ADC angles degrees
    # ADC1 is controlled by its CBin rolls (lines 12 & 16)
    # ADC2 is controlled by its CBout rolls (lines 17 and 21)
    # Call this only after .OPT table has been loaded & parsed.
    Oarray[12, OROLL] = adc1
    Oarray[16, OROLL] = -adc1
    Oarray[17, OROLL] = -adc2
    Oarray[21, OROLL] = adc2
    setEulers()


def runOneRay(iray):
    # uses Rarray[iray, 0] to get starting coordinates, so be sure to prep Rarray[] before calling this method.
    howfar = 0 
    # print('Starting runOneRay() with iray = ', iray)
    code = OK
    for jtarget in range(1, Nsurfs+1): 
        isCBout = OCBOUTACTION == int(Oarray[jtarget, OACTIONTYPE])
        if isCBout:
            vxtovx(iray, jtarget) 
            vxtolab(iray, jtarget)
            howfar = jtarget
            continue
        labtovx(iray, jtarget)
        code = intercept(iray, jtarget)
        # showLocals(iray, jtarget)
        if code == OK:
            code = validate(iray, jtarget)
        if code == OK:
            code = redirect(iray, jtarget)
        # showLocals(iray, jtarget)
        if code == OK:
            howfar = jtarget
        vxtolab(iray, jtarget) 
        if code != OK:
            break
        testUVW(iray, jtarget)
    # print('runOneRay() is returning howfar, code = ', howfar, failures[code])
    return howfar
    
#------set up global buffers for one plot's outputs-----------

xpoints = list() # raw labframe X values of individual rays in group; mm; base=0
ypoints = list() # raw labframe Y values of individual rays in group; mm, base=0
cpoints = list() # corresponding ray colors in this group; base=0
xAVE = 0.        # global labframe, mm
yAVE = 0.        # global labframe, mm
xRMS = 0.        # global mm
yRMS = 0.        # global mm
rRMS = 0.        # global mm
xRMC = 0.        # global mm
yRMC = 0.        # global mm

        
def prepMatrixRefractedRay(M, uv):
    global Rarray
    # Requires precalculated matrix  Mtel(fieldcenterAzEl)
    # Steers and refracts each table ray, writes Rarray[iray, 0, property]
    # za is zenith angle in degrees; u,v are ray direction cosines
    # Sign: downward is increasing za is increasing u
    # Sign: refraction is always upward, decreasing zenith angle slightly
    # Be sure to precompute Mtel matrix based on CenterAZEL
    wavel = 0.588  # better to use getFraunhofer(Rwavename[iray])

    for iray in range(1, Nrays+1): 
        # first get each ray's gross direction uvw
        Rarray[iray, 0, RWAVE]  = Raystarts[iray, RWAVE]           
        x = Rarray[iray, 0, RX] = Raystarts[iray, RX]
        y = Rarray[iray, 0, RY] = Raystarts[iray, RY]
        z = Rarray[iray, 0, RZ] = Raystarts[iray, RZ]
        u = Rarray[iray, 0, RU] = uv[0]                    # force u
        v = Rarray[iray, 0, RV] = uv[1]                    # force v
        w = Rarray[iray, 0, RW] = np.sqrt(1. - u*u - v*v)  # force w
        uvw = np.array([u,v,w])                            # telescope frame
        newUVW = rayUVW2bentUVW(M, uvw)                    # getAtmosphereDeg(elev)
        # now deposit the refracted rays into Rarray[][]
        nu = Rarray[iray, 0, RU] = newUVW[0]
        nv = Rarray[iray, 0, RV] = newUVW[1]
        nw = Rarray[iray, 0, RW] = np.sqrt(1. - nu*nu - nv*nv)

    
def runCentroid():  # get centroid of all good rays in current table; no output lists.
    ngood = 0
    xvals = list()
    yvals = list()
    for iray in range(1, Nrays+1):                     # base=1
        howfar  = runOneRay(iray)
        if howfar == Nsurfs:                           # this is a "good" ray
            ngood += 1
            xvals.append(Rarray[iray, Nsurfs, RX])
            yvals.append(Rarray[iray, Nsurfs, RY])
    if ngood == 0:
       return 0, 0., 0., 0., 0.
    xave = np.average(xvals)
    yave = np.average(yvals)
    xrms = np.std(xvals)
    yrms = np.std(yvals)
    return ngood, xave, yave, xrms, yrms
    
    
    
#------RAY PATTERN GENERATOR--------------
#------RAY PATTERN GENERATOR--------------
#------RAY PATTERN GENERATOR--------------
#------RAY PATTERN GENERATOR--------------

def isOdd(n):
    return n % 2 != 0
    
def getHexagonalNumbers(nrings, rmax):
    # produces an array of {xy} pairs of hexagonal numbers
    droll = 0.
    nrays = 1 + 3*nrings + 3*nrings**2
    dr = rmax/nrings
    xypairs = list()
    xypairs.append([0., 0.])
    for ring in range(1, nrings+1):
        r = ring * dr
        dangle = 60./ring
        offset = droll+dangle/2 if isOdd(ring) else droll
        for ray in range(6*ring):
            rad = np.radians(offset + dangle*ray)
            x = r*np.cos(rad)
            y = r*np.sin(rad)
            xypairs.append([x, y])
    return np.array(xypairs)






#-------GRAPHICS OUTPUT ROUTINES------------       
#-------GRAPHICS OUTPUT ROUTINES------------       
#-------GRAPHICS OUTPUT ROUTINES------------       
#-------GRAPHICS OUTPUT ROUTINES------------       
#-------GRAPHICS OUTPUT ROUTINES------------            
    
def plotTableRayGroup(title, jsurf):    # Plot one ray group
    global xpoints, ypoints, cpoints, xDEVs, yDEVs, xAVE, yAVE, xRMS, yRMS, rRMS
    xAVEstring  = 'xAVE={:12.6f}'.format(xAVE)
    yAVEstring  = 'yAVE={:12.6f}'.format(yAVE)
    ngood = len(xpoints)
    ngoodstring = 'Ngood=  '+str(ngood) 
    xRMSstring  = 'xRMS={:12.6f}'.format(xRMS)
    yRMSstring  = 'yRMS={:12.6f}'.format(yRMS)

    fig = plt.figure(figsize = (6.0, 6.0))           # inches; 6.5=max LaTeX.
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])      # about right.  75x75mm on screen.
    ax.tick_params(direction='in')                   # required for Python3
    for i in range(len(xpoints)):                    # xpoints are good rays, base=0
        color = "".join(cpoints[i] + 'o')            # cpoints has base=0
        x = xpoints[i]
        y = ypoints[i]
        # print 'plotTableRayGroup() has x, y, color = ', x, y, color
        ax.plot(x, y, color, markersize=2)
        
    xleft, xright = plt.xlim()
    xspan = xright - xleft
    xcenter = 0.5*(xleft + xright) 
    ybot, ytop = plt.ylim()
    yspan = ytop - ybot
    ycenter = 0.5*(ybot + ytop)
    span = max(xspan, yspan)
    # now double the spans, original centers, equal spans
    plt.xlim((xcenter-span, xcenter+span))
    plt.ylim((ycenter-span, ycenter+span))

    ax.set_xlabel('Xfp, mm, eastward')
    ax.set_ylabel('Yfp, mm, zenithward')
    ngood = len(xpoints)

    ax.text(0.02, 0.96, ngoodstring, transform=ax.transAxes, fontsize=8) 
    ax.text(0.02, 0.93, xAVEstring, transform=ax.transAxes, fontsize=8)  
    ax.text(0.02, 0.90, yAVEstring, transform=ax.transAxes, fontsize=8) 
    ax.text(0.02, 0.87, xRMSstring, transform=ax.transAxes, fontsize=8)  
    ax.text(0.02, 0.84, yRMSstring, transform=ax.transAxes, fontsize=8) 
    ax.text(0.02, 0.08, myOptFileName, transform=ax.transAxes, fontsize=8)  
    ax.text(0.02, 0.05, myRayFileName, transform=ax.transAxes, fontsize=8)      
    ax.text(0.02, 0.02, myMedFileName, transform=ax.transAxes, fontsize=8)    
    ax.text(0.86, 0.02, PROGNAME, transform=ax.transAxes, fontsize=8)  
    plt.title(title) 
    fig.savefig(title+'.eps', format='eps', dpi=600)
    plt.show()
    

def plotOneSpot(u, v):                          # Single spot diagram at angles u,v
    prepForcedRayGroup(u, v)                    # set up forced ray direction.
    ngoodrays = runTableRayGroup(Nsurfs)  
    print("ngoodrays = ", ngoodrays )
    title = 'U0 eastward={:8.4f},    V0 southward={:8.4f}'.format(u,v)
    plotTableRayGroup(title, Nsurfs, PLOTBOXMM)  # default is at final surface

def plotOnePupil(u, v):                         # draw all valid pupil rays at angle u,v
    prepForcedRayGroup(u, v)                    # set up forced ray direction.
    ngoodrays = runTableRayGroup(Ojmirror)      # request plotted rays at mirror surface
    print("ngoodrays = ", ngoodrays)
    title = 'U0 eastward={:8.4f},    V0 southward={:8.4f}'.format(u,v)
    plotTableRayGroup(title, Ojmirror, 4000)    # at mirror surface Ojmirror=7


def plotFullField(adc1, adc2):  
    fig = plt.figure(figsize = (6.0, 6.0))           # inches; 6.5=max LaTeX.
    ax = fig.add_axes([0.10, 0.10, 0.80, 0.80])      # about right.  75x75mm on screen.
    ax.tick_params(direction='in')                  
    
    nrings= 3
    rmax = 0.025  # radians in field of view
    uvpairs = getHexagonalNumbers(nrings, rmax)
    xvals = list()
    yvals = list()
    npts = len(uvpairs)    
    for i in range(npts):   
        u = uvpairs[i,0]
        v = uvpairs[i,1]
        prepForcedRayGroup(u,v)
        ngood, xave, yave, xrms, yrms = runCentroid()
        print('{:6d}{:6d}{:10.3f}{:10.3f}{:10.3f}{:10.3f}'.format(i,ngood,xave,yave,xrms,yrms))
        xvals.append(xave)
        yvals.append(yave)
    ax.scatter(xvals, yvals)
    plt.show()







#------sky coordinate converters and refractors------------
#------sky coordinate converters and refractors------------
#------sky coordinate converters and refractors------------
#------sky coordinate converters and refractors------------


DEG = np.degrees(1.)
RAD = np.radians(1.)

def printVector(name, vector):
    print(name, end='')
    for item in vector:
        print('{:12.6f}'.format(item), end='')
    print()   


def sind(a):
    return np.sin(RAD*a)

def cosd(a):
    return np.cos(RAD*a)

def tand(a):
    return np.tan(RAD*a)
    
def arcsind(x):
    return DEG*np.arcsin(x)

def arccosd(x):
    return DEG*np.arccos(x)

def arctan2d(y, x):
    return put360(DEG*np.arctan2(y, x))

def getUXYZ(lonlat):  # Convert spherical angles degrees into its xyz unit triplet
    return np.array([cosd(lonlat[0]) * cosd(lonlat[1]), 
                     sind(lonlat[0]) * cosd(lonlat[1]), 
                     sind(lonlat[1])])    
                     
def getLONLAT(xyz): # Convertunit sphere  xyz triplet into its spherical angles, degrees
    xyz = getNormalized(xyz)  # usually unnecessary
    return np.array([arctan2d(xyz[1],xyz[0]), arcsind(xyz[2])])
    
def getNorm(xyz):
    return np.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)
    
def getNormalized(xyz):
    return xyz/getNorm(xyz)  
    
def refX(xdeg):  # Rolls reference frame clockwise about +X
    c=cosd(xdeg); s=sind(xdeg)
    return np.array([[1,0,0], [0,c,+s], [0,-s,c]])

def refY(ydeg):  # rotates reference frame clockwise about +Y
    c=cosd(ydeg); s=sind(ydeg)
    return np.array([[c,0,-s], [0,1,0], [+s,0,c]])
    
def refZ(zdeg):  # Rotates reference frame clockwise about +Z
    c=cosd(zdeg); s=sind(zdeg)
    return np.array([[c,+s,0], [-s,c,0], [0,0,1]])

def put360(degrees): # Puts an angle into range 0 to 360.
    return np.fmod(720.+degrees, 360)

def getAtmosphereDeg(elev):
    # evaluates the increase in target elevation due to atmosphere
    # COEF here is evaluated in Sellmeier-MAXADC.xlsx spreadsheet
    COEF = 0.00022     # radians @45deg upwards, 800mb, 283K, 0.5878um    
    refraction = np.degrees(COEF*np.tan(np.radians(90.-elev)))
    return refraction  # in degrees; typically 0.01
    
                     
#-----Matrix converters AzEl to/from Telescope View----------

def getSky2TelMatrix(FieldCenterAZEL): 
    # Matrix generator: evaluate this ONCE per field center
    # Returns the 3x3 unitary transform from sky directions to telescope ray {u,v,w}
    # Its inverse, ray {uvw} to sky {xyz}, is its transpose
    Caz  = FieldCenterAZEL[0]
    Cel  = FieldCenterAZEL[1]
    Maz  = refZ(Caz)                    # brings field center to azimuth zero
    Mel  = np.dot(refY(90-Cel), Maz)    # brings Z axis to field center elevation
    Mtel = np.dot(refZ(-90), Mel)       # rolls eastward V to vertical, for ADCs
    # print(Mtel)
    return Mtel
    
def azel2rayUVW(Mtel, targetAZEL):
    # Matrix user
    # Convert sky coords in degrees to incoming telescope frame {u,v,w}
    targetXYZ = getUXYZ(targetAZEL)    
    return np.dot(Mtel, targetXYZ)
    
def rayUVW2azel(Mtel, rayUVW):
    # Matrix user
    # Converts telescope incoming ray direction  {u,v,w} to sky AzEl, degrees
    Msky = np.transpose(Mtel)
    xyz = np.dot(Msky, rayUVW)
    return getLONLAT(xyz)   
    
def rayUVW2bentUVW(Mtel, rayUVW):
    # Matrix user
    # Given an unrefracted telescope frame rayUVW, returns the refracted rayUVW
    azel = rayUVW2azel(Mtel, rayUVW)               # go to geographic coords
    lift = getAtmosphereDeg(azel[1])               # get the refraction lift angle
    newazel = np.array([azel[0], azel[1]+lift])    # get the new AzEl
    newUVW = azel2rayUVW(Mtel, newazel)            # return to telescope coordinates
    return newUVW 
    
    
    


#---Dispersion and ADC correction--------
#  From RT165 with ADCs={90,-90} "i" line vs "t" line, measured centroid image separation 200um.
#  From RT165 DESI-9 files determined FL = 13.93 meters on axis.
#  Infer max correction is 0.200/13.93 = 14.36 microradians.
#  so then CORRECTION(phi) = 14.36 * sin(phi)
#  Then, from Sellmeier-AtmoDispersion.xlsx look at KPNO=warm case:
#  p=790mbar, T=293K, totalDev(ZA=45deg)=212urad, diffDev=8.3urad.
#  making DISPERSION(za) = 8.3urad*tan(za)
#  set DISPERSION = CORRECTION to get sin(phi) = 0.577 tan(za)
#   ZA,deg    Phi,deg
#     0.000     0.000
#    30.000    19.459
#    45.000    35.240
#    60.000    88.004
def getPhiDeg(zaDeg):
    return np.degrees(np.arcsin(0.577*np.tan(np.radians(zaDeg))))
    
def getZAdeg(phiDeg):
    return np.degrees(np.arctan(np.sin(np.radians(phiDeg))/0.577))
    







print('\n--------------STARTING THE MAIN PROGRAM---------------------')

#--------set up the input filenames--------------
myOptFileName = 'DESI-9.OPT.CSV'
myRayFileName = 'DESI-9.RAY.CSV'  # 408 on-axis pupil rays; 0.588 microns monochromatic
myMedFileName = 'DESI-9.MED.CSV'

print("Loading files: " + myOptFileName + '  ' + myRayFileName + '  ' + myMedFileName)
getOpticsCSV(myOptFileName)
getRaysCSV(myRayFileName)

#----Set up global output Rarray now that we know Nrays and Nsurfs---------------
Rarray = np.zeros([Nrays+1, Nsurfs+1, RFINAL])
Rarray.fill(-0.0)

#---See if we need a .MED lookup table. If so, load and verify it-------
if OneedsMedia:
    getMedia(myMedFileName)
print('\nInput files have now been processed without error. \n')        

          
#---------------------ZA----- ADC1------ADC2-----
tasks = np.array([[   0.0,     0.0,     0.0],   
                  [  45.0,   36.74,  -36.74], 
                  [  60.0,    90.0,  -90.00]])             
                  
                  
tasks = np.array([[   0.0,     0.0,    0.00],   
                  [  60.0,     0.0,    0.00], 
                  [  60.0,    90.0,  -90.00]])         
                  
#---The .RAY table populates the on-axis pupil-----------------------
#---Here we generate the 217 sky directions out to 28 milliradians---     
#---One of these is on axis; the other 216 are off axis--------------             
NRINGS = 8
RMAX = 28
uvpairs = getHexagonalNumbers(NRINGS, 0.001*RMAX)


#------simpler: five targets per ZA-----
uvpairs = np.array([[ 0.00,  0.00],
                    [ 0.01,  0.00],
                    [ 0.00,  0.01],
                    [-0.01,  0.00],
                    [ 0.00, -0.01]])
                    

Mtel = np.eye(3)   # 3x3 matrix converting sky to telescope ray {uvw} directions

noext = PROGNAME.rsplit('.',1)[0]        
outfilename = noext + '_Nstars_5' + '.txt'

print('   ZA   ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms')

with open(outfilename, 'w') as outfile:
    for t in tasks:
        xpoints = list()
        ypoints = list()
        cpoints = list()
        modifyADC(t[1], t[2])
        CenterAZEL = np.array([0., 90.-t[0]])
        Mtel = getSky2TelMatrix(CenterAZEL)
        for uv in uvpairs:
            prepMatrixRefractedRay(Mtel, uv)
            ru = Rarray[1,0,RU]
            rv = Rarray[1,0,RV]
            ngood, xave, yave, xrms, yrms = runCentroid()
            result = '{:6.1f}{:6.1f}{:6.1f}{:10.6f}{:10.6f}{:6d}{:10.3f}{:10.3f}{:10.3f}{:10.3f}'  \
                .format(t[0],t[1],t[2],ru,rv, ngood, xave,yave,xrms,yrms)
            print(result)
            outfile.writelines(result + '\n')  # yes plural even for one line
            
            

""" 
RT171 calling prepMatrixRefractedRay()------
   ZA   ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms
   0.0   0.0   0.0  0.000000  0.000000   408     0.007     0.000     0.007     0.007
   0.0   0.0   0.0  0.009998 -0.000000   363   140.060    -0.000     0.013     0.009
   0.0   0.0   0.0  0.000000  0.009998   364     0.007   140.052     0.009     0.012
   0.0   0.0   0.0 -0.009998  0.000000   363  -140.044    -0.000     0.012     0.009
   0.0   0.0   0.0 -0.000000 -0.009998   364     0.007  -140.052     0.009     0.012
  60.0   0.0   0.0 -0.000000 -0.000381   408     0.007    -5.309     0.007     0.007
  60.0   0.0   0.0  0.009998 -0.000381   364   140.061    -5.338     0.013     0.009
  60.0   0.0   0.0  0.000000  0.009610   364     0.007   134.564     0.009     0.012
  60.0   0.0   0.0 -0.009998 -0.000381   364  -140.045    -5.338     0.012     0.009
  60.0   0.0   0.0 -0.000000 -0.010372   364     0.007  -145.359     0.009     0.012
  60.0  90.0 -90.0 -0.000000 -0.000381   408     0.000    -0.336     0.008     0.008
  60.0  90.0 -90.0  0.009998 -0.000381   364   140.052    -0.326     0.012     0.009
  60.0  90.0 -90.0  0.000000  0.009610   364     0.000   139.604     0.009     0.012
  60.0  90.0 -90.0 -0.009998 -0.000381   364  -140.052    -0.326     0.012     0.009
  60.0  90.0 -90.0 -0.000000 -0.010372   364     0.000  -140.314     0.010     0.013


RT171 calling on prepRefractedRayGroup()-----
   ZA   ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms
   0.0   0.0   0.0  0.000000  0.000000   408     0.007     0.000     0.007     0.007
   0.0   0.0   0.0  0.009998 -0.000000   363   140.060    -0.000     0.013     0.009
   0.0   0.0   0.0  0.000000  0.009998   364     0.007   140.052     0.009     0.012
   0.0   0.0   0.0 -0.009998  0.000000   363  -140.044    -0.000     0.012     0.009
   0.0   0.0   0.0 -0.000000 -0.009998   364     0.007  -140.052     0.009     0.012
  60.0   0.0   0.0 -0.000000 -0.000381   408     0.007    -5.308     0.007     0.007
  60.0   0.0   0.0  0.009998 -0.000381   364   140.061    -5.337     0.013     0.009
  60.0   0.0   0.0  0.000000  0.009610   364     0.007   134.565     0.009     0.012
  60.0   0.0   0.0 -0.009998 -0.000381   364  -140.045    -5.337     0.012     0.009
  60.0   0.0   0.0 -0.000000 -0.010372   364     0.007  -145.359     0.009     0.012
  60.0  90.0 -90.0 -0.000000 -0.000381   408     0.000    -0.335     0.008     0.008
  60.0  90.0 -90.0  0.009998 -0.000381   364   140.052    -0.325     0.012     0.009
  60.0  90.0 -90.0  0.000000  0.009610   364     0.000   139.605     0.009     0.012
  60.0  90.0 -90.0 -0.009998 -0.000381   364  -140.052    -0.325     0.012     0.009
  60.0  90.0 -90.0 -0.000000 -0.010372   364     0.000  -140.313     0.010     0.013


  RT165a.py:  result for NRINGS=8, RMAX=28, DESI-9.RAY.CSV using 408 rays monochromatic wavel="d"
   ZA   ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms
   0.0   0.0   0.0  0.000000  0.000000   408     0.007     0.000     0.007     0.007
   0.0   0.0   0.0  0.003030  0.001750   378    42.255    24.392     0.009     0.009
   0.0   0.0   0.0  0.000000  0.003499   378     0.007    48.784     0.008     0.009
   0.0   0.0   0.0 -0.003030  0.001750   378   -42.241    24.392     0.009     0.008
   0.0   0.0   0.0 -0.003030 -0.001750   378   -42.241   -24.392     0.009     0.008
   0.0   0.0   0.0 -0.000000 -0.003499   378     0.007   -48.784     0.008     0.009
   0.0   0.0   0.0  0.003030 -0.001750   378    42.255   -24.392     0.009     0.009
   .... total 217*3 = 651 records 

  RT165.py:  result for NRINGS=1, RMAX=25, DESI-9.RAY.CSV using 408 rays monochromatic wavel="d"
  
   ZA   ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms
   5.0   0.0   0.0 -0.000000 -0.000019   408     0.007    -0.268     0.007     0.007
   5.0   0.0   0.0  0.021646  0.012478   328   312.128   179.925     0.010     0.012
   5.0   0.0   0.0  0.000000  0.024975   328     0.008   360.112     0.013     0.008
   5.0   0.0   0.0 -0.021646  0.012478   328  -312.104   179.920     0.010     0.012
   5.0   0.0   0.0 -0.021646 -0.012516   328  -312.121  -180.485     0.010     0.012
   5.0   0.0   0.0 -0.000000 -0.025014   328     0.008  -360.706     0.013     0.008
   5.0   0.0   0.0  0.021646 -0.012516   328   312.144  -180.490     0.010     0.012
   
   5.0  36.7 -36.7 -0.000000 -0.000019   408     0.006     2.707     0.008     0.008
   5.0  36.7 -36.7  0.021646  0.012478   328   312.176   183.089     0.010     0.013
   5.0  36.7 -36.7  0.000000  0.024975   328     0.007   363.363     0.014     0.007
   5.0  36.7 -36.7 -0.021646  0.012478   328  -312.156   183.081     0.009     0.012
   5.0  36.7 -36.7 -0.021646 -0.012516   328  -312.074  -177.326     0.010     0.012
   5.0  36.7 -36.7 -0.000000 -0.025014   328     0.005  -357.460     0.011     0.009
   5.0  36.7 -36.7  0.021646 -0.012516   328   312.092  -177.326     0.010     0.011
   
  45.0   0.0   0.0 -0.000000 -0.000220   408     0.007    -3.065     0.007     0.007
  45.0   0.0   0.0  0.021646  0.012274   329   312.040   176.940     0.009     0.012
  45.0   0.0   0.0  0.000000  0.024769   328     0.008   356.932     0.013     0.008
  45.0   0.0   0.0 -0.021646  0.012274   329  -312.016   176.935     0.010     0.012
  45.0   0.0   0.0 -0.021646 -0.012715   329  -312.208  -183.392     0.010     0.012
  45.0   0.0   0.0 -0.000000 -0.025209   328     0.008  -363.724     0.013     0.009
  45.0   0.0   0.0  0.021646 -0.012715   329   312.231  -183.397     0.010     0.012
  
  45.0  36.7 -36.7 -0.000000 -0.000220   408     0.006    -0.090     0.008     0.008
  45.0  36.7 -36.7  0.021646  0.012274   329   312.088   180.102     0.010     0.013
  45.0  36.7 -36.7  0.000000  0.024769   328     0.007   360.178     0.014     0.007
  45.0  36.7 -36.7 -0.021646  0.012274   329  -312.068   180.094     0.009     0.012
  45.0  36.7 -36.7 -0.021646 -0.012715   329  -312.160  -180.231     0.010     0.012
  45.0  36.7 -36.7 -0.000000 -0.025209   328     0.005  -360.473     0.012     0.010
  45.0  36.7 -36.7  0.021646 -0.012715   329   312.178  -180.231     0.010     0.011
  
  RT164.py:  result for NRINGS=1, RMAX=25, DESI-9.RAY.CSV using 408 rays monochromatic wavel="d"  

   ZA   ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms
   
   5.0   0.0   0.0 -0.000000 -0.000019   408     0.007    -0.268     0.007     0.007
   5.0   0.0   0.0  0.021646  0.012478   328   312.128   179.925     0.010     0.012
   5.0   0.0   0.0  0.000000  0.024975   328     0.008   360.112     0.013     0.008
   5.0   0.0   0.0 -0.021646  0.012478   328  -312.104   179.920     0.010     0.012
   5.0   0.0   0.0 -0.021646 -0.012516   328  -312.121  -180.485     0.010     0.012
   5.0   0.0   0.0 -0.000000 -0.025014   328     0.008  -360.706     0.013     0.008
   5.0   0.0   0.0  0.021646 -0.012516   328   312.144  -180.490     0.010     0.012
   
   5.0  36.7 -36.7 -0.000000 -0.000019   408     0.006     2.707     0.008     0.008
   5.0  36.7 -36.7  0.021646  0.012478   328   312.176   183.089     0.010     0.013
   5.0  36.7 -36.7  0.000000  0.024975   328     0.007   363.363     0.014     0.007
   5.0  36.7 -36.7 -0.021646  0.012478   328  -312.156   183.081     0.009     0.012
   5.0  36.7 -36.7 -0.021646 -0.012516   328  -312.074  -177.326     0.010     0.012
   5.0  36.7 -36.7 -0.000000 -0.025014   328     0.005  -357.460     0.011     0.009
   5.0  36.7 -36.7  0.021646 -0.012516   328   312.092  -177.326     0.010     0.011
   
  45.0   0.0   0.0 -0.000000 -0.000220   408     0.007    -3.065     0.007     0.007
  45.0   0.0   0.0  0.021646  0.012274   329   312.040   176.940     0.009     0.012
  45.0   0.0   0.0  0.000000  0.024769   328     0.008   356.932     0.013     0.008
  45.0   0.0   0.0 -0.021646  0.012274   329  -312.016   176.935     0.010     0.012
  45.0   0.0   0.0 -0.021646 -0.012715   329  -312.208  -183.392     0.010     0.012
  45.0   0.0   0.0 -0.000000 -0.025209   328     0.008  -363.724     0.013     0.009
  45.0   0.0   0.0  0.021646 -0.012715   329   312.231  -183.397     0.010     0.012
  
  45.0  36.7 -36.7 -0.000000 -0.000220   408     0.006    -0.090     0.008     0.008
  45.0  36.7 -36.7  0.021646  0.012274   329   312.088   180.102     0.010     0.013
  45.0  36.7 -36.7  0.000000  0.024769   328     0.007   360.178     0.014     0.007
  45.0  36.7 -36.7 -0.021646  0.012274   329  -312.068   180.094     0.009     0.012
  45.0  36.7 -36.7 -0.021646 -0.012715   329  -312.160  -180.231     0.010     0.012
  45.0  36.7 -36.7 -0.000000 -0.025209   328     0.005  -360.473     0.012     0.010
  45.0  36.7 -36.7  0.021646 -0.012715   329   312.178  -180.231     0.010     0.011



"""