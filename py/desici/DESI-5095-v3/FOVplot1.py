
# FOVplot1.py 
# showing monochromatic differential refraction over DESI field, 45deg elevation
# Includes monochromatic atmospheric differential refraction in field coordinates
#     getDifferentialRefractionXY(centerAZEL, xy, centerRefraction)
# Includes monochromatic ray trace centroid locations from RT164 output tables
#     RT164-Nrings8_Rmax25.txt  408 sky locations, ZA=0 and ZA=45deg
#
# M.Lampton UCB SSL 2019, mlampton@berkeley.edu

# Three alternative ray direction descriptions:
# 1. Cartesian on unit sphere: {u,v,w} are three components Rx,Ry,Rz that
#    satisfy u^2+v^2+w^2=1, and usually +w is the incoming optical axis;
# 2. Tangent plane description {x,y,z} with x=u/w, y=v/w, z=1=optical axis.
# 3. Angles theta & phi: theta=off-axis angle, phi=roll from North. 
# Note that {u,v} are sub-linear in theta, varying like sin(theta)
# Note that {x,y} are super-linear in theta, varying like tan(theta)
#
# Conversion from sky AzEl angles to {u,v,w}:
#    Define CenterAzEl which is the AzEl of the FOV center axis on the sky;
#    Define TargetAzEl which is the AzEl of any target in the field; 
#    Define the telescope Z axis to be opposite of CenterAzEl because
#       Z axis is incoming direction from the given sky FOV center. 




import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'FOVplot1.py'

#---trig, vector, output helpers------------  

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

def getXYZ(lonlat):  
    # Convert spherical angles (degrees) into unit equatorial xyz
    return np.array([cosd(lonlat[0])*cosd(lonlat[1]), 
                     sind(lonlat[0])*cosd(lonlat[1]), 
                     sind(lonlat[1])])   
def getLONLAT(xyz): 
    # Convert equatorial xyz into its spherical angles
    xyz = getNormalized(xyz)  # often unnecessary
    return np.array([arctan2d(xyz[1],xyz[0]), arcsind(xyz[2])])
                         
def getPolarUVW(tf):  
    # Convert spherical angles {theta,phi} degrees into unit polar uvw
    return np.array([sind(tf[0])*cosd(tf[1]), 
                     sind(tf[0])*sind(tf[1]),
                     cosd(tf[0])])
        
def getPolarTF(uvw): 
    # Convert polar uvw into polar angles {theta,phi}   
    uvw = getNormalized(uvw)  # usually unnecessary
    return np.array([arccosd(uvw[2]), arctan2d(uvw[1],uvw[0])]) 
                     
def swapEW(lonlat):   # reverses east-west longitude sign
    return np.array([put360(-lonlat[0]), lonlat[1]])

def swapELZA(lonlat): # reverses zenith angle vs. elevation
    return np.array([lonlat[0], 90.-lonlat[1]])
    
def getNorm(xyz):
    return np.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)
    
def getNormalized(xyz):
    return xyz/getNorm(xyz)  
    
def vecX(xdeg):  # For positive xdeg=cwRoll: +y=>+z; +z=>-y.
    c=np.cosd(xdeg); s=np.sind(xdeg)
    return np.array([[1,0,0], [0,c,-s], [0,+s,c]])

def vecY(ydeg):  # For positive ydeg=-elev: +z=>+x; +x=>-z.
    # do not overlook this minus sign: positive ydeg pitches downward.
    c=np.cosd(ydeg); s=np.sind(ydeg)
    return np.array([[c,0,+s], [0,1,0], [-s,0,c]])
    
def vecZ(zdeg):  # For positive zdeg=+az: +x=>+y; +y=>-x.
    c=np.cosd(zdeg); s=np.sind(zdeg)
    return np.array([[c,-s,0], [+s,c,0], [0,0,1]])

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
    
def put180(degrees): # Put an angle into -180 to +180.
    d360 = put360(degrees) 
    if d360 <= 180.:
        return d360
    else:
        return d360-360.
       
#---------------------coordinate converters-----------------------------
#--Azimuths are eastward from south, to agree with equatorial {xyz}-----

def azel2hadec(azel, lat):
    # azimuth = eastward from south
    xyz = getXYZ(azel)
    colat = 90.-lat
    abc = np.dot(refY(-colat), xyz)  # pitch frame upward by colatitude
    return swapEW(getLONLAT(abc))    # ha=west, az=east

def eclip2radec(eclip):  # same epoch
    xyz = getXYZ(eclip)
    equ = np.dot(refX(-OBLIQ), xyz)  # roll frame counterclockwise by obliq
    return getLONLAT(equ)

def hadec2azel(hadec, lat):  
    azdec = swapEW(hadec)            # ha=west, az=eastward from south
    xyz = getXYZ(azdec)
    colat = 90.-lat
    abc = np.dot(refY(colat), xyz)  # pitch frame downward by colatitude
    return getLONLAT(abc)  

def hadec2radec(hadec, lst):         # all in degrees
    return np.array([put360(lst - hadec[0]), hadec[1]])
    
def radec2eclip(radec):  # same epoch
    xyz = getXYZ(radec)
    ecl = np.dot(refX(OBLIQ), xyz)  # roll frame clockwise
    return getLONLAT(ecl)
    
def radec2hadec(radec, lst):
    return np.array([put360(lst-radec[0]), radec[1]])
    
def radec2azel(radec, lst, lat):
    # North equatorial pole stays at {az,el} = {180,lat}
    # Equator passes through {az,el} = {0., colat}
    # South equatorial pole stays at {az,el} = {0,-lat}
    # azimuth = eastward from south
    hadec = radec2hadec(radec, lst)
    return hadec2azel(hadec, lat)

def azel2radec(azel, lst, lat):
    # zenith is always at {ra,dec} = {lst, lat}
    # azimuth = eastward from south
    hadec = azel2hadec(azel, lat)
    return hadec2radec(hadec, lst)
    




#-----coordinate converters AzEl to/from Telescope View----------

def azel2rayUVW(CenterAZEL, targetAZEL):
    # Convert sky coords in degrees to incoming Cartesian ray {u,v,w}
    Caz = CenterAZEL[0]                   # field center degrees eastward from south
    Cel = CenterAZEL[1]                   # field center degrees upward from horizon
    targetXYZ = getXYZ(targetAZEL)        # +X to S horizon, +Y to E horizon, +Z to zenith 
    uvw = np.dot(refZ(Caz), targetXYZ)    # reduces field center azimuth to zero
    uvw = np.dot(refY(90-Cel), uvw)       # rotates Z to field center
    uvw = np.dot(refZ(-90), uvw)          # rolls eastward V to vertical for ADC
    return uvw
    
def rayUVW2azel(CenterAZEL, rayUVW):
    # Converts telescope incoming ray direction to sky AzEl, degrees
    Caz = CenterAZEL[0]                   # field center degrees eastward from south
    Cel = CenterAZEL[1]                   # field center degrees upward from horizon
    uvw = np.dot(refZ(90), rayUVW)        # rolls V vertical back to eastward
    uvw = np.dot(refY(Cel-90), uvw)       # move the field center Z axis to the zenith
    uvw = np.dot(refZ(-Caz), uvw)         # rotate the X axis back to the south horizon
    targetAZEL = getLONLAT(uvw)           # convert to horizon angle system. 
    return targetAZEL


""" 
Test azel2rayUVW() and rayUVW2azel(), target West of center, positive U = eastward ray
center     =    63.000000   25.000000
target     =    60.000000   25.000000
teleUVW    =     0.047432   -0.000525    0.998874
targetAZEL =    60.000000   25.000000

Test azel2rayUVW() and rayUVW2azel(), target South of center, positive V = ascending ray
center     =    63.000000   25.000000
target     =    63.000000   22.000000
teleUVW    =     0.000000    0.052336    0.998630
targetAZEL =    63.000000   22.000000

Test azel2rayUVW() and rayUVW2azel(), target SouthWest: U, V, W all positive
center     =    63.000000   25.000000
target     =    60.000000   22.000000
teleUVW    =     0.048525    0.051799    0.997478
targetAZEL =    60.000000   22.000000

"""



#-------sky evaluators----------


def getViewXY(CenterAZEL, targetAZEL):
    # Convert sky coords into telescope view cube coordinates
    print('getViewXY is deprecated')
    Caz = CenterAZEL[0]   # degrees; eastward from south
    Cel = CenterAZEL[1]   # degrees
    targetXYZ = getXYZ(targetAZEL)
    Vxyz = np.dot(refY(-Cel), np.dot(refZ(Caz), targetXYZ))
    Vxyz =  DEG*Vxyz
    return np.array([Vxyz[1], Vxyz[2]])
    
def getHorizAZEL(CenterAZEL, ViewXY):
    # Convert telescope view cube coords into sky angles
    print('getHorizAZEL is deprecated')
    Caz = CenterAZEL[0]           # degrees, eastward from south
    Cel = CenterAZEL[1]           # degrees
    y = RAD*ViewXY[0]             # build a normalized triplet
    z = RAD*ViewXY[1]             # build a normalized triplet
    x = np.sqrt(1 - y**2 -z**2)   # build a normalized triplet
    Vxyz = np.array([x, y, z])
    Hxyz = np.dot(refZ(-Caz), np.dot(refY(Cel), Vxyz))
    return getLONLAT(Hxyz)
    
def getUV(CenterAZEL, targetAZEL):
    # Convert sky coords into telescope view angles {u,v}
    print('getUV is deprecated')
    Caz = CenterAZEL[0]     # degrees; eastward from south
    Cel = CenterAZEL[1]     # degrees
    targetXYZ = getXYZ(targetAZEL)
    Vxyz = np.dot(refY(-Cel), np.dot(refZ(Caz), targetXYZ))
    return getLONLAT(Vxyz)

def getAZEL(CenterAZEL, teleUV):
    # Convert telescope view angles {u,v} into sky AZel in degrees
    print('getAZEL is deprecated')
    Caz = CenterAZEL[0]     # degrees, eastward from south
    Cel = CenterAZEL[1]     # degrees
    Vxyz = getXYZ(teleUV)   # X=toward u=v=0; Y toward v=1; Z toward XxY 
    Hxyz = np.dot(refZ(-Caz), np.dot(refY(Cel), Vxyz))
    return getLONLAT(Hxyz)
"""
Test getViewXY() and getHorizAZEL()...deprecated...
center     =    63.000000   45.000000
target     =    66.000000   40.000000
teleXY     =     2.297083   -4.951123
targetAZEL =    66.000000   40.000000

Test getUV() and getAZEL()......also deprecated...
center     =    63.000000   45.000000
target     =    66.000000   40.000000
teleUV     =     2.306331   -4.957306
targetAZEL =    66.000000   40.000000

"""


#-------sky evaluators----------

def getRefraction(wavel, elev):
    # Sellmeier approximation, Peck & Reeder JOSA v.62 (1972)
    # 0.365000   45.000000  226.138006
    # 0.588000   45.000000  219.961045
    # 1.014000   45.000000  217.544915
    A = 57918.27
    B = 238.0185
    C = 1679.09
    D = 57.362
    P = 790  # mbar
    T = 283  # kelvin
    wm2 = wavel**(-2)
    sell = A/(B-wm2) + C/(D-wm2)   # microradians
    return (P*288)/(T*1013) * sell * tand(elev)
    
def getRefraction588(elev):
    # COEF here is evaluated in Sellmeier-MAXADC.xlsx spreadsheet
    COEF = 0.00022000     # radians towards zenith, 790mb, 283K, 0.588um    
    coel = 90.0-elev      # degrees
    refraction = COEF*tand(coel)
    return refraction     # in radians; typically 0.000222
    
def getDifferentialRefractionUV(centerAZEL, xy, centerRefraction):
    azelTail = getHorizAZEL(CenterAZEL, xy)
    dElev = getRefraction588(azelTail[1]) - centerRefraction
    azelHead = azelTail + np.array([0., dElev]) 
    return getViewXY(CenterAZEL, azelHead)  
        

#-----refraction divergence evaluator---------------

def getDivergence(theAZEL):
    DELTA = 1E-3  # equal steps in tangent plane
    dx = np.array([DELTA, 0.])
    dy = np.array([0., DELTA])
    eastAZEL = getHorizAZEL(theAZEL, dx)
    eastAZEL += np.array([0., getRefractionElev(eastAZEL[1])])
    eastXY   = getViewXY(theAZEL, eastAZEL) - dx
    westAZEL = getHorizAZEL(theAZEL, -dx)
    westAZEL += np.array([0., getRefractionElev(westAZEL[1])])
    westXY   = getViewXY(theAZEL, westAZEL) + dx
    upAZEL   = getHorizAZEL(theAZEL, dy)
    upAZEL   += np.array([0., getRefractionElev(upAZEL[1])])
    upXY     = getViewXY(theAZEL, upAZEL) - dy
    downAZEL = getHorizAZEL(theAZEL, -dy)
    downAZEL += np.array([0., getRefractionElev(downAZEL[1])])
    downXY   = getViewXY(theAZEL, downAZEL) + dy
    divx = (eastXY - westXY)/(2*DELTA)
    divy = (upXY - downXY)/(2*DELTA)
    #---now prepare the resulting record: elev, divx[0], divy[1]
    result = np.array([theAZEL[1], divx[0], divy[1]])
    printVector('elev, divx, divy = ', result)

def getRefractionCurl(theAZEL): 
    print()
    printVector('Starting getRefractionCurl() with AZEL = ', theAZEL)
    DELTA = 1E-3  # equal steps in tangent plane
    # sample refraction in four nearby locations
    dx = np.array([DELTA, 0.])
    dy = np.array([0., DELTA])
    eastAZEL = getHorizAZEL(theAZEL, dx)
    eastAZEL += np.array([0., getRefractionElev(eastAZEL[1])])
    eastXY   = getViewXY(theAZEL, eastAZEL) - dx
    eastCurl = eastXY[1]  # upward component
    westAZEL = getHorizAZEL(theAZEL, -dx)
    westAZEL += np.array([0., getRefractionElev(westAZEL[1])])
    westXY   = getViewXY(theAZEL, westAZEL) + dx
    westCurl = -westXY[1]  # downward component
    upAZEL   = getHorizAZEL(theAZEL, dy)
    upAZEL   += np.array([0., getRefractionElev(upAZEL[1])])
    upXY     = getViewXY(theAZEL, upAZEL) - dy
    upCurl   = -upXY[0]   # rightward component
    downAZEL = getHorizAZEL(theAZEL, -dy)
    downAZEL += np.array([0., getRefractionElev(downAZEL[1])])
    downXY   = getViewXY(theAZEL, downAZEL) + dy
    downCurl = downXY[0]  # leftward component 
    curl = eastCurl + westCurl + upCurl + downCurl
    print('Curl = {:12.6f}'.format(curl))

LAT = 32.0
LST = 0.

def curlField(anyAZEL):  # rotates 1 degree of LST about north equatorial pole
    DLST = 1.0
    hadec = azel2hadec(anyAZEL, LAT)
    radec = hadec2radec(hadec, LST+DLST)
    hadec = radec2hadec(radec, LST)
    azel = hadec2azel(hadec, LAT)
    return azel

def getRotationCurl(theAZEL):
    print()
    DELTA = 1E-3  # equal steps in tangent plane
    # sample distortion at four nearby locations
    dx = np.array([DELTA, 0.])
    dy = np.array([0., DELTA])
    eastAZEL = curlField(getHorizAZEL(theAZEL, dx))
    eastXY   = getViewXY(theAZEL, eastAZEL) - dx
    eastCurl = eastXY[1]  # upward component
    westAZEL = curlField(getHorizAZEL(theAZEL, -dx))
    westXY   = getViewXY(theAZEL, westAZEL) + dx
    westCurl = -westXY[1]  # downward component
    upAZEL   = curlField(getHorizAZEL(theAZEL, dy))
    upXY     = getViewXY(theAZEL, upAZEL) - dy
    upCurl   = -upXY[0]   # rightward component
    downAZEL = curlField(getHorizAZEL(theAZEL, -dy))
    downXY   = getViewXY(theAZEL, downAZEL) + dy
    downCurl = downXY[0]  # leftward component 
    curl = (eastCurl + westCurl + upCurl + downCurl)/4.0
    return curl

#-----------file unpacker---------------------

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
    # Make nums2D rectangular so asarray() returns an array not a list
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
    
    
    
#=================PROGRAM STARTS HERE=================


# first get the data table from RT164
# 10 fields:  0=ZA, 1=ADC1, 2=ADC2, 3=U0, 4=V0, 5=ngood, 6=Xave, 7=Yave, 8=Xrms, 9=Yrms

print()
filename = 'RT164_Nrings_12_Rmax_25'
filenametxt = filename + '.txt'
print(filenametxt)
myArray = getFileArray(filenametxt)
nrows = len(myArray)
ncols = len(myArray[0])
print('nrows, ncols = ' + str(nrows) + ' ' + str(ncols))
# for row in myArray:
#     print(row)
nADC = 2  # how many ADC settings are in the table
nZA  = 2  # how many ZAs are in the table
nGroups = nADC*nZA
nEach = nrows//nGroups
print('nEach = ' + str(nEach))
my0000 = myArray[0 : nEach]
my0037 = myArray[nEach : 2*nEach]
my4500 = myArray[2*nEach : 3*nEach]
my4537 = myArray[3*nEach : 4*nEach]

# build the graphic. units are mm in focal plane.
# idea is to compare the sky, adc, and sky+adc distortions 


fig = plt.figure(figsize=(7.,7.))  # inches
ax = fig.add_subplot(1,1,1)
ax.set_xlim(-450.,+450.)  # mm
ax.set_ylim(-450.,+450.)  # mm
ax.set_aspect(1.0)

ArrowMag = 100

plt.title('Monochromatic Distortion at 0.589um')

headstr = 'Heads: ZA=45, ADC=37'
tailstr = 'Tails: ZA=0, ADC=00'
for i in range(nEach):
    xTail = my0000[i,6]
    yTail = my0000[i,7]
    xHead = my4537[i,6]
    yHead = my4537[i,7]
    width = ArrowMag * (xHead - xTail)
    height = ArrowMag * (yHead - yTail)
    ax.arrow(xTail, yTail, width, height, head_width=5, head_length=8, fc='black')

ax.tick_params(direction='in')
ax.set_xlabel('field location x,mm')
ax.set_ylabel('field location y,mm')
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.95, headstr, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.91, tailstr, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.87, "ArrowMag = "+str(ArrowMag), transform=ax.transAxes, fontsize=10)
noext = PROGNAME.rsplit('.',1)[0]
plotname = noext + '_' + filename + '.png'             
plt.savefig(plotname) 
plt.show()           

