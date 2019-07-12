
# TangentPlane28.py  see also FOVplot1.py
# showing monochromatic differential refraction over DESI field, 45deg elevation
# M.Lampton UCB SSL 2019, mlampton@berkeley.edu
#
# Exploring improved coordinate changers sky <=> telescope {u,v,w}
# Three alternative ray direction descriptions:
# 1. Cartesian on unit sphere: {u,v,w} are three components Rx,Ry,Rz that
#    satisfy u^2+v^2+w^2=1, and usually +w is the incoming optical axis;
# 2. Tangent plane description {x,y,z} with x=u/w, y=v/w, z=1=optical axis.
# 3. Angles theta & phi: theta=off-axis angle, phi=roll from North. 
# Note that {u,v} are sub-linear in theta, varying like sin(theta)
# Note that {x,y} are super-linear in theta, varying like tan(theta)
#
# Conversion from sky AzEl angles to {u,v,w}:  azel2rayUVW(CenterAZEL, targetAZEL)
#    Define CenterAzEl which is the AzEl of the FOV center axis on the sky;
#    Define TargetAzEl which is the AzEl of any target in the field; 
#    Define the telescope Z axis to be opposite of CenterAzEl because
#       Z axis is incoming direction from the given sky FOV center. 




import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'TangentPlane27.py'

OBSLAT   = +31.9634     # latitude of KPNO
OBSELON  = -111.6003    # east longitude of KPNO
OBLIQ    = 23.4393      # obliquity of the ecliptic; so dec(NEP)=66.5607
YEARS    = 5



#---trig, vector, output helpers------------  

DEG = np.degrees(1.)
RAD = np.radians(1.)

def printVector(name, vector):
    print(name, end='')
    for item in vector:
        print('{:12.6f}'.format(item), end='')
    print()   

def sind(a):
    return np.sin(np.radians(a))

def cosd(a):
    return np.cos(np.radians(a))

def tand(a):
    return np.tan(np.radians(a))
    
def arcsind(x):
    return np.degrees(np.arcsin(x))

def arccosd(x):
    return np.degrees(np.arccos(x))

def arctan2d(y, x):
    return put360(np.degrees(np.arctan2(y, x)))

def getXYZ(lonlat):  # Convert spherical angles into equatorial xyz triplet
    return np.array([cosd(lonlat[0])*cosd(lonlat[1]), 
                     sind(lonlat[0])*cosd(lonlat[1]), 
                     sind(lonlat[1])])    
                     
def getLONLAT(xyz): # Convert equatorial xyz into its spherical angles
    xyz = getNormalized(xyz)  # usually unnecessary
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
    c=np.cos(np.radians(xdeg)); s=np.sin(np.radians(xdeg))
    return np.array([[1,0,0], [0,c,-s], [0,+s,c]])

def vecY(ydeg):  # For positive ydeg=-elev: +z=>+x; +x=>-z.
    # do not overlook this minus sign: positive ydeg pitches downward.
    c=np.cos(np.radians(ydeg)); s=np.sin(np.radians(ydeg))
    return np.array([[c,0,+s], [0,1,0], [-s,0,c]])
    
def vecZ(zdeg):  # For positive zdeg=+az: +x=>+y; +y=>-x.
    c=np.cos(np.radians(zdeg)); s=np.sin(np.radians(zdeg))
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


#--------other coordinate converters----------

def getViewXY(CenterAZEL, targetAZEL):
    # Convert sky coords into telescope view cube coordinates
    # print('getViewXY is deprecated')
    Caz = CenterAZEL[0]   # degrees; eastward from south
    Cel = CenterAZEL[1]   # degrees
    targetXYZ = getXYZ(targetAZEL)
    Vxyz = np.dot(refY(-Cel), np.dot(refZ(Caz), targetXYZ))
    Vxyz =  DEG*Vxyz
    return np.array([Vxyz[1], Vxyz[2]])
    
def getHorizAZEL(CenterAZEL, ViewXY):
    # Convert telescope view cube coords into sky angles
    # print('getHorizAZEL is deprecated')
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
    # print('getUV is deprecated')
    Caz = CenterAZEL[0]     # degrees; eastward from south
    Cel = CenterAZEL[1]     # degrees
    targetXYZ = getXYZ(targetAZEL)
    Vxyz = np.dot(refY(-Cel), np.dot(refZ(Caz), targetXYZ))
    return getLONLAT(Vxyz)

def getAZEL(CenterAZEL, teleUV):
    # Convert telescope view angles {u,v} into sky AZel in degrees
    # print('getAZEL is deprecated')
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

  
#----timekeeping------

def getGSTdeg(jd):  # Converts JD+frac to Greenwich Siderial Time in degrees
    return put360(360.9856469*(jd-2451545.22034))
    
def getDOY(jd):  # Converts JD.fraction to DayOfYear + fractional part
    date = getGreenwichDate(jd)
    year = (int) (date/10000)
    dateBOY = 10000*year + 0101.0
    jdBOY = getjulian(dateBOY)
    return jd - jdBOY

def getJD(yyyymmddf): # converts yyyymmdd.ffff to JD.fraction
    # T.C.van Flandern, Sky & Tel April 1981 p.312
    lymd = int(yyyymmddf)
    y = lymd // 10000
    m = np.fmod(lymd, 10000) // 100
    d = np.fmod(yyyymmddf, 100)  # float; includes the fraction
    return 367*y - (7*(y+((m+9)//12)))//4 + (275*m)//9 + d + 1721013.5
    
def getGDT(jd):  # Converts JD.fraction to yyyymmdd.fffff
    # T.C.van Flandern, Sky & Tel April 1981 p.312
    jd += 0.5
    f, z = np.modf(jd)  # frac and integer parts of jd
    if (z < 2299161):   
        a = z
    else:
       alpha = int((z - 1867216.25) / 36524.25)
       a = z + 1 + alpha - (alpha // 4)
    b = a + 1524 
    c = int((b - 122.1)/365.25)
    k = int(365.25*c)
    e = int((b-k)/30.6001)
    day = int(30.6001*e)
    day = b - k - day + f
    if e <= 13:
        month = e-1
    else:
        month = e-13
    if month > 2:
        year = c - 4716
    else:
        year = c - 4715
    if year < 0:
        return 10000*year - 100*month - day
    else:
        return 10000*year + 100*month + day
              
def getLSTdeg(jd, elon):  # Converts JD.fraction to local siderial time, degrees
    return put360(getGSTdeg(jd) + elon)
    
def getGMTdeg(jd):  # Converts JD to Greenwich Mean time, degrees
    return 360*(jd+0.5 % 1)  # JD starts at noon; GMT starts at midnight
    
def getMSThours(jd):   # Hours; Arizona is 7 hours behind Greenwich or 17 hours ahead.
    return (getGMTdeg(jd)/15.0 + 17) % 24.

def separation(lon1, lat1, lon2, lat2): # Angular separation of two lonlats, all in degrees
    # a1, a2 = RA or Lon in degrees; d1, d2 = dec or lat in degrees.
    # Returns the angular separation of these two directions.
    return arccosd(cosd(lon1-lon2)*cosd(lat1)*cosd(lat2) + sind(lat1)*sind(lat2))

def getSunLon(jd):  # Given JD.fff, returns sun ecliptic longitude degrees
    # https://en.wikipedia.org/wiki/Position_of_the_Sun
    days = jd - 2451545.0  # days from Greenwich noon 2000 Jan 01
    mean = 280.460 + 0.9856474*days
    anom = 357.528 + 0.9856003*days
    elon = mean + 1.915*sind(anom) + 0.020*sind(2*anom)
    return put360(elon);

def getUnitEarthVelocity(jd):  # ecliptic frame {x,y,z} unit vector
    # assumes circular orbit
    earthLon = getSunLon(jd) + 180
    return np.array([sind(earthLon), cosd(earthLon), 0.])







    
#-------sky evaluators----------
    
def getRefractionElev(elev):
    # COEF here is evaluated in Sellmeier-MAXADC.xlsx spreadsheet
    COEF = 0.00022275     # radians towards zenith, 800mb, 283K, 0.5878um    
    refraction = np.degrees(COEF*np.tan(np.radians(90.-elev)))
    return refraction  # in degrees; typically 0.01
    
def getDifferentialRefractionXY(centerAZEL, xy, centerRefraction):
    azelTail = getHorizAZEL(CenterAZEL, xy)
    dElev = getRefractionElev(azelTail[1]) - centerRefraction
    azelHead = azelTail + np.array([0., dElev]) 
    return getViewXY(CenterAZEL, azelHead)  
        
def getAberration(target, jd, magnif):
    # Given target ecliptic lonlat, and JD,
    # returns the fully aberrated & magnified LonLat.
    speed     = 0.99365E-4 # Meeus p.151 in radians
    apexlon   = getSunLon(jd) - 90.
    apexXYZ   = getXYZ(np.array([apexlon, 0.]))
    targetXYZ = getXYZ(target)
    VxT       = np.cross(apexXYZ, targetXYZ)
    TxVxT     = np.cross(targetXYZ, VxT)
    aberXYZ   = speed * TxVxT
    plotXYZ   = targetXYZ + magnif * aberXYZ
    return getLONLAT(getNormalized(plotXYZ))

def getPrecess(radec, years): 
    deltaELON = 0.013972*years
    ecl = getXYZ(radec2eclip(radec))
    pre = np.dot(vecZ(deltaELON), ecl) 
    eclip = getLONLAT(pre)
    return eclip2radec(eclip)   
    
def getParallactic(hadec, lat):
    # returns -180 < parallactic < +180
    azel  = hadec2azel(hadec, lat) 
    ha    = hadec[0]
    za    = 90. - azel[1]
    if za == 0.:
        return 0.
    codec = 90. - hadec[1] 
    colat = 90. - lat
    sinP  = sind(ha)*sind(colat)/sind(za)
    cosP  = (cosd(colat)-cosd(codec)*cosd(za))/(sind(codec)*sind(za))
    P     = np.degrees(np.arctan2(sinP, cosP))
    return P
    
    

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


"""
#=================PROGRAM STARTS HERE=================
#===Polar plot coordinates are AZ(deg) and ZA(deg)====
#==some Greenwich dates + times are listed here=======

# GDT = 20250320.29  # Midnight March20
# GDT = 20250320.79  # Noon March20  
# GDT = 20250620.29  # Midnight June20
# GDT = 20250920.29  # Midnight Sept20
# GDT = 20251220.29  # Midnight Sept20
# GDT = 20250320.90  # 2PM; nice view of AntiApex and Sun
# GDT = 20250320.50  # nice 2PM view of Apex and ecliptic pole
GDT = 20250320.815 # KPNO LST=0deg; view March noon sun & ecliptic pole
# GDT = 20250320.690 # KPNO LST=315deg like precession all sky plot
# GDT = 20250320.40  # 2AM view; neither apex nor antiapex
# GDT = 20250320.50  # 5AM view showing ecliptic pole and apex


def isOdd(n):
    return n % 2 != 0
    
def getHexagonalNumbers(nrings, rmax):
    # produces an array of {xy} pairs of hexagonal numbers
    # 1 degree angular offset avoids matplotlib "ValueError: Given lines do not intersect."
    nrays = 1 + 3*nrings + 3*nrings**2
    dr = rmax/nrings
    xypairs = list()
    xypairs.append([0., 0.])
    for ring in range(1, nrings+1):
        r = ring * dr
        dangle = 60./ring
        offset = 1.+dangle/2 if isOdd(ring) else 1.0
        for ray in range(6*ring):
            rad = np.radians(offset + dangle*ray)
            x = r*np.cos(rad)
            y = r*np.sin(rad)
            xypairs.append([x, y])
    return np.array(xypairs)

    
#----list the starting data--------------
print('GDT = {:15.4f}'.format(GDT))
jd  = getJD(GDT)
print('JD  = {:15.4f}'.format(jd))
lst = getLSTdeg(jd, OBSELON)
print('LST,deg = {:11.4f}'.format(lst))
print('MST,hours = {:9.4f}'.format(getMSThours(jd)))

#---locate the ecliptic pole in azel----
eclip = np.array([0., 90.])
radec = eclip2radec(eclip)
azel = radec2azel(radec, lst, OBSLAT)
printVector('eclip pole = ', eclip)
printVector('its radec  = ', radec)
printVector('its azel   = ', azel)

#-----------set up for our elevated tangent plane plot----------------

CenterAZ = 180.
CenterEL = 45.0
CenterAZEL = np.array([CenterAZ, CenterEL])

fig = plt.figure(figsize=(6.,6.))  # inches
ax = fig.add_subplot(1,1,1)
ax.set_xlim(-2.,+2.)
ax.set_ylim(-2.,+2.)
ax.set_aspect(1.0)
ArrowMag = 500 
plt.suptitle('Differential Atmospheric Refraction')
plt.title('0.589um   Elev='+str(int(CenterEL))+'  ArrowMag='+str(ArrowMag))

# draw constant azim contours for given view direction Vaz,Vel:
azims = np.arange(CenterAZ-25, CenterAZ+25, 1)
evals = np.arange(CenterEL-25, CenterEL+25, 1)

hue = 'blue'
for a in azims:
    aztemp = np.full_like(evals, a)
    azels = np.array([aztemp, evals])
    xys = getViewXY(CenterAZEL, azels)
    plt.plot(xys[0], xys[1], color=hue, linewidth=0.3)
    
# draw constant elevation contours

for e in evals:
    eltemp = np.full_like(azims, e)
    azels = np.array([azims, eltemp])
    xys = getViewXY(CenterAZEL, azels)
    plt.plot(xys[0], xys[1], color=hue, linewidth=0.3)      
    
  
centerRefraction = getRefractionElev(CenterEL)
for trial in range(200):
    xTail = np.random.uniform(-2, 2)    # plotbox coordinates deg on sky
    yTail = np.random.uniform(-2, 2)
    xyTail = np.array([xTail, yTail])
    
    
    # azelTail = getHorizAZEL(CenterAZEL, xyTail)
    # arrowLength = getRefractionElev(azelTail[1]) - centerRefraction
    # azelHead = azelTail + np.array([0., arrowLength]) 
    # xyHead = getViewXY(CenterAZEL, azelHead)
    
    xyHead = getDifferentialRefractionXY(CenterAZEL, xyTail, centerRefraction)
    width = ArrowMag *(xyHead[0] - xyTail[0])
    height = ArrowMag*(xyHead[1] - xyTail[1])
    ax.arrow(xTail, yTail, width, height, head_width=0.05, head_length=0.05, fc='black')

ax.tick_params(direction='in')
ax.set_xlabel('field angle x degrees')
ax.set_ylabel('field angle y degrees')
ax.text(0.01, 0.01, PROGNAME, transform=ax.transAxes, fontsize=8)

noext = PROGNAME.rsplit('.',1)[0]
filename = noext + '-' + '{:11.2f}'.format(GDT) + '.png'             
plt.savefig(filename) 
plt.show()           
"""