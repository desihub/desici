
# TangentPlane8.py   
# Exploring divergence and curl operators
# Includes an allsky plot routine from Allsky19.py
# Explanatory writeup in DESI-DOC-4957
# M.Lampton UCB SSL 2019, mlampton@berkeley.edu

import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'TangentPlane8.py'

KPNOLAT  = +31.9634     # latitude of KPNO
KPNOELON = -111.6003    # east longitude of KPNO
OBLIQ    = 23.4393      # obliquity of the ecliptic; so dec(NEP)=66.5607
OBSLAT   = 32.0
OBSLON   = -111.6
YEARS    = 25

ZAMAX    = 70.0   # max zenith angle to calculate
RMAX     = ZAMAX  # max polar radius to plot
NRINGS   = 18     # how many rings of little arrows
MAGNIF   = 1000.  # arrow magnification




#---trig, vector, output helpers------------  

def printVector(name, vector):
    print(name, end='')
    for item in vector:
        print('{:12.6f}'.format(item), end='')
    print()   

def sind(xa):
    return np.sin(np.radians(xa))

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

def getXYZ(lonlat):  # Convert spherical angles into xyz triplet
    return np.array([cosd(lonlat[0])*cosd(lonlat[1]), 
                     sind(lonlat[0])*cosd(lonlat[1]), 
                     sind(lonlat[1])])    
                     
def getHorizonXYZ(az,el):
    return getXYZ(np.array([az,el]))
    
def getLONLAT(xyz): # Convert xyz into its spherical angles
    xyz = getNormalized(xyz)  # usually unnecessary
    return np.array([arctan2d(xyz[1],xyz[0]), arcsind(xyz[2])])
                     
def swapEW(lonlat):   # reverses east-west longitude sign
    return np.array([put360(-lonlat[0]), lonlat[1]])

def swapELZA(lonlat): # reverses zenith angle vs. elevation
    return np.array([lonlat[0], 90.-lonlat[1]])
    
def getNorm(xyz):
    # return np.sqrt((xyz**2).sum())  # fails for multiple xyz arrays
    return np.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)  # OK
    
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

def refX(xdeg):  # Rolls reference frame clockwise about X
    return vecX(-xdeg)

def refY(elev):  # Elevates reference frame about Y
    return vecY(+elev)
    
def refZ(azim):  # Rotates reference frame to new azimuth
    return vecZ(-azim)

def put360(degrees): # Puts an angle into range 0 to 360.
    return np.fmod(720.+degrees, 360)
    

       
#---------------coordinate converters--------------

def azel2hadec(azel, lat):
    xyz = getXYZ(azel)
    abc = np.dot(refY(90-lat), xyz) # pitch frame upward by colatitude
    return swapEW(getLONLAT(abc))   # ha=west, az=east

def eclip2radec(eclip):  # same epoch
    xyz = getXYZ(eclip)
    equ = np.dot(refX(-OBLIQ), xyz)  # roll frame counterclockwise by obliq
    return getLONLAT(equ)

def hadec2azel(hadec, lat):  
    azdec = swapEW(hadec)        # ha=west, az=east
    xyz = getXYZ(azdec)
    abc = np.dot(refY(lat-90), xyz)  # pitch frame downward by colatitude
    return getLONLAT(abc)  

def hadec2radec(hadec, lst):   # all in degrees
    return np.array([put360(lst - hadec[0]), hadec[1]])
    
def precess(radec, years): 
    deltaELON = 0.013972*years
    ecl = getXYZ(radec2eclip(radec))
    pre = np.dot(vecZ(deltaELON), ecl) 
    eclip = getLONLAT(pre)
    return eclip2radec(eclip)   
    
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
    hadec = radec2hadec(radec, lst)
    return hadec2azel(hadec, lat)

def azel2radec(azel, lst, lat):
    # zenith is always at {ra,dec} = {lst, lat}
    hadec = azel2hadec(azel, lat)
    return hadec2radec(hadec, lst)
    

  
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

#-------sky coordinate modifiers----------
    
def getRefractionAngle(elev):
    COEF = 0.00022     # radians, for 600mmHg, 800mb, 0.7um    
    refraction = np.degrees(COEF*np.tan(np.radians(90.-elev)))
    return refraction  # in degrees; typically 0.01
        
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
    # DANGER with getNormalized(): use only one triplet per call.
    return getLONLAT(getNormalized(plotXYZ))
    

    
    

#-----coordinate converters AzEl to/from Tangent Plane View----------

def getViewXY(CenterAZEL, otherAZEL):
    Caz = CenterAZEL[0]   # degrees
    Cel = CenterAZEL[1]   # degrees
    otherXYZ = getXYZ(otherAZEL)
    Vxyz = np.dot(refY(Cel), np.dot(refZ(Caz), otherXYZ))
    Vxyz =  np.degrees(Vxyz)
    return np.array([Vxyz[1], Vxyz[2]])

def getHorizAZEL(CenterAZEL, ViewXY):
    Caz = CenterAZEL[0]           # degrees
    Cel = CenterAZEL[1]           # degrees
    y = np.radians(ViewXY[0])     # build a normalized triplet
    z = np.radians(ViewXY[1])     # build a normalized triplet
    x = np.sqrt(1 - y**2 -z**2)   # build a normalized triplet
    Vxyz = np.array([x, y, z])
    Hxyz = np.dot(refZ(-Caz), np.dot(refY(-Cel), Vxyz))
    return getLONLAT(Hxyz)

#-----refraction divergence evaluator---------------

def getDivergence(theAZEL):
    DELTA = 1E-3  # equal steps in tangent plane
    dx = np.array([DELTA, 0.])
    dy = np.array([0., DELTA])
    eastAZEL = getHorizAZEL(theAZEL, dx)
    eastAZEL += np.array([0., getRefractionAngle(eastAZEL[1])])
    eastXY   = getViewXY(theAZEL, eastAZEL) - dx
    westAZEL = getHorizAZEL(theAZEL, -dx)
    westAZEL += np.array([0., getRefractionAngle(westAZEL[1])])
    westXY   = getViewXY(theAZEL, westAZEL) + dx
    upAZEL   = getHorizAZEL(theAZEL, dy)
    upAZEL   += np.array([0., getRefractionAngle(upAZEL[1])])
    upXY     = getViewXY(theAZEL, upAZEL) - dy
    downAZEL = getHorizAZEL(theAZEL, -dy)
    downAZEL += np.array([0., getRefractionAngle(downAZEL[1])])
    downXY   = getViewXY(theAZEL, downAZEL) + dy
    divx = (eastXY - westXY)/(2*DELTA)
    divy = (upXY - downXY)/(2*DELTA)
    #---now prepare the resulting record: elev, divx[0], divy[1]
    result = np.array([theAZEL[1], divx[0], divy[1]])
    printVector('elev, divx, divy = ', result)

      
#-------getDivergence() checks out OK--------
print()
print('testing getDivergence')
testAZELs =  np.array([[0., 30.],
                       [0., 45.],
                       [0., 60.],
                       [0., 75.],
                       [0., 90.]])    
for item in testAZELs:
    getDivergence(item)
    
# elev, divx, divy =    30.000000   -0.000220   -0.000880
# elev, divx, divy =    45.000000   -0.000220   -0.000440
# elev, divx, divy =    60.000000   -0.000220   -0.000293
# elev, divx, divy =    75.000000   -0.000220   -0.000236
# elev, divx, divy =    90.000000   -0.000220   -0.000220
#
# Looks like divx = -0.000220 for all elevations & azimuths
# Looks like divy = -0.000220*cosecant^2(elev)
  
        
def getRefractionCurl(theAZEL): 
    print()
    printVector('Starting getRefractionCurl() with AZEL = ', theAZEL)
    DELTA = 1E-3  # equal steps in tangent plane
    # sample refraction in four nearby locations
    dx = np.array([DELTA, 0.])
    dy = np.array([0., DELTA])
    eastAZEL = getHorizAZEL(theAZEL, dx)
    eastAZEL += np.array([0., getRefractionAngle(eastAZEL[1])])
    eastXY   = getViewXY(theAZEL, eastAZEL) - dx
    eastCurl = eastXY[1]  # upward component
    westAZEL = getHorizAZEL(theAZEL, -dx)
    westAZEL += np.array([0., getRefractionAngle(westAZEL[1])])
    westXY   = getViewXY(theAZEL, westAZEL) + dx
    westCurl = -westXY[1]  # downward component
    upAZEL   = getHorizAZEL(theAZEL, dy)
    upAZEL   += np.array([0., getRefractionAngle(upAZEL[1])])
    upXY     = getViewXY(theAZEL, upAZEL) - dy
    upCurl   = -upXY[0]   # rightward component
    downAZEL = getHorizAZEL(theAZEL, -dy)
    downAZEL += np.array([0., getRefractionAngle(downAZEL[1])])
    downXY   = getViewXY(theAZEL, downAZEL) + dy
    downCurl = downXY[0]  # leftward component 
    curl = eastCurl + westCurl + upCurl + downCurl
    print('Curl = {:12.6f}'.format(curl))
"""#----getRefractionCurl() is always zero: checks out OK-------    
print()
testAZELs =  np.array([[0.,  45.],
                       [45., 45.],
                       [45., 60.],
                       [45., 90.]])    
for item in testAZELs:
    getRefractionCurl(item)
"""

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

"""-----curl results------------
# Expect max curl at the north equatorial pole
# Expect positive curl at midnorthern declinations
# Expect zero curl along the celestial equator
# Expect negative curl at midsouthern declinations
# Expect max neg curl at the south equatorial pole.
#
# North equatorial pole stays at {az,el} = {180,lat}
# Equator passes through {az,el} = {0., colat}
# South equatorial pole stays at {az,el} = {0,-lat}

testAZELs =  np.array([[180., 32.],  # NEP:   curl = +17.4ppm
                       [180., 77.],  # +45:   curl = +12.3ppm
                       [0.,   58.],  # Equat: curl = 0
                       [0.,   13.],  # -45:   curl = -12.3ppm
                       [0.,  -32.]]) # SEP:   curl = -17.4ppm
            
print()
for item in testAZELs:
    curl = getRotationCurl(item)
    printVector('testAZELs = ',item)
    print('curl = {:12.9f}'.format(curl))
"""        





"""
#---test analytical FdotH() that moves azel into purely elevated Vxyz-------
#-----this is a reference frame shift to elevated view----
def FdotH(elea):     # Correct: if e=a=0 then x=1, all elevations
    el = elea[0]     # elevation of the target direction keeping az=0
    e  = elea[1]     # small change in elevation
    a  = elea[2]     # small change in azimuth
    x = cosd(el)*cosd(el+e)*cosd(a) + sind(el)*sind(el+e)
    y = cosd(el+e)*sind(a)
    z = -sind(el)*cosd(el+e)*cosd(a) + cosd(el)*sind(el+e)
    return np.array([x,y,z])
#--------checks out OK-----------    
print()
print('Running the FdotH test....')
vals = np.array([[45.,  0.,  0.],  # expect & get 1st order zero: {1., 0., 0.}
                 [45.,  1,   0.],  # expect & get 1st order in Z: {1., 0., +0.017}
                 [45.,  0.,  1.],  # expect & get 1st order in Y: {1., +0.012, 0.}
                 [45.,  1.,  1.]]) # expect & get 1st order both: {1., +0.012, +0.017}
for item in vals:
    xyz = FdotH(item)
    printVector('el,e,a = ', item)
    printVector('xyz    = ', xyz)
    print('norm =   {:12.6f}'.format(getNorm(xyz)))
    print()
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
LST = getLSTdeg(jd, OBSLON)
print('LST,deg = {:11.4f}'.format(LST))
print('MST,hours = {:9.4f}'.format(getMSThours(jd)))

fig = plt.figure(figsize=(8,7))          # inches
ax = fig.add_subplot(111, projection='polar')
ax.set_xticks([])
ax.set_theta_zero_location('S')          # thetaZero = South = bottom
ax.set_theta_direction(-1)               # increase clockwise = east from south
azticks = np.radians(30.)*np.arange(12)  # every 30 degrees
azlabels = ['S=0', '', '', 'E=90', '', '', 'N=180', '', '',  'W=270', '', '']
ax.set_xticks(azticks)
ax.set_xticklabels(azlabels)
ax.set_ylim(0, ZAMAX)

# plt.grid(linestyle='dotted')
ax.grid(False)

ax.tick_params(axis='x', pad=-4, labelsize=8, labelcolor='black') # 'x' is azimuthal
ax.tick_params(axis='y', labelsize=6, labelcolor='black')         # 'y' is radial.
plt.title("Annual Aberration, KPNO, {AZ,ZA}")

#------plot north equatorial pole--- 
ha = 0.
dec = 90.
azel = hadec2azel(np.array([ha, dec]), OBSLAT)
azza = swapELZA(azel)
plt.scatter(np.radians(azza[0]), azza[1], s=100,  c='blue', marker='*')

#-----plot north ecliptic pole------- 
ra = 270.
dec = 90. - OBLIQ
radecFafnir = np.array([ra,dec])
hadecFafnir = radec2hadec(radecFafnir, LST)
azel = hadec2azel(hadecFafnir, OBSLAT)
azza = swapELZA(azel)
plt.scatter(np.radians(azza[0]), azza[1], s=100,  c='red', marker='*')

#------red+blue dot at the vernal equinox for this map's LST-------
ra = 0
dec = 0
ha = LST - ra
azel = hadec2azel(np.array([ha, dec]), OBSLAT)
azza = swapELZA(azel)
plt.scatter(np.radians(azza[0]), azza[1], s=50,  facecolors='red', edgecolors='blue', marker='o')

#---green triangle at Earth-Velocity-Apex using ecliptic coordinates---------
EVAxyz    = getUnitEarthVelocity(jd)         # ecliptic xyz
EVAlonlat = getLONLAT(EVAxyz)                # ecliptic lonlat
EVAradec  = eclip2radec(EVAlonlat)
EVAhadec  = radec2hadec(EVAradec, LST)  
EVAazel   = hadec2azel(EVAhadec, OBSLAT)
EVAazza   = swapELZA(EVAazel)
plt.scatter(np.radians(EVAazza[0]), EVAazza[1], s=50,  facecolors='green', edgecolors='green', marker='^')

#---red triangle at AntiApex---------------
EVAxyz    = getUnitEarthVelocity(jd)
EVAxyz    *= -1.0    
EVAlonlat = getLONLAT(EVAxyz)
EVAradec  = eclip2radec(EVAlonlat)
EVAhadec  = radec2hadec(EVAradec, LST)  
EVAazel   = hadec2azel(EVAhadec, OBSLAT)
EVAazza   = swapELZA(EVAazel)
plt.scatter(np.radians(EVAazza[0]), EVAazza[1], s=50,  facecolors='red', edgecolors='red', marker='v')

#---draw some RA gridlines-------
hue = 'blue'
decs = np.arange(-45., 87., 1.)         # fine range of decs
for ra in range(0, 355, 15):            # coarse stepsize=15 is OK here
    ras = np.full_like(decs, ra)        # match up pairs for plotting
    radec = np.vstack((ras, decs))
    hadec = radec2hadec(radec, LST)
    azel = hadec2azel(hadec, OBSLAT)
    azza = swapELZA(azel)
    lw = 1.5  if (ra==0) else 0.5
    plt.plot(np.radians(azza[0]), azza[1], hue, lw=lw)
    
#--draw some dec gridlines--------
hue = 'blue'
has = np.arange(0., 360.)                 # fine range of hour angles
for dec in range(-15, 90, 15):            # coarse steps in declination
    decs = np.full_like(has, dec)         # match up pairs for plotting
    azel = hadec2azel(np.vstack((has, decs)), OBSLAT)
    azza = swapELZA(azel)
    lw = 1.5  if (dec==0) else 0.5
    plt.plot(np.radians(azza[0]), azza[1], hue, lw=lw)

#---draw some ecliptic longitude gridlines---------
#---Unwelcome extra gridline around elon = 195 if step=15 and LST=0
hue = 'red'
elats = np.arange(-55., 87., 1.)         # fine range in latitudes
for elon in range(0, 355, 15):           # stepsize=15 sometimes fails here
    elons = np.full_like(elats, elon)    # match up pairs for plotting
    radec = eclip2radec(np.vstack((elons, elats)))
    hadec = radec2hadec(radec, LST)
    azel = hadec2azel(hadec, OBSLAT)
    azza = swapELZA(azel)
    lw = 1.5 if (elon==0) else 0.5
    plt.plot(np.radians(azza[0]), azza[1], hue, lw=lw, linestyle='dashed')        

#----draw some ecliptic latitude gridlines------
elons = np.arange(0., 360.)              # fine range of longitudes
for elat in range(-75, 90, 15):          # steps in latitudes
    elats = np.full_like(elons, elat)    # match up pairs for plotting
    radec = eclip2radec(np.vstack((elons, elats)))
    hadec = radec2hadec(radec, LST)
    azel = hadec2azel(hadec, OBSLAT)
    azza = swapELZA(azel)
    lw = 1.5 if (elat==0) else 0.5
    plt.plot(np.radians(azza[0]), azza[1], hue, lw=lw, linestyle='dashed')
    
#--big yellow dot at Sun coordinates----------
SUNlonlat = np.array([getSunLon(jd), 0.])
SUNradec  = eclip2radec(SUNlonlat)
SUNhadec  = radec2hadec(SUNradec, LST)
SUNazel   = hadec2azel(SUNhadec, OBSLAT)
SUNazza   = swapELZA(SUNazel)
plt.scatter(np.radians(SUNazza[0]), SUNazza[1], s=250,  facecolors='yellow', edgecolors='black')

xyvalues = getHexagonalNumbers(NRINGS, RMAX)
npts = len(xyvalues)
for p in range(npts):
    x         = xyvalues[p,0]    # Cartesian plot coords   
    y         = xyvalues[p,1]    # Cartesian plot coords
    azzaTail  = np.array([arctan2d(y,x), np.sqrt(x**2 + y**2)])  
    azelTail  = swapELZA(azzaTail)
    hadecTail = azel2hadec(azelTail, KPNOLAT)
    radecTail = hadec2radec(hadecTail, LST)
    eclipTail = radec2eclip(radecTail)
    eclipLat  = eclipTail[1]
    eclipHead = getAberration(eclipTail, jd, MAGNIF) 
    radecHead = eclip2radec(eclipHead) 
    hadecHead = radec2hadec(radecHead, LST)
    azelHead  = hadec2azel(hadecHead, KPNOLAT)
    azzaHead  = swapELZA(azelHead) 
    ax.annotate('', xy=(np.radians(azzaHead[0]), azzaHead[1]), 
                xytext=(np.radians(azzaTail[0]), azzaTail[1]),
                arrowprops=dict(width=0.1, headwidth=2, headlength=3))

gdtstring = 'Greenwich Date Time ={:14.3f}'.format(GDT)
plt.text(-0.15, 1.0, gdtstring, transform=ax.transAxes, fontsize=8) 
lststring = 'KPNO Siderial Time, deg = ' + str(int(LST))
plt.text(-0.15, 0.97, lststring, transform=ax.transAxes, fontsize=8) 
latstring = 'KPNO Longitude, deg = {:6.1f}'.format(OBSLON)
plt.text(-0.15, 0.94, latstring, transform=ax.transAxes, fontsize=8)
latstring = 'KPNO Latitude, deg = {:6.1f}'.format(OBSLAT)
plt.text(-0.15, 0.91, latstring, transform=ax.transAxes, fontsize=8)
mst = getMSThours(jd)
mststring = 'MST, hours = {:6.2f}'.format(mst)
plt.text(-0.15, 0.88, mststring, transform=ax.transAxes, fontsize=8)

magstring = 'Arrow magnif = ' + str(MAGNIF)
plt.text(-0.15, 0.85, magstring, transform=ax.transAxes, fontsize=8)
plt.text(0.1, 0.0, PROGNAME, transform=ax.transAxes, fontsize=8)
noext = PROGNAME.rsplit('.',1)[0]
noext += '_allsky'
filename = noext + '-' + '{:11.2f}'.format(GDT) + '.png'             
plt.savefig(filename)
plt.show()









#-----------set up for our elevated plot----------------

CenterAZEL = np.array([0., 45.0])

fig, ax = plt.subplots()
ax.set_xlim(-10., +10.)
ax.set_ylim(-10., +10.)

# draw constant azim contours for given view direction Vaz,Vel:
evals = np.arange(20, 60, 1)
azims = range(-20, 20, 1)
for a in azims:
    aztemp = np.full_like(evals, a)
    azels = np.array([aztemp, evals])
    xys = getViewXY(CenterAZEL, azels)
    plt.plot(xys[0], xys[1], color='blue', linewidth=0.3)
    
# draw constant elevation contours
evals = range(20, 60, 1)
azims = np.arange(-20., 20., 1.)
for e in evals:
    eltemp = np.full_like(azims, e)
    azels = np.array([azims, eltemp])
    xys = getViewXY(CenterAZEL, azels)
    plt.plot(xys[0], xys[1], color='red', linewidth=0.3)    

# draw 100 arrows but CAUTION tangent approximation is no good near horizon
ArrowMag = 100.0 
for trial in range(100):
    xTail = np.random.uniform(-10, 10)                        # plotbox coordinates
    yTail = np.random.uniform(-10, 10)
    xyTail = np.array([xTail, yTail])
    azelTail = getHorizAZEL(CenterAZEL, xyTail)               # get AzEl coordinates
    arrowLength = ArrowMag * getRefractionAngle(azelTail[1])  # degrees in elevation
    azelHead = azelTail + np.array([0., arrowLength])         # same azim, raised elev
    xyHead = getViewXY(CenterAZEL, azelHead)                  # get plotbox head coords
    xHead = xyHead[0]
    yHead = xyHead[1]
    width = xHead - xTail
    height = yHead - yTail
    ax.arrow(xTail, yTail, width, height, head_width=0.25, head_length=0.5, fc='black')

ax.tick_params(direction='in')
ax.set_xlabel('field angle x degrees')
ax.set_ylabel('field angle y degrees')
ax.text(0.01, 0.01, PROGNAME, transform=ax.transAxes, fontsize=8)
plt.title('Az-El grid and refractions in field at center elevation 45deg')

noext = PROGNAME.rsplit('.',1)[0]
noext += '_elevated'
filename = noext + '-' + '{:11.2f}'.format(GDT) + '.png'             
plt.savefig(filename) 
plt.show()           
