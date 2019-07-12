
# MountError9.py  
# Now exploring how to track center of FOV using
# four points per arrow iteration: North & Center on sky, Nom and Dev on Earth.
#
# M.Lampton UCB SSL 8 July 2019
# mlampton@berkeley.edu

import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'MountError9.py'

OBSLAT   = +31.9634     # latitude of KPNO
OBSELON  = -111.6003    # east longitude of KPNO
OBLIQ    = 23.4393      # obliquity of the ecliptic; so dec(NEP)=66.5607
YEARS    = 25

#  allsky plot constants
ZAMAX    = 70.0   # max zenith angle to calculate
RMAX     = ZAMAX  # max polar radius to plot
NRINGS   = 18     # how many rings of little arrows


MAGNIF   = 1000.  # arrow magnification

"""  Google definition:
"Azimuth: the direction of a celestial object from the observer, 
expressed as the angular distance from the north or south point 
of the horizon to the point at which a vertical circle passing 
through the object intersects the horizon."

Azimuth: four possible definitions.
  Eastward from North: Left handed {xyz}, does not converge to RaDec
  Westward from North: does not converge to RaDec.
  Eastward from South: OK; adopted here
  Westward from South: Left handed {xyz}, does not converge to RaDec.

"""


#---trig, vector, output helpers------------  

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
    # return put360(np.degrees(np.arctan2(y, x)))
    return np.degrees(np.arctan2(y,x))

def getXYZ(lonlat):  # Convert spherical angles into xyz triplet
    return np.array([cosd(lonlat[0])*cosd(lonlat[1]), 
                     sind(lonlat[0])*cosd(lonlat[1]), 
                     sind(lonlat[1])])    
                     
def getLONLAT(xyz): # Convert xyz into its spherical angles in degrees
    xyz = getNormalized(xyz)  # usually unnecessary
    return np.array([arctan2d(xyz[1],xyz[0]), arcsind(xyz[2])])
                     
def swapEW(lonlat):   # reverses east-west longitude sign
    return np.array([put360(-lonlat[0]), lonlat[1]])

def swapELZA(lonlat): # reverses zenith angle vs. elevation
    return np.array([lonlat[0], 90.-lonlat[1]])
    
def getNorm(xyz):
    return np.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)
    
def getNormalized(xyz):
    return xyz/getNorm(xyz)  
    
#---unitary matrices for vector rotations and reference frame rotations------

def vecX(xdeg):  # For positive xdeg clockwise around +X: +y=>+z; +z=>-y.
    c=cosd(xdeg); s=sind(xdeg)
    return np.array([[1,0,0], [0,c,-s], [0,+s,c]])

def vecY(ydeg):  # For positive ydeg clockwise around +Y: +z=>+x; +x=>-z.
    c=cosd(ydeg); s=sind(ydeg)
    return np.array([[c,0,+s], [0,1,0], [-s,0,c]])
    
def vecZ(zdeg):  # For positive zdeg=azim clockwise around +Z: +x=>+y; +y=>-x.
    c=cosd(zdeg); s=sind(zdeg)
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
    abc = np.dot(refY(-colat), xyz)   # pitch frame upward by colatitude
    return swapEW(getLONLAT(abc))     # ha=west, az=east

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
    # azimuth = eastward from south
    hadec = radec2hadec(radec, lst)
    return hadec2azel(hadec, lat)

def azel2radec(azel, lst, lat):
    # zenith is always at {ra,dec} = {lst, lat}
    # azimuth = eastward from south
    hadec = azel2hadec(azel, lat)
    return hadec2radec(hadec, lst)
    
def azza2radec(azza, lst, lat):
    azel = swapELZA(azza)
    return azel2radec(azel, lst, lat)
    
def radec2azza(radec, lst, lat):
    azel = radec2azel(radec, lst, lat)
    return swapELZA(azel)
    
    

"""
def radec2telexy(radec, lst, lat, dazmech, elmech, hangle, dangle, roll, flag):
    # Misaligned telescope, given radec gives image direction in telescope. Everything in degrees.
    azel = radec2azel(radec, lst, lat)                 # convert radec to azel
    azelxyz = getXYZ(azel)                             # convert azel to 3D Cartesian vector
    # opportunity here for atmospheric refraction
    dazelxyz = np.dot(refZ(dazmech), azelxyz)          # mechanical axis azimuth west of north 
    dazdelxyz = np.dot(refY(elmech-90.), dazelxyz)     # mech axis elevation
    clockedxyz = np.dot(refZ(-hangle), dazdelxyz)      # drive hour angle about this mech axis
    clockeddec = np.dot(refY(-dangle), clockedxyz)     # drive declination axis by dangle
    rolledxyz = np.dot(refX(roll), clockeddec)         # allow for telescope/FP roll
    telexy =   getLONLAT(rolledxyz)                    # get image xy in telescope
    return telexy
    
def telexy2radec(telexy, roll, dangle, hangle, elmech, dazmech, lat, lst, flag):
    # Misaligned telescope, given image direction in telescope, gives radec. Everything in degrees.
    rolledxyz = getXYZ(telexy)                         # get 3D Cartesian telescope image vector
    clockeddec = np.dot(refX(-roll), rolledxyz)        # remove telescope optical axis roll
    clockedxyz = np.dot(refY(dangle), clockeddec)      # remove dec drive angle setting
    dazdelxyz = np.dot(refZ(hangle), clockedxyz)       # remove hour angle setting
    dazelxyz = np.dot(refY(90.-elmech), dazdelxyz)     # remove mount elevation
    azelxyz = np.dot(refZ(-dazmech), dazelxyz)         # remove mount azimuth error
    azel = getLONLAT(azelxyz)                          # convert 3D Cartesian to azel angles
    radec = azel2radec(azel, lst, lat)                 # convert azel to radec
    return radec


"""
    
def radec2telexy(radec, lst, lat, dazmech, elmech, hangle, dangle, roll, flag):
    # Misaligned telescope, given radec gives image direction in telescope. Everything in degrees.
    if flag: printVector('radec      = ', radec)
    azel = radec2azel(radec, lst, lat)                 # convert radec to azel
    if flag: printVector('azel       = ', azel)
    azelxyz = getXYZ(azel)                             # convert azel to 3D Cartesian vector
    if flag: printVector('azelxyz    = ', azelxyz)
    # opportunity here for atmospheric refraction
    dazelxyz = np.dot(refZ(dazmech), azelxyz)          # mechanical axis azimuth west of north
    if flag: printVector('dazelxyz   = ', dazelxyz)    
    dazdelxyz = np.dot(refY(elmech-90.), dazelxyz)     # mech axis elevation
    if flag: printVector('dazdelxyz  = ', dazdelxyz)
    clockedxyz = np.dot(refZ(-hangle), dazdelxyz)      # drive hour angle about this mech axis
    if flag: printVector('clockedxyz = ', clockedxyz)
    clockeddec = np.dot(refY(-dangle), clockedxyz)     # drive declination axis by dangle
    if flag: printVector('clockeddec = ', clockeddec)
    rolledxyz = np.dot(refX(roll), clockeddec)         # allow for telescope/FP roll
    if flag: printVector('rolledxyz  = ', rolledxyz)
    telexy =   getLONLAT(rolledxyz)                    # get image xy in telescope
    if flag: printVector('telexy     = ', telexy)
    # opportunity here for corrector distortion
    return telexy
    
def telexy2radec(telexy, roll, dangle, hangle, elmech, dazmech, lat, lst, flag):
    # Misaligned telescope, given image direction in telescope, gives radec. Everything in degrees.
    if flag: printVector('telexy     = ', telexy)
    rolledxyz = getXYZ(telexy)                         # get 3D Cartesian telescope image vector
    if flag: printVector('rolledxyz  = ', rolledxyz)
    clockeddec = np.dot(refX(-roll), rolledxyz)        # remove telescope optical axis roll
    if flag: printVector('clockddec  = ', clockeddec)
    clockedxyz = np.dot(refY(dangle), clockeddec)      # remove dec drive angle setting
    if flag: printVector('clockedxyz = ', clockedxyz)
    dazdelxyz = np.dot(refZ(hangle), clockedxyz)       # remove hour angle setting
    if flag: printVector('dazdelxyz  = ', dazdelxyz)
    dazelxyz = np.dot(refY(90.-elmech), dazdelxyz)     # remove mount elevation
    if flag: printVector('dazelxyz   = ', dazelxyz)
    azelxyz = np.dot(refZ(-dazmech), dazelxyz)         # remove mount azimuth error
    if flag: printVector('azelxyz    = ', azelxyz)
    azel = getLONLAT(azelxyz)                          # convert 3D Cartesian to azel angles
    if flag: printVector('azel       = ', azel)
    radec = azel2radec(azel, lst, lat)                 # convert azel to radec
    if flag: printVector('radec      = ', radec)
    return radec

    
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
    
def getRefractionAngle(elev):
    COEF = 0.00022     # radians, for 600mmHg, 800mb, 0.7um  DESI-0656  
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
    
    
#------graphics helpers--------------

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



#=================PROGRAM STARTS HERE=================

GDT = 20250320.815 # KPNO LST=0deg; view March noon sun & ecliptic pole

#------test misalignment routines-------------------
print()
print('testing elmech=31 deg, lat=32deg, target in southern sky')
radec = np.array([0.0, 0.0])
lst = 0.
lat = 32.
dazmech = 0.
elmech = 31.  
hangle = 0.
dangle = 0.
roll   = 0.
telexy = radec2telexy(radec, lst, lat, dazmech, elmech, hangle, dangle, roll, True)
# teley =  elmech-lat = -1.000 IS CORRECT.
# If elmech is less elevated than the pole elevation = lat,
# then the telescope southern dec=0 plane will be elevated above equator,
# and a star on the equator will appear beneath the telescope axis.
# (Unless dangle is aimed downwards (towards south) to correct this.) 
print('inverse:')
testradec = telexy2radec(telexy, roll, dangle, hangle, elmech, dazmech, lat, lst, True)

print()
print('testing dazmech=+1deg WfromN, target in eastern sky')
radec = np.array([90.0, 0.0])
lst = 0.
lat = 32.
dazmech = +1.  # west of north or east of south
elmech = 32.  
hangle = -90.
dangle = 0.
roll   = 0.
telexy = radec2telexy(radec, lst, lat, dazmech, elmech, hangle, dangle, roll, True)
print('inverse:')
testradec = telexy2radec(telexy, roll, dangle, hangle, elmech, dazmech, lat, lst, True)

print('testing dazmech=+1deg WfromN, target in western sky')
radec = np.array([-90.0, 0.0])
lst = 0.
lat = 32.
dazmech = +1.  # west of north or east of south
elmech = 32.  
hangle = +90.
dangle = 0.
roll   = 0.
telexy = radec2telexy(radec, lst, lat, dazmech, elmech, hangle, dangle, roll, True)
print('inverse:')
testradec = telexy2radec(telexy, roll, dangle, hangle, elmech, dazmech, lat, lst, True)
# CORRECT: Positive dazmech drives eastern targets south so negative fpy
# and also drives western targets north so positive fpy.
#
# RESULTS
"""
testing elmech=31 deg, lat=32deg, target in southern sky
radec      =     0.000000    0.000000
azel       =     0.000000   58.000000
azelxyz    =     0.529919    0.000000    0.848048
dazelxyz   =     0.529919    0.000000    0.848048
dazdelxyz  =     0.999848    0.000000   -0.017452
clockedxyz =     0.999848    0.000000   -0.017452
clockeddec =     0.999848    0.000000   -0.017452
rolledxyz  =     0.999848    0.000000   -0.017452
telexy     =     0.000000   -1.000000
inverse:
telexy     =     0.000000   -1.000000
rolledxyz  =     0.999848    0.000000   -0.017452
clockddec  =     0.999848    0.000000   -0.017452
clockedxyz =     0.999848    0.000000   -0.017452
dazdelxyz  =     0.999848    0.000000   -0.017452
dazelxyz   =     0.529919    0.000000    0.848048
azelxyz    =     0.529919    0.000000    0.848048
azel       =     0.000000   58.000000
radec      =     0.000000    0.000000

testing dazmech=+1deg WfromN, target in eastern sky
radec      =    90.000000    0.000000
azel       =    90.000000    0.000000
azelxyz    =     0.000000    1.000000    0.000000
dazelxyz   =     0.017452    0.999848    0.000000
dazdelxyz  =     0.009248    0.999848   -0.014800
clockedxyz =     0.999848   -0.009248   -0.014800
clockeddec =     0.999848   -0.009248   -0.014800
rolledxyz  =     0.999848   -0.009248   -0.014800
telexy     =    -0.529958   -0.848036
inverse:
telexy     =    -0.529958   -0.848036
rolledxyz  =     0.999848   -0.009248   -0.014800
clockddec  =     0.999848   -0.009248   -0.014800
clockedxyz =     0.999848   -0.009248   -0.014800
dazdelxyz  =     0.009248    0.999848   -0.014800
dazelxyz   =     0.017452    0.999848    0.000000
azelxyz    =     0.000000    1.000000    0.000000
azel       =    90.000000    0.000000
radec      =    90.000000   -0.000000
testing dazmech=+1deg WfromN, target in western sky
radec      =   -90.000000    0.000000
azel       =   -90.000000   -0.000000
azelxyz    =     0.000000   -1.000000   -0.000000
dazelxyz   =    -0.017452   -0.999848   -0.000000
dazdelxyz  =    -0.009248   -0.999848    0.014800
clockedxyz =     0.999848   -0.009248    0.014800
clockeddec =     0.999848   -0.009248    0.014800
rolledxyz  =     0.999848   -0.009248    0.014800
telexy     =    -0.529958    0.848036
inverse:
telexy     =    -0.529958    0.848036
rolledxyz  =     0.999848   -0.009248    0.014800
clockddec  =     0.999848   -0.009248    0.014800
clockedxyz =     0.999848   -0.009248    0.014800
dazdelxyz  =    -0.009248   -0.999848    0.014800
dazelxyz   =    -0.017452   -0.999848   -0.000000
azelxyz    =     0.000000   -1.000000   -0.000000
azel       =   -90.000000   -0.000000
radec      =   270.000000   -0.000000



"""
#-------concludes testing radec2telexy and telexy2radec----------




#----list the starting data--------------
print()
print('GDT = {:15.4f}'.format(GDT))
jd  = getJD(GDT)
print('JD  = {:15.4f}'.format(jd))
lst = getLSTdeg(jd, OBSELON)
print('LST,deg = {:11.4f}'.format(lst))
print('MST,hours = {:9.4f}'.format(getMSThours(jd)))

RMAX     = ZAMAX = 70.0 # degrees
NRINGS   = 18
ARROWMAG = 50000.
DAZIM    = 0.01 # degrees west of north
DELEV    = 0.01 # degrees above latitude

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

plt.title('Guided Roll Error, Misaligned Telescope, DEL={:6.3f},  DAZ={:6.3f}'.format(DELEV,DAZIM))

#---show the polar coordinate grid----------
plt.grid(linestyle='solid')
# ax.grid(False) # this kills the coordinate grid
ax.tick_params(axis='x', pad=-4, labelsize=8, labelcolor='black') # 'x' is azimuth radians
ax.tick_params(axis='y', labelsize=6, labelcolor='black')         # 'y' is radial degrees


#------red+blue dot at the vernal equinox for this map's LST-------
ra = 0
dec = 0
ha = lst - ra
azel = hadec2azel(np.array([ha, dec]), OBSLAT)
azza = swapELZA(azel)
plt.scatter(np.radians(azza[0]), azza[1], s=50,  facecolors='red', edgecolors='blue', marker='o')

#------plot north equatorial pole blue star--- 
ha = 0.
dec = 90.
azel = hadec2azel(np.array([ha, dec]), OBSLAT)
azza = swapELZA(azel)
plt.scatter(np.radians(azza[0]), azza[1], s=150,  c='blue', marker='*', zorder=5)

#---draw some RA gridlines-------
hue = 'blue'
decs = np.arange(-45., 87., 1.)         # fine range of decs
for ra in range(0, 355, 15):            # coarse stepsize=15 is OK here
    ras = np.full_like(decs, ra)        # match up pairs for plotting
    radec = np.vstack((ras, decs))
    hadec = radec2hadec(radec, lst)
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

#--big yellow dot at Sun coordinates----------
SUNlonlat = np.array([getSunLon(jd), 0.])
SUNradec  = eclip2radec(SUNlonlat)
SUNhadec  = radec2hadec(SUNradec, lst)
SUNazel   = hadec2azel(SUNhadec, OBSLAT)
SUNazza   = swapELZA(SUNazel)
plt.scatter(np.radians(SUNazza[0]), SUNazza[1], s=250,  c='yellow', zorder=5)


xyvalues = getHexagonalNumbers(NRINGS, RMAX)
npts = len(xyvalues)
for p in range(npts):
    x                = xyvalues[p,0]    # Cartesian plot coords   
    y                = xyvalues[p,1]    # Cartesian plot coords
    azzaNomCenter    = np.array([arctan2d(y,x), np.sqrt(x**2 + y**2)])  # polar plot coordinates
    azelNomCenter    = swapELZA(azzaNomCenter)
    hadecNomCenter   = azel2hadec(azelNomCenter, OBSLAT)
    radecNomCenter   = hadec2radec(hadecNomCenter, lst)
    radecNomNorth    = radecNomCenter + np.array([0., 1.])  # one degree north
    hangle           = hadecNomCenter[0]
    dangle           = hadecNomCenter[1]
    roll             = 0.
    
    # now evaluate the four telescope {x,y} locations: Center & North sky, Nom and Deviated Mount
    
    telexyNomCenter  = radec2telexy(radecNomCenter, lst, OBSLAT, 0., OBSLAT, hangle, dangle, roll, False)
    telexyNomNorth   = radec2telexy(radecNomNorth,  lst, OBSLAT, 0., OBSLAT, hangle, dangle, roll, False)
    telexyDevCenter  = radec2telexy(radecNomCenter, lst, OBSLAT, DAZIM, OBSLAT+DELEV, hangle, dangle, roll, False)
    telexyDevNorth   = radec2telexy(radecNomNorth,  lst, OBSLAT, DAZIM, OBSLAT+DELEV, hangle, dangle, roll, False)    
    
    # now form the interesting difference into a tiny arrow
        
    arrowxy          = telexyDevCenter - telexyNomCenter  # simple difference
    arrowxy   = telexyDevNorth - telexyDevCenter - telexyNomNorth + telexyNomCenter
    
    telexyDevCenter = telexyNomCenter + ARROWMAG*arrowxy
    
    # now put this back into plot coordinates
    
    radecDevCenter  = telexy2radec(telexyDevCenter, roll, dangle, hangle, OBSLAT, 0., OBSLAT, lst, False)
    hadecDevCenter  = radec2hadec(radecDevCenter, lst)
    azelDevCenter   = hadec2azel(hadecDevCenter, OBSLAT)
    azzaDevCenter   = swapELZA(azelDevCenter) 
    ax.annotate('', xy=(np.radians(azzaDevCenter[0]), azzaDevCenter[1]), 
                xytext=(np.radians(azzaNomCenter[0]), azzaNomCenter[1]),
                arrowprops=dict(width=0.1, headwidth=2, headlength=3))
    
gdtstring = 'Greenwich Date Time ={:14.3f}'.format(GDT)
plt.text(-0.19, 1.0, gdtstring, transform=ax.transAxes, fontsize=8) 
latstring = 'KPNO E.Longitude, deg = {:9.4f}'.format(OBSELON)
plt.text(-0.19, 0.97, latstring, transform=ax.transAxes, fontsize=8)
latstring = 'KPNO Latitude, deg = {:8.4f}'.format(OBSLAT)
plt.text(-0.19, 0.94, latstring, transform=ax.transAxes, fontsize=8)
lststring = 'KPNO Siderial Time, deg = {:3.0f}'.format(lst)
plt.text(-0.19, 0.91, lststring, transform=ax.transAxes, fontsize=8) 
mst = getMSThours(jd)
mststring = 'MST, hours = {:6.2f}'.format(mst)
plt.text(-0.19, 0.88, mststring, transform=ax.transAxes, fontsize=8)
dazstring = 'DAZ deg W of North = {:6.3f}'.format(DAZIM)
plt.text(-0.19, 0.85, dazstring, transform=ax.transAxes, fontsize=8)
delstring = 'DEL deg above LAT = {:6.3f}'.format(DELEV)
plt.text(-0.19, 0.82, delstring, transform=ax.transAxes, fontsize=8)
arrowstr  = 'Arrow magnif = {:4.0f}'.format(ARROWMAG)
plt.text(-0.19, 0.79, arrowstr, transform=ax.transAxes, fontsize=8)

plt.text(0.1, 0.0, PROGNAME, transform=ax.transAxes, fontsize=8)
noext = PROGNAME.rsplit('.',1)[0]
noext += '_allsky'
filename = noext + '-' + '{:11.2f}'.format(GDT) + '.png'             
plt.savefig(filename)
plt.show()
