from math import *
import numpy as np

'''
characteristics common to each CCD
'''

CCDsize1 = 3072
CCDsize2 = 2048
CRPIX1 = CCDsize1/2.+.5
CRPIX2 = CCDsize2/2.+.5
 
'''
these are in SKY coordinates
compare to DESI 3347
CIW is C5
CIS is C4
CIC is C3
CIN is C2
CIE is C1
''' 

inddic = {'CIW': 0,'CIS': 1, 'CIC': 2, 'CIN': 3, 'CIE': 4} #dictionary for indexing

'''
orientations of CCDs
using 
ci_cameras = Table.read('/project/projectdirs/desi/cmx/ci/tiles/v3/ci-corners.ecsv',format='ascii.ecsv' )
ci = desimodel.focalplane.gfa.GFALocations(ci_cameras)
to determine plate scale
'''

rad = 1.563*pi/180. #average number of radians from center for outer cameras, all cameras are within 8 arcsec

platescalec = 0.0041 #deg/mm
platescaleor = 0.00365 #deg/mm
platescaleot = 0.00395 #deg/mm

pixscale = .009 #mm/pix
psc = pixscale*platescalec #deg/pix
psr = pixscale*platescaleor
pst = pixscale*platescaleot

print(psc,psr,pst)
print('plate scale, arcseconds/pix')
print('center, outer radial, outer transverse')
print(psc*3600,psr*3600,pst*3600)


#below assumes pixel x,y oriented orthogonally to CS5 X,Y; measured values cause < 4 arcsec change in position 

#1 west on sky C5, positive y is E (radial, increasing ra), x is N (transverse, increasing dec)
#Using 
cam1CD = np.zeros((2,2))
cam1CD[0][1] = 1.*psr #in units deg./pix, x is transverse
cam1CD[1][0] = 1.*pst


#2 south on sky C4, positive y is N (radial, increasing Dec), x is w (transverse, decreasing ra)
cam2CD = np.zeros((2,2))
cam2CD[0][0] = -1.*pst #in units deg./pix, x is transverse
cam2CD[1][1] = 1.*psr

#3 center, positive y is S (radial, decreasing dec), x is E (transverse, increasing ra)
cam3CD = np.zeros((2,2))
cam3CD[0][0] = 1.*psc #in units deg./pix
cam3CD[1][1] = -1.*psc

#4 north on sky, positive y is S (radial, decreasing dec), x is E (transverse, increasing ra)
cam4CD = np.zeros((2,2))
cam4CD[0][0] = 1.*pst
cam4CD[1][1] = -1.*psr

#5 east on sky, positive y is W (radial, decreasing ra), x is S (transverse, decreasing dec)
cam5CD = np.zeros((2,2))
cam5CD[0][1] = -1.*psr
cam5CD[1][0] = -1.*pst

camCDs = [cam1CD,cam2CD,cam3CD,cam4CD,cam5CD]


def getCD(ind,dec,stretch=False): #ra stretch depends on declination, I think
	'''
	get the CD matrix for the camera with index ind
	if stretch = True, declination stretch factor is added based on the declination given by dec
	'''
	camCD = camCDs[ind]
	if stretch:
		if ind == 0:
		#1 west, positive y is E (radial), x is N (transverse), index 0
			camCD[1][0] *= 1./cos(dec*pi/180.)
			return camCD
	
		if ind == 1:
		#2 south, positive y is N (radial), x is w (transverse)
			camCD[0][0] *= 1./cos(dec*pi/180.) #in units deg./pix, x is transverse
			return camCD

		if ind == 2:
		#3 center, positive y is N, x is w
			camCD[0][0] *= 1./cos(dec*pi/180.) #in units deg./pix
			return camCD

		if ind == 3:
		#4 north, positive y is S (radial), x is E (transverse)
			camCD[0][0] *= 1./cos(dec*pi/180.)
			return camCD

		if ind == 4:
		#5 east, positive y is W (radial), x is S (transverse
			camCD[1][0] *= 1./cos(dec*pi/180.)
			return camCD
	else:
		return camCD

def calcCRVALme(ra,dec):
	'''
	pass ra dec in degrees, return E/W dec and ra separation
	'''
	delta_0 = dec*pi/180.
	alpha_0 = ra*pi/180.
	alphap = 0
	if alpha_0 > 180:
		alpha_0 = alpha_0-360.
	if abs(dec + rad*180/pi) < 90.:
		dec2 = (dec)*pi/180.+rad
	else:
		dec2 = (dec)*pi/180.-rad
		alphap = 180.
		#decn = dec + rad*180/pi
		#return 'NEED to fix pole crossing!'
		#dec2 = 	(dec)*pi/180.
	sinDECEW = (cos(rad)*cos(dec2)/cos(delta_0)-cos(rad)**2.)/(tan(delta_0)*cos(dec2)-sin(dec2))
	DECEW = asin(sinDECEW)
	cosdra = (cos(rad)**2.-sin(delta_0+rad)*sin(DECEW))/(cos(delta_0+rad)*cos(DECEW))
	#print cosdra
	dra = acos(cosdra)*180/pi 
	
	return DECEW*180/pi,dra


def calcCRVAL_allcam(ra,dec):
	CRVALS = []
	#cam 1 is west
	decew,raoff = calcCRVALme(ra,dec)
	CRVALS.append([(ra-raoff),decew])
	#cam 2 is south
	dec2 = dec-rad*180./pi
	if dec2 < -90:
		dec2 = -90.-(dec2+90.)
	CRVALS.append((ra,dec2))
	#cam 3 is center
	CRVALS.append((ra,dec))
	#cam 4 is north
	dec4 = dec+rad*180./pi
	if dec4 > 90:
		dec4 = 90.-(dec4-90.)
	CRVALS.append((ra,dec4))
	#cam 5 is east
	CRVALS.append([(ra+raoff),decew])
	return CRVALS


def get_wcs(cam,tel_ra,tel_dec):
	'''
	input 'CIx', with x = W,S,C,N, or E and the telescope pointing tel_ra,tel_dec in degrees
	return dictionary with WCS info that can be added to header
	'''
	ind = inddic[cam]
	CRVALS = calcCRVAL_allcam(tel_ra,tel_dec)[ind]
	CDm = getCD(ind,CRVALS[1])
	outdic ={
	'CRVAL1' : CRVALS[0],
	'CRVAL2' : CRVALS[1],
	'CRPIX1' : CRPIX1,
	'CRPIX2' : CRPIX2,
	'CD1_1' : CDm[0][0],
	'CD1_2' : CDm[0][1],
	'CD2_1' : CDm[1][0],
	'CD2_2' : CDm[1][1]
	}
	return outdic
