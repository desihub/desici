import fitsio
from astropy.table import Table
import desimodel.focalplane
import numpy as np
from matplotlib import pyplot as plt
import numpy as np
ci_tile_min = 58002
ci_tile_max = 58995

#- Load a targets file that has columns RA, DEC
dirci = '/project/projectdirs/desi/cmx/ci/tiles/v3/'
dirout = '/project/projectdirs/desi/cmx/ci/tiles/v3/tileinfo/'
#dirout = ''
targets = fitsio.read(dirci+'ci-targets-v3.fits')

#- Load the CI camera location definitions
ci_cameras = Table.read('/project/projectdirs/desi/cmx/ci/tiles/v3/ci-corners.ecsv',format='ascii.ecsv' )

'adding some text to make changes for new commit'

def get_Donut_tiles(gmax=15,scale=0.5,ncammin=3,fout='Donut_tiles.txt'):
	'''
	find tiles that should be good for donut analysis
	gmax is maximum gaia mag to accept for star
	scale is how much of the detector to search over
	ncammin is the minimum number of cameras to accept
	file fout gets written out in dirout directory with tile number, total number of stars < gmax in the scaled area, the number of cameras with at least on star
	'''
	fo = open(dirout+fout,'w')
	fo.write('# written with options \n')
	fo.write('# gmax ='+str(gmax)+'\n')
	fo.write('# scale ='+str(scale)+'\n')
	fo.write('# ncammin ='+str(ncammin)+'\n')
	fo.write('# tile_num num_stars num_cam\n')
	nl = []
	ci = desimodel.focalplane.gfa.GFALocations(ci_cameras,scale=scale)
	for i in range(58002,58995):
		try:
		#tile = fitsio.read(dirci+'citile-0'+str(i)+'.fits')
			tileheader = fitsio.read_header(dirci+'citile-0'+str(i)+'.fits')
			ra,dec = tileheader['TILERA'],tileheader['TILEDEC']	
			print('testing tile '+str(i)+' with center '+str(ra)+','+str(dec))		
			#print('need way to get ra/dec of tile center')
			#ra,dec = tileheader(ra,dec)
			citargets = ci.targets_on_gfa(telra=ra, teldec=dec, targets=targets)
			#print(len(citargets))
			sel = (citargets['GAIA_PHOT_G_MEAN_MAG'] < gmax)
			#print(len(citargets[sel]))
			locs = np.unique(citargets[sel]['GFA_LOC'])
		
			if len(locs) >= ncammin:
				print('#targets over #cameras')
				print(len(citargets[sel]),len(locs))
				fo.write(str(i)+' '+str(len(citargets[sel]))+' '+str(len(locs))+'\n')
			else:
				print('failed')	
		except:
			print('tile '+str(i) +' not found')
	fo.close()
	return True
	
def get_brightstar_info(magmaxcen=7,scale=1,fout='brightarginfo.txt',mkplot=False):
	'''
	select all ci targets brighter than magmaxcen
	write out info about what is on other CCDs
	'''
	ci = desimodel.focalplane.gfa.GFALocations(ci_cameras,scale=scale)
	sel = (targets['GAIA_PHOT_G_MEAN_MAG'] < magmaxcen)
	brighttarg = targets[sel]
	caml = [3,2,1,4,5] #order to match viewer top to bottom
	camdir = ['center','north','east','south','west']

	fo = open(dirout+fout)
	fo.write('#info about points with star brighter than '+str(magmaxcen)+' at the center\n')
	fo.write('#RA DEC Nstar_CIC Brightest_star_CIC Nstar_CIN Brightest_star_CIN  Nstar_CIE Brightest_star_CIE  Nstar_CIS Brightest_star_CIS  Nstar_CIW Brightest_star_CIW \n')
	print(str(len(brightarg))+' stars brighter than '+str(magmaxcen))
	print('finding info using them as the pointing center')
	for i in range(0,len(brighttarg)):
		telra,teldec = brighttarg['RA'],brighttarg['DEC']
		fo.write(str(brighttarg['RA']+' '+str(brighttarg['DEC']))+' ')
		citargets = ci.targets_on_gfa(telra, teldec, targets=brighttarget)
		print(str(i)+' get target info for each camera for pointing '+str(telra)+','+str(teldec)+' and scale '+str(scale))
		for i in range (0,len(caml)):
			cam = caml[i]
			sel = (citargets['GFA_LOC'] == cam)
			CIC = citargets[sel]
			if mkplot:
				plotcam(cam,CIC)
			fo.write(str(len(CIC))+' '+str(np.min(CIC['GAIA_PHOT_G_MEAN_MAG']))+' ')	
		fo.write('\n')		
	fo.close()
def plotcam(cam,CIC,winfac=0.01):
	sizes = (18.1-CIC['GAIA_PHOT_G_MEAN_MAG'])*10
	#draw_camera(130.36, 24.525,3)
	ra,dec = CIC['RA'],CIC['DEC']
	disra = np.max(ra) - np.min(ra)
	disdec = np.max(dec) - np.min(dec)
	xdec = winfac*disdec+np.max(dec)
	mdec = winfac*disdec+np.min(dec)
	xra = winfac*disra+np.max(ra)
	mra = winfac*disra+np.min(ra)
	if cam == 3 or cam == 2: #these go from positive to negative dec with positive Y and from negative to postive RA with positive X
		#3 is Center, 2, is North
		plt.ylim(xdec,mdec)
		plt.xlim(mra,xra)
		plt.scatter(CIC['RA'],CIC['DEC'],s=sizes)
		plt.xlabel('RA')
		plt.ylabel('DEC')
		if cam == 3:
			plt.title('Center')
		if cam == 2:
			plt.title('North')
	if cam == 4: #these go from negative to positivedec with positive Y and from postive to negative RA with positive X
		#Sky South
		plt.ylim(mdec,xdec)
		plt.xlim(xra,mra)
		plt.scatter(CIC['RA'],CIC['DEC'],s=sizes)
		plt.xlabel('RA')
		plt.ylabel('DEC')
		plt.title('South')

	if cam == 1: #Sky East, X/Y flipped RA/DEC both from positive to negative
		plt.xlim(xdec,mdec)
		plt.ylim(xra,mra)
		plt.scatter(CIC['DEC'],CIC['RA'],s=sizes)
		plt.xlabel('DEC')
		plt.ylabel('RA')
		plt.title('East')

	if cam == 5: #Sky West, X/Y flipped RA/DEC both from negative to positive 
		plt.xlim(mdec,xdec)
		plt.ylim(mra,xra)
		plt.scatter(CIC['DEC'],CIC['RA'],s=sizes)
		plt.xlabel('DEC')
		plt.ylabel('RA')
		plt.title('West')

	plt.show()	

	
def get_info(telra,teldec,scale=1.,mkplot=True):
	'''
	get summary of CI targets for given ra,dec pointing
	scale allows one to look in smaller or larger area than the detector area
	prints some basic info
	plots stars with sizes scaled to magnitude, with orientation that matches viewer, if mkplot is set to True
	'''
	caml = [3,2,1,4,5] #order to match viewer top to bottom
	camdir = ['center','north','east','south','west']
	ci = desimodel.focalplane.gfa.GFALocations(ci_cameras,scale=scale)
	citargets = ci.targets_on_gfa(telra, teldec, targets=targets)
	print('get target info for each camera for pointing '+str(telra)+','+str(teldec)+' and scale '+str(scale))
	for i in range (0,len(caml)):
		cam = caml[i]
		sel = (citargets['GFA_LOC'] == cam)
		CIC = citargets[sel]
		if mkplot:
			plotcam(cam,CIC)
		print(str(len(CIC))+' stars on '+camdir[i]+' camera')
		print('brightest is '+str(np.min(CIC['GAIA_PHOT_G_MEAN_MAG'])))
		print('that is it for now, easy to add more info, histograms, etc...')
		#plt.plot(ra, dec)
	
	