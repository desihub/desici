import fitsio
from astropy.table import Table
from astropy.io import fits
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

def isolate_guidestars(isocrit=10.):
	'''
	add +2 to the GUIDE_FLAG if the angular distance to nearest star is < isocrit, in arcsec
	distance is calculated with flat-sky approximation taking cos(dec) factor into account for delta ra
	tile files are written back out into new version directory
	'''
	from astropy.table import Table
	oldversion = 'v3/'
	newversion = 'v4/'
	tiledir = '/project/projectdirs/desi/cmx/ci/tiles/'
	indir = tiledir+oldversion
	outdir = tiledir+newversion
	isocrit = isocrit/3600.
	isocritsq = isocrit**2.

	for tl in range(58002,58995):
		try:
			f = fitsio.read(indir+'citile-0'+str(tl)+'.fits')
			gf = f['GUIDE_FLAG']
			dl = np.zeros(len(gf)) #need to flag star when close pair is found so it doesn't get flagged multiple times
			inmax = np.max(gf)
			for i in range(0,len(f)):
				for j in range(i+1,len(f)):
					rad = f[i]['TARGET_RA']-f[j]['TARGET_RA']
					decd = f[i]['TARGET_DEC']-f[j]['TARGET_DEC']
					if abs(decd) < isocrit and dl[i] == 0 and dl[j] == 0:
						cosd = np.cos(f[i]['TARGET_DEC']*np.pi/180.) #cosdec factor
						td = (decd**2.+(rad*cosd)**2.)
						if td < isocritsq:
							#print(i,j,td)
							gf[i] += 2
							gf[j] += 2		
							dl[i] = 1
							dl[j] = 1
			tab = Table.read(indir+'citile-0'+str(tl)+'.fits')
			tab['GUIDE_FLAG'] = gf
			print(tl,np.max(tab['GUIDE_FLAG']),inmax)
			tab.write(outdir+'citile-0'+str(tl)+'.fits',overwrite=True)	
		except:
			print('no '+str(tl)+'?')	
	return True

def remove_tycho():
	'''
	remove objects that are Tycho stars
	'''
	from astropy.table import Table
	oldversion = 'v4/'
	newversion = 'v4/'
	tiledir = '/project/projectdirs/desi/cmx/ci/tiles/'
	indir = tiledir+oldversion
	outdir = tiledir+newversion

	for tl in range(58002,58995):
		try:
			tab = Table.read(indir+'citile-0'+str(tl)+'.fits')
			tych = (0 < tab['REF_ID']) 
			tych &= ( tab['REF_ID'] < 1e10)
			fo = tab[~tych]
			print(len(tab),len(fo))
			fo.write(outdir+'citile-0'+str(tl)+'.fits',overwrite=True)
			print(len(fitsio.read(outdir+'citile-0'+str(tl)+'.fits')))	
		except:
			print('no '+str(tl)+'?')	
	return True


def fixheader():
	'''
	get primary header from v3 into v4
	'''
	from astropy.table import Table
	oldversion = 'v3/'
	newversion = 'v4/'
	tiledir = '/project/projectdirs/desi/cmx/ci/tiles/'
	indir = tiledir+oldversion
	outdir = tiledir+newversion

	for tl in range(58002,58995):
		try:
			fold = fits.open(indir+'citile-0'+str(tl)+'.fits')
			fnew = fits.open(outdir+'citile-0'+str(tl)+'.fits')
			fnew[0].header = fold[0].header
			fnew.writeto(outdir+'citile-0'+str(tl)+'.fits',overwrite=True)
			print(len(fitsio.read_header(outdir+'citile-0'+str(tl)+'.fits')))	
		except:
			print('no '+str(tl)+'?')	
	return True

	
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

def get_tile_info(ramin=135,ramax=180,decmin=10,decmax =70,magtest=10,nstartest=1,mkplot=False,searchrange=5):
	'''
	get information about CI targets for each tile within selection
	searchrange cuts the target list to +/- around the center RA,DEC; RA is corrected by /cos(bt['DEC']*pi/180.)
	information is printed in caps if all five cameras have >= nstartest with mag < magtest
	9-12 hours first half, 9-16 full
	135 to 180 in ra 240 total
	10-70
	'''
	from math import cos,pi
	caml = [3,2,1,4,5] #order to match viewer top to bottom
	camdir = ['center','north','east','south','west']
	dirv4 = '/project/projectdirs/desi/cmx/ci/tiles/v4/'
	gtile = []
	for i in range(58002,58995):
		try:
		#tile = fitsio.read(dirci+'citile-0'+str(i)+'.fits')
			tileheader = fitsio.read_header(dirv4+'citile-0'+str(i)+'.fits',ext=1)
			ra,dec = tileheader['TILERA'],tileheader['TILEDEC']	
			if ra > ramin and ra < ramax and dec > decmin and dec < decmax:
				targets = fitsio.read(dirv4+'citile-0'+str(i)+'.fits')
				nmag = 0
				#print(str(i)+' get target info for each camera for pointing '+str(telra)+','+str(teldec)+' and scale '+str(scale))
				log = []
				for j in range (0,len(caml)):
					cam = caml[j]
					sel = (targets['GFA_LOC'] == cam)
					CIC = targets[sel]
					if mkplot:
						plotcam(cam,CIC)
					if len(CIC) > nstartest:
						ming = np.min(CIC['GAIA_PHOT_G_MEAN_MAG'])
						if ming < magtest:
							nmag += 1 
							log.append(str(len(CIC))+' stars on '+camdir[j]+' camera')
							log.append('brightest is '+str(ming))

					else:
						if len(CIC) > 0:
							ming = np.min(CIC['GAIA_PHOT_G_MEAN_MAG'])
						else:
							ming = 'NaN'		
					#fo.write(str(len(CIC))+' '+str(ming)+' ')	
				if nmag == 5:
					gtile.append(i)
					print(str(i)+' got target info for each camera for tile '+str(i)+' centered on '+str(ra) +' '+str(dec))
					print(str(nstartest)+' STARS ON ALL 5 CCDS PASSING MAG TEST for TILE '+str(i)+' centered on '+str(ra) +' '+str(dec))
					for ln in log:
						print(ln)
			else:
				pass
		except:
			print('no '+str(i)+'?')	
		#fo.write('\n')		
	#fo.close()
	print('good tiles are:')
	print(gtile)
	return True

	
def get_brightstar_info(magmaxcen=7,scale=1,ramin=135,ramax=180,decmin=10,decmax =70,fout='brightarginfo.txt',magtest=10,nstartest=1,mkplot=False,searchrange=5):
	'''
	select all ci targets brighter than magmaxcen and within ramin,ramx and decmin,decmax
	write out info about what is on other CCDs
	searchrange cuts the target list to +/- around the center RA,DEC; RA is corrected by /cos(bt['DEC']*pi/180.)
	9-12 hours first half, 9-16 full
	135 to 180 in ra 240 total
	10-70
	'''
	from math import cos,pi
	ci = desimodel.focalplane.gfa.GFALocations(ci_cameras,scale=scale)
	sel = (targets['GAIA_PHOT_G_MEAN_MAG'] < magmaxcen) & (targets['RA'] < ramax) & (targets['RA'] > ramin) & (targets['DEC'] < decmax) & (targets['DEC'] > decmin)
	brighttarg = targets[sel]
	caml = [3,2,1,4,5] #order to match viewer top to bottom
	camdir = ['center','north','east','south','west']

	#fo = open(dirout+fout,'w')
	#fo.write('#info about points with star brighter than '+str(magmaxcen)+' and with RA between '+str(ramin)+','+str(ramax)+' and DEC between '+str(decmin)+' '+str(decmax)+' at the center\n')
	#fo.write('#RA DEC Nstar_CIC Brightest_star_CIC Nstar_CIN Brightest_star_CIN  Nstar_CIE Brightest_star_CIE  Nstar_CIS Brightest_star_CIS  Nstar_CIW Brightest_star_CIW \n')
	print(str(len(brighttarg))+' stars brighter than '+str(magmaxcen)+ ' and with RA between '+str(ramin)+','+str(ramax)+' and DEC between '+str(decmin)+' '+str(decmax))
	print('finding info using them as the pointing center')
	for i in range(0,len(brighttarg)):
		bt = brighttarg[i]
		telra,teldec = bt['RA'],bt['DEC']
		#fo.write(str(bt['RA'])+' '+str(bt['DEC'])+' ')
		tar_sel = (targets['RA']>telra-searchrange/cos(bt['DEC']*pi/180.)) & (targets['RA']<telra+searchrange/cos(bt['DEC']*pi/180.)) & (targets['DEC']>teldec-searchrange) & (targets['DEC']<teldec+searchrange)
		ptargets = targets[tar_sel]
		citargets = ci.targets_on_gfa(telra, teldec, targets=ptargets)
		nmag = 0
		#print(str(i)+' get target info for each camera for pointing '+str(telra)+','+str(teldec)+' and scale '+str(scale))
		log = []
		for j in range (0,len(caml)):
			cam = caml[j]
			sel = (citargets['GFA_LOC'] == cam)
			CIC = citargets[sel]
			if mkplot:
				plotcam(cam,CIC)
			if len(CIC) > nstartest:
				ming = np.min(CIC['GAIA_PHOT_G_MEAN_MAG'])
				if ming < magtest:
					nmag += 1 
					log.append(str(len(CIC))+' stars on '+camdir[j]+' camera')
					log.append('brightest is '+str(ming))

			else:
				if len(CIC) > 0:
					ming = np.min(CIC['GAIA_PHOT_G_MEAN_MAG'])
				else:
					ming = 'NaN'		
			#fo.write(str(len(CIC))+' '+str(ming)+' ')	
		if nmag == 5:
			print(str(i)+' got target info for each camera for pointing '+str(telra)+','+str(teldec)+' and scale '+str(scale))
			print(str(nstartest)+' STARS ON ALL 5 CCDS PASSING MAG TEST AT POINTING '+str(telra)+' '+str(teldec))
			for ln in log:
				print(ln)
			
		#fo.write('\n')		
	#fo.close()
	return True
	
def plotcam(cam,CIC,rap,decp,winfac=0.01,tar=''):
	sizes = (18.1-CIC['GAIA_PHOT_G_MEAN_MAG'])*10
	#draw_camera(130.36, 24.525,3)
	#ra,dec = CIC['RA'],CIC['DEC']
	ii = ci_cameras['GFA_LOC'] == cam
	xx = ci_cameras['X'][ii]
	yy = ci_cameras['Y'][ii]
	ra, dec = desimodel.focalplane.xy2radec(rap,decp, xx, yy)
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
		plt.scatter(CIC[tar+'RA'],CIC[tar+'DEC'],s=sizes)
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
		plt.scatter(CIC[tar+'RA'],CIC[tar+'DEC'],s=sizes)
		plt.xlabel('RA')
		plt.ylabel('DEC')
		plt.title('South')

	if cam == 1: #Sky East, X/Y flipped RA/DEC both from positive to negative
		plt.xlim(xdec,mdec)
		plt.ylim(xra,mra)
		plt.scatter(CIC[tar+'DEC'],CIC[tar+'RA'],s=sizes)
		plt.xlabel('DEC')
		plt.ylabel('RA')
		plt.title('East')

	if cam == 5: #Sky West, X/Y flipped RA/DEC both from negative to positive 
		plt.xlim(mdec,xdec)
		plt.ylim(mra,xra)
		plt.scatter(CIC[tar+'DEC'],CIC[tar+'RA'],s=sizes)
		plt.xlabel('DEC')
		plt.ylabel('RA')
		plt.title('West')
	#if cam != 1 and cam != 5:
	#	plt.scatter(CIC[tar+'RA'],CIC[tar+'DEC'],s=sizes)
	#else:
	#	plt.scatter(CIC[tar+'DEC'],CIC[tar+'RA'],s=sizes)	
	plt.show()	

	
def get_info(telra,teldec,scale=1.,mkplot=True,searchrange=5):
	'''
	get summary of CI targets for given ra,dec pointing
	scale allows one to look in smaller or larger area than the detector area
	searchrange cuts the targets to ra,dec +/- within it (might need to be changed as function of dec)
	prints some basic info
	plots stars with sizes scaled to magnitude, with orientation that matches viewer, if mkplot is set to True
	'''
	caml = [3,2,1,4,5] #order to match viewer top to bottom
	camdir = ['center','north','east','south','west']
	ci = desimodel.focalplane.gfa.GFALocations(ci_cameras,scale=scale)
	tar_sel = (targets['RA']>telra-searchrange) & (targets['RA']<telra+searchrange) & (targets['DEC']>teldec-searchrange) & (targets['DEC']<teldec+searchrange)
	ptargets = targets[tar_sel]
	citargets = ci.targets_on_gfa(telra, teldec, targets=ptargets)
	print('get target info for each camera for pointing '+str(telra)+','+str(teldec)+' and scale '+str(scale))
	for i in range (0,len(caml)):
		cam = caml[i]
		sel = (citargets['GFA_LOC'] == cam)
		CIC = citargets[sel]
		if mkplot:
			plotcam(cam,CIC,telra, teldec)
		print(str(len(CIC))+' stars on '+camdir[i]+' camera')
		print('brightest is '+str(np.min(CIC['GAIA_PHOT_G_MEAN_MAG'])))
		print('that is it for now, easy to add more info, histograms, etc...')
		#plt.plot(ra, dec)
	
	