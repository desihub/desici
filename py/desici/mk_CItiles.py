'''
This mainly copies Stephen's CITiles notebook
'''


import fitsio
from astropy.table import Table
from astropy.io import fits
import desimodel.focalplane
import numpy as np
from matplotlib import pyplot as plt
import numpy as np
import os

version = 'v5'

def mktiles(scale=2,faintlim=18,mint=0,maxt=1e6,searchrange=5):
	'''
	This mainly copies Stephen's CITiles notebook
	scale allows the footprint for each CCD to be that factor greater in each dimension
	objects with Gaia G > faintlim are assigned the bit 3 too faint flag
	'''
	outdir = '/project/projectdirs/desi/cmx/ci/tiles/'+version+'/'
	alltiles = desimodel.io.load_tiles(onlydesi=False, extra=True)
	tiles = alltiles[alltiles['PASS'] == 0]
	keep = np.zeros(len(tiles), dtype=bool)
	keep[0] = True
	for i in range(1, len(tiles)):
		ra = tiles['RA'][keep]
		dec = tiles['DEC'][keep]
		r = angdist(tiles['RA'][i], tiles['DEC'][i], ra, dec)
		if np.min(r) > 5:
			keep[i] = True

	keep &= tiles['DEC'] > -30     
	citiles = Table(tiles[keep])	
	citiles.meta['EXTNAME'] = 'TILES'
	citiles.write('{}/ci-tiles-{}.fits'.format(outdir, version), overwrite=True)

	#- Give CI tiles new TILEIDs, making sure they don't overlap with the standard tiles
	citileid_min = ((np.max(alltiles['TILEID'])+1001)//1000)*1000
	citiles['CENTERID'] = citiles['TILEID']
	citiles['TILEID'] = np.arange(citileid_min, citileid_min+len(citiles), dtype=alltiles['TILEID'].dtype)
	print('Assigning CI TILEIDs {}-{}'.format(np.min(citiles['TILEID']), np.max(citiles['TILEID'])))
	assert np.all(~np.in1d(citiles['TILEID'], alltiles['TILEID']))

	#- Update PROGRAM
	citiles['PROGRAM'] = 'CI'


	#That prepared the tiles, now the targets and info need to be written to them
	targfile = '{}/ci-targets-{}.fits'.format(outdir, version)
	targets, tgthdr = fitsio.read(targfile, 1, header=True)

	onlyGaia = (targets['MORPHTYPE'] == b'GGAL') | (targets['MORPHTYPE'] == b'GPSF')
	print(np.count_nonzero(onlyGaia), 'targets with only Gaia morphologies')
	isPSF = (targets['MORPHTYPE'] == b'PSF ')
	gaia_g = targets['GAIA_PHOT_G_MEAN_MAG']
	isGPSF = (gaia_g<=19.0) & (targets['GAIA_ASTROMETRIC_EXCESS_NOISE'] < 10**0.5)
	isGPSF |= (gaia_g>=19.0) & (targets['GAIA_ASTROMETRIC_EXCESS_NOISE'] < 10**(0.5 + 0.2*(gaia_g - 19.0)))
	print('Targets with both LS and Gaia info', np.count_nonzero(~onlyGaia))
	ci_corners = Table.read(os.getenv('DESIMODEL')+'/data/focalplane/ci-corners.ecsv', format='ascii.ecsv')
	ciloc = desimodel.focalplane.GFALocations(ci_corners,scale=scale)
	#tile_indices = desimodel.footprint.find_points_in_tiles(citiles, targets['RA'], targets['DEC'], radius=1.7)
	#print('tile got indices')
	#above takes a long time, reducing time per tile (at risk of increasing total time) by pre-selecting per tile
	notargets = list()
	badtileids = list()
	#if mint == 0:
	#	mint = mint(citiles)
	if maxt == 1e6:
		maxt = len(citiles)+1	
	for i in range(mint,maxt):
	#     if citiles['IN_DESI'][i] == 0:
	#         continue
		print('tile '+str(i))
		tileid = citiles['TILEID'][i]
		#ii = tile_indices[i]
		#if len(ii) == 0:
		#	# print('Skipping tile ID {} with no input targets'.format(tileid))
		#	notargets.append(tileid)
		#	continue
		#else:
		#	# print(tileid)
		#	pass

		telra, teldec = citiles['RA'][i], citiles['DEC'][i]
		tar_sel = (targets['RA']>telra-searchrange/np.cos(teldec*np.pi/180.)) & (targets['RA']<telra+searchrange/np.cos(teldec*np.pi/180.)) & (targets['DEC']>teldec-searchrange) & (targets['DEC']<teldec+searchrange)
		ptargets = targets[tar_sel]
		if len(ptargets) == 0:
		# print('Skipping tile ID {} with no input targets'.format(tileid))
			notargets.append(tileid)
			continue
		#else:
		# print(tileid)
		#	pass
		
		if teldec > 85:
			print('might need different selection for this tile centered at dec '+str(teldec))
		#citargets = ciloc.targets_on_gfa(telra, teldec, targets[tile_indices[i]])
		citargets = ciloc.targets_on_gfa(telra, teldec, ptargets)
		if len(citargets) == 0:
			print("ERROR: no targets cover CI cameras on tile {}".format(tileid))
			badtileids.append(tileid)
			continue

		#- Copied & modified from fiberassign; should refactor that
		flag = np.ones(len(citargets), dtype="i2")
		ii = (citargets["MORPHTYPE"] == "PSF ") | (citargets["MORPHTYPE"] == "GPSF")
		if np.count_nonzero(ii) == 0:
			print("ERROR: no good GFA targets for "
					  "ETC/GUIDE/FOCUS on tile {}".format(tileid))
			badtileids.append(tileid)
			#- but keep going and write file even without guide stars
		tych = (0 < citargets['REF_ID']) 
		tych &= ( citargets['REF_ID'] < 1e10)

		flag[ii] = 0
		citargets["ETC_FLAG"] = flag
		citargets["GUIDE_FLAG"] = flag
		citargets[tych]["GUIDE_FLAG"] += 2**2
		faint = citargets['GAIA_PHOT_G_MEAN_MAG'] > 18
		citargets[faint]["GUIDE_FLAG"] += 2**3
		citargets["FOCUS_FLAG"] = flag
		ig =  np.zeros(len(citargets), dtype="bool") #this will be 0 if not a guide star, 1 if one
		ig[ii] =1

		ig[tych] = 0
		ig[faint] = 0
		citargets["IS_GUIDE"] = ig #this will be 0 if not a guide star, 1 if one
		
		#- Rename some columns to prevent ambiguity for ICS
		citargets.rename_column('RA', 'TARGET_RA')
		citargets.rename_column('DEC', 'TARGET_DEC')
		citargets.rename_column('RA_IVAR', 'TARGET_RA_IVAR')
		citargets.rename_column('DEC_IVAR', 'TARGET_DEC_IVAR')

		#- ADD KEYWORDS
		citargets.meta['EXTNAME'] = 'GFA_TARGETS'
		citargets.meta['TILEID'] = tileid
		citargets.meta['FIELDNUM'] = 0
		citargets.meta['TILERA'] = citiles['RA'][i]
		citargets.meta['TILEDEC'] = citiles['DEC'][i]
		citargets.meta['REQRA'] = citiles['RA'][i]
		citargets.meta['REQDEC'] = citiles['DEC'][i]
		citargets.meta['REFEPOCH'] = tgthdr['REFEPOCH']
	
		hdulist = fits.HDUList()
	
		#- Primary HDU with header keywords
		hdr = fits.Header()
		hdr['EXTNAME'] = 'PRIMARY'
		for key in ['TILEID', 'FIELDNUM', 'TILERA', 'TILEDEC', 'REQRA', 'REQDEC', 'REFEPOCH']:
			hdr[key] = citargets.meta[key]
		
		hdulist.append(fits.PrimaryHDU(data=None, header=hdr))
		hdulist.append(fits.convenience.table_to_hdu(citargets))
	
		outfile = os.path.join(outdir, 'citile-{:06d}.fits'.format(tileid))
		hdulist.writeto(outfile, overwrite=True)
		print('citile-{:06d}.fits'.format(tileid))
	 
def isolate_guidestars(isocrit=5.):
	'''
	add +2 to the GUIDE_FLAG if the angular distance to nearest star is < isocrit, in arcsec
	distance is calculated with flat-sky approximation taking cos(dec) factor into account for delta ra
	tile files are written back out into new version directory
	also, cut columns to those specified as needed
	'''
	from astropy.table import Table
	oldversion = 'v5/'
	newversion = 'v6/'
	tiledir = '/project/projectdirs/desi/cmx/ci/tiles/'
	indir = tiledir+oldversion
	outdir = tiledir+newversion
	isocrit = isocrit/3600.
	isocritsq = isocrit**2.

	for tl in range(58001,58995):
		try:
			f = fitsio.read(indir+'citile-0'+str(tl)+'.fits')
			w = f['GAIA_PHOT_G_MEAN_MAG'] < 18
			fw = f[w]
			gf = fw['GUIDE_FLAG']
			dl = np.zeros(len(gf)) #need to flag star when close pair is found so it doesn't get flagged multiple times
			inmax = np.max(gf)
			for i in range(0,len(fw)):
				for j in range(i+1,len(fw)):
					rad = fw[i]['TARGET_RA']-fw[j]['TARGET_RA']
					decd = fw[i]['TARGET_DEC']-fw[j]['TARGET_DEC']
					if abs(decd) < isocrit and dl[i] == 0 and dl[j] == 0:
						cosd = np.cos(fw[i]['TARGET_DEC']*np.pi/180.) #cosdec factor
						td = (decd**2.+(rad*cosd)**2.)
						if td < isocritsq:
							#print(i,j,td)
							gf[i] += 2
							gf[j] += 2		
							dl[i] = 1
							dl[j] = 1
			tab = Table.read(indir+'citile-0'+str(tl)+'.fits')
			tab[w]['GUIDE_FLAG'] = gf
			print(tl,np.max(tab['GUIDE_FLAG']),inmax)
			tab.keep_columns(['RELEASE', 'TARGETID', 'BRICKID', 'BRICK_OBJID', 'TARGET_RA', 'TARGET_DEC', 'TARGET_RA_IVAR', 'TARGET_DEC_IVAR', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'REF_ID', 'REF_CAT', 'PMRA', 'PMDEC', 'PMRA_IVAR', 'PMDEC_IVAR', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_FLUX_OVER_ERROR', 'GAIA_PHOT_RP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_FLUX_OVER_ERROR', 'GAIA_ASTROMETRIC_EXCESS_NOISE', 'HPXPIXEL', 'GFA_LOC', 'ETC_FLAG', 'GUIDE_FLAG', 'FOCUS_FLAG'])
			tab.write(outdir+'citile-0'+str(tl)+'.fits',overwrite=True)	
			#put the headers on
			fold = fits.open(indir+'citile-0'+str(tl)+'.fits')
			fnew = fits.open(outdir+'citile-0'+str(tl)+'.fits')
			fnew[0].header = fold[0].header
			fnew.writeto(outdir+'citile-0'+str(tl)+'.fits',overwrite=True)

		except:
			print('no '+str(tl)+'?')	


	return True
	
	

def angdist(ra1, dec1, ra2, dec2):
    '''
    Returns angular distance in degrees between (ra1,dec1) and (ra2,dec2),
    both measured in degrees.
    '''
    dec1 = np.radians(dec1)
    dec2 = np.radians(dec2)
    dd = dec1 - dec2
    dr = np.radians(ra1 - ra2)
    #- Haversine formula https://en.wikipedia.org/wiki/Haversine_formula
    r = 2 * np.arcsin(np.sqrt(
        np.sin(dd/2)**2 + (np.cos(dec1) * np.cos(dec2) * np.sin(dr/2)**2)
        ))
    return np.degrees(r)	