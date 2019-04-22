#just get info from pre-made tiles, so desimodel is not needed

import fitsio
import numpy as np
from matplotlib import pyplot as plt
import numpy as np
ci_tile_min = 58002
ci_tile_max = 58995

def get_tile_info(ramin=135,ramax=180,decmin=10,decmax =70,magtest=10,nstartest=1,mkplot=False,searchrange=5):
	'''
	get information about CI targets for each tile within selection
	searchrange cuts the target list to +/- around the center RA,DEC
	information is printed in caps if all five cameras have >= nstartest with mag < magtest
	9-12 hours first half, 9-16 full
	135 to 180 in ra 240 total
	10-70
	'''
	from math import cos,pi
	caml = [3,2,1,4,5] #order to match viewer top to bottom
	camdir = ['center','north','east','south','west']
	dirv4 = '/project/projectdirs/desi/cmx/ci/tiles/v4/'
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
					print(str(i)+' got target info for each camera for tile '+str(i)+' centered on '+str(ra) +' '+str(dec))
					print(str(nstartest)+' STARS ON ALL 5 CCDS PASSING MAG TEST for TILE '+str(i)+' centered on '+str(ra) +' '+str(dec))
					for ln in log:
						print(ln)
		except:
			print('no '+str(i)+'?')	
		#fo.write('\n')		
	#fo.close()
	return True
