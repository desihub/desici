import argparse
import time
import logging
import numpy as np
import numpy
import sys
import os
import scipy
from scipy.interpolate import interp1d
from operator import itemgetter
from math import *
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
from photutils import DAOStarFinder
import fitsio
max_star_num=100
cic_count=0
cis_count=0
cin_count=0
ciw_count=0
cie_count=0


def get_parser():
	'''return parser object, tells it what options to look for
	options can come from command line'''
	parser = argparse.ArgumentParser(formatter_class=argparse.
								 ArgumentDefaultsHelpFormatter,
								 description='DECaLS Focus, seeing, alignment analysis')
	parser.add_argument('-co', '--count', type=int, default=1, required = False, help = 'number of images used to perform the FOCUS analysis, required in FOCUS mode and alignment mode')
	parser.add_argument('--expid', type=int, default=None, required = True, help = 'starting expose id used for analysis, 4 digit number')
	parser.add_argument('-d', '--obsday', type=str, default=None, required = True, help = 'date of observation, example: 20190306')
	parser.add_argument('-ca' , '--camera', type=str, default='ALL', required = False, choices = ['ALL','CIS','CIW','CIE','CIN','CIC'], help = 'choose which camera to use')
	parser.add_argument('-m', '--mode', type=str, default='FOCUS', required = False, choices=['FOCUS', 'alignment','seeing','star_sizes','star_sizes_plot'], help = 'choose the mode of the script.\n FOCUS: comput the FOCUS using input image.\n seeing: comput fwhm for input image.\n alignment: compute tip tilt value')
	parser.add_argument('-np', '--noplot', action='store_true', default = False, required = False, help='do not show the Focus curves plot')
	parser.add_argument('--skipid', nargs='+', help='skip certain expids, 4 digit number', required = False)
	parser.add_argument('-v','--verbose', action='store_true', default = False, required = False, help = 'verbose level of the script')
	return parser


class filesystem():
	def __init__(self,count,expid,date,skipid,method = 'kpno',imgdir = './plots/'):
		log = logging.getLogger('MK_ULTRA')
		if method == 'cosmos':
			# for my own test on cosmos
			filesets = self.test(count, expid, date, skipid)
		if method == 'kpno':
			#using data stored on kpno
			filesets = self.official(count, expid, date, skipid)
		self.filelist = np.array([filename_i for filename_i in filesets],dtype=np.str).flatten()
		log.debug('input filenames dirs are:%s' %self.filelist)

	def test(self, count,expid,date, skipid):
		log = logging.getLogger('MK_ULTRA')
		if skipid is None:
			filelist = ['expid%08ddate%s.fits' %(i+expid, date) for i in range(count)]
		else:
			filelist = []
			for expid_i in [i+expid for i in range(count)]:
				if expid_i not in np.array(skipid).astype(np.int):
					filelist.append('expid%08ddate%s.fits' %(expid_i, date))
		log.debug('input filenames are:%s' %filelist)
		return filelist

	def official(self, count,expid, date, skipid):
		#still need information of how the official one works
		log = logging.getLogger('MK_ULTRA')
		if skipid is None:
			filelist = ['/exposures/desi/{0}/{1:08d}/ci-{2:08d}.fits.fz'.format(date, expid+i, expid+i) for i in range(count)]
		else: 
			filelist = []
			for expid_i in [i+expid for i in range(count)]:
				if expid_i not in np.array(skipid).astype(np.int):
					filelist.append('/exposures/desi/{0}/{1:08d}/ci-{2:08d}.fits.fz'.format(date, expid_i, expid_i))
				log.debug('input filenames are:%s' %filelist)
		return filelist


class star_sizes(object):
	'''
	AJR editing this to get sizes of stars for particular
	copying the seeing class and re-purposing it
	'''
	def __init__(self, pix_scale=9,expid=0):
		self.pix_scale = pix_scale
		self.micron_to_seeing = {'CIC': 0.01418, 'CIS': 0.013665, 'CIW':0.013665, 'CIE': 0.013665, 'CIN': 0.013665}
		self.expid = expid    		

	def plot_sources(self, data, camera, maxp=50000.,nplot=10,background=None, noise=None, sources=None, model='simp', det_sigma=5.0, base_fwhm=8.0,plot_stuff=False,expid=None):
		'''
		plot light distribution for top ten sources
		'''
		log = logging.getLogger('star_sizes')
		from astropy.modeling import models, fitting
		if background == None or noise == None:
			log.debug("Finding background level...")
			mean, median, std = sigma_clipped_stats(data, sigma=5.0, iters=2)
			background = median
			noise = std
			print(background, noise)
		global cic_count
		global   cis_count
		global   cin_count
		global   ciw_count
		global   cie_count
		# Detect sources
		#print('base fwhm'+str(base_fwhm))
		if sources == None:
			log.debug("Detecting sources..")
			daofind = DAOStarFinder(fwhm=base_fwhm, threshold=det_sigma * noise)
			sources = daofind(data - background)
			print('found '+str(len(sources))+' sources')			
		# Sort by SN
		w = sources['peak'] < maxp
		sources = sources[w]
		sources.sort('flux')
		sources.reverse()

		fwhm_array = []
		fwhm_weight_array = []
		x_coordinate=[]
		y_coordinate=[]
		final_func_val=[]
		if len(sources)>nplot:
			sources = sources[:nplot]
		print(len(sources))
		model_init = models.Gaussian1D(amplitude=1.0, mean=0, stddev=base_fwhm)
		for isource in sources:

			flux_counts = []
			pixel_distance = []
			npix = 0

			x_cen = int(isource['xcentroid'])
			y_cen = int(isource['ycentroid'])
			##print(x_cen,y_cen)
			analysis_radius = 20
			print(isource['peak'])
			try:
				for x in range(x_cen - analysis_radius, x_cen + analysis_radius):
					for y in range(y_cen - analysis_radius, y_cen + analysis_radius):
						#flux_counts.append((data[y][x] - 0.) / isource['peak'])
						val = (data[y][x] - background) / isource['peak']
						
						flux_counts.append(val)
						dis = numpy.linalg.norm((isource['xcentroid'] - x, isource['ycentroid'] - y))
						pixel_distance.append(dis)
				flux_counts = np.array(flux_counts)
				pixel_distance = np.array(pixel_distance)
				w = (flux_counts > 0.5) #& (flux_counts < 0.65)
				psel = pixel_distance[w]
				fsel = flux_counts[w]
				mpsel = np.mean(psel)
				wn = psel < 2.*mpsel
				psel = psel[wn]
				fsel = fsel[wn]
				A = np.vstack([fsel, np.ones(len(psel))]).T
				m, c = np.linalg.lstsq(A, psel, rcond=None)[0]
				ssq = sum((psel-m*fsel+c)**2.)/float(len(psel))
				fitter = fitting.SimplexLSQFitter()
				fitted_model = fitter(model_init, psel, fsel) 

				iFWHM = 2.355 * fitted_model.stddev * self.pix_scale * self.micron_to_seeing[camera]
				print(np.mean(psel),m*.5+c,ssq,iFWHM,fitter.fit_info['final_func_val'])
				plt.plot(pixel_distance,flux_counts,'k-')
				plt.plot(psel,fsel,'ro')
				plt.show()		
			except:
				print('hit edge?')

	def write_sources(self, data, camera, maxp=50000.,nplot=100,background=None, noise=None, sources=None, model='simp', det_sigma=5.0, base_fwhm=8.0,plot_stuff=False,expid=None):
		'''
		plot light distribution for top ten sources
		'''
		log = logging.getLogger('star_sizes')
		from astropy.modeling import models, fitting
		if background == None or noise == None:
			log.debug("Finding background level...")
			mean, median, std = sigma_clipped_stats(data, sigma=5.0, iters=2)
			background = median
			noise = std
			print(background, noise)
		global cic_count
		global   cis_count
		global   cin_count
		global   ciw_count
		global   cie_count
		# Detect sources
		#print('base fwhm'+str(base_fwhm))
		if sources == None:
			log.debug("Detecting sources..")
			daofind = DAOStarFinder(fwhm=base_fwhm, threshold=det_sigma * noise)
			sources = daofind(data - background)
			print('found '+str(len(sources))+' sources')			
		# Sort by SN
		w = sources['peak'] < maxp
		sources = sources[w]
		sources.sort('flux')
		sources.reverse()

		fwhm_array = []
		fwhm_weight_array = []
		x_coordinate=[]
		y_coordinate=[]
		final_func_val=[]
		if len(sources)>nplot:
			sources = sources[:nplot]
		print(len(sources))
		model_init = models.Gaussian1D(amplitude=1.0, mean=0, stddev=base_fwhm)
		fo = open(camera.lower()+'_star_halffwhm%s.txt' %self.expid,'w')
		fo.write('#X_pix Y_pix fwhm ssq peak_flux\n')

		for isource in sources:

			flux_counts = []
			pixel_distance = []
			npix = 0

			x_cen = int(isource['xcentroid'])
			y_cen = int(isource['ycentroid'])
			##print(x_cen,y_cen)
			analysis_radius = 20
			print(isource['peak'])
			try:
				for x in range(x_cen - analysis_radius, x_cen + analysis_radius):
					for y in range(y_cen - analysis_radius, y_cen + analysis_radius):
						#flux_counts.append((data[y][x] - 0.) / isource['peak'])
						val = (data[y][x] - background) / isource['peak']
						
						flux_counts.append(val)
						dis = numpy.linalg.norm((isource['xcentroid'] - x, isource['ycentroid'] - y))
						pixel_distance.append(dis)

				flux_counts = np.array(flux_counts)
				pixel_distance = np.array(pixel_distance)

				w = (flux_counts > 0.5) #& (flux_counts < 0.65)
				psel = pixel_distance[w]
				fsel = flux_counts[w]
				mpsel = np.mean(psel)
				wn = psel < 2.*mpsel
				psel = psel[wn]
				fsel = fsel[wn]
				A = np.vstack([fsel, np.ones(len(psel))]).T
				m, c = np.linalg.lstsq(A, psel, rcond=None)[0]
				ssq = sum((psel-m*fsel+c)**2.)/float(len(psel))
				fitter = fitting.SimplexLSQFitter()
				fitted_model = fitter(model_init, psel, fsel) 

				iFWHM = 2.355 * fitted_model.stddev * self.pix_scale * self.micron_to_seeing[camera]
				print(np.mean(psel),m*.5+c,ssq,iFWHM,fitter.fit_info['final_func_val'])

				fo.write(str(x_cen)+' '+str(y_cen)+' '+str(iFWHM)+ ' '+str(fitter.fit_info['final_func_val'])+' '+str(isource['peak'])+'\n')
				#plt.plot(pixel_distance,flux_counts,'k-')
				#plt.plot(psel,fsel,'ro')
				#plt.show()		
			except:
				print('hit edge?')

		fo.close()
		return True

	def return_sources(self, data, camera, maxp=50000.,nplot=1000,background=None, noise=None, sources=None, model='simp', det_sigma=5.0, base_fwhm=8.0,plot_stuff=False,expid=None):
		'''
		plot light distribution for top ten sources
		'''
		log = logging.getLogger('star_sizes')
		from astropy.modeling import models, fitting
		if background == None or noise == None:
			log.debug("Finding background level...")
			mean, median, std = sigma_clipped_stats(data, sigma=5.0, iters=2)
			background = median
			noise = std
			#print(background, noise)
		global cic_count
		global   cis_count
		global   cin_count
		global   ciw_count
		global   cie_count
		# Detect sources
		#print('base fwhm'+str(base_fwhm))
		if sources == None:
			log.debug("Detecting sources..")
			daofind = DAOStarFinder(fwhm=base_fwhm, threshold=det_sigma * noise)
			sources = daofind(data - background)
			#print('found '+str(len(sources))+' sources')			
		# Sort by SN
		w = sources['peak'] < maxp
		sources = sources[w]
		sources.sort('flux')
		sources.reverse()

		fwhm_array = []
		fwhm_weight_array = []
		x_coord=[]
		y_coord=[]
		final_func_val=[]
		if len(sources)>nplot:
			sources = sources[:nplot]
		#print(len(sources))
		model_init = models.Gaussian1D(amplitude=1.0, mean=0, stddev=base_fwhm)

		for isource in sources:

			flux_counts = []
			pixel_distance = []
			npix = 0

			x_cen = int(isource['xcentroid'])
			y_cen = int(isource['ycentroid'])
			##print(x_cen,y_cen)
			analysis_radius = 20
			#print(isource['peak'])
			try:
				for x in range(x_cen - analysis_radius, x_cen + analysis_radius):
					for y in range(y_cen - analysis_radius, y_cen + analysis_radius):
						#flux_counts.append((data[y][x] - 0.) / isource['peak'])
						val = (data[y][x] - background) / isource['peak']
						
						flux_counts.append(val)
						dis = numpy.linalg.norm((isource['xcentroid'] - x, isource['ycentroid'] - y))
						pixel_distance.append(dis)

				flux_counts = np.array(flux_counts)
				pixel_distance = np.array(pixel_distance)

				w = (flux_counts > 0.5) #& (flux_counts < 0.65)
				psel = pixel_distance[w]
				fsel = flux_counts[w]
				mpsel = np.mean(psel)
				wn = psel < 2.*mpsel
				psel = psel[wn]
				fsel = fsel[wn]
				#A = np.vstack([fsel, np.ones(len(psel))]).T
				#m, c = np.linalg.lstsq(A, psel, rcond=None)[0]
				#ssq = sum((psel-m*fsel+c)**2.)/float(len(psel))
				fitter = fitting.SimplexLSQFitter()
				fitted_model = fitter(model_init, psel, fsel) 

				iFWHM = 2.355 * fitted_model.stddev * self.pix_scale * self.micron_to_seeing[camera]
				x_coord.append(x_cen)
				y_coord.append(y_cen)
				fwhm_array.append(iFWHM)
				final_func_val.append(fitter.fit_info['final_func_val'])
				#print(np.mean(psel),m*.5+c,ssq,iFWHM,fitter.fit_info['final_func_val'])

				#fo.write(str(x_cen)+' '+str(y_cen)+' '+str(iFWHM)+ ' '+str(fitter.fit_info['final_func_val'])+' '+str(isource['peak'])+'\n')
				#plt.plot(pixel_distance,flux_counts,'k-')
				#plt.plot(psel,fsel,'ro')
				#plt.show()		
			except:
				#print('hit edge?')
				pass

		#fo.close()
		return np.array(x_coord),np.array(y_coord),np.array(fwhm_array),np.array(final_func_val)

	def goodimage(self, data, camera,psfx=1.3, fitmax=.2,xrange=[1024,2048],yrange=[512,1536],maxp=50000.,nplot=100,background=None, noise=None, sources=None, model='simp', det_sigma=5.0, base_fwhm=8.0,plot_stuff=False,expid=None):
		'''
		test if an image has a star with a good psf in its center
		'''
		log = logging.getLogger('star_sizes')
		from astropy.modeling import models, fitting
		if background == None or noise == None:
			log.debug("Finding background level...")
			mean, median, std = sigma_clipped_stats(data, sigma=5.0, iters=2)
			background = median
			noise = std
			print(background, noise)
		global cic_count
		global   cis_count
		global   cin_count
		global   ciw_count
		global   cie_count
		# Detect sources
		#print('base fwhm'+str(base_fwhm))
		if sources == None:
			log.debug("Detecting sources..")
			daofind = DAOStarFinder(fwhm=base_fwhm, threshold=det_sigma * noise)
			sources = daofind(data - background)
			#print('found '+str(len(sources))+' sources')			
		# Sort by SN
		w = sources['peak'] < maxp
		sources = sources[w]
		sources.sort('flux')
		sources.reverse()

		fwhm_array = []
		fwhm_weight_array = []
		x_coordinate=[]
		y_coordinate=[]
		final_func_val=[]
		if len(sources)>nplot:
			sources = sources[:nplot]
		#print(len(sources))
		model_init = models.Gaussian1D(amplitude=1.0, mean=0, stddev=base_fwhm)

		for isource in sources:

			flux_counts = []
			pixel_distance = []
			npix = 0

			x_cen = int(isource['xcentroid'])
			y_cen = int(isource['ycentroid'])
			##print(x_cen,y_cen)
			if x_cen >= xrange[0] and x_cen < xrange[1] and y_cen >= yrange[0] and y_cen < yrange[1]:
				analysis_radius = 20
				#print(isource['peak'])
				try:
					for x in range(x_cen - analysis_radius, x_cen + analysis_radius):
						for y in range(y_cen - analysis_radius, y_cen + analysis_radius):
							#flux_counts.append((data[y][x] - 0.) / isource['peak'])
							val = (data[y][x] - background) / isource['peak']
						
							flux_counts.append(val)
							dis = numpy.linalg.norm((isource['xcentroid'] - x, isource['ycentroid'] - y))
							pixel_distance.append(dis)

					flux_counts = np.array(flux_counts)
					pixel_distance = np.array(pixel_distance)

					w = (flux_counts > 0.5) #& (flux_counts < 0.65)
					psel = pixel_distance[w]
					fsel = flux_counts[w]
					mpsel = np.mean(psel)
					wn = psel < 2.*mpsel
					psel = psel[wn]
					fsel = fsel[wn]
					A = np.vstack([fsel, np.ones(len(psel))]).T
					m, c = np.linalg.lstsq(A, psel, rcond=None)[0]
					ssq = sum((psel-m*fsel+c)**2.)/float(len(psel))
					fitter = fitting.SimplexLSQFitter()
					fitted_model = fitter(model_init, psel, fsel) 

					iFWHM = 2.355 * fitted_model.stddev * self.pix_scale * self.micron_to_seeing[camera]
					#print(np.mean(psel),m*.5+c,ssq,iFWHM,fitter.fit_info['final_func_val'])
					if iFWHM < psfx and fitter.fit_info['final_func_val'] < fitmax:
						return True
				except:
					pass
					#print('hit edge?')
					
		return False


	def estimate_fwhm(self,data, camera, nstar=10, background=None, noise=None, sources=None, model='simp', maxfitv=1.,det_sigma=5.0, base_fwhm=8.0,plot_stuff=False,expid=None):
		'''estimate fwhm for a single image
		plot_stuff:plot the fit for each single source
		return:fwhm of this image'''
		log = logging.getLogger('star_sizes')
		
		if background == None or noise == None:
			log.debug("Finding background level...")
			mean, median, std = sigma_clipped_stats(data, sigma=5.0, iters=2)
			background = median
			noise = std
			print(background, noise)
		global cic_count
		global   cis_count
		global   cin_count
		global   ciw_count
		global   cie_count
		# Detect sources
		#print('base fwhm'+str(base_fwhm))
		if sources == None:
			log.debug("Detecting sources..")
			daofind = DAOStarFinder(fwhm=base_fwhm, threshold=det_sigma * noise)
			sources = daofind(data - background)
			print('found '+str(len(sources))+' sources')			
		# Sort by SN
		sources.sort('flux')
		sources.reverse()

		fwhm_array = []
		fwhm_weight_array = []
		x_coordinate=[]
		y_coordinate=[]
		final_func_val=[]
		if len(sources)>nstar:
			sources = sources[:nstar]
		model_init = models.Gaussian1D(amplitude=1.0, mean=0, stddev=base_fwhm)
		for isource in sources:
			flux_counts = []
			pixel_distance = []
			npix = 0

			x_cen = int(isource['xcentroid'])
			y_cen = int(isource['ycentroid'])
			##print(x_cen,y_cen)
			analysis_radius = 20
			#print(isource['peak'])
			try:
				for x in range(x_cen - analysis_radius, x_cen + analysis_radius):
					for y in range(y_cen - analysis_radius, y_cen + analysis_radius):
						#flux_counts.append((data[y][x] - 0.) / isource['peak'])
						val = (data[y][x] - background) / isource['peak']
						
						flux_counts.append(val)
						dis = numpy.linalg.norm((isource['xcentroid'] - x, isource['ycentroid'] - y))
						pixel_distance.append(dis)

				flux_counts = np.array(flux_counts)
				pixel_distance = np.array(pixel_distance)

				w = (flux_counts > 0.5) #& (flux_counts < 0.65)
				psel = pixel_distance[w]
				fsel = flux_counts[w]
				mpsel = np.mean(psel)
				wn = psel < 3.*mpsel
				psel = psel[wn]
				fsel = fsel[wn]
				#A = np.vstack([fsel, np.ones(len(psel))]).T
				#m, c = np.linalg.lstsq(A, psel, rcond=None)[0]
				#ssq = sum((psel-m*fsel+c)**2.)/float(len(psel))
				fitter = fitting.SimplexLSQFitter()
				fitted_model = fitter(model_init, psel, fsel) 
				fitv = fitter.fit_info['final_func_val']
				iFWHM = 2.355 * fitted_model.stddev * self.pix_scale * self.micron_to_seeing[camera]
				#if iFWHM/(m*.5+c) > 1:
				#	fitv = 999
				#print(np.mean(psel),m*.5+c,ssq,iFWHM,fitv)

				#fo.write(str(x_cen)+' '+str(y_cen)+' '+str(iFWHM)+ ' '+str(fitter.fit_info['final_func_val'])+' '+str(isource['peak'])+'\n')
				#plt.plot(pixel_distance,flux_counts,'k-')
				#plt.plot(psel,fsel,'ro')
				#plt.show()		
						#if val > 0.5:
						#	print(dis)
				log.debug("\n>>>>> New source")
			
				log.debug(("Fit value:",  fitter.fit_info['final_func_val']))
				log.debug(("SN:", isource['flux'] * det_sigma))
				log.debug(("single FWHM estimated: ", iFWHM))

				#if fitv < maxfitv:
					#print('final value less than 100, info going into arrays')
				x_coordinate.append(x_cen)
				y_coordinate.append(y_cen)
				fwhm_array.append(iFWHM)
				final_func_val.append(fitv)
				#print(iFWHM,fitter.fit_info['final_func_val'])
				fwhm_weight_array.append(isource['flux'] * det_sigma / fitter.fit_info['final_func_val'])
				#else:
					#print('final value greater than 100, no info for you')
					##print(fitter.fit_info['final_func_val'])
					#if fitter.fit_info['final_func_val'] < 5.0:
					#   fwhm_array.append(iFWHM)
					#   fwhm_weight_array.append(isource['flux'] * det_sigma / fitter.fit_info['final_func_val'])
				#	print(fitv)
				#	color = 'red'
			except:
				#print('exception')
				pass

		if np.min(final_func_val) > maxfitv:
			maxfitv = 1.01*np.min(final_func_val) #ensure at least one point gets fit										
		if len(fwhm_array) > 0:
			print(str(len(fwhm_array))+' stars found with some fwhm')
			#thredshold = np.array(fwhm_weight_array).max()*0.1
			fwhm_array_new = np.array(fwhm_array)[np.array(final_func_val)<maxfitv]
			print(str(len(fwhm_array_new))+' stars found with accepted fwhm')
			fwhm_weight_array_new = np.array(fwhm_weight_array)[np.array(final_func_val)<maxfitv]
			x_coordinate_new = np.array(x_coordinate)[np.array(final_func_val)<maxfitv]
			y_coordinate_new = np.array(y_coordinate)[np.array(final_func_val)<maxfitv]

			final_func_val_new = np.array(final_func_val)[np.array(final_func_val)<maxfitv]
			#FWHM_estimation = numpy.average(fwhm_array_new, weights=fwhm_weight_array_new)
			FWHM_estimation = numpy.median(fwhm_array_new)
			log.info('Final FWHM estimated is %.2f"' % FWHM_estimation)
 
		if len(fwhm_array) > 0:     
		  return FWHM_estimation,len(fwhm_array_new)
		else:
		  return FWHM_estimation,None




def plotpsf(exp,cam,fitmax=.1):
	from matplotlib import pyplot as plt
	import numpy as np
	plt.xlim(0,3072)
	plt.ylim(0,2048)
	d = np.loadtxt(cam+'_star_halffwhm'+exp+'.txt').transpose()
	xl = []
	yl = []
	sl = []
	for i in range(0,len(d[0])):
		#print(i,d[-2][i])
		if d[-2][i] < fitmax:
	 		xl.append(d[0][i])
	 		yl.append(d[1][i])
	 		sl.append(d[2][i])
	print(xl)
	plt.scatter(xl,yl,c=sl,s=10)
	plt.colorbar()
	plt.title(cam)
	plt.show()
	

def main(args=None):
    time0 = time.time()

    if args is None:
        # Read from cmd line
        parser= get_parser()
        args = parser.parse_args(args=args)
    if args is None:
        ##print('ERROR! you need parameters to make this work!')
        return None

    if args.verbose:
        lvl = logging.DEBUG
    else:
        lvl = logging.INFO
    logging.basicConfig(level=lvl, stream=sys.stdout)
    log = logging.getLogger('MK_ULTRA')

    if args.mode == 'FOCUS':
        if args.count<3:
            log.error('number of counts MUST be greater than 3!')
            return None
        log = logging.getLogger('FOCUS')
    elif args.mode == 'alignment':
        if args.count<3:
             log.error('number of counts MUST be greater than 3!')
             return None
        log = logging.getLogger('alignment')
    elif args.mode == 'seeing':
        if args.count!=1:
           log.error('seeing mode accept ONE image at a time. set count=1')
           return None
        log = logging.getLogger('seeing')
    elif args.mode == 'star_sizes':
        if args.count!=1:
           log.error('star_sizesw mode accept ONE image at a time. set count=1')
           return None
        log = logging.getLogger('star_sizes')

    fns = filesystem(args.count,args.expid,args.obsday,args.skipid)  
    if args.mode == 'seeing':
         target = star_sizes(pix_scale=9)
         if args.camera == 'ALL':
             cameras = ['CIN','CIW','CIC','CIE','CIS']
         else:
             cameras = [args.camera]
         fwhm_dict = dict()
         for camera in cameras:
             log.info('Measuring FWHM on %s' %camera)
             fwhm_estimation,nsource = target.estimate_fwhm(fits.getdata(fns.filelist[0],extname = camera), camera = camera)
             fwhm_dict.update({camera:fwhm_estimation})
         log.info(fwhm_dict)
         return None
    if args.mode == 'star_sizes':
         target = star_sizes(pix_scale=9,expid=args.expid)
         if args.camera == 'ALL':
             cameras = ['CIN','CIW','CIC','CIE','CIS']
         else:
             cameras = [args.camera]
         fwhm_dict = dict()
         for camera in cameras:
             log.info('Measuring FWHM on %s' %camera)
             target.plot_sources(fits.getdata(fns.filelist[0],extname = camera), camera = camera)
             #fwhm_dict.update({camera:fwhm_estimation})
         #log.info(fwhm_dict)
         return None
 
    if args.mode == 'star_sizes_plot':
         target = star_sizes(pix_scale=9,expid=args.expid)
         if args.camera == 'ALL':
             cameras = ['CIN','CIW','CIC','CIE','CIS']
         else:
             cameras = [args.camera]
         for camera in cameras:
             log.info('Plotting stars %s' %camera)
             target.write_sources(fits.getdata(fns.filelist[0],extname = camera), camera = camera)

         return None

 
    if args.mode == 'FOCUS':
         if args.noplot == False:
              target = Focus(pix_scale=9, cameras = args.camera, plot = True)
         else:
              target = Focus(pix_scale=9, cameras = args.camera, plot = False)
         target.process_FOCUS_set(fns.filelist)
         #target.sum_plot()
    if args.mode == 'alignment':
           if args.noplot == False:
                target = alignment(pix_scale=9,plot=True)
           else:
                target = alignment(pix_scale=9,plot=False)
           target.alignment(fns.filelist)
  
def findgoodimages(date,filemin,filemax,camera='CIC'):
	target = star_sizes(pix_scale=9)
	fo = open('goodCICimages'+date+'.dat','w')
	for i in range(filemin,filemax):
		filename ='/exposures/desi/'+date+'/0000'+str(i)+'/ci-0000'+str(i)+'.fits.fz'
		target.expid=i
		try:			
			h = fitsio.read_header(filename,ext='CIC')
			prog = h['PROGRAM'].split()[0]
			if prog != 'focus' and prog != 'alignment' and prog != 'Alignment':		
				if target.goodimage(fits.getdata(filename,extname = camera), camera = camera):    			
					print(i)
					fo.write(str(i)+'\n')
			else:
				print(h['PROGRAM'])		
		except:		 
			pass

def mkimagequalitymap_onenight(date,fitmax=0.2):
	imagel = np.loadtxt('goodCICimages'+date+'.dat')
	cameras = ['CIN','CIW','CIC','CIE','CIS']
	target = star_sizes(pix_scale=9)
	rl = np.zeros((11,13))
	nl = np.zeros((11,13))
	for image in imagel:
		filen = str(int(image))
		print(filen)
		filename ='/exposures/desi/'+date+'/0000'+str(filen)+'/ci-0000'+str(filen)+'.fits.fz'
		#first get norm
		camera = 'CIC'
		xl,yl,psfl,fitl = target.return_sources(fits.getdata(filename,extname = camera), camera = camera)
		w = (xl > 1024) & (xl <= 2048) & (yl > 512) & (yl <= 1536) & (fitl < fitmax)
		norm = np.median(psfl[w])
		print(norm,np.mean(psfl[w]))
		try:
			for camera in cameras:
				if camera == 'CIN':
					xp = 5
					yp = 8
				if camera == 'CIW':
					xp = 0
					yp = 4
				if camera == 'CIC':
					xp = 5
					yp = 4
				if camera == 'CIE':
					xp = 10
					yp = 4
				if camera == 'CIS':
					xp = 5
					yp = 0
				xl,yl,psfl,fitl = target.return_sources(fits.getdata(filename,extname = camera), camera = camera)
				wfit = fitl < fitmax
				xl = xl[wfit]
				yl = yl[wfit]
				psfl = psfl[wfit]/norm
				xind = (xl/1024).astype(int)
				yind = (yl/1024).astype(int)
				for i in range(0,3):
					for j in range(0,2):
						w = (xind == i) & (yind == j)
						mpsf = np.median(psfl[w])
						if 0*mpsf == 0 and len(psfl[w]) > 2:
							rl[j+yp][i+xp] += mpsf
							nl[j+yp][i+xp] += 1.
							print(mpsf,i+xp,j+yp,len(psfl[w]))
		except:
			print(camera)
	fo = open('imqual'+date+'.dat','w')
	
	print(rl)
	print(nl)
	print(rl/nl)
	ql = rl/nl
	for i in range(0,len(ql)):
		for j in range(0,len(ql[i])):
			fo.write(str(ql[i][j])+' ')
		fo.write('\n')	
	fo.close()
	plt.imshow(rl/nl,origin='lower')
	plt.colorbar()
	plt.show()
	return True

def mkimagequalitymap_oneimage(date,image,fitmax=0.2):
	imagel = np.loadtxt('goodCICimages'+date+'.dat')
	cameras = ['CIN','CIW','CIC','CIE','CIS']
	target = star_sizes(pix_scale=9)
	rl = np.zeros((11,13))
	nl = np.zeros((11,13))
	filen = str(image)
	print(filen)
	filename ='/exposures/desi/'+date+'/0000'+str(filen)+'/ci-0000'+str(filen)+'.fits.fz'
	#first get norm
	camera = 'CIC'
	xl,yl,psfl,fitl = target.return_sources(fits.getdata(filename,extname = camera), camera = camera)
	w = (xl > 1024) & (xl <= 2048) & (yl > 512) & (yl <= 1536) & (fitl < fitmax)
	norm = np.median(psfl[w])
	print(norm,np.mean(psfl[w]))
	for camera in cameras:
		if camera == 'CIN':
			xp = 5
			yp = 8
		if camera == 'CIW':
			xp = 0
			yp = 4
		if camera == 'CIC':
			xp = 5
			yp = 4
		if camera == 'CIE':
			xp = 10
			yp = 4
		if camera == 'CIS':
			xp = 5
			yp = 0
		xl,yl,psfl,fitl = target.return_sources(fits.getdata(filename,extname = camera), camera = camera)
		wfit = fitl < fitmax
		xl = xl[wfit]
		yl = yl[wfit]
		psfl = psfl[wfit]/norm
		xind = (xl/1024).astype(int)
		yind = (yl/1024).astype(int)
		for i in range(0,3):
			for j in range(0,2):
				w = (xind == i) & (yind == j)
				mpsf = np.median(psfl[w])
				if 0*mpsf == 0 and len(psfl[w]) > 2:
					rl[j+yp][i+xp] += mpsf
					nl[j+yp][i+xp] += 1.
					print(mpsf,i+xp,j+yp,len(psfl[w]))
	fo = open('imqual'+date+str(image)+'.dat','w')
	
	print(rl)
	print(nl)
	print(rl/nl)
	ql = rl/nl
	for i in range(0,len(ql)):
		for j in range(0,len(ql[i])):
			fo.write(str(ql[i][j])+' ')
		fo.write('\n')	
	fo.close()
	plt.imshow(rl/nl,origin='lower')
	plt.colorbar()
	plt.show()
	return True


if __name__ == '__main__':
	print('analysis started at %s' % time.strftime("%Y-%m-%d %H:%M:%S"))
	import sys
	
	from matplotlib import pyplot as plt
	date = '20190409'

	filen = int(sys.argv[3])
	
	if sys.argv[1] == 'seeing':
		target = star_sizes(pix_scale=9)
		target.expid=filen

	if sys.argv[2] == 'ALL':
		cameras = ['CIN','CIW','CIC','CIE','CIS']
	else:
		cameras = sys.argv[2]

	for camera in cameras:
		focusl = []
		fwl = []
		filename ='/exposures/desi/'+date+'/0000'+str(filen)+'/ci-0000'+str(filen)+'.fits.fz'
		h = fitsio.read_header(filename,ext=1)
		focus = h['FOCUS'][2]
		target.write_sources(fits.getdata(filename,extname = camera), camera = camera)
		plotpsf(str(filen),camera.lower(),fitmax=.2)
		#plt.title(camera)	
		#plt.show()
	#plt.ylim(0.2,3)
	#plt.legend(cameras)
	#plt.show()	
     
     
	print('analysis finished at %s' % time.strftime("%Y-%m-%d %H:%M:%S"))
