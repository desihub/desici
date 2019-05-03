import numpy as np
from matplotlib import pyplot as plt
import fitsio
from operator import itemgetter
from math import *
import scipy
from scipy.interpolate import interp1d
from astropy.table import Table
focusdat = np.loadtxt('focus_temp_AJR.txt').transpose()

def mktempfits():
	focusdatn = np.loadtxt('focus_temp_AJR_new.txt').transpose()
	#make table from txt file
	focustab = Table([focusdatn[0],focusdatn[1],focusdatn[2],focusdatn[3],focusdatn[4],focusdatn[6],focusdatn[7],focusdatn[8],focusdatn[9],focusdatn[10],focusdatn[5].astype(int),focusdatn[11].astype(int)],names=('focus_N','focus_W', 'focus_C', 'focus_E', 'focus_S', 'fwhm_N', 'fwhm_W', 'fwhm_C', 'fwhm_E', 'fwhm_S', 'expmin','date'))
	#make empty arrays for everything to extract from headers, using array of quantities
	hquant = ['TARGTRA','TARGTDEC','MOUNTHA','MOUNTAZ','MOUNTEL','AIRMASS','CI-T1','CI-T2','CI-T3','CI-T4','CI-T5','TDEWPNT','TAIRFLOW','TAIRITMP','TAIROTMP','TAIRTEMP','TCASITMP','TCASOTMP','TCSITEMP','TCSOTEMP','TDBTEMP','TPMNIBT','TPMEOBT','TTRSTEMP','TTRWTEMP','TTRUETBT','TTRUETTT','TTRUNTBT','TTRUNTTT','TTRUSTBT','TTRUSTST','TTRUSTTT','TTRUTSBT','TTRUTSMT','TTRUTSTT','TTRUWTBT','TTRUWTTT','AMNIENTN','AMBIENTS','OUTTEMP','TELBASE']
	for q in hquant:
		print(q)
		focustab[q] = np.zeros(len(focusdatn[0]))
		for i in range(0,len(focusdatn[0])):
			night = focustab['date'][i]
			emin = focustab['expmin'][i]
			if emin < 10000:
				zer = '0000'
			else:
				zer = 	'000'
			fl = '/project/projectdirs/desi/spectro/data/'+str(night)+'/'+zer+str(emin)+'/ci-'+zer+str(emin)+'.fits.fz'
			f = fitsio.read_header(fl,ext=1)
			try:
				focustab[q][i] = f[q]
			except:
				print(f[q])	
	focustab.write('focusdata.fits', format='fits', overwrite=True)  
	return True

def writeDT():
	f = fitsio.read('focusdata.fits')
	fo = open('dTdt.dat','w')
	for i in range(0,len(f)):
		dt = getdTdt(f[i]['date'],f[i]['expmin'])
		print(dt)
		fo.write(str(dt)+'\n')
	fo.close()
	return True
			
def dTdt(night,npoints=1001,minexp=1000,maxexp=10000):
	'''
	get dTdt relationship for a given night
	translated from Arjun Dey's IDL

	pro calcdt,date,mjd,dtemp,plot=plot,npts=npts

	; calculates the rate of change of truss temperature  with time 
	; (in units of deg celcius per min) for a given night. Truss temps 
	; and MJD are taken from image headers and then fit using a polynomial. 
	; Derivative is computed from the polynomaial.
	;
	; Inputs: date - string (e.g., '20190904')
	; Keywords : npts - number of points for which polynomial and derivative are computed (default 1000)
	;	plot - make a plot (default - don't)
	; Outputs: mjd - double array - MJD-OBS for the night
	; 	   dtemp - float array - derivative of temperature (C/min)
	;
	; - Arjun
	'''
	time = []
	temp = []
	dir = '/project/projectdirs/desi/spectro/data/'
	#dir = '/exposures/desi/'
	for i in range(minexp,maxexp):
		#print(fl)
		try:
			fl = dir+str(night)+'/0000'+str(i)+'/ci-0000'+str(i)+'.fits.fz'
			h = fitsio.read_header(fl,ext=1)
			#print(i)
			ti = np.mean([h['TTRUNTTT'],h['TTRUETBT'],h['TTRUETTT'],h['TTRUNTBT'],h['TTRUSTBT'],h['TTRUSTST'],h['TTRUSTTT'],h['TTRUTSBT'],h['TTRUTSMT'],h['TTRUTSTT'],h['TTRUWTBT'],h['TTRUWTTT']])
			temp.append(ti)
			#print(i)
			time.append(h['MJD-OBS'])
		except:
			pass	
	
	nf=len(time)
	if nf == 0:
		print('no files found!')
		return('ERROR, no files found!')
	time = np.array(time)
	t0 = np.min(time)
	time = (time-t0)*24.*60. # convert to minutes
	cc = np.polyfit(time,temp,5) # fit the variation as a 5th order polynomial
	f_fit = np.poly1d(cc)
	# calculate new x's and y's
	xx = np.linspace(time[0], time[-1], npoints)
	yy = f_fit(xx)

	dyy=np.gradient(yy,xx)

	return(xx,dyy,t0)
	
def getdTdt(night,exp):	
	xx,dyy,t0 = dTdt(night)
	dir = '/project/projectdirs/desi/spectro/data/'
	#dir = '/exposures/desi/'

	fl = dir+str(night)+'/0000'+str(exp)+'/ci-0000'+str(exp)+'.fits.fz'
	h = fitsio.read_header(fl,ext=1)
	time = (h['MJD-OBS']-t0)*24.*60. # convert to minutes
	dT = scipy.interpolate.interp1d(xx, dyy, kind='linear')(time)
	return dT

def linfit(cam,bad_date=20190409):
	fv = focusdat[cam]
	print(fv)
	temp = focusdat[-4]
	w = focusdat[-3] != bad_date
	az = focusdat[-2]
	el = focusdat[-1]
	try:
		dtl = np.loadtxt('dTdtfocus.dat').transpose()
	except:
		dtl = []
		fo = open('dTdtfocus.dat','w')
		for i in range(0,len(temp)):
			date = int(focusdat[-3][i])
			exp = int(focusdat[6][i])
			print(date,exp)
			dt = getdTdt(date,exp)
			dtl.append(dt)
			print(dt)
			fo.write(str(dt)+'\n')
		fo.close()
	print(len(fv),len(temp))
	A = np.vstack([temp[w], np.ones(len(temp[w]))]).T
	m, c = np.linalg.lstsq(A, fv[w], rcond=None)[0]
	B = np.vstack([dtl, np.ones(len(dtl))]).T
	mr, cr = np.linalg.lstsq(B, fv-(m*temp+c), rcond=None)[0]

	if cam == 1:
		camt = 'CIN'
	if cam == 2:
		camt = 'CIW'
	if cam == 3:
		camt = 'CIC'
	if cam == 4:
		camt = 'CIE'
	if cam == 5:
		camt = 'CIS'
		
	print(m,c)
	plt.clf()
	plt.plot(temp,m*temp+c,'k--')
	plt.plot(temp, -8400 + (7.0 - temp) * 110,'r:') 
	
	plt.plot(temp[w],fv[w]-(mr*dtl[w]+cr),'ko')
	plt.plot(temp[~w],fv[~w]-(mr*dtl[~w]+cr),'ro')
	plt.title(camt)
	plt.xlabel('mean truss temperature (C)')
	plt.ylabel('focus value (um)-('+str(np.round(mr,2))+'dT/dt+'+str(np.round(cr,0))+')')
	plt.text(10,-8200,'least squares fit '+str(np.round(m,1))+'T '+str(np.round(c,0)))
	plt.text(10,-8300,'value at 7 deg. '+str(np.round(m*7+c,1)))	
	plt.savefig(camt+'bestfit.png')
	plt.show()
	plt.clf()
	plt.plot(temp,np.ones(len(temp)),'k--')
	plt.plot(temp,fv-(m*temp+c)-(mr*dtl+cr),'ko')
	rms = np.std(fv-(m*temp+c)-(mr*dtl+cr))
	print(rms)
	plt.plot(temp[~w],fv[~w]-(m*temp[~w]+c)-(mr*dtl[~w]+cr),'ro')
	diff = fv- (m*temp+c)
	plt.title(camt+' RMS ='+str(np.round(rms,1)))
	plt.xlabel('mean truss temperature (C)')
	plt.ylabel('residual focus value from linear focus+dT/dt fits (um)')
	plt.savefig(camt+'subfit.png')
	plt.show()
	plt.clf()
	dtl = np.array(dtl)
	plt.plot(dtl,dtl*mr+cr,'k--')
	plt.plot(dtl,diff,'ko')
	plt.plot(dtl[~w],diff[~w],'ro')
	plt.xlabel('dT/dt')
	plt.ylabel('residual focus value from linear fit (um)')
	plt.title(camt)
	plt.savefig(camt+'dTdtfit.png')
	plt.show()
	diff = fv-(m*temp+c)-(mr*dtl+cr)
	plt.plot(az,np.ones(len(temp)),'k--')
	plt.plot(az,diff,'ko')
	#plt.plot(az[w],diff[w],'ro')
	plt.xlabel('mount azimuth (deg)')
	plt.ylabel('residual focus value from linear fit (um)')
	plt.title(camt)
	plt.savefig(camt+'subvsaz.png')
	plt.show()
	plt.clf()
	plt.plot(el,np.ones(len(temp)),'k--')
	plt.plot(el,diff,'ko')
	#plt.plot(el[w],diff[w],'ro')
	plt.xlabel('mount elevation (deg)')
	plt.ylabel('residual focus value from linear fit (um)')
	plt.title(camt)
	plt.savefig(camt+'subvsel.png')
	plt.show()

def linfit_all(bad_date=20190409):
	
	#print(fv)
	temp = focusdat[-2]
	w = focusdat[-1] != bad_date
	fv = np.hstack([focusdat[1][w],focusdat[2][w],focusdat[3][w],focusdat[4][w],focusdat[5][w]])
	temps = np.hstack([temp[w],temp[w],temp[w],temp[w],temp[w]])
	print(temps)
	A = np.vstack([temps, np.ones(len(temps))]).T
	m, c = np.linalg.lstsq(A, fv, rcond=None)[0]
	camt= 'All'		
	print(m,c)
	plt.clf()
	plt.plot(temps,m*temps+c,'k--')
	# -8400 + (7.0 – truss_temp) * 110
	plt.plot(temps, -8400 + (7.0 - temps) * 110,'r:') 
	plt.plot(temps,fv,'ko')
	#plt.plot(temp[~w],fv[~w],'ro')
	plt.title(camt)
	plt.xlabel('mean truss temperature (C)')
	plt.ylabel('focus value (um)')
	plt.text(10,-8200,'least squares fit '+str(np.round(m,1))+'T '+str(np.round(c,0)))
	plt.text(10,-8300,'value at 7 deg. '+str(np.round(m*7+c,1)))
	#plt.show()
	plt.savefig(camt+'bestfit.png')
	plt.clf()
	plt.plot(temps,np.ones(len(temps)),'k--')
	plt.plot(temps,fv-(m*temps+c),'ko')
	#plt.plot(temp[~w],fv[~w]-(m*temp[~w]+c),'ro')
	plt.title(camt)
	plt.xlabel('mean truss temperature (C)')
	plt.ylabel('residual focus value from linear fit (um)')
	plt.savefig(camt+'subfit.png')

def fitstar(star):
	foc = star['focus']
	fwhm = star['fwhm']
	new_length = len(fwhm) * 10

	# Interpolation
	foc_values_interp = np.linspace(min(foc), max(foc), new_length)

	z_fit,residuals, rank, singular_values, rcond = np.polyfit(foc, fwhm, 2,full=True)
	f_fit = np.poly1d(z_fit)
	# calculate new x's and y's
	foc_fit = np.linspace(foc_values_interp[0], foc_values_interp[-1], new_length)
	ccd_fwhms_fit = f_fit(foc_fit)

	# Find minimum of POLYNOMIAL FIT
	min_index = min(enumerate(ccd_fwhms_fit), key=itemgetter(1))[0]
	fwhm_at_foc_fit = ccd_fwhms_fit[min_index]
	foc_pos_fit = foc_fit[min_index]
	return(foc_pos_fit,residuals,foc_values_interp[0],foc_values_interp[-1])
	

def starinfo(bad_date=20190409,resmax=1,rf=0.15):
	exp = focusdat[6]
	date = focusdat[-3]
	temp = focusdat[-4]
	mlin = -111.4 #fit for CIC from data through April 15th
	blin = -7548.
	caml = ['CIN','CIW','CIC','CIE','CIS']
	fo = open('camxyfocusv.dat','w')
	for i in range(0,len(exp)):
		datei = int(date[i])
		if datei != bad_date:
#			cfit = temp[i]*mlin+blin
			f = fitsio.read('thunderstruck/'+str(int(date[i]))+'/'+str(int(exp[i]))+'/ci-star-catalogue-'+str(int(exp[i]))+'-CIC.fits')
			starl = np.unique(f['star_id'])
			cl = []
			for star in starl:
				sel = f['star_id'] == star
				fv,res,minf,maxf = fitstar(f[sel])
				if len(res) > 0 and (fv-minf)/(maxf-minf) > rf and (maxf-fv)/(maxf-minf) < 1-rf:
					if res[0] < resmax:
						cl.append(fv)
			cfit = np.median(cl)
			print(cfit)
			if cfit*0 == 0:
				for cam in caml:
					f = fitsio.read('thunderstruck/'+str(int(date[i]))+'/'+str(int(exp[i]))+'/ci-star-catalogue-'+str(int(exp[i]))+'-'+cam+'.fits')
					starl = np.unique(f['star_id'])
					if cam == 'CIN':
						ind = 0
					if cam == 'CIW':
						ind = 1
					if cam == 'CIC':
						ind = 2
					if cam == 'CIE':
						ind = 3
					if cam == 'CIS':
						ind = 4
				
					for star in starl:
						sel = f['star_id'] == star
						fv,res,minf,maxf = fitstar(f[sel])
						if len(res) > 0 and abs(fv-minf)/abs(maxf-minf) > rf and abs(fv-minf)/abs(maxf-minf)  < 1-rf and abs(maxf-fv)/abs(maxf-minf) < 1-rf and abs(maxf-fv)/abs(maxf-minf) > rf:
							fo.write(str(ind)+' '+str(np.median(f[sel]['x']))+' '+str(np.median(f[sel]['y']))+' '+str(fv-cfit)+' '+str(res[0])+'\n')
							if fv-cfit == 0:
								print(fv,minf,maxf,abs(fv-minf)/abs(maxf-minf),abs(maxf-fv)/abs(maxf-minf))
						#print(fv-cfit)
					#if np.max(f['peak']) > 50000:
					#	print(i,cam,np.max(f['peak'])) 
	fo.close()
	return True
	
def plotstarinfo(cam,pixres=1,resmax=1):
	d = np.loadtxt('camxyfocusv.dat').transpose()
	w = d[0] == cam
	l0 = []
	l1 = []
	l2 = []
	l3 = []
	l4 = []
	l5 = []
	div = int(1024/pixres)
	pixlx = []
	pixly = []
	focl = []
	for i in range(0,len(d[1][w])):
		x = d[1][w][i]
		y = d[2][w][i]
		res = d[4][w][i]
		if res < resmax:
			if x < 1024 and y < 1024:
				l0.append(d[3][w][i])
			if x < 1024 and y >= 1024:
				l1.append(d[3][w][i])
			if x >= 1024 and x < 2048 and y < 1024:
				l2.append(d[3][w][i])
			if x >= 1024 and x < 2048 and y >= 1024:
				l3.append(d[3][w][i])
			if x >= 2048 and y < 1024:
				l4.append(d[3][w][i])
			if  x >= 2048 and y >= 1024:
				l5.append(d[3][w][i])
			pixx = int(x/div)
			pixlx.append(pixx)
			pixy = int(y/div)
			pixly.append(pixy)
			focl.append(d[3][w][i]*1.7)	
	focl = np.array(focl)
	pixlx = np.array(pixlx).astype(int)
	pixly = np.array(pixly).astype(int)
	#print(np.unique(pixlx))
	print(np.median(l0),np.std(l0)/sqrt(float(len(l0))),len(l0))
	plt.hist(l0)
	plt.show()
	print(np.median(l1),np.std(l1)/sqrt(float(len(l1))),len(l1))
	plt.hist(l1)
	plt.show()
	print(np.median(l2),np.std(l2)/sqrt(float(len(l2))),len(l2))
	plt.hist(l2)
	plt.show()
	print(np.median(l3),np.std(l3)/sqrt(float(len(l3))),len(l3))
	plt.hist(l3)
	plt.show()
	print(np.median(l4),np.std(l4)/sqrt(float(len(l4))),len(l4))
	plt.hist(l4)
	plt.show()
	print(np.median(l5),np.std(l5)/sqrt(float(len(l5))),len(l5))
	plt.hist(l5)
	plt.show()

	v0 = np.median(l0)	
	v1 = np.median(l1)	
	v2 = np.median(l2)	
	v3 = np.median(l3)	
	v4 = np.median(l4)	
	v5 = np.median(l5)	
	reg = np.zeros((2*pixres,3*pixres))
	for i in range(0,len(reg)):
		for j in range(0,len(reg[0])):
			w = (pixlx == j)# and (pixly = j)
			w &= (pixly == i)
			print(i,j,len(focl[w]),np.median(focl[w]))
			reg[i][j] = np.median(focl[w])
	#reg[0][0] = v0
	#reg[1][0] = v1
	#reg[0][1] = v2
	#reg[1][1] = v3
	#reg[0][2] = v4
	#reg[1][2] = v5
	plt.imshow(reg,origin='lower',extent=[0,3072,0,2048],cmap='viridis')
	plt.xlabel('x pixel')
	plt.ylabel('y pixel')
	cbar = plt.colorbar()
	cbar.ax.set_ylabel(r'1.7$\Delta$ microns relative to CIC for given focus sequence')
	if cam == 0:
		titl = 'CIN'
	if cam == 1:
		titl = 'CIW'
	if cam == 2:
		titl = 'CIC'
	if cam == 3:
		titl = 'CIE'
	if cam == 4:
		titl = 'CIS'
	plt.title(titl)
	plt.savefig(titl+'focuschange.png')
	plt.show()	

if __name__ == '__main__':
# 	linfit_all()
# 	for i in range(1,6):
# 		linfit(i)
	starinfo()