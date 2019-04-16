import numpy as np
from matplotlib import pyplot as plt
focusdat = np.loadtxt('focus_temp_AJR.txt').transpose()

def linfit(cam,bad_date=20190409):
	fv = focusdat[cam]
	print(fv)
	temp = focusdat[-2]
	w = focusdat[-1] != bad_date
	A = np.vstack([temp[w], np.ones(len(temp[w]))]).T
	m, c = np.linalg.lstsq(A, fv[w], rcond=None)[0]
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
	plt.plot(temp[w],fv[w],'ko')
	plt.plot(temp[~w],fv[~w],'ro')
	plt.title(camt)
	plt.xlabel('mean truss temperature (C)')
	plt.ylabel('focus value (um)')
	plt.text(10,-8200,'least squares fit '+str(np.round(m,1))+'T '+str(np.round(c,0)))
	plt.text(10,-8300,'value at 7 deg. '+str(np.round(m*7+c,1)))
	#plt.show()
	plt.savefig(camt+'bestfit.png')
	plt.clf()
	plt.plot(temp,np.ones(len(temp)),'k--')
	plt.plot(temp[w],fv[w]-(m*temp[w]+c),'ko')
	rms = np.std(fv[w]-(m*temp[w]+c))
	print(rms)
	plt.plot(temp[~w],fv[~w]-(m*temp[~w]+c),'ro')
	plt.title(camt+' RMS ='+str(np.round(rms,1)))
	plt.xlabel('mean truss temperature (C)')
	plt.ylabel('residual focus value from linear fit (um)')
	plt.savefig(camt+'subfit.png')
	#plt.show()

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

	
if __name__ == '__main__':
	linfit_all()
	for i in range(1,6):
		linfit(i)