'''
functions to go from CI pixel to focal plane x,y
plug into desimodel.focalplane to get to sky position
'''

import numpy as np

#coefficients for transform between pixel and CS5 x,y
#X_0	Y_0	Z_0	a	b	c 	d	e	f
#would have been better to do this as a dictionary
Center = [ 14.053,8.85,-0.004,-1.,-0.0085,	0.0085,	-1.,-0.0002,	-0.0001]
North= [-13.457,-404.841,-18.319,1.	-0.0018,0.002,0.9958,-0.0024,0.0914]
South = [12.64,404.366,-18.401,-1.,0.005,-0.005,-0.9956,0.0001,0.0934]
West = [-405.399,12.16,	-18.455,0.0038,	0.9959,	-1.,0.0038,0,0.091]
East =	[405.298,-13.192,-18.43,-0.0056,-0.9957,1.,	-0.0056,0.0008,	0.0924]
pixsize = 0.009

def pixtoCS5(x,y,cam):
	if cam == 'CIC':
		cf = Center
	if cam == 'CIW':
		cf = West
	if cam == 'CIN':
		cf = North
	if cam == 'CIE':
		cf = East
	if cam == 'CIS':
		cf = South
	X_0 = cf[0]
	Y_0 = cf[1]
	Z_0 = cf[2]
	a = cf[3]*pixsize
	b = cf[4]*pixsize
	c = cf[5]*pixsize
	d = cf[6]*pixsize
	e = cf[7]*pixsize
	f = cf[8]*pixsize	
	cx = X_0+a*x+b*y
	cy = Y_0+c*x+d*y
	cz = Z_0+e*x+f*y
	return cx,cy,cz
	
def CS5topix(x,y,cam):
	if cam == 'CIC':
		cf = Center
	if cam == 'CIW':
		cf = West
	if cam == 'CIN':
		cf = North
	if cam == 'CIE':
		cf = East
	if cam == 'CIS':
		cf = South
	X_0 = cf[0]
	Y_0 = cf[1]
	Z_0 = cf[2]
	a = cf[3]*pixsize
	b = cf[4]*pixsize
	c = cf[5]*pixsize
	d = cf[6]*pixsize
	e = cf[7]*pixsize
	f = cf[8]*pixsize	
	xp = 1./(c-d*a/b)*(y-Y_0-d/b*x+d/b*X_0)
	yp = 1./b*(x-X_0-a*xp)
	return xp,yp
