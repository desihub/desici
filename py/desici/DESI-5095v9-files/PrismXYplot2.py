
# PrismXYplot2.py 
# showing how a field-center spot roams about the focal surface
# and now with ACTUAL SKY DATA to compare against theoretical RT170a.py
#
# The 25 central spots are 0, 5, 10, 15, ...
#
# Usual environment:  Python/RT170


import numpy as np
import matplotlib.pyplot as plt


PROGNAME = 'PrismXYplot2.py'

infile = 'RT170a_Nstars_5_SkyXYinferredUV.txt'  
# theory predicted spot locations
# 125 rows = 5 locations * 25 ADC settings
# Nine columns:  ADC1, ADC2, skyU, skyV, Ngood, Xfpmm, Yfpmm, sigmaX, sigmaY

plotfile = 'PrismXYplot2.png'

# These labels & star IDs are from Python/FITS/on_sky_adc/AutoGuideStars.py
camlabels = ['CIN: 714791001685998720 16.79',
             'CIW: 714955241235254656 16.94',
             'CIC: 714073845226040576 14.83', 
             'CIE: 702345732650570624 14.25',
             'CIS: 701824289260562304 15.90']
             
# These focal plane coordinates (mm) are from Python/FITS/on_sky_adc/PlotGuideStars2.py
nsky = np.array([
    [    -2.25,   394.25],
    [    -2.56,   395.64],
    [    -3.49,   396.62],
    [    -4.83,   397.00],
    [    -6.13,   396.62],
    [    -7.09,   395.61],
    [    -7.42,   394.23],
    [    -7.02,   392.83],
    [    -6.05,   391.83],
    [    -4.68,   391.46],
    [    -3.33,   391.81],
    [    -2.37,   392.83],
    [    -2.02,   394.20],
    [    -1.64,   392.87],
    [    -0.63,   391.89],
    [     0.69,   391.54],
    [     2.01,   391.90],
    [     2.96,   392.91],
    [     3.34,   394.26],
    [     2.98,   395.63],
    [     2.05,   396.59],
    [     0.75,   396.95],
    [    -0.52,   396.56],
    [    -1.42,   395.58],
    [    -1.76,   394.22]])
    
wsky = np.array([
    [   397.06,    -8.79],
    [   396.74,    -7.48],
    [   395.75,    -6.53],
    [   394.37,    -6.16],
    [   393.00,    -6.51],
    [   392.00,    -7.50],
    [   391.67,    -8.80],
    [   392.06,   -10.14],
    [   393.07,   -11.12],
    [   394.53,   -11.47],
    [   395.92,   -11.12],
    [   396.93,   -10.16],
    [   397.31,    -8.86],
    [   397.69,   -10.18],
    [   398.77,   -11.08],
    [   400.14,   -11.40],
    [   401.51,   -11.06],
    [   402.49,   -10.08],
    [   402.89,    -8.80],
    [   402.52,    -7.48],
    [   401.55,    -6.54],
    [   400.20,    -6.20],
    [   398.89,    -6.56],
    [   397.93,    -7.51],
    [   397.57,    -8.83]])
    
csky = np.array([
    [    -7.60,    -2.99],
    [    -7.88,    -1.75],
    [    -8.76,    -0.86],
    [   -10.01,    -0.52],
    [   -11.25,    -0.85],
    [   -12.11,    -1.78],
    [   -12.45,    -3.01],
    [   -12.08,    -4.27],
    [   -11.17,    -5.19],
    [    -9.86,    -5.50],
    [    -8.62,    -5.18],
    [    -7.70,    -4.27],
    [    -7.37,    -3.03],
    [    -7.02,    -4.24],
    [    -6.04,    -5.13],
    [    -4.82,    -5.45],
    [    -3.59,    -5.11],
    [    -2.69,    -4.21],
    [    -2.33,    -2.98],
    [    -2.66,    -1.76],
    [    -3.54,    -0.88],
    [    -4.76,    -0.55],
    [    -5.94,    -0.89],
    [    -6.81,    -1.79],
    [    -7.12,    -3.01]])
    
esky = np.array([
    [  -396.28,    -2.93],
    [  -396.58,    -1.62],
    [  -397.55,    -0.67],
    [  -398.93,    -0.31],
    [  -400.30,    -0.66],
    [  -401.30,    -1.64],
    [  -401.66,    -2.96],
    [  -401.27,    -4.30],
    [  -400.22,    -5.28],
    [  -398.80,    -5.61],
    [  -397.39,    -5.28],
    [  -396.39,    -4.31],
    [  -396.01,    -2.99],
    [  -395.61,    -4.27],
    [  -394.56,    -5.21],
    [  -393.19,    -5.54],
    [  -391.82,    -5.19],
    [  -390.82,    -4.23],
    [  -390.44,    -2.93],
    [  -390.79,    -1.62],
    [  -391.77,    -0.69],
    [  -393.12,    -0.35],
    [  -394.43,    -0.70],
    [  -395.39,    -1.66],
    [  -395.73,    -2.97]])
    
ssky = np.array([
    [     3.90,  -396.51],
    [     3.58,  -395.14],
    [     2.64,  -394.13],
    [     1.31,  -393.78],
    [     0.00,  -394.13],
    [    -0.92,  -395.17],
    [    -1.27,  -396.53],
    [    -0.88,  -397.93],
    [     0.07,  -398.96],
    [     1.46,  -399.32],
    [     2.79,  -398.95],
    [     3.77,  -397.95],
    [     4.13,  -396.56],
    [     4.50,  -397.90],
    [     5.52,  -398.89],
    [     6.83,  -399.22],
    [     8.15,  -398.87],
    [     9.11,  -397.85],
    [     9.49,  -396.50],
    [     9.13,  -395.15],
    [     8.19,  -394.17],
    [     6.89,  -393.81],
    [     5.63,  -394.18],
    [     4.70,  -395.18],
    [     4.38,  -396.52]])

allsky = np.array((nsky, wsky, csky, esky, ssky))

#-----------file unpacker---------------------

def getFileArray(filename):
    nums2D = list()    
    with open(filename) as f:
        data = f.readlines()
    for row in data:
        numrow = list()
        words = row.split()
        for word in words:
            try:
                x = float(word)
            except:
                x = -0.0
            numrow.append(x)
        nums2D.append(numrow)
    # Make nums2D rectangular so asarray() returns an array not a list
    nrows = len(nums2D)
    ncols = 0
    for i in range(nrows):
        ncols = max(ncols, len(nums2D[i]))
    print('ncols = ' + str(ncols))
    for i in range(nrows):
        while len(nums2D[i]) < ncols:
            nums2D[i].append(-0.0)
    myArray = np.asarray(nums2D)
    return myArray
    
    
    
#=================PROGRAM STARTS HERE=================

# The 8 diagnostics are:  CW=01, CN=02, CE=03, CS=04, NW=21, SW=41, NE=23, SE=43 
# Nine columns:  ADC1, ADC2, skyU, skyV, Ngood, Xfpmm, Yfpmm, sigmaX, sigmaY
#                  0     1     2     3     4      5      6      7       8

labels = ['CW',  'CN',  'CE',  'CS',  'NW',  'SW',  'NE',  'SE']
pairs  = [[0,1], [0,2], [0,3], [0,4], [2,1], [4,1], [2,3], [4,3]]
colors = ['red', 'blue', 'green', 'black']
styles = ['-', '--']


print()

myArray = getFileArray(infile)
nrows = len(myArray)      # 125 rows
ncols = len(myArray[0])   # 9 parameters
print('nrows, ncols = ' + str(nrows) + ' ' + str(ncols))

fig = plt.figure(figsize=(7.,7.))  # inches
ax = fig.add_subplot(1,1,1)
plt.title('Monochromatic Shift of Central Field Star ')

center = 2 

hlist = list()
theoxlist = list()
theoylist = list()
theox0 = myArray[center, 5]
theoy0 = myArray[center, 6]

npts = 25

measxlist = list()
measylist = list()
measx0 = allsky[center, 0, 0]
measy0 = allsky[center, 0, 1]

for iadc in range(npts):
    hlist.append(iadc)
    theoxlist.append(myArray[iadc*5+center, 5] - theox0)
    theoylist.append(myArray[iadc*5+center, 6] - theoy0)
    measxlist.append(allsky[center, iadc, 0] - measx0)
    measylist.append(allsky[center, iadc, 1] - measy0)
    
plt.plot(hlist, theoxlist, 'b-')
plt.plot(hlist, theoylist, 'r-')
plt.errorbar(hlist, measxlist, yerr=0.2, fmt='o', color='b')
plt.errorbar(hlist, measylist, yerr=0.2, fmt='o', color='r') 
    
#---get the raw statistics--------

xsos = 0.
ysos = 0.
for i in range(npts):
    xsos += (theoxlist[i] - measxlist[i])**2
    ysos += (theoylist[i] - measylist[i])**2
xms = xsos/npts
yms = ysos/npts
xrms = np.sqrt(xms)
yrms = np.sqrt(yms)
print('RMS theory minus data in X and Y, mm = {:12.3f}{:12.3f}'.format(xrms,yrms))  # 0.271 and 0.039


trends = np.zeros((5,2,3))   # five cameras,  X and Y, adc=0,12,24
for icam in range(5):
    trends[icam,0,1] = allsky[icam, 12, 0] - allsky[icam, 0, 0]  # x,12
    trends[icam,0,2] = allsky[icam, 24, 0] - allsky[icam, 0, 0]  # x,24
    trends[icam,1,1] = allsky[icam, 12, 1] - allsky[icam, 0, 1]  # y,12
    trends[icam,1,2] = allsky[icam, 24, 1] - allsky[icam, 0, 1]  # y,24  

print('\nX trend                             ADC=0        12          24')
for icam in range(5):
    print(camlabels[icam] + '{:12.3f}{:12.3f}{:12.3f}'.format(trends[icam,0,0], trends[icam,0,1], trends[icam,0,2]))
    
print('\nY trend                             ADC=0        12          24')
for icam in range(5):
    print(camlabels[icam] + '{:12.3f}{:12.3f}{:12.3f}'.format(trends[icam,1,0], trends[icam,1,1], trends[icam,1,2]))    
        
""" Results
X trend                             ADC=0        12          24
CIN: 714791001685998720 16.79       0.000       0.230       0.490
CIW: 714955241235254656 16.94       0.000       0.250       0.510
CIC: 714073845226040576 14.83       0.000       0.230       0.480
CIE: 702345732650570624 14.25       0.000       0.270       0.550
CIS: 701824289260562304 15.90       0.000       0.230       0.480

Y trend                             ADC=0        12          24
CIN: 714791001685998720 16.79       0.000      -0.050      -0.030
CIW: 714955241235254656 16.94       0.000      -0.070      -0.040
CIC: 714073845226040576 14.83       0.000      -0.040      -0.020
CIE: 702345732650570624 14.25       0.000      -0.060      -0.040
CIS: 701824289260562304 15.90       0.000      -0.050      -0.010
"""

x0 = allsky[center, 0, 0] - measx0
x12 = allsky[center, 12, 0] - measx0
x24 = allsky[center, 24, 0] - measx0
print('x0, x12, x24 = {:12.3f}{:12.3f}{:12.3f}'.format(x0, x12, x24))  # 0.000, 0.0230, 0.480
# suggests subtracting 0.02mm per step

measzlist = list()
for iadc in range(npts):
    measzlist.append(measxlist[iadc] - 0.02*iadc)
    
plt.errorbar(hlist, measzlist, yerr=0.1, fmt='o', color='k')    

zsos = 0.
for i in range(len(theoxlist)):
    zsos += (theoxlist[i] - measzlist[i])**2
zms = zsos/npts
zrms = np.sqrt(zms)
print('RMS theory minus data in X, Y, Z, mm = {:12.3f}{:12.3f}{:12.3f}'.format(xrms,yrms,zrms))  # 0.271, 0.039, 0.042


ax.tick_params(direction='in')
ax.set_xlabel('ADC setting')
ax.set_ylabel('Centroid location shift, mm: X(blue), Y(red), X-0.02*iadc(black)')
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.95, infile, transform=ax.transAxes, fontsize=10)

ax.text(0.6, 0.10, 'X difference RMS=0.271mm', transform=ax.transAxes, fontsize=10, color='b')
ax.text(0.6, 0.06, 'Y difference RMS=0.039mm', transform=ax.transAxes, fontsize=10, color='r')
ax.text(0.6, 0.02, 'X de-trended RMS=0.042mm', transform=ax.transAxes, fontsize=10, color='k')



plt.savefig(plotfile) 
plt.show()           

