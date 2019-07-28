
# PlotCombined.py
# Reads the RT171a.py output file, 217 * 3 ADC combinations,
# and plots their differences as distortion grid maps.
#
# see also TangentPlane6.py
#
# M.Lampton July 2019
# 
#
#  RT171a.py output file RT171a_Nstars_217.txt  has ten columns...
#     ZA   ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms
#      0     1     2      3         4        5       6         7          8         9



import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'PlotCombined.py'

infile = 'RT171a_Nstars_217.txt'


def printArray(name, array):
    print(name, end='')
    for item in array:
        print('{:12.2f}'.format(item), end='')
    print()   
    
    
    
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
    #---Make nums2D rectangular---
    nrows = len(nums2D)
    ncols = 0
    for i in range(nrows):
        ncols = max(ncols, len(nums2D[i]))
    # print(' ncols = ' + str(ncols))
    for i in range(nrows):
        while len(nums2D[i]) < ncols:
            nums2D[i].append(-0.0)
    myArray = np.asarray(nums2D)
    return myArray
    
    
#----------MAIN PROGRAM STARTS HERE-----------------

raw2D = getFileArray(infile)
print('raw2D.shape = ' + str(raw2D.shape))  # (217+217+217) x 10

group0 = raw2D[:217,  :]  # ZA=0, ADC=0,0
group1 = raw2D[217:434,:] # ZA=60, ADC=0,0
group2 = raw2D[434:,  :]  # ZA=60, ADC=90,-90

print('group0.shape = ' + str(group0.shape))
print('group1.shape = ' + str(group1.shape))    
print('group2.shape = ' + str(group2.shape))    


    

#-----do the plot: case 1 vs case 0----------

fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(-400., +400.)
ax.set_ylim(-400., +400.)
plt.axis('equal')

sos = 0.0

nrays = 217
ArrowMag = 10.0 
for iray in range(nrays):
    x0 = group0[iray, 6]
    y0 = group0[iray, 7]
    x1 = group1[iray, 6]
    y1 = group1[iray, 7]
    dx = x1 - x0
    dy = y1 - y0
    
    sos += dx**2 + dy**2
    
    width = ArrowMag * dx
    height = ArrowMag * dy
    ax.arrow(x0, y0, width, height, head_width=5., head_length=12., fc='black')
    
ms = sos/(2*nrays)
rms = np.sqrt(ms)
rmsstr = 'RMSmm={:6.3f}'.format(rms)
ax.text(0.02, 0.90, rmsstr, transform=ax.transAxes, fontsize=10)

ax.tick_params(direction='in')
ax.set_xlabel('focal plane X mm')
ax.set_ylabel('focal plane Y mm')
arrowstr = 'ArrowMag=' + str(int(ArrowMag))
ax.text(0.02, 0.94, arrowstr, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=8)
plt.title('Difference Field: ZA=60 w/no ADC vs ZA=0')

noext = PROGNAME.rsplit('.',1)[0]
plt.savefig(noext + '_case_1_vs_0.png') 
plt.show()          

#-----do the plot: Case 2 vs case 0-----------

fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(-400., +400.)
ax.set_ylim(-400., +400.)
plt.axis('equal')

sos = 0.0

nrays = 217
ArrowMag = 100.0 
for iray in range(nrays):
    x0 = group0[iray, 6]
    y0 = group0[iray, 7]
    x1 = group2[iray, 6]
    y1 = group2[iray, 7]
    dx = x1 - x0
    dy = y1 - y0
    
    sos += dx**2 + dy**2
    
    width = ArrowMag * dx
    height = ArrowMag * dy
    ax.arrow(x0, y0, width, height, head_width=5., head_length=12., fc='black')
    
ms = sos/(2*nrays)
rms = np.sqrt(ms)
rmsstr = 'RMSmm={:6.3f}'.format(rms)
ax.text(0.02, 0.90, rmsstr, transform=ax.transAxes, fontsize=10)

ax.tick_params(direction='in')
ax.set_xlabel('focal plane X mm')
ax.set_ylabel('focal plane Y mm')
arrowstr = 'ArrowMag=' + str(int(ArrowMag))
ax.text(0.02, 0.94, arrowstr, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=8)
plt.title('Difference Field: ZA=60 with ADC vs ZA=0')

noext = PROGNAME.rsplit('.',1)[0]
plt.savefig(noext + '_case_2_vs_0.png') 
plt.show()          



#-----do the plot: Case 2 vs case 0 Nulled at Center-----------

fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(-400., +400.)
ax.set_ylim(-400., +400.)
plt.axis('equal')

sos = 0.

nrays = 217
ArrowMag = 100.0 
dxcenter = group2[0,6] - group0[0,6]
dycenter = group2[0,7] - group0[0,7] 
for iray in range(nrays):
    x0 = group0[iray, 6]
    y0 = group0[iray, 7]
    x1 = group2[iray, 6]
    y1 = group2[iray, 7]
    dx = x1 - x0 - dxcenter
    dy = y1 - y0 - dycenter
    
    sos += dx*dx + dy*dy
    
    width = ArrowMag * dx
    height = ArrowMag * dy
    ax.arrow(x0, y0, width, height, head_width=5., head_length=12., fc='black')

ms = sos/(2*nrays)
rms = np.sqrt(ms)
rmsstr = 'RMSmm={:6.3f}'.format(rms)
ax.text(0.02, 0.90, rmsstr, transform=ax.transAxes, fontsize=10)

ax.tick_params(direction='in')
ax.set_xlabel('focal plane X mm')
ax.set_ylabel('focal plane Y mm')
arrowstr = 'ArrowMag=' + str(int(ArrowMag))
ax.text(0.02, 0.94, arrowstr, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=8)
plt.title('Difference Field: ZA=60 with ADC vs ZA=0, Matched')

noext = PROGNAME.rsplit('.',1)[0]
plt.savefig(noext + '_case_2_vs_0_Matched.png') 
plt.show()          








 
