
# PlotCombo.py
# Reads the RT171a.py output file, 217 * 3 ADC combinations,
# and plots their differences as distortion grid maps.
#
# derived from TangentPlane6.py
#
# M.Lampton 26 July 2019
# 
#
#  RT171a.py output file RT171a_Nstars_217.txt  has 868 records and ten columns...
#     ZA   ADC1  ADC2    RU        RV      Ngood    Xave      Yave       Xrms      Yrms
#      0     1     2      3         4        5       6         7          8         9
#
#-----------------------ZA----- ADC1------ADC2-----
# tasks = np.array([[   0.0,     0.0,    0.00],   records 0:217
#                   [   0.0,    90.0,  -90.00],   records 217:434
#                   [  60.0,     0.0,    0.00],   records 434:651
#                   [  60.0,    90.0,  -90.00]])  records 651:868
#
# Interesting plot 'Pure ADC' keeping ZA=0, subtracting 0000 from 9090 and zeroing center
# Interesting plot 'Pure ZA' keeping ADC=none, subtracting ZA=0 from ZA=60 and zeroing center
# Interesting plot "Combo" vary both, subtract 0,00 FROM 60,9090 and zeroing center. 

import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'PlotCombo.py'

infile = 'RT171a_Nstars_217.txt'


def printArray(name, array):
    print(name, end='')
    for item in array:
        print('{:12.3f}'.format(item), end='')
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
    
    
    
def plotArrows(headgroup, tailgroup, label, ArrowMag):
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_xlim(-400., +400.)
    ax.set_ylim(-400., +400.)
    plt.axis('equal')
    sos = 0.
    nrays = 217
    dxcenter = headgroup[0,6] - tailgroup[0,6]
    dycenter = headgroup[0,7] - tailgroup[0,7] 
    for iray in range(nrays):
        x0 = tailgroup[iray, 6]
        y0 = tailgroup[iray, 7]
        x1 = headgroup[iray, 6]
        y1 = headgroup[iray, 7]
        dx = x1 - x0 - dxcenter
        dy = y1 - y0 - dycenter
        sos += dx*dx + dy*dy
        width = ArrowMag * dx
        height = ArrowMag * dy
        ax.arrow(x0, y0, width, height, head_width=5., head_length=12., fc='black')
    ms = sos/nrays
    rms = np.sqrt(ms)
    rmsstr = 'RMSmm={:6.3f}'.format(rms)

    ax.tick_params(direction='in')
    ax.set_xlabel('focal plane X mm')
    ax.set_ylabel('focal plane Y mm')
    arrowstr = 'ArrowMag=' + str(int(ArrowMag)) 
    ax.text(0.02, 0.95, arrowstr, transform=ax.transAxes, fontsize=8)
    ax.text(0.02, 0.91, rmsstr, transform=ax.transAxes, fontsize=8)
    ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=8)
    plt.title(label)
    plt.savefig('PlotCombo_' + label + '.png') 
    plt.show()          

def plotSpots(group):
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_xlim(-400., +400.)
    ax.set_ylim(-400., +400.)
    plt.axis('equal')
    nrays = 217
    xpos = 0.
    ypos = 0.
    xneg = 0.
    yneg = 0.
    for iray in range(nrays):
        x = group[iray, 6]
        y = group[iray, 7]
        color = 'black'
        if iray==216: 
            color='red'     # +X
            xpos = x
        if iray==204: 
            color='green'   # -Y
            yneg = y
        if iray==192: 
            color='blue'    # -X
            xneg = x
        if iray==180: 
            color='yellow'  # +Y
            ypos = y
        ax.scatter(x, y, color=color)
    ax.tick_params(direction='in')
    ax.set_xlabel('focal plane X mm')
    ax.set_ylabel('focal plane Y mm')
    ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=8)
    # plt.savefig('PlotCombo_' + label + '.png') 
    plt.show()          
    dx = xpos - xneg
    dy = ypos - yneg
    return dx, dy
     


#----------MAIN PROGRAM STARTS HERE-----------------

raw2D = getFileArray(infile)
print('raw2D.shape = ' + str(raw2D.shape))  # (217+217+217+217) x 10

group0 = raw2D[0:217,  :]  # ZA=0, ADC=0,0
group1 = raw2D[217:434,:]  # ZA=0, ADC=+90-90
group2 = raw2D[434:651,:]  # ZA=60, ADC=0,0
group3 = raw2D[651:868,:]  # ZA=60, ADC=+90-90

groups = np.array([group0, group1, group2, group3])

"""
print('group0.shape = ' + str(group0.shape))
print('group1.shape = ' + str(group1.shape))    
print('group2.shape = ' + str(group2.shape))    
print('group3.shape = ' + str(group3.shape))    
"""

"""#-----Case 1 minus Case 0:  Pure ADC, keeping ZA=0, nulled at center-----------
headgroup = group1
tailgroup = group0
label = 'ZA-zero-ADC-Motion'
ArrowMag = 100.0
plotArrows(headgroup, tailgroup, label, ArrowMag)

#-----Case 2 minus case 0: pure ZA, keeping ADC=0, nulled at center----
headgroup = group2
tailgroup = group0
label = 'ZA-Motion-ADC-zero'
ArrowMag = 100.0
plotArrows(headgroup, tailgroup, label, ArrowMag)

#----Case 3 minus case 0:  Full ZA and full ADC, nulled at center----
headgroup = group3
tailgroup = group0
label = 'ZA-Motion-ADC-compensating'
ArrowMag = 100.0
plotArrows(headgroup, tailgroup, label, ArrowMag)
"""

for group in groups:
    diamx, diamy = plotSpots(group)
    both = np.array([diamx, diamy])
    printArray('DiamX, Diamy = ', both)

""" Result:
DiamX, Diamy =      814.489     814.490
DiamX, Diamy =      814.492     814.504
DiamX, Diamy =      814.496     813.924
DiamX, Diamy =      814.492     813.905

Are these consistent with the claims of my DESI-4957v11 chart 5:
Horizontal compression = 220 ppm
Vertical compression = 220ppm * sec^2(ZA)
"""
 
