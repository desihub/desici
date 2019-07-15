
# PrismDevPlot.py 
# showing how a field-center spot roams about the focal surface
#
# The 24 central spots are 0, 217, 434, ....


import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'PrismDevPlot.py'


infile = 'RT170_Nrings_8_Rmax_28.txt' 
# Nine columns:  ADC1, ADC2, skyU, skyV, Ngood, Xfpmm, Yfpmm, sigmaX, sigmaY

plotfile = 'PrismDevPlot.png'



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
    print(' ncols = ' + str(ncols))
    for i in range(nrows):
        while len(nums2D[i]) < ncols:
            nums2D[i].append(-0.0)
    myArray = np.asarray(nums2D)
    return myArray
    
    
    
#=================PROGRAM STARTS HERE=================


print()


myArray = getFileArray(infile)
nrows = len(myArray)
ncols = len(myArray[0])
print('nrows, ncols = ' + str(nrows) + ' ' + str(ncols))

fig = plt.figure(figsize=(7.,7.))  # inches
ax = fig.add_subplot(1,1,1)
plt.title('Monochromatic ADC Prism Deviation')

angles = range(25)
xdevs = list()
ydevs = list()
for i in angles:
    j = i % 24
    xdevs.append(myArray[j*217, 5])
    ydevs.append(myArray[j*217, 6])


plt.plot(angles, xdevs, 'blue')
plt.plot(angles, ydevs, 'red')

ax.tick_params(direction='in')
ax.set_xlabel('ADC setting')
ax.set_ylabel('X (blue) or Y (red) deviation, mm')
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.95, infile, transform=ax.transAxes, fontsize=10)

plt.savefig(plotfile) 
plt.show()           

