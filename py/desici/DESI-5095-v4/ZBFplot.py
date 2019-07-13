
# ZBFplot.py 
# showing how the eleven Zhao-Burge coefficients vary with ADC! and ADC2 roll angles


import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'ZBFplot.py'

infile = 'ZBF_Nrings_8_Rmax_27.txt'
plotfile = 'ZBFplot.png'



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

LUT = [3,4,6,7,8,9,10]  #ignore translations, magnification, and InOut
nplots = len(LUT)

names = ['0: S2 = Xtranslate', 
         '1: S3 = Ytranslate',
         '2: S4 = Magnify',
         '3: S7 = UpDownUp',
         '4: S8 = RightLeftRight',
         '5: S11 = InOut',
         '6: S22 = OutInOut',
         '7: T4 = Roll',
         '8: T7 = LeftRightCurl',
         '9: T8 = UpDownCurl',
         '10: T11 = RollAntiRoll']
         
colors = ['r', 'b', 'k', 'g', 'c', 'm', 'y', 'r', 'b', 'k', 'g']
         
myArray = getFileArray(infile)
nrows = len(myArray)
ncols = len(myArray[0])
print('nrows, ncols = ' + str(nrows) + ' ' + str(ncols))

fig = plt.figure(figsize=(7.,7.))  # inches
ax = fig.add_subplot(1,1,1)
plt.title('Monochromatic Zhao-Burge Coefficients vs ADC1, ADC2')

angles = range(24)
for i in range(nplots): 
# for icoef in icoefs:
    icoef = LUT[i]
    coefs = myArray[0:24, icoef]
    plt.plot(angles, coefs, color=colors[icoef])
    ax.text(0.02, 0.90-0.03*i, names[icoef], color=colors[icoef], \
        transform=ax.transAxes, fontsize=10)

ax.tick_params(direction='in')
ax.set_xlabel('ADC setting')
ax.set_ylabel('Z-B coefficient')
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.95, infile, transform=ax.transAxes, fontsize=10)

plt.savefig(plotfile) 
plt.show()           

