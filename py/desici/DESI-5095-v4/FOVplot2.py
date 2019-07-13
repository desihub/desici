
# FOVplot2.py 
# Shows field of view target locations or focal plane locations


import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'FOVplot2.py'

infile = 'RT166_Nrings_8_Rmax_27.txt'
plotfile = 'FOVplot2.png'



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

# ADC1  is column 0
# ADC2  is column 1
# Sky U is column 2
# Sky V is column 3
# Ngood is column 4
# FP  X is column 5
# FP  Y is column 6

# First 217 rows are ADC=0,0


         
myArray = getFileArray(infile)
nrows = len(myArray)
ncols = len(myArray[0])
print('nrows, ncols = ' + str(nrows) + ' ' + str(ncols))

fig = plt.figure(figsize=(7.,7.))  # inches
ax = fig.add_subplot(1,1,1)
plt.title('')

u = myArray[0:218, 2]
v = myArray[0:218, 3]
plt.scatter(u,v)

ax.tick_params(direction='in')
ax.set_xlabel('Sky coordinate U')
ax.set_ylabel('Sky coordinate V')
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.95, infile, transform=ax.transAxes, fontsize=10)

plt.savefig(plotfile) 
plt.show()           

