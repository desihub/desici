
# EightDistancesPlot.py 
# showing how a field-center spot roams about the focal surface
#
# The 25 central spots are 0, 5, 10, 15, ...
# The five fields are:  C=0, W=1, N=2, E=3, S=4
# The 8 diagnostics are:  CW=01, CN=02, CE=03, CS=04, NW=21, SW=41, NE=23, SE=43 

import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'EightDistancesPlot.py'


infile = 'RT167_Nfields_5_Rmax_27.txt' 
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

# The 8 diagnostics are:  CW=01, CN=02, CE=03, CS=04, NW=21, SW=41, NE=23, SE=43 
# Nine columns:  ADC1, ADC2, skyU, skyV, Ngood, Xfpmm, Yfpmm, sigmaX, sigmaY
#                  0     1     2     3     4      5      6      7       8

labels = ['CW',  'CN',  'CE',  'CS',  'NW',  'SW',  'NE',  'SE']
pairs  = [[0,1], [0,2], [0,3], [0,4], [2,1], [4,1], [2,3], [4,3]]
colors = ['red', 'blue', 'green', 'black', 'cyan', 'magenta', 'gray', 'yellow']

print()

myArray = getFileArray(infile)
nrows = len(myArray)
ncols = len(myArray[0])
print('nrows, ncols = ' + str(nrows) + ' ' + str(ncols))

fig = plt.figure(figsize=(7.,7.))  # inches
ax = fig.add_subplot(1,1,1)
plt.title('Variation in Diagnostic Distances')

angles = range(25)
cw = list()
cn = list()
ce = list()
cs = list()
nw = list()
sw = list()
ne = list()
se = list()

for j in range(8):  # eight diagnostic distances to plot
    dist = list()
    f = pairs[j][0]  # one of the two fields
    g = pairs[j][1]  # other of the two fields
    dx0 = myArray[f, 5] - myArray[g, 5]
    dy0 = myArray[f, 6] - myArray[g, 6]    
    
    for i in angles:  # 25 angles of ADC rotation states
        dx = myArray[i*5+f, 5] - myArray[i*5+g, 5] - dx0
        dy = myArray[i*5+f, 6] - myArray[i*5+g, 6] - dy0
        dist.append(np.sqrt(dx**2 + dy**2))
    plt.plot(angles, dist, colors[j])
    ax.text(0.02, 0.90-0.04*j, labels[j], color=colors[j], transform=ax.transAxes, fontsize=12)

ax.set_ylim(-0.04, 0.65)
ax.tick_params(direction='in')
ax.set_xlabel('ADC setting')
ax.set_ylabel('change in image distance, mm')
ax.text(0.02, 0.02, PROGNAME, transform=ax.transAxes, fontsize=10)
ax.text(0.02, 0.95, infile, transform=ax.transAxes, fontsize=10)



plt.savefig(plotfile) 
plt.show()           

