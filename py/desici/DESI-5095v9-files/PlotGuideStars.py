
# PlotGuideStars.py
# similar to GetGuideStars.py but automated selection of guide stars
# based on magnitude sorted list seen in ADC=0 frames. 
#
# M.Lampton July 2019
# 
#
# To match stars in a five-sensor-image-sequence starlists, I have 3D input:
# First dimension is ADC sequence number, 0 to 26.
# Second dimension is a varying length list of stars each with its sensor name
# Third dimension is the attribute number 0 to 4: x,y,mag,flag,gaia.
# I will want sorted magnitudes by sensor axis "mag" item. 
#
# Immediate task is that each ADC sequence number has a jumble of stars & sensors
# so these will have to be reorganized to make stars rankable by sensor.
#
# I will want to rearrange some axes:
#     Level 0: which sensor
#     Level 1: which ADC setting
#     Level 2: list of stars in image, each with attributes x,y,mag,flags,gaia
# so, FOR EACH SENSOR:
#    sort the stars in its ADC=0 frame
#    Starting with the brightest, see if it is present in all other ADC images
#       if so, keep it and list its 5 parms in each frame
#       if not, go back to the next brightest and search.
#
# Some potentially useful column identifiers:
#  XCENTROID
#  YCENTROID
#  DQ_FLAGS
#  RA
#  DEC
#  DET_ID   strings like 'ci-00004576o000010eCIN' CIN=North, CIW=west...
#  CAMERA   strings like 'CIN', etc
#  SOURCE_ID
#  PHOT_G_MEAN_MAG
#  PHOT_G_MEAN_FLUX_OVER_ERROR
#
# Result from GetGuideStars.py manual method:
#  CIN good choice is 714791001685998720   16.79mag
#  CIW good choice is 714955241235254656   16.94mag
#  CIC good choice is 714074669859822208   16.29mag
#  CIE good choice is 702439775254198272   14.59mag
#  CIS good choice is 701824289260562304   15.90mag

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt

PROGNAME = 'PlotGuideStars'

def printArray(name, array):
    print(name, end='')
    for item in array:
        print('{:12.2f}'.format(item), end='')
    print()   
    
def printStar(name, p):
    print(name + '{:9.2f} {:9.2f} {:9.2f} {:20d} {:4d}'.format(p[0],p[1],p[2],p[3],p[4]))
    

# Level zero: the sensors

cin = list()
ciw = list()
cic = list()
cie = list()
cis = list()
allcams = [cin,   ciw,   cic,   cie,   cis]
camIDs = ['CIN', 'CIW', 'CIC', 'CIE', 'CIS']
ncams = len(allcams)
styles = ['ro-', 'bo-', 'go-', 'mo-', 'ko-']
dots   = ['r.',  'b.',  'g.',  'm.',  'k.']
colors = ['r',   'b',   'g',   'm',   'k']

# Level 1: the tasks, one task per ADC setting

tasks =[[4576, 'ci-00004576_gaia-a.fits',   0,   0], 
        [4577, 'ci-00004577_gaia-a.fits',   0,   0],
        [4578, 'ci-00004578_gaia-a.fits',   0,   0],
        [4579, 'ci-00004579_gaia-a.fits',  29,   0],
        [4580, 'ci-00004580_gaia-a.fits',  60,   0],
        [4581, 'ci-00004581_gaia-a.fits',  89, 359],
        [4582, 'ci-00004582_gaia-a.fits', 119,   0], 
        [4583, 'ci-00004583_gaia-a.fits', 150,   0],
        [4584, 'ci-00004584_gaia-a.fits', 179,   0], 
        [4585, 'ci-00004585_gaia-a.fits', 210,   0],
        [4586, 'ci-00004586_gaia-a.fits', 240, 359],
        [4587, 'ci-00004587_gaia-a.fits', 269,   0],
        [4588, 'ci-00004588_gaia-a.fits', 299,   0],
        [4589, 'ci-00004589_gaia-a.fits', 329,   0],
        [4590, 'ci-00004590_gaia-a.fits', 359,   0],
        [4591, 'ci-00004591_gaia-a.fits',   0,  29],
        [4592, 'ci-00004592_gaia-a.fits',   0,  60],
        [4593, 'ci-00004593_gaia-a.fits',   0,  90],
        [4594, 'ci-00004594_gaia-a.fits', 359, 120], 
        [4595, 'ci-00004595_gaia-a.fits', 359, 150], 
        [4596, 'ci-00004596_gaia-a.fits', 359, 180],
        [4597, 'ci-00004597_gaia-a.fits', 359, 210],
        [4598, 'ci-00004598_gaia-a.fits', 359, 240],  # 250 was a typo
        [4599, 'ci-00004599_gaia-a.fits', 359, 270],
        [4600, 'ci-00004600_gaia-a.fits', 359, 300],
        [4601, 'ci-00004601_gaia-a.fits', 359, 330],
        [4602, 'ci-00004602_gaia-a.fits', 359,   0]]

ntasks = len(tasks)
print('\n\n\nNtasks = ' + str(ntasks))  # expect 27 tasks

for i in range(ntasks):  # cruise through tasks, creating unsorted "per sensor" lists
    task = tasks[i]
    expno = task[0]
    file = task[1]
    
    # set up the empty star lists for this task, for the five sensors
    cintask = list()
    cin.append(cintask)
    ciwtask = list()
    ciw.append(ciwtask)
    cictask = list()
    cic.append(cictask)
    cietask = list()
    cie.append(cietask)
    cistask = list()
    cis.append(cistask) 
    
    # extract the lists of all stars in this ADC setting, all five sensors
    tab = fits.getdata(file)
    x_ci = tab['XCENTROID']   
    y_ci = tab['YCENTROID'] 
    mag  = tab['PHOT_G_MEAN_MAG']
    gaia = tab['SOURCE_ID']  
    flag = tab['DQ_FLAGS'] 
    cams = tab['CAMERA']    
    nstars = len(x_ci)     # total stars per ADC setting
    # print('Task, Nstars = ' + str(i) + '  ' + str(nstars)) # all five imagers

    for istar in range(nstars):
        mylist = list()
        # for each star, gather its five parms
        mylist.append(x_ci[istar]) 
        mylist.append(y_ci[istar]) 
        mylist.append(mag[istar]) 
        mylist.append(gaia[istar])
        mylist.append(flag[istar]) 
        # now figure out where to append it
        cam = cams[istar]
        if cam == 'CIN': cintask.append(mylist)
        if cam == 'CIW': ciwtask.append(mylist)
        if cam == 'CIC': cictask.append(mylist)
        if cam == 'CIE': cietask.append(mylist)
        if cam == 'CIS': cistask.append(mylist)
        
        
""" Quality control: plenty stars/camera?  yes.         
for i in range(ntasks): 
    n = len(cin[i])
    w = len(ciw[i])
    c = len(cic[i])
    e = len(cie[i])
    s = len(cis[i])
    t = n+w+c+e+s 
    print('Task, Stars/cam[n,w,c,e,s], total = {:4d}  {:4d}{:4d}{:4d}{:4d}{:4d}  {:4d}'.format(i,n,w,c,e,s,t))

Task, Stars/cam[n,w,c,e,s], total =    0    12  13   9  13  14    61
Task, Stars/cam[n,w,c,e,s], total =    1    15   9  13  14  13    64
Task, Stars/cam[n,w,c,e,s], total =    2    14   8  10  14  12    58
Task, Stars/cam[n,w,c,e,s], total =    3    13   8  14  13  12    60
Task, Stars/cam[n,w,c,e,s], total =    4    12  12  12  12  13    61
Task, Stars/cam[n,w,c,e,s], total =    5    11  13  10  14  14    62
Task, Stars/cam[n,w,c,e,s], total =    6    12  12  12  13  14    63
Task, Stars/cam[n,w,c,e,s], total =    7    10  11  10  12  13    56
Task, Stars/cam[n,w,c,e,s], total =    8    15   9  13  15  14    66
Task, Stars/cam[n,w,c,e,s], total =    9    14   8  12  15  14    63
Task, Stars/cam[n,w,c,e,s], total =   10    13  10  12  16  13    64
Task, Stars/cam[n,w,c,e,s], total =   11    14  10  13  15  14    66
Task, Stars/cam[n,w,c,e,s], total =   12    13   9  13  14  12    61
Task, Stars/cam[n,w,c,e,s], total =   13    15   7  11  15  15    63
Task, Stars/cam[n,w,c,e,s], total =   14    16  10  14  12  13    65
Task, Stars/cam[n,w,c,e,s], total =   15    15   8  15  12  14    64
Task, Stars/cam[n,w,c,e,s], total =   16    14   9  18  12  13    66
Task, Stars/cam[n,w,c,e,s], total =   17    15   7  22  11  11    66
Task, Stars/cam[n,w,c,e,s], total =   18    13   9  23  10  14    69
Task, Stars/cam[n,w,c,e,s], total =   19    15  13  20  10  13    71
Task, Stars/cam[n,w,c,e,s], total =   20    17  14  19   9  13    72
Task, Stars/cam[n,w,c,e,s], total =   21    15  13  19   9  14    70
Task, Stars/cam[n,w,c,e,s], total =   22    14  14  20   9  13    70
Task, Stars/cam[n,w,c,e,s], total =   23    16  10  19  10  11    66
Task, Stars/cam[n,w,c,e,s], total =   24    14  10  15   9  12    60
Task, Stars/cam[n,w,c,e,s], total =   25    14  12  14   9  15    64
Task, Stars/cam[n,w,c,e,s], total =   26    15  10  15  12  13    65
"""



# Construct the allcams[] array, organized by camera number 0...4
# And for each camera, sort the magnitudes in its ADC=0
# And save the best star for each camera

bestID = list()

ADC = 0

def getMag(starparms):
    return starparms[2]

def getID(starparms):
    return starparms[3]

ncams = len(allcams)
for icam in range(5):
    myCamID = camIDs[icam]
    # print('\nNow starting camera = '+ myCamID)
    # print('len(allcams[icam][ADC=0]) = nrows =    '+str(len(allcams[icam][ADC])))
    # print('allcams[icam][ADC]: ')
    # for star in allcams[icam][ADC]:
    #     printStar('raw   ', star)
    
    sortedstars = sorted(allcams[icam][ADC], key=getMag)

    # now replace the raw star listing with the sorted star listing
    
    allcams[icam][ADC] = sortedstars
    # print('len(allcams[icam][ADC=0]) = nrows =    '+str(len(allcams[icam][ADC])))
    # print('allcams[icam][ADC]: ')
    # for star in allcams[icam][ADC]:
    #     printStar('sorted', star)
        
    # next, starting with the brightest, see if EVERY frame has this star...
    # print('...continuing with Camera = '+ myCamID)    
    trialID = None
    ntrials = len(sortedstars)
    for istar in range(ntrials):                # start with the ADC=0 list
        starparms = sortedstars[istar]          # shorthand 
        trialID = getID(starparms)              # search for this ID
        # print('trying  '+str(trialID))
        trialOK = False
        for adc in range(ntasks):               # examine each image
            trialOK = False                     # assume the worst each image       
            availables = allcams[icam][adc]
            for i in range(len(availables)):    # maybe this star?
                thisID = getID(availables[i])
                if trialID == thisID:
                    trialOK = True              # HA! found it
                    # print('matched '+str(thisID)+'  at ADC= '+str(adc))
                    break                       # proceed to the next ADC setting
            if not trialOK:                     # no stars matched this trialID
                # print('failed at ADC= '+str(adc))
                # print('so, trying next trialID star.')
                trialID = None                
                break                           # abandon adc loop: try another trialID
        if trialOK:                             # found our good star
            break     
    # print(myCamID +'........goodID '+ str(trialID))      
    bestID.append(trialID)
    
print('Best choices: ')
for icam in range(ncams):
    for starparms in allcams[icam][0]:
        starID = getID(starparms)
        if starID == bestID[icam]:
            printStar('cam='+str(icam), starparms)
            
""" Results:
Best choices: 
how many bestID findings =  5
cam=0  2033.37    886.34     16.79   714791001685998720    0
cam=1   705.50   1330.54     16.94   714955241235254656    0
cam=2  2744.73   1129.64     14.83   714073845226040576    0
cam=3  1459.55    599.77     14.25   702345732650570624    0
cam=4  1541.94   1160.01     15.90   701824289260562304    0
"""
camlabels = ['CIN: 714791001685998720 16.79',
             'CIW: 714955241235254656 16.94',
             'CIC: 714073845226040576 14.83', 
             'CIE: 702345732650570624 14.25',
             'CIS: 701824289260562304 15.90']
             
#------Next we create a list of the {x,y} centroids for each ADC frame, for each sensor


nxy = list()
wxy = list()
cxy = list()
exy = list()
sxy = list()
allxy = [nxy,   wxy,   cxy,   exy,   sxy]

for icam in range(ncams):
    for iadc in range(ntasks):
        nstars = len(allcams[icam][iadc])
        for istar in range(nstars):
            starparms = allcams[icam][iadc][istar]
            if getID(starparms) == bestID[icam]:
                xy = '{:12.2f},{:12.2f}'.format(starparms[0], starparms[1])
                allxy[icam].append(xy)  # strings
"""
print('\nAllXY :')
for k in range(len(allxy)):
    cam = allxy[k]
    print(camIDs[k])
    for j in range(len(cam)):
        xy = cam[j]
        print('    [{:3d}'.format(j) +', ' + cam[j] +'],')
"""
        
#-------Results, as numpy arrays, to allow slicing-------

CINXY=np.array([
    [  0,      2033.37,      886.34],
    [  1,      1713.58,     1125.66],
    [  2,      1659.75,     1120.09],
    [  3,      1693.44,      965.77],
    [  4,      1796.27,      855.44],
    [  5,      1944.99,      812.03],
    [  6,      2089.83,      853.44],
    [  7,      2196.57,      966.37],
    [  8,      2234.48,     1119.55],
    [  9,      2190.77,     1276.29],
    [ 10,      2083.11,     1388.78],
    [ 11,      1931.62,     1430.35],
    [ 12,      1781.96,     1392.10],
    [ 13,      1673.71,     1278.96],
    [ 14,      1634.28,     1126.62],
    [ 15,      1592.83,     1274.64],
    [ 16,      1481.37,     1385.08],
    [ 17,      1335.15,     1424.92],
    [ 18,      1187.55,     1385.65],
    [ 19,      1082.15,     1272.93],
    [ 20,      1039.13,     1122.53],
    [ 21,      1078.64,      969.95],
    [ 22,      1181.33,      861.70],
    [ 23,      1324.80,      820.99],
    [ 24,      1466.03,      863.27],
    [ 25,      1567.53,      972.29],
    [ 26,      1606.07,     1123.74]])
    
CIWXY=np.array([
    [  0,       705.50,     1330.54],
    [  1,       488.68,      974.27],
    [  2,       494.62,      917.05],
    [  3,       640.10,      951.93],
    [  4,       746.32,     1061.33],
    [  5,       787.73,     1215.09],
    [  6,       750.11,     1368.64],
    [  7,       640.63,     1480.70],
    [  8,       496.57,     1518.45],
    [  9,       347.23,     1475.65],
    [ 10,       237.77,     1363.25],
    [ 11,       197.83,     1200.80],
    [ 12,       235.96,     1045.15],
    [ 13,       341.68,      931.94],
    [ 14,       486.54,      888.57],
    [ 15,       339.19,      847.10],
    [ 16,       239.29,      727.57],
    [ 17,       202.25,      573.98],
    [ 18,       239.56,      421.84],
    [ 19,       347.23,      310.89],
    [ 20,       489.81,      265.83],
    [ 21,       636.26,      306.41],
    [ 22,       741.22,      414.09],
    [ 23,       780.12,      564.76],
    [ 24,       741.38,      710.52],
    [ 25,       635.41,      818.88],
    [ 26,       489.14,      859.44]])
    
CICXY=np.array([
    [  0,      2744.73,     1129.64],
    [  1,      2428.67,     1337.44],
    [  2,      2378.18,     1332.82],
    [  3,      2410.71,     1195.31],
    [  4,      2508.88,     1096.97],
    [  5,      2648.26,     1060.91],
    [  6,      2785.81,     1098.88],
    [  7,      2880.76,     1202.18],
    [  8,      2917.03,     1339.22],
    [  9,      2875.42,     1479.65],
    [ 10,      2772.90,     1580.47],
    [ 11,      2627.31,     1614.19],
    [ 12,      2490.03,     1577.38],
    [ 13,      2388.88,     1475.06],
    [ 14,      2352.50,     1337.64],
    [ 15,      2312.47,     1471.37],
    [ 16,      2203.12,     1569.43],
    [ 17,      2067.09,     1603.07],
    [ 18,      1930.60,     1564.92],
    [ 19,      1832.24,     1463.73],
    [ 20,      1792.98,     1326.69],
    [ 21,      1831.07,     1191.53],
    [ 22,      1928.90,     1094.29],
    [ 23,      2065.63,     1059.00],
    [ 24,      2196.20,     1097.50],
    [ 25,      2291.58,     1199.15],
    [ 26,      2324.84,     1335.00]])
    
CIEXY=np.array([
    [  0,      1459.55,      599.77],
    [  1,      1685.73,      954.73],
    [  2,      1680.98,     1010.51],
    [  3,      1534.37,      977.83],
    [  4,      1428.31,      870.21],
    [  5,      1388.51,      716.51],
    [  6,      1426.52,      563.84],
    [  7,      1535.28,      451.37],
    [  8,      1681.30,      410.96],
    [  9,      1830.67,      453.71],
    [ 10,      1939.82,      569.87],
    [ 11,      1977.68,      728.63],
    [ 12,      1940.61,      886.16],
    [ 13,      1833.76,      998.64],
    [ 14,      1686.86,     1040.69],
    [ 15,      1829.58,     1085.45],
    [ 16,      1934.37,     1202.03],
    [ 17,      1972.06,     1354.67],
    [ 18,      1933.36,     1507.96],
    [ 19,      1827.29,     1619.51],
    [ 20,      1682.52,     1662.29],
    [ 21,      1537.40,     1623.96],
    [ 22,      1434.02,     1514.81],
    [ 23,      1394.67,     1364.39],
    [ 24,      1433.86,     1218.03],
    [ 25,      1539.55,     1110.58],
    [ 26,      1685.34,     1072.15]])
    
CISXY=np.array([
    [  0,      1541.94,     1160.01],
    [  1,      1876.46,      921.51],
    [  2,      1930.03,      925.41],
    [  3,      1894.72,     1078.79],
    [  4,      1790.59,     1191.23],
    [  5,      1643.23,     1230.68],
    [  6,      1497.85,     1191.68],
    [  7,      1394.47,     1076.82],
    [  8,      1355.33,      924.76],
    [  9,      1398.39,      768.73],
    [ 10,      1504.36,      653.40],
    [ 11,      1659.03,      612.73],
    [ 12,      1806.93,      653.55],
    [ 13,      1915.60,      765.53],
    [ 14,      1955.88,      920.22],
    [ 15,      1996.41,      770.22],
    [ 16,      2110.18,      659.31],
    [ 17,      2255.60,      622.54],
    [ 18,      2401.89,      661.66],
    [ 19,      2509.33,      774.78],
    [ 20,      2550.93,      925.47],
    [ 21,      2511.27,     1076.05],
    [ 22,      2407.40,     1185.66],
    [ 23,      2262.63,     1226.59],
    [ 24,      2122.99,     1185.40],
    [ 25,      2019.26,     1073.96],
    [ 26,      1983.69,      924.03]])

ALLXY = [CINXY, CIWXY, CICXY, CIEXY, CISXY]


""" #---plot these five in a single diagram----

fig, ax = plt.subplots(figsize=(6,5))
for icam in range(ncams):                 # labels go into the upper left corner
    ax.text(0.03, 0.95-0.04*icam, camlabels[icam], transform=ax.transAxes, color=colors[icam], fontsize=8)    
for icam in range(ncams): 
    xvals = ALLXY[icam][:,1]              # all rows, second column
    yvals = ALLXY[icam][:,2]              # all rows, third column
    plt.plot(xvals, yvals, styles[icam])  # scatterplot with connected dots
ax.text(0.03, 0.02, PROGNAME+'.py', transform=ax.transAxes, color='black', fontsize=10)
ax.set_ylim(0, 2400)
ax.axis('equal')
ax.tick_params(direction='in')
plt.savefig('Single' + PROGNAME + '.png')
plt.show()
"""



""" #---plot these five in five diagrams------------

camboxes = [1, 3, 4, 5, 7]
fig, axarr = plt.subplots(nrows=3, ncols=3, figsize=(6,6))
fig.tight_layout()   # gives a looser layout: more room
axarr[0,0].axis('off')
icam = 0
for ibox in range(9):
    row = ibox//3
    col = ibox % 3
    if ibox in camboxes:
        xyarray = ALLXY[icam]
        xvals = xyarray[0:,1]
        yvals = xyarray[0:,2]
        axarr[row,col].plot(xvals, yvals, styles[icam])
        axarr[row,col].tick_params(direction='in', labelsize='6')
        axarr[row,col].axis('equal')
        axarr[row,col].set_title(camlabels[icam], fontsize=8)
        icam += 1
    else:   
        axarr[row,col].axis('off')
        
axarr[2,0].text(0.03, 0.02, PROGNAME+'.py', transform=ax.transAxes, color='black', fontsize=10)    
plt.savefig('Multi'+ PROGNAME + '.png')
plt.show()
"""


#-------Now I introduce the Coordinate Transform Pixels to FPxy---------
#
# Functions to go from CI pixel to focal plane x,y
# plug into desimodel.focalplane to get to sky position
# from Ashley ross via Aaron Meisner 15 July 2019
#
#  Coefficients for transform between pixel and CS5 x,y
#  X_0	Y_0	Z_0	a	b	c 	d	e	f
#  would have been better to do this as a dictionary
#  These are sky
Center = [ 14.053,  8.85,    -0.004,  -1.,    -0.0085,  0.0085,  -1.,   -0.0002, -0.0001]
South  = [-13.457, -404.841, -18.319,  1.,    -0.0018,  0.002,  0.9958, -0.0024,0.0914]
North  = [12.64,    404.366, -18.401, -1.,     0.005,  -0.005, -0.9956,  0.0001,0.0934]
East   = [-405.399, 12.16,   -18.455, 0.0038,  0.9959,	-1.,    0.0038,   0,     0.091]
West   = [405.298,  -13.192, -18.43, -0.0056,  -0.9957,  1.,	-0.0056, 0.0008, 0.0924]
pixsize = 0.009

def pixtoCS5(x,y,camstring):
	xp = 0
	yp = 0
	if camstring == 'CIC':
		cf = Center
		xp = 16
		yp = 3
	if camstring == 'CIW':
		cf = West
	if camstring == 'CIN':
		cf = North
	if camstring == 'CIE':
		cf = East
	if camstring == 'CIS':
		cf = South
	x += xp #corrects for overscan issue
	y += yp
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
	return cx, cy, cz
	
def CS5topix(x,y,camstring):
	if camstring == 'CIC':
		cf = Center
	if camstring == 'CIW':
		cf = West
	if camstring == 'CIN':
		cf = North
	if camstring == 'CIE':
		cf = East
	if camstring == 'CIS':
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


# Using ALLXY pixel table, build the arrays of focal CSS {x,y} coordinates
print('\nFocal Plane Coordinates...')

fpxy = np.empty((ncams,ntasks,2))

for icam in range(ncams):
    camID = camIDs[icam]
    for iadc in range(ntasks):
        xpix = ALLXY[icam][iadc][1]
        ypix = ALLXY[icam][iadc][2]
        xfp, yfp, zfp = pixtoCS5(xpix, ypix, camID)
        # print('    [{:9.2f},{:9.2f}],'.format(xfp,yfp))
        fpxy[icam][iadc][0] = xfp
        fpxy[icam][iadc][1] = yfp
        
for icam in range(ncams):
    print(camIDs[icam])
    for iadc in range(2,ntasks):
        print('    [{:9.2f},{:9.2f}],'.format(fpxy[icam][0], fpxy[icam][1]))
""" Results lookin good...
Focal Plane Coordinates...

"""

fig, ax = plt.subplots(figsize=(6,6))
ax.text(0.03, 0.95, camlabels[0], transform=ax.transAxes, color='red', fontsize=8)
ax.text(0.03, 0.91, camlabels[1], transform=ax.transAxes, color='blue', fontsize=8)
ax.text(0.03, 0.87, camlabels[2], transform=ax.transAxes, color='green', fontsize=8)
ax.text(0.03, 0.83, camlabels[3], transform=ax.transAxes, color='magenta', fontsize=8)
ax.text(0.03, 0.79, camlabels[4], transform=ax.transAxes, color='black', fontsize=8)
for icam in range(ncams):
    xvals = fpxy[icam,:,0]
    yvals = fpxy[icam,:,1]
    plt.plot(xvals, yvals, dots[icam])
ax.text(0.03, 0.02, PROGNAME+'.py', transform=ax.transAxes, color='black', fontsize=10)
ax.axis('equal')
ax.tick_params(direction='in')
plt.savefig('FocalPlane' + PROGNAME + '.png')
plt.show()