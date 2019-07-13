
# CSVtoB4.py
#
#
# Translates Excel CSV file format to B4 file format.
# For Sellmeier interpolation, assumes that most Schott Product Data columns have been deleted,
#    leaving only Glass, Nd, Ne, perhaps SellSix, wavelengths with refraction.   
# Output files will be name same as input file but with internal extension appended.
# Adaptive field widths and custom ruler fabrication.
# Accepts messy native .xlsx files with trailing empty fields, dangling rows; Eliminates that mess.
#
#  To get 3-digit wavelengths, explicitly set Excel to 3-digits.
#  Otherwise it will truncate enough to eliminate trailing zeros. 
#  Or: use ray table wavelengths having no trailing zeros.
#  
#  USAGE AT COMMAND LINE PROMPT:
#    python CSVtoB4.py  myfile1  myfile2  myfile3 ...
#    
#  FILE INPUT FORMAT:
#     title row at top
#     any number of blank rows such as those found in Schott glass spreadsheets
#     header row next
#     any number of blank rows
#     Finally the data records.
#     any number of blank rows
#
#  FILE OUTPUT FORMAT:
#     guideNumber number and title at top
#     header row next, spaced by colons
#     ruler row, colons, justified for max width each field
#     finally the data, justified to appropriate width each field.
#
#
# M.Lampton UCB SSL 2017, 2018; rewritten for Py2 + Py3 compatibility

from __future__ import print_function   # allows Python2 to use print(end=' ') from Python3
import sys
import csv
import warnings                   # suppress annoying FutureWarning
warnings.filterwarnings("ignore")

#--------conversion routine starts here--------------------

def isEmpty(anyList):
    if len(anyList) < 1:
        return True
    width = 0
    for i in range(0, len(anyList)):
        width = max(width, len(anyList[i]))
    return True if width==0 else False
    
def suckInt(s):
    # extracts a positive integer from a messy string
    slen = len(s)
    if slen<1:
        return 0
    numstring = ""
    i = 0
    while i<slen and not s[i].isdigit():
        i += 1
    while i<slen and s[i].isdigit():
        numstring += s[i]
        i += 1
    if len(numstring) ==0:
        return 0
    try:
        result = int(numstring)
        return result
    except ValueError:
        return 0
        
        
def getGuideNumber(s):
    # extracts a leading number from string "s"
    slen = len(s)
    if slen<1:
        return 0
    i = 0
    numstr = ""
    while i<slen and s[i]==' ': 
        i += 1
    while i<slen and s[i].isdigit():
        numstr += s[i]
        i += 1
    if len(numstr) < 1:
        return 0
    try:
        result = int(numstr)
        return result
    except ValueError:
        return 0
    
    
def receiveCSV(fname):
    data = list()
    print('\nTrying: ', fname)
    try:
        data = list(csv.reader(open(fname, 'rU')))  # 2D list of snippets
    except IOError:
        print("Could not open that file. Quitting this file.")
        return data
    if len(data) < 3:
        print("Fewer than three CSV records are found. Quitting this file.")
        return data
    for irow in range(0, len(data)):    
        for jcol in range(0, len(data[irow])):
            data[irow][jcol] = data[irow][jcol].strip()  # unnecessary from Excel
    print("Initial nrows = ", len(data))

    # delete all empty rows
    initialLength = len(data)
    for irow in range(initialLength-1, 0, -1):  # don't eliminate title row even if empty
        if isEmpty(data[irow]):
            del data[irow]
    print("After removing empties, nrows = ", len(data))

    # get number of fields from longest row that holds any actual data
    nfields = 0
    for irow in range(0, len(data)):
       for jcol in range(0, len(data[irow])):
           if len(data[irow][jcol]) > 0:
               nfields = max(nfields, jcol+1)
               
    print("Nfields = ", nfields)
    if nfields < 2:
        print("Nfields is < 2.  Quitting this file.")
        del data[:]  # empty it of everything
        return data

    # make all rows have nfields by appending empty fields where needed.
    for irow in range(len(data)):
        data[irow] = data[irow][:nfields]  # truncate beyond nfields
        while len(data[irow]) < nfields:   # append empty fields
            data[irow].append("")

    print("Retaining the original title, nrows = ", len(data))
       
    # Cleanup is now complete.
    # data[0] is the original title
    # data[1] is the header row
    # data[2].... are the actual data
    
    return data
    
    
    
    
def createB4(data, fname):
    linelist = list()
    if len(data) < 1:
        print("Got zilch from receiveCSV().  Quitting this file.")
        return linelist   # empty
        
    # For CSVtoB4 reformatting, need uniform widths, colons, ruler, title. 
    # Get and adjust lengths of snippets in each column while building ruler.
    
    guide = suckInt(data[0][0])
    del data[0]  # abandon old title row
    print("After removing initial title, nrows = ", len(data))    
    
    nfields = len(data[0])
    ruler = list()
    for jcol in range(0, min(nfields, len(data[0]))):    # header row
        width = 0
        for irow in range(0, len(data)):
            width = max(width, len(data[irow][jcol]))
        width += 1
        for irow in range(0, len(data)):
            while len(data[irow][jcol]) < width:
               data[irow][jcol] = data[irow][jcol] + ' '  # inflate to equal width
        for iseg in range(0, width):
            ruler.append('-')                             
        ruler.append(':')
    rulerstring = "".join(ruler)           
    
    # build the output linelist[] using colon separators
    for row in range(0, len(data)):
        mychars = list()
        for col in range(0, len(data[row])):
            snippet = data[row][col]        # extract each snippet
            mychars.append(snippet)         # append it to mychars
            mychars.append(':')
        linelist.append("".join(mychars))   # join them into a single string.
    headerstring = linelist[0]
    
    # Next: build a new title line.
    # What extension is appropriate?
    upfname = fname
    upfname.upper()
    # print("now testing upfname = ", upfname)
    ext = ".B4"
    if ".OPT" in upfname:
        ext = ".OPT"
    if ".RAY" in upfname:
        ext = ".RAY"
    if ".MED" in upfname:
        ext = ".MED"
    outname = fname + ext
    nrecords = len(linelist) -1
    print("nrecords = ", nrecords)

    print("Found guideNumber = ", guide)
    if guide>0:
        nrecords = min(guide, nrecords)
    newtitle = "   " +str(nrecords) + '  glasses   ' + outname
    linelist.insert(0, newtitle)    # headerstring becomes linelist[1]
    linelist.insert(2, rulerstring) # rulerstring becomes linelist[2]
    return linelist, outname


def convertCSVtoB4(fname):
    myData = receiveCSV(fname);
    linelist, outname = createB4(myData, fname)
    if len(linelist) < 2:
        print("No linelist was generated. Quitting this file.")
        return

    # show the results.
    print('\n')
    for row in range(0, len(linelist)):
        print(linelist[row])
        
    # save the file.
    outfile = open(outname, 'w')
    for row in range(0, len(linelist)):
        outfile.write(linelist[row])
        outfile.write('\n')
    outfile.close()

#---------Main program starts here-------------

nfiles = len(sys.argv) - 1
print("Number of files specified : ", nfiles)
if nfiles < 1:
   print("Please specify one or more filenames on command line.  Quitting.")
   quit()

del sys.argv[0]    # delete self name
for name in sys.argv:
    convertCSVtoB4(name)
    
print("All done.")













