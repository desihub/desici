

# newB4toCSV.py
# Translates B4 file format to CSV file format.
# Input files may have their usual B4 filenames.
# Output files have same names but with ".CSV" appended. 
# Each record gets commas at the ruler-colon positions.
# Two practical issues: in B4, some colons may be missing;
# Also some commas may be present within headers or data fields.
# This revision fixes both. 
# May 2018: revised to be both Py2 and Py3 compatible.
#
#  USAGE AT COMMAND LINE PROMPT:
#    python B4toCSV.py  myfile1  myfile2  myfile3 ...
#    
#  FILE INPUT FORMAT:
#     guide number and title row at top
#     header row next
#     ruler row
#     Finally the data records.
#
#  FILE OUTPUT FORMAT:
#     guideNumber number and title at top
#     header row next, spaced by commas
#     finally the data, comma separated.
#
#
# M.Lampton UCB SSL 2017

from __future__ import print_function   # allows Python2 to use print(end='') from Python3
import sys
import warnings                   # suppress annoying FutureWarning
warnings.filterwarnings("ignore")
             
#---helper function---------------

def getGuideEnd(string):
    # finds location of the end of the guide number
    slen = len(string)
    if slen<1:
        return -0
    i = 0
    while (i<slen) and not string[i].isdigit():
        i += 1
    while (i<slen) and string[i].isdigit():
        i += 1
    return i


#--------conversion routine starts here--------------------

def convertB4toCSV(fname):
    print('trying ', fname)
    try:
        infile = open(fname)
    except IOError:
        print("Could not open that file. Quitting.")
        return
    
    linelist = infile.readlines()  # these retain their EOLs
    for i in range(0, len(linelist)):
        print(linelist[i], end=" ")         # space not newline
    
    # package the guide number, if any...
    end = getGuideEnd(linelist[0])
    titlechars = list(linelist[0])
    titlechars.insert(end, ' ')
    titlechars.insert(end+1, ' ')
    titlechars.insert(end+2, ',')
    linelist[0] = "".join(titlechars)

    # Takes a bit of explicit work to get a 2D "list" ....
    allchars = [] 
    for row in range(0, len(linelist)):
        thesechars = list(linelist[row])  # retains EOLs
        allchars.append(thesechars)
        
    rulerstr = "".join(allchars[2])       # sequester the ruler. 
    
    for row in range(1, len(allchars)):     # clean up allchars
        for i in range(0, len(allchars[row])):  
            thischar = allchars[row][i]
            if thischar == ',' or thischar == ':':
                allchars[row][i] = '-'

    # Lengthen the header if the ruler needs the spaces...
    while (len(allchars[1]) < len(rulerstr)):
        allchars[1].remove('\n')
        allchars[1].append(' ')
        allchars[1].append('\n')

    # Put commas at all ruler colon locations except the title....
    for row in range(1, len(allchars)):
        lesser = min(len(rulerstr), len(allchars[row]))
        for i in range(0, lesser):
            if rulerstr[i]==':':                   # search for colon in the ruler
                allchars[row][i] = ','               # insert a colon into the headers
        linelist[row] = "".join(allchars[row])      # substitute the new line.
        
    # Discard the ruler.
    linelist.remove(linelist[2])  
    
    # Show the revised CSV file....
    print(' ')
    for i in range(0, len(linelist)):
        print(linelist[i], end=" ")  # space but no newline

    # Save the CSV file...
    outfile = open(fname + ".CSV", "w")
    outfile.writelines(linelist)
    outfile.close()
    print(' ')

#---------Main program starts here-------------


nfiles = len(sys.argv) - 1
print("Number of files specified : ", nfiles)
if nfiles < 1:
   print("Please specify one or more filenames on command line. Wildcards ok.  Quitting.")
   quit()

del sys.argv[0]    # delete self name
for name in sys.argv:
    convertB4toCSV(name)
    
print("All done.")













