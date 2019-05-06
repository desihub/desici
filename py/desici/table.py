#!/usr/bin/python
"""
Provides a lookup table.
ROOT removed 20120223
ROOT for display only readded 20120409
Now uses matplotlib for delaunay interpolation
Provides nearest point (no interpolation) in case of failure
"""
#######################
import sys
import math
import copy
import random
import os
import time
import math
keeplist={}

try:
    import multiprocessing
    use_threading=False
except:
    import threading
    use_threading=True

ROOT=None
try:
    import ROOT
except:
    pass

triang=None
import pylab
import numpy
try:
    import matplotlib.tri as triang
    #import matplotlib.delaunay.interpolate as inter
    #import matplotlib.delaunay as triang
except:
    print("LUT_table has no matplotlib! Will use nearest points")
    triang=None

# KAR debugging
def fei(maxTBlevel=15):
    """ Prints formatted exception information """
    import traceback
    cla, exc, trbk = sys.exc_info()
    excName = cla.__name__
    try:
        excArgs = exc.__dict__["args"]
    except KeyError:
        excArgs = "<no args>"
    excTb = traceback.format_tb(trbk, maxTBlevel)
    return (excName, excArgs, excTb)

class LUT_table:
    # list of lists of len 8
    table=[]
    # Triangulation
    tri=None
    # Arrays of 6 parameters
    zval={}
    # X and Y
    az=None
    el=None
    # Interpolator class
    interp=None

    # Parameters and options
    par={}
    # Store the current table version number
    par['version']=0
    # When writing, allow the clobber of an existing file/version
    par['clobber']=False
    # Verbosity
    par['verbose']=-1
    # When loading, create wrap for (-360..0) and (360..720)
    par['pm360']=False
    # Do we need to recalculate the triangulating
    par['update']=True
    # A thread for displaying
    par['show']=None
    # A list of figures
    par['figs']=[]
    # Tells which method was used to interpolate
    par['method']=['',]*6
    # Information for the user about the last lookup
    par['info']=''
    # Default location
    # User must allow display (since it can "hang" on exit until windows close
    par['allow_display']=False
    

    def __init__(self,fname=None):
        """
        Initialize
        """
        self.minvalue=None
        self.maxvalue=None
        self.method=["None","None","None",
                     "None","None","None",]
        self.info="None"
        try:
            if not fname==None and os.path.exists(fname):
                self.read(fname)
                print("Loaded %s"%fname)
                return
        except:
            pass
        #print "Table not loaded (but initialized)"

    #######################################################

    def get(self,key):
        """
        Returns the value of a self.par parameter key
        Use 'get keys' to get a list of keys.
        """
        if not type(key) is str: key=repr(key)
        # Special case empty/keys
        if key in ['None', None, ''] or key.lower()=='keys':
            rval='get [keys, all, table, '
            for k in list(self.par.keys()):
                rval+=k+', '
            rval+=']'
        # Special case 'all'
        elif key.lower()=='all':
            return self.par
        # Special case table
        elif key.lower()=='table':
            return self.table
        # Key or fail
        else:
            try:
                rval=self.par[key]
            except:
                rval="FAILED: No parameter %s (try get keys)"%repr(key)

        return rval


    #######################################################

    def set(self,key,val):
        """
        Set the value for one of the self.par parameter keys
        entries. Will cast val to maintain types.
        """
        rval=""
        if not type(key) is str: key=repr(key)
        if not key in list(self.par.keys()):
            rval="FAILED: (set) You cannot set %s"%repr(key)
        else:
            if type(self.par[key]) is bool:
                try:
                    self.par[key]=bool(val)
                    rval="SUCCESS: %s now %s"%(key,repr(val))
                except:
                    rval="FAILED: Cannot set %s to %s"%(key,repr(val))
            if type(self.par[key]) is int:
                try:
                    self.par[key]=int(val)
                    rval="SUCCESS: %s now %s"%(key,repr(val))
                except:
                    rval="FAILED: Cannot set %s to %s"%(key,repr(val))
            if type(self.par[key]) is float:
                try:
                    self.par[key]=float(val)
                    rval="SUCCESS: %s now %s"%(key,repr(val))
                except:
                    rval="FAILED: Cannot set %s to %s"%(key,repr(val))
            if type(self.par[key]) is str:
                try:
                    self.par[key]=str(val)
                    rval="SUCCESS: %s now %s"%(key,repr(val))
                except:
                    rval="FAILED: Cannot set %s to %s"%(key,repr(val))
            if rval=="":
                self.par[key]=val
                rval="SUCCESS: set to type %s"%(repr(type(val)))

        if self.get('verbose')>=2:
            print("   "*2+"%s"%rval)
        return rval

    #######################################################
    def interpolate_delaunay(self,az,el,name="lut"):
        """
        Handles bad az,el,name before calling
        interpolate_delaunay_t(self,az,el,name)
        """
        try:
            az=float(az)
            el=float(el)
            name=str(name)
        except:
            self.set('method',["FAILED",]*6)
            rval=[None,]*6
            rval.append("FAILED")
            return rval
            
        return self.interpolate_delaunay_t(az,el,name)

    #######################################################

    def in_triangle(self,az,el,t1,t2,t3):
        """
        (deprecated)
        Determines if point (x,y) or here (az,el) lies
        within the triangle defined by points t1,t2,t3
        """
        p1=numpy.array((t1[0],t1[1]))
        p2=numpy.array((t2[0],t2[1]))
        p3=numpy.array((t3[0],t3[1]))

        p=numpy.array((az,el))
        
        v0=p2-p1
        v1=p3-p1
        v2=p-p1

        #print v0,v1,v2

        dot00=numpy.dot(v0,v0)
        dot01=numpy.dot(v0,v1)
        dot02=numpy.dot(v0,v2)
        dot11=numpy.dot(v1,v1)
        dot12=numpy.dot(v1,v2)

        try:
            invDenom=1.0/((dot00 * dot11) - (dot01 *dot01))
        except:
            return False
        u = ((dot11 * dot12) - (dot01 * dot02)) * invDenom
        v = ((dot00 * dot12) - (dot01 * dot02)) * invDenom
        
        return (u>=v) and (v>=0) and (u+v<1.0)

    #######################################################

    def interpolate_delaunay_t(self,az,el,name="lut"):
        """
        Use an delaunay triangulation from the matplotlib library
        Does the interpolation on self.table
        """

        self.set('method',['delaunay',]*6)
        self.set('info','')


        if self.tri==None or self.par['update']:
            self.update()

        try:
            tri=self.tri
        except:
            if self.get('verbose')>=0:
                print('   '*1+"Failed to find element %s"%repr(element))
            self.set('method',["FAILED",]*6)
            rval=[None,]*6
            rval.append("FAILED")
            return rval

        rval=[None,]*6
        rval.append("SUCCESS")
        for element in range(6):
            cint=self.interp[element]
            #y1=el-0.01
            #y2=el+0.01
            #x1=az-0.01
            #x2=az-0.02
            # Must be a 3x3 grid (and then pick center)
            #agrid=cint[y1:y2:1*3j,x1:x2:1*3j]
            #midline=grid[1]
            val=cint(az,el)
            #print "grid",len(grid),'x',len(grid[0])
            if repr(val)=='nan' or math.isnan(val):
                val=self.interpolate_nearest(az,el)[element]
                print("Replace nan with nearest",val, end=' ')
                print("for element",element)
            rval[element]=float(val)

        return rval

    #######################################################

    def interpolate_nearest(self,az,el,name="lut"):
        """
        Finds the point in the lookup table nearest (angular) to the requested
        loaction and returns that value
        """

        table=self.table

        # Find the table entry with the smallest distance
        neardist=None
        neari=0
        bestv=None
        for vi in range(len(table)):
            v=table[vi]
            dz=math.fabs(el-v[1])
            dist=self.get_dist_alt(az,el,v[0],v[1])
            if neardist==None or dist<neardist:
                neardist=dist
                bestv=v

        v=bestv
        self.set('method',['nearest',]*6)
        self.set('info','nearest point az=%f el=%f'%(v[0],v[1]))
        return v[2:]

    #######################################################

    def interpolate_inverse_d(self,az,el,
                              maxdist=30.0,maxdz=5,
                              power=3.0,depth=0,
                              usetable=None):
        """
        Uses an inverse distance weight to find best values for all
        six lookup table entries.

        finds values for azimuth (az) and elevation (el) (degrees)
        ignores (zero weight) for points at a distance greater than maxdist
        (degrees)

        Has a separate distance cut (maxdz) for zenith angle under the
        assumption that values will be more stable at constant zenith angle
        cuts off

        The inverse weight is set via power (w=1/d^power) where power should
        be greater than 2 for a two dimensional data set (or distant values
        will dominate do to increasing density)

        Finally, if there are insufficient (<4) points within the distance
        range cuts then the cuts will be loosened (20% per iteration) until 
        sufficient points are found (or max depth of 5 iterations)
        
        uses Shepard's method

        If no good interpolation is found within 5 iterations then
        interpolation_nearest(az,el) is called instead (simply returns the
        values from the nearest point in the lookup table.

        
        """
        #return self.interpolate_nearest(az,el,element,name)

        
        self.set('method',['inverse_d',]*6)
        self.set('info','')

        if usetable==None:
            table=self.table
        else:
            table=usetable


        
        weight=[]
        sum1=0.0
        nterms=0
        # Do sum over all points of distance from point vi to given az,el
        # Drop terms outside of max zenith and max distance
        # Count the number of "good" terms
        neardist=1E9
        neari=1
        for vi in range(len(table)):
            v=table[vi]

            dz=math.fabs(el-v[1])
            dist=self.get_dist_alt(az,el,v[0],v[1])
            if dist<neardist:
                neardist=dist
                neari=vi

            # If we hit a point closer than 0.01 degrees simply return
            # the 6-tuple of that point
            if dist<1E-2:
                if self.get('verbose')>=2:
                    print('   '*2*"Distance is %s"%repr(dist))
                return v[2:]
            w=math.pow(dist,-1.0*power)
            if dist>maxdist and maxdist>0: w=0.0
            if dz>maxdz and maxdz>0: w=0.0
            if w>0: nterms+=1
            weight.append(w)
            sum1+=w
        if neardist<0.01:
            v=table[neari]
            self.set('method',['nearest',]*6)
            return v[2:]

        # iterate if less than 4 points used in average
        if nterms<4:
            if depth<5:
                return self.interpolate_inverse_d(az,el,
                                                  maxdist=maxdist*1.2,
                                                  maxdz=maxdz*1.2,
                                                  power=power,
                                                  depth=depth+1)
            else:
                print("Interpolate_inverse_d failed.", end=' ')
                print("Returning interpolate_nearest() instead")
                return self.interpolate_nearest(az,el)

        # build sum over all
        sum2=[0.0,0.0,0.0,
              0.0,0.0,0.0,]
        for i in range(len(weight)):
            vv=self.table[i]
            w=weight[i]
            for j in range(6):
                sum2[j]+=vv[2+j]*w
        for j in range(6):
            sum2[j]/=sum1
            #print j, sum2[j],self.minvalue[j]
            if sum2[j]>self.maxvalue[j+2]:
                print(j,"change to max",self.maxvalue[j+2])
                sum2[j]=self.maxvalue[j+2]
            if sum2[j]<self.minvalue[j+2]:
                print(j,"change to min",self.minvalue[j+2])
                sum2[j]=self.minvalue[j+2]

        nterms=0
        for w in weight:
            if w>0: nterms+=1
        return sum2

    #######################################################

    def unshow(self):
        """
        (deprecated) Using ROOT to display table for now.
        """
        try:
            p=self.get('show')
            p.join()
            #p.stop()
        except:
            return "FAILED: could not unshow"

    #######################################################

    def show(self):
        """
        (deprecated) Using ROOT to display table for now.
        """
        if not self.get('allow_display'):
            rval="FAILED: (show) allow_display not set True"
            print(rval)
            return rval
        try:
            p=self.get('show')
        except:
            pass
        if not p==None:
            print("Cannot show a second time")
            return
        if not use_threading:
            p=multiprocessing.Process(target=self.detach_graph,args=[])
        else:
            p=threading.Thread(target=self.detach_graph,args=[])
        p.start()
        self.set('show',p)
        print("WARNING:  You must close figure window or thread will not exit!")
        return "SUCCESS: You must close figure window or thread will not exit!"

    #######################################################

    def detach_graph(self,*args):
        """
        (deprecated) Using ROOT to display table for now.
        """
        print("Closing the graphs will exit this thread")
        pylab.show()

    #######################################################

    def update(self):
        """
        Builds delaunay triangles and interpolators from table
        """
        # Grab the table
        table=self.table
        # Empty list of x=az, y=el and z=table element
        x=[]
        y=[]
        zz=[[],[],[],[],[],[]]
        for v in table:
            x.append(v[0])
            y.append(v[1])
            for i in range(6):
                zz[i].append(v[i+2])
        # Convert to numpy array
        x=numpy.array(x)
        y=numpy.array(y)
        # Store the z values for 6 elements
        for i in range(6):
            self.zval[i]=numpy.array(zz[i])
        # Store the x,y (or rather az,el)
        self.az=x
        self.el=y
        # Build a triangulation (needs to be converted)
        self.tri=triang.Triangulation(x,y)
        # Build an interpolator for each element
        self.interp=[]
        #import pdb
        #pdb.set_trace()
        for i in range(6):
            ni=triang.LinearTriInterpolator(self.tri, self.zval[i])  # needs to be converted
            self.interp.append(ni)

        # Store the triangles for drawing
        # For drawing (new matplotlib.tri module)
        #tri = triang(x, y).triangles
        #cens,edg,tri,neig = triang.delaunay(x,y)
        #self.tri2=tri
        # Turn off 'forced' update
        self.set('update',False)
        return "SUCCESS"

    #######################################################
    def display_root(self,element=0,name="nt"):
        """
        Uses ROOT library to display self.table
        LUT.png is created.
        """
        if ROOT==None:
            rval="FAILED: ROOT not imported"
            print(rval)
            return rval

        if 'f' in list(keeplist.keys()):
            try:
                keeplist['f'].Close()
                del keeplist['f']
            except:
                pass
        f=ROOT.TFile("LUT.root","RECREATE")
        keeplist['f']=f
        try:
            c1=ROOT.FindObject('c1')
        except:
            c1=ROOT.TCanvas('c1','LUT',1000,500)
        c1.SetFillColor(0)
        c1.Clear()
        c1.Divide(2,1,0.00,0.0)
        c1.Modified()
        c1.Update()
        keeplist['c1']=c1
        #print "Drawing element",repr(element)
        tg=ROOT.TGraph2D()
        tg.SetTitle("")
        keeplist['tg']=tg
        i=0
        for v in self.table:
            try:
                x=v[0]
                y=v[1]
                z=v[2+element]
                if element==2: 
                    z=v[2+element]-2500
                    print(x,y,z)
                if x<-215 or x>215: continue
                tg.SetPoint(i,x,y,z)
                i+=1
            except:
                print(fei())     
        #print "Loaded",tg.GetN(),"entries from table len",len(self.table)
        c1.cd(2)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.15)
        #tg.Draw("TRIW P0")
        tg.GetXaxis().SetTitle("az")
        tg.GetXaxis().SetTitleOffset(1.6)
        tg.GetYaxis().SetTitle("el")
        tg.GetYaxis().SetTitleOffset(1.6)
        if element==2:
            tg.GetZaxis().SetTitle("Hexapod %d [-2500]"%(element+1))
        else:
            tg.GetZaxis().SetTitle("Hexapod %d"%(element+1))
        tg.GetZaxis().SetTitleOffset(1.6)
        ROOT.gPad.SetPhi(90+60)
        ROOT.gPad.SetTheta(30)
        #c1.cd(2)
        tg.Draw("surf2z")
        tg.Draw("TRIW P0 same")
        ROOT.gPad.SetPhi(90+60)
        ROOT.gPad.SetTheta(30)
        #c1.cd(3) 
        #tg.Draw("TRIW P0")
        #ROOT.gPad.SetPhi(-30-180)
        #ROOT.gPad.SetTheta(90)
        c1.cd(1)
        ROOT.gPad.SetRightMargin(0.10)
        tg.Draw("surf2")
        tg.Draw("TRIW P0 same")
        ROOT.gPad.SetPhi(-30-180)
        ROOT.gPad.SetTheta(90)

        if 1==2:
            c1.Clear()
            c1.cd()
            tg.Draw("surf2z")
            tg.Draw("TRIW P0 same")
            ROOT.gPad.SetPhi(60+90)
            ROOT.gPad.SetTheta(30)
            
        
        c1.Modified()
        c1.Update()
        c1.Print("LUT%d.png"%(element+1))
        pass

    #######################################################

    def display(self,element=0,name="nt"):
        """
        (deprecated) Using ROOT to display table for now.
        Creates a picture of the llokup table for element
        call show() for it to appear on screen.
        also need to set('display_enabled',True)
        WARNING:
        python cannot exit while the figure is still displayed on the screen
        """
        if not triang:
            print("Cannot display without matplotlib")
            return

        if self.par['update']: self.update()
        if len(self.table)<1:
            return "FAILED: Table is not filled"

        if self.get('update'): self.update()
        
        figs=self.get('figs')
        n=len(figs)+1
        f=pylab.figure(n)
        figs.append(f)

        try:
            tri=self.tri2
            z=self.zval[element]
            x=self.az
            y=self.el
        except:
            print(fei())
            return "FAILED: Could not get triangulation"

        for t in tri:
            t_i=[t[0],t[1],t[2],t[0]]
            pylab.plot(x[t_i],y[t_i])

        pylab.plot(x,y,'o')
        pylab.draw()
        rval="SUCCESS"
        return rval

    #######################################################

    def get_dist(self,az1,el1,az2,el2):
        """
        Angular distance beween to points on sphere.
        Units for both input and output are degrees.
        """
        s1=math.sin(az1*math.pi/180.0)
        s2=math.sin(az2*math.pi/180.0)
        c1=math.cos(az1*math.pi/180.0)
        c2=math.cos(az2*math.pi/180.0)
        dele=el1-el2
        sdele=math.sin(dele*math.pi/180.0)
        cdele=math.cos(dele*math.pi/180.0)
        #Vincenty formula 
        t1=math.pow(c2*sdele,2.0)
        t2=math.pow(c1*s2-s1*c2*cdele,2.0)
        top=math.sqrt(t1+t2)
        bot=s1*s2+c1*c2*cdele
        dist=math.atan2(top,bot)*180.0/math.pi
        return dist

    def get_dist_alt(self,az1,el1,az2,el2):
        """
        Angular distance beween to points on sphere.
        Units for both input and output are degrees.
        AJR using formulae used for angular clustering
        """
        
        sa1=math.sin(az1*math.pi/180.0)
        sa2=math.sin(az2*math.pi/180.0)
        ca1=math.cos(az1*math.pi/180.0)
        ca2=math.cos(az2*math.pi/180.0)
        se1=math.sin(el1*math.pi/180.0)
        se2=math.sin(el2*math.pi/180.0)
        ce1=math.cos(el1*math.pi/180.0)
        ce2=math.cos(el2*math.pi/180.0)
        ac = ce1*ce2*(ca1*ca2 + sa1*sa2) + se1*se2
        dist=math.acos(ac)*180.0/math.pi
        return dist


    #######################################################

    def fake(self,az,el,term=0):
        """
        For testing only - used in create(writetype='fake')
        """
        rval=0.0
        falloff=math.cos(el*math.pi/180.0)
        rval+=0.05*math.cos(el*math.pi/180.0)*falloff
        rval-=0.05*float(term)*falloff
        rval-=0.003*math.cos(2*az*math.pi/180.0)*falloff
        return rval*100
        
    #######################################################

    def read(self,name=None):
        """
        Reads in a text file. Must be 8 columns.
        # Comments allowed
        # Version on first line
        """
        self.par['update']=True
        try:
            f=open(name)
            if self.get('verbose')>=1:
                print("   "*1+"Reading",name)
        except:
            print("FAILED: (read) Could not open %s"%repr(name))
            return -1

        self.table=[]
        self.minvalue=None
        self.maxvalue=None
        for line in f:
            line=line.strip()
            if line.startswith("# Version"):
                try:
                    self.set('version',line.split()[-1])
                except:
                    pass
            if line[0]=="#": continue
                
            l2=line.split()
            if not len(l2)==8:
                print("Bad line",line)
                continue
            
            l=[]
            ok=True
            for v in l2:
                try:
                    l.append(round(float(v),2))
                except:
                    ok=False
                    break
            if not ok: continue
            if self.minvalue==None or self.maxvalue==None:
                # Must copy!
                self.minvalue=copy.copy(l)
                self.maxvalue=copy.copy(l)
            for i in range(len(l)):
                val=l[i]
                cmin=self.minvalue[i]
                cmax=self.maxvalue[i]
                if val>cmax:
                    self.maxvalue[i]=val
                if val<cmin:
                    self.minvalue[i]=val
            self.add(copy.copy(l))
            if self.get('pm360'):
                l[0]-=360.
                if l[0]<0:
                    self.add(copy.copy(l))
                l[0]+=360.
                l[0]+=360.
                if l[0]>360:
                    self.add(copy.copy(l))
            
        if self.get('verbose')>=1:
            print("   "*1+"Loaded %d values from %s"%(len(self.table),repr(name)))
        return len(self.table)

    #######################################################

    def dump(self):
        """
        Dumps the table
        """
        for v in self.table:
            print(v)
        
    #######################################################

    def write(self,name,version=0):
        """
        Writes table to .txt file. name is the file stub
        version written to first line and appends file name
        IE /data/lut/LUT,1 becomes /data/lut/LUT.00000001.txt
        /data/lut/LUT,20120313 becomes /data/lut/LUT.20120313.txt
        """
        try:
            ver=int(version)
            name+=".%08d.txt"%(ver)
            if os.path.exists(name) and not self.par['clobber']:
                print("FAILED: (write) %s exists and I will not clobber."%name)
                return "FAILED"
            if self.get('verbose')>=1:
                print("   "*1+"Open",name)
            fout=open(name,"w")
        except:
            print("FAILED: Could not open %s for writing"%name)
            return "FAILED"
        
        if self.get('verbose')>=1:
            print("   "*1+"Writing %s "%(name))
        fout.write("# Version %08d\n"%ver)
                   
        vlist=['AZ','EL','H1','H2','H3','H4','H5','H6',]
        fout.write("# ")
        for v in vlist:
            fout.write("%6s "%v)
        fout.write("\n")

        for v in self.table:
            if not len(v)==8: continue
            for v2 in v:
                fout.write("%6.2f "%v2)
            fout.write("\n")
        fout.write("# End of file\n")
        return name
            
    def create(self,name,version=0,writetype="zero",number=100):
        """
        Creates a zeroed, faked or random table
        writetype=['zero', 'fake', 'rand', or 'fake2']
        number=(deprecated) now returns
        15 degrees steps in elevation
        30 degrees steps in azimuth
        """

        self.par['update']=True

        elrange=list(range(90,-1,-15))
        azrange=list(range(0,361,30))

        self.table=[]
        nwritten=0
        for el in elrange:
            if writetype=="fake2":break
            for az in azrange:
                nl=[]
                nwritten+=1
                nl.append(round(float(az),2))
                nl.append(round(float(el),2))
                for hexele in range(6):
                    if writetype=="zero":
                        nl.append(0.0)
                    if writetype=="rand":
                        nl.append(round(random.gauss(1.0,0.1),2))
                    if writetype=="fake" or writetype=="fake2":
                        nl.append(round(self.fake(az,el,hexele),2))
                        continue
                if len(nl)==8:
                    self.add(nl)
                    #print nwritten,nl[0],nl[1],nl[2]
        if self.get('verbose')>=0:
            print("   "*0+"Write (LUT_table)",name,version)
        nname=self.write(name,version)
        if self.get('verbose')>=1:
            print("   "*1+"Wrote",nname)
        if not nname=="FAILED" and os.path.exists(nname):
            return self.read(nname)
        else:
            return -1

    #######################################################

    def add(self,input):
        """ Add a new point into lookup table
        input must be convertable to list(input) of len 8
        add('[22,33,1,2,3,4,5,6]')
        add('(22,33,1,2,3,4,5,6)')
        add((22,33,1,2,3,4,5,6))
        add([22,33,1,2,3,4,5,6])
        should all work
        """
        # Handle string input
        if type(input) is str:
            try:
                input=eval(input)
            except:
                return "FAILED: Could not evaluate %s"%repr(input)
        # Can only add lists of len 8 (az,el,v1,v2,v3,v4,v5,v6)
        try:
            input=list(input)
            if not len(input)==8:
                return "FAILED: Need 8 terms az,el,v1,...v6"
        except:
            return "FAILED: Bad input."
        # Convert each element to float as needed
        for v in input:
            if not type(v) is float:
                if self.get('verbose')>=1:
                    print('   '*1+"Converting add element to float")
                try:
                    v=float(v)
                except:
                    return "FAILED: Bad element in input"
        # Tell self that an update to triangulation etc will be needed
        self.set('update',True)
        # See if this add is an overwrite
        caz=round(input[0],2)
        cel=round(input[1],2)
        for i in range(len(self.table)):
            az=round(self.table[i][0],2)
            el=round(self.table[i][1],2)
            if az==caz and el==cel:
                print("Replacing existing element at",az,el)
                self.table[i]=copy.copy(input)
                return "SUCCESS"
        # If not overwrite, append
        self.table.append(copy.copy(input))
        return "SUCCESS"

    #######################################################

    def example(self):
        """
        Simple example using 2 different 'fake' grids of differing density
        """
        table=self
        print("-"*80)
        print("Testing hexapod lookup table")
        table.set('clobber',True)
        table.set('verbose',-1)
        table.set('version',0)
        print('     %s'%repr(table.get('')))
        print('     %s'%repr(table.get('keys')))
        print('     %s'%repr(table.get('all')))
        
        
        for testi in range(2):
            print('     '+"-"*(80-5))
            n=[10,200]
            cver=n[testi]
            nn=n[testi]
            n=table.create("LUT",cver,"fake",nn)
            print("   %d entries loaded"%n)
            if triang:
                print("      Display via pylab")
                rval=table.display(name="nt%d"%(testi))
                print("      NOTE: Close the figure window(s) to exit",rval)

            val1=table.interpolate_nearest(22,22)
            val3=table.interpolate_delaunay(22,22,name="del%d"%testi)
            print("--------------------")
            print("- Comparing interpolation method results")
            print("--------------------")
            for i in range(len(val1)):
                tval=table.fake(22,22,term=i)
                print('  ', end=' ')
                print("%02d near=%8.3f delaunay=%8.3f   true=%8.3f"%(i,val1[i],val3[i],tval))

        if triang:
            table.show()
    def example2(self):
        """
        Print some example calls
        """
        print() 
        print("Interactive testing of table")
        print()
        print("Try:")
        print()
        print("# table.set('verbose',2)")
        print("table.set('clobber',True)")
        print("table.create('test',20120101,'fake',100)")
        print("# Increase 100 for more grig points")
        print("az=44.54")
        print("el=46.32")
        print("table.interpolate_delaunay(az,el)[0]")
        print("table.get('info')")
        print("table.interpolate_nearest(az,el)[0]")
        print("table.get('info')")
        print("# Compare to functional value")
        print("print table.fake(az,el,0)")
        print() 
        print("# table.display()")
        print("# table.show()")
    #######################################################

def main():
    print("use 'python -m doctest test.txt'")
    print(" or")
    print("t=table.LUT_table()")
    print("t.example2()")
    if 1==1:
        t=LUT_table()
        keeplist['tt']=t
        n=t.read("/usr/remote/user/sispi/decam/HexapodLUT/TelescopeLUT.txt")
        print("Delaunay 64.2,71.0",t.interpolate_delaunay(64.2,71.9)[0])
        t.display_root(5)
        t.display_root(4)
        t.display_root(3)
        t.display_root(2)
        t.display_root(1)
        t.display_root(0)
        try:
            os.remove("LUT.20120101.txt")
        except:
            pass
    return

if __name__ == "__main__":
    main()
