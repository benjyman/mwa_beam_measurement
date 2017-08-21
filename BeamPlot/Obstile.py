
##########################################################
################## MWA OBS AND PLOTTING ##################
##########################################################


import numpy as np, matplotlib.pyplot as plt
from glob                import glob
from subprocess          import Popen, PIPE
from BeamPlot.Settings     import obs_range
from BeamPlot.TimeMethods        import tm
from collections         import namedtuple


tilelist = (['0'+str(i)+'XX' for i in range(51,59)]+['rf0XX','rf1XX']+
            ['0'+str(i)+'YY' for i in range(51,59)]+['rf0YY','rf1YY'])

Tile = namedtuple('Tile', 'desig ID pol cls f_file t_file')


class Obstile(object):

    def __init__(self,tiles=tilelist, data_dir='./Obs/f-Rows',freqmin=136800000, freqmax=138200000):
        self.data_dir = data_dir
        self.tiles = {}
        self.data = {}
        self.freqmin = freqmin
        self.freqmax = freqmax
        self.freqspan = (freqmax-freqmin)/112
        
        for q in tiles:
            cls = 'mwa' if q[:2]!='rf' else 'ref'                 # type ref/mwa
            #frows = glob('./Obs/f-Rows/%s*'%q)[0]                 # file loc.
            frow_file_loc='%s/%s*'%(self.data_dir,q)
            print 'searching for f-row data in %s' % frow_file_loc
            frows = glob(frow_file_loc)[0]   
            #trows = glob('./Obs/t-Rows/%s*'%q)[0]                 # file loc.
            self.tiles[q] = Tile(q,q[:3],q[3:],cls,frows,'')   # tile dict.
    
    def printfreq(self):
        oc_chan0 = 137000000
        oc_span = 2500
        print '\n\n'
        fmin = self.freqmin
        for c in range(112):
            fmax = fmin+self.freqspan
            oc_chan = (fmin-oc_chan0)/oc_span
            print 'RFEx: %3d  f-min: %d  f-max: %d'%(c,fmin,fmax),
            if oc_chan>=0:
                print '  OC_chan: %d'%oc_chan
            else:
                print
            fmin = fmax
        return
    
    def getdata(self, tmin=None, tmax=None, c=range(1,113)):
        for tile in self.tiles:
            with open(self.tiles[tile].f_file) as dfile:
                d = dfile.readline().strip()
                stamps = [float(d[i:i+13]) for i in range(0,len(d),13)]
                ti, tf = self.restrict(stamps, tmin, tmax)
                chan = [[-1*ord(d)/2. for d in dfile.readline().strip()[ti:tf]]
                        if i in c else dfile.readline()[0:0] for i in range(1,113)]
            self.data[tile] = np.array([stamps[ti:tf]]+chan)
        return self

    def restrict(self, stamps, tmin=None, tmax=None):
        # Execute restriction before loading data set
        if tmax and tmax<stamps[0]:
            return 0,0
        if tmin and tmin>stamps[-1]:
            return -1,-1
        points = len(stamps)
        if tmin:
            for j in range(points):
                ti=j
                if stamps[j]>tmin:
                    break
        else:
            ti = 0
        if tmax:
            for j in range(ti,points):
                tf=j
                if stamps[j]>tmax:
                    break
        else:
            tf = -1
        return ti,tf
        
    def rdata(self, tile, tmin=None, tmax=None):
        # Execute restriction on full data set
        #print "rdata tmin tmax"
        #print tmin,tmax
        index = self.restrict(self.data[tile][0],tmin,tmax)
        return np.array([d[index[0]:index[1]] for d in self.data[tile]])

    def plot_tiles(tmin, tmax, tstep, chan, restrictions=None):

        if restrictions:
            e = Sateph(restrictions).ephemeris(tmin, tmax, tstep, quiet=False, restrict='full')
        else:
            e = Sateph().ephemeris(tmin, tmax, tstep, quiet=False, restrict='full')

        d = self.getdata(tmin.utc,tmax.utc,[chan])

        itp = {}
        for q in d.tiles:
            data = np.copy(d.rdata(q))
            if not data[0]:
                continue
            itp[q] = interp.interp1d(data[0],data[chan])

        for p in e:
            f, ax = plt.subplots(6, figsize=(16,10))
            plt.subplots_adjust(hspace=0.25, left=.06, right=.96)
            ax[0]=plt.subplot2grid((3,2), (0,0))
            ax[1]=plt.subplot2grid((3,2), (0,1), sharex=ax[0])
            ax[2]=plt.subplot2grid((3,2), (1,0), sharex=ax[0])
            ax[3]=plt.subplot2grid((3,2), (1,1), sharex=ax[0])
            ax[4]=plt.subplot2grid((3,2), (2,0), sharex=ax[0])
            ax[5]=plt.subplot2grid((3,2), (2,1), sharex=ax[0])

            pdata,eph = p[0],np.copy(p[1])
            sat, tiles = pdata[0], pdata[4]
            tmin, tmax = eph[0][0], eph[0][-1]
            eph[0]=eph[0]-tmin
            
            Satplot().plotpass([ax[0]],['alt'],sat,eph)
            ax[0].set_ylim([0,90])
            Satplot().plotpass([ax[1]],['alt'],sat,eph)
            ax[1].set_ylim([0,90])

            for q in tiles:
                anum = 2 if q['type']=='ref' else 4
                data = np.copy(d.rdata(q['desig'],tmin,tmax))
                data[0]=data[0]-tmin
                ax[anum].plot(data[0],data[chan])

                #change here to use sattimes
                xnew=np.arange(tmin,tmax,10)
                ynew=itp[q['desig']](xnew)-itp['rf0XX'](xnew)
                xnew=xnew-tmin
                ax[anum+1].plot(xnew,ynew)

            Satplot().setplot([ax[0]],['alt'], tm(0), tm(tmax-tmin))
            plt.show()


    def obs_exist(self, passlist, rmode='full'):
        """Create new passlist removing unobserved passes"""

        print 'Restricting passes to observation periods...\n'
        obs_passes = []
        for p in passlist:
            # Determine which tiles observed pass
            tlist = []
            for q in self.tiles:
                for period in obs_range[q]:
                    if p.rise.utc > period[0] and p.fall.utc < period[1]:
                        tlist.append(self.tiles[q])
            # If rmode is full, only include cases with both rec&mwa observations
            if rmode=='partial':
                if tlist:
                    p.tiles=' '.join([t.ID for t in tlist])
                    obs_passes.append(p)
                continue
            for pol in ['XX','YY']:
                if [q for q in tlist if q.cls=='ref' and q.pol==pol]:
                    if [q for q in tlist if q.cls=='mwa' and q.pol==pol]:
                        print p.tiles
                        print tlist
                        print ' '.join([t.ID for t in tlist])
                        p.tiles=' '.join([t.ID for t in tlist])
                        obs_passes.append(p)
                        break
        return obs_passes

            
###END###
