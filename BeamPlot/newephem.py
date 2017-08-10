
class ephemeris(object):

    

    def ephemeris(self, sat, timestamps, quiet=False, showplot=True, restrict=None):
        """Calculates ephemeris for passes within trange using step size tstep (excluding
        right bound if tmax-tmin%tstep!=0) to produce a list of numpy arrays in t,alt,azi"""
        
        def compute(u, seek=False):
            # Keep constant tle for ephemeris calculation over single pass
            self.obs.date = ephem.Date(self.obs.date+u)
            self.target.compute(self.obs)
            return [tm(self.obs.date).utc, self.target.alt, self.target.az]
        
        def ephemeris(tstep, satpass):
            # Compute eph for single pass- operates in the form: target.compute(obs)
            self.target = ephem.readtle('TLE0',*satpass.tle.split(','))
            tstep = tstep*ephem.second
            rise = tmin.eph() + tstep*((satpass['rise'].eph()-tmin.eph())//tstep)
            fall = tmin.eph() + tstep*((satpass['fall'].eph()-tmin.eph())//tstep)
            steps = int((fall-rise)//tstep)
            self.obs.date = rise
            return np.swapaxes([compute(tstep) for i in range(steps+1)],0,1)

        # Compute ephemeris for passes over requested range
        print '\n\nComputing ephemeris...',
        eph = [
        #eph = [[p, ephemeris(tstep, p)] for p in plist]
        if not quiet:
            # Show plot of computed passes unless showplot=False
            ref = Satplot().ephemref(tmin, tmax, eph, self.sats, show=showplot)
            print "%d passes calculated, reference plot:\n>'%s'\n"%(len(eph),ref)+lb
        return eph

