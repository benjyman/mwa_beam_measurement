##########################################################
#################### TIME METHODS ########################
##########################################################

import ephem, numpy
from datetime import datetime as dt, timedelta as tdelta
from BeamPlot.Settings import fault

class tm(object):
    
    """Converts times to Unix format, dform gives cal/TLE format, zone in hours
          *Accepted time formats:    -str        : (calender, TLE, JD/Unix)
                                     -int/float  : (JD/Unix)
                                     -class      : (self, datetime, ephem.Date)"""
    
    def __init__(self, d=None, dform='%Y/%m/%d %H:%M:%S.%f', zone=0):
        try:
            d = float(d)                              #detect JD/Unix/numpy strings
            if d<100000:
	        d = (d-25567.5)*86400                 #convert Julian dates to Unix
        except:
            if dform=='tle':                          #extract TLE time to datetime
                d = dt.strptime(d[18:23],'%y%j') + tdelta(float(d[23:32]))
            elif d and type(d)==str:                  #convert dform to datetime
                d = dt.strptime(d,dform)            
        if isinstance(d, dt):                         #datetime class
            d = (d-dt(1970,01,01)).total_seconds()
        elif isinstance(d, tm):                       #time of this class
            d = d.utc
        elif type(d)!=float:                          #eg. null strings
            raise AttributeError(fault[4])
        self.utc = d-zone*3600                        #UTC time in Unix format
        return

    #The following output a time of this class in the listed manner

    def dtime(self):                                  #datetime format
        return dt.utcfromtimestamp(self.utc)

    def eph(self):                                    #pyEphem date format
        return ephem.Date(self.utc/86400.+25567.5)    #dt(1970,01,01)-dt(1899,12,31,12)=25567.5

    def cal(self, dform='%Y/%m/%d %H:%M:%S.%f'):      #string format of type dform
        return self.dtime().strftime(dform)

    #The following describe operations between either class/time objects

    def __str__(self):                                #when using str() or %s on class object
        return self.cal('%Y-%m-%d %H:%M:%S')

    def __add__(self,t):                              #Addition
        if type(t) in (int,float):                                          
            return tm(self.utc+t)                         #outside class (+delta seconds)
        else:
            raise AttributeError(fault[5])                #b/w class/others (non-sensical)

    def __sub__(self,t):                              #Subtraction
        if type(t) in (int,float):                                          
            return tm(self.utc-t)                         #outside class (-delta seconds)
        elif type(t)==type(self):                        
            return self.utc-t.utc                         #b/w class (return delta seconds)
        else:
            raise AttributeError(fault[5])                #other objects
        
    def __cmp__(self, cls):                           #Comparisons
        if type(cls)==type(self):
            return (0 if self.utc==cls.utc else (1 if self.utc>cls.utc else -1))
        else:
            raise AttributeError(fault[5])

###END###
