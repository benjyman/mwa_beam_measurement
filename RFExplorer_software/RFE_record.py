#!/usr/bin/python

###### RFEXPLORER RECORDING SCRIPT ######

import serial, time, sys, os#, socket
import RFE_config_tiles as RFE
# Set maximum filesize for recording (10 hours)
max_file = 0

import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--time_obs', type=int, default=10,
                    help='Lendgth of observation in seconds')
parser.add_argument('--tile_index', type=int,
                    help='Index of tile, from 0 to 7 inclusive')
parser.add_argument('--date', type=str,
                    help='Date/time of the obs in the following format: hh:mm-dd-MM-YY e.g 16:15-11-Aug-2017')

args = parser.parse_args()

time_obs = args.time_obs
tile_index = args.tile_index
date = args.date

month_dict = {'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,
              'Sep':9,'Oct':10,'Nov':11,'Dec':12}

this_time = date.split('-')[0]
day = int(date.split('-')[1])
month = month_dict[date.split('-')[2]]
year = date.split('-')[3]

obs_date = "%s-%02d-%02d-%s" %(year,month,day,this_time)

class obs(object):
    #change argument to portlist eg. [(0,tile51),(1,tile52),........]
    #def __init__(self,obsname,portlist)
    def __init__(self,tile=None,obs_date=None):
        """Establish connection to RFE and hold all RFE output"""
        #self.obsname = obsname
        self.tile = tile
        #self.port = get_serial_port()
        self.port = int(tile)
        self.portname = '/dev/ttyUSB%d' %int(tile)
        self.ser = serial.Serial(self.portname, 500000, timeout=2)
        self.hold()
        
    def hold(self):
        """Hold RFE output and clear input (from RFE) buffer"""
        if self.ser.readline():
            self.ser.write(RFE.send["Hold"]) #result is RFExplorer will stop dumping data in input buffer
            #time.sleep(1)          #2->1
            self.ser.flushInput()
           # print "Mohit hold function is right"
    
    def configure(self):
        """Set RFE sweep configuration"""
        self.ser.write(RFE.send["Config"])  #writing config to RFExplorer via serial port
        if not self.ser.readline():
            error_msg = "%s: error> failed to write config to RFE"%host
            raise serial.SerialException(error_mesg)
        confirmation = self.ser.readline()
        #print confirmation
        self.hold()
        return confirmation
        
    def record(self,obslen=None):#,data='trialrecord.txt',log='triallog.txt'):
        
        # Observe for given seconds saving to filehandles
        if obslen:
            end = time.time()+obslen
        
        self.ser.write(RFE.send["Resume"])
        max_record = time.time() + obslen 
        
        #datafile = 'RFEdata-%s-%s-%s.txt'%(self.obsname,self.port,self.tile)
        #logfile = 'RFElog-%s-%s-%s.txt'%(self.obsname,self.port,self.tile)
        
        datafile = '%s%s_%s.txt'%(RFE.tiles[self.tile],RFE.pol,obs_date)
        logfile = 'RFElog_%s%s_%s.txt'%(RFE.tiles[self.tile],RFE.pol,obs_date)
        
        data=open(datafile,'w+')
        log=open(logfile,'w+') 
        
        header = 'tile%s-ttyUSB%s-%s-pol-%s'%(self.tile,self.tile,obs_date,RFE.pol)
        
        data.write(header+'\n')
        log.write(header+'\n')

        print 'COMMENCING RECORDING'
                 
        while(time.time()<max_record):
            
            if obslen and time.time()>end:
                break
            print 'sleeping',
            #time.sleep() 
            print 'wakeup'
            sweep = self.ser.readline() #your function is wrong! "readlines"
            print 'sweep-data'
            #print sweep
            sweep1 = sweep[-1] 
	    #print len(sweep)
 
            # Send initial pars and inconsistencies to log
            if len(sweep)==117:                #how that works?
                # Save valid data
                print 'valid data'
                data.write(('%s'%time.time())+sweep)
            else:
                print 'invalid data'
                log.write(('%s'%time.time())+sweep)
            # Check for early termination signal
            print 'looping'
            if os.path.exists("terminate"):
                os.remove("terminate")
                break
        self.hold()
        
        
if __name__=="__main__":
# 
    # run in terminal as python RFE_record.py Feb2017 0 test
    #if os.path.exists('locked-%s'%sys.argv[3]):
    #    print "%s:obs already in progress!"
    #    sys.exit()
    #open('locked-%s'%sys.argv[3],'a') 
    
    RFEnew=obs(tile=tile_index,obs_date=obs_date)
    RFEnew.configure()
    
    RFEnew.record(obslen=time_obs)
    sys.exit()

    
    ### OBSERVATION COMPLETE ###
    
