#!/usr/bin/python

###### RFEXPLORER RECORDING SCRIPT ######

import serial, time, sys, os#, socket
import RFE_config as RFE
time_obs = 60.0*25
# Set maximum filesize for recording (10 hours)
max_file = 0


def get_serial_port():
    """Identify serial port by locating port to which cp210x is attached"""   
#the use of docstrings (triple inverted commas) is to provide documentation when viewed using command help or man
    output = os.popen("dmesg | grep 'cp210x converter now attached'").read() # Open a pipe to/from a command returning a file object.
    try:
        port = output.split("ttyUSB")[1][0] #here split first split the output in two  ## I removed '/dev/ttyUSB%s'%
    except IndexError:
        error_msg = "%s: error> cannot find cp210x driver attached"%obsname
        raise serial.SerialException(error_msg)
    return port


class obs(object):
    #change argument to portlist eg. [(0,tile51),(1,tile52),........]
    #def __init__(self,obsname,portlist)
    def __init__(self,obsname,tile,port):
        """Establish connection to RFE and hold all RFE output"""
        self.obsname = obsname
        self.tile = tile
        self.port = get_serial_port() 
        self.portname = '/dev/ttyUSB%s'%self.port
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
        
    def record(self,obslen):#,data='trialrecord.txt',log='triallog.txt'):
        
        # Observe for given seconds saving to filehandles
        if obslen:
            end = time.time()+obslen
        
        self.ser.write(RFE.send["Resume"])
        max_record = time.time() + obslen 
        
        datafile = 'RFEdata-%s-%s-%s.txt'%(self.obsname,self.port,self.tile)
        logfile = 'RFElog-%s-%s-%s.txt'%(self.obsname,self.port,self.tile)
        
        data=open(datafile,'w+')
        log=open(logfile,'w+') 
        
        header = 'tile%s-ttyUSB%s-%s-pol-%s'%(sys.argv[2],sys.argv[3],sys.argv[1],RFE.pol)
        
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
    if os.path.exists('locked-%s'%sys.argv[3]):
        print "%s:obs already in progress!"
        sys.exit()
    #open('locked-%s'%sys.argv[3],'a') 
    #os.system('mkdir '+sys.argv[1])
    RFEnew=obs(sys.argv[1],sys.argv[2],sys.argv[3])
    RFEnew.configure()
    time_obs = sys.argv[4]
    #open(recording-%s'%RFEnew.port,'a')
    RFEnew.record(float(time_obs))
    
    #open('complete-%s'%RFEnew.port,'a')
    sys.exit()
    
    
    
    # Open connection and write RFE config to file 
    #obs = RFE("/dev/ttyUSB%s"%sys.argv[1])
    #open("./Log/status_%s.log"%host,'w').write(obs.configure())
    
    # Run observations and mark completion
    #filename = '_'.join(sys.argv[1],sys.argv[2],host)
    #datafile = open("./Obs/%s.dat"%filename,'a')   #########print initial info to file
    #logfile = open("./Log/log_%s.log"%host,'a')
    #obs.record(sys.argv[2],datafile,logfile)
    #open('complete','a')
    
    ### OBSERVATION COMPLETE ###
    
