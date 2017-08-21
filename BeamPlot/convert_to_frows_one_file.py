### CONVERSION OF RFEXPLORER RAW DATA TO F-ROWS FORMAT ###

"""To use, either import the function into a separate script or
program such as IDLE, or just edit the filenames using their
relative location into the function at the bottom of this script"""

import os, sys
from subprocess import Popen, PIPE

def distribute(obs_files, target='Converted/', converted_file_name=None):
    """Distribute raw obs data into 'f-rows' format, that is,
    a file with a first line containing the observation timestamps
    (10 leading digits plus 2 digits following a decimal point)
    without separating characters, then a further 112 lines
    (corresponding to 112 RFExplorer frequency channels) of data
    presented in a similar manner"""

    #Expect list, if single file given, convert to list
    if type(obs_files) == str:
        obs_files = [obs_files]

    #Make target directory for converted files if absent
    if not os.path.isdir(target):
        os.makedirs(target)

    #new file name
    new_file_name=target+converted_file_name
    if os.path.exists(new_file_name):
       cmd="rm -rf %s " % new_file_name
       os.system(cmd)
       

    #Cycle through each file performing conversion
    #Do the stamps first
    stamp_list=[]
    for obsfile in obs_files:
        #print "obsfile is %s obsfile" % obsfile
        for line in open(obsfile):
            if len(line.split('$Sp'))!=2:
                continue
            else:
                #open(tmpfile,'a').write(line)
                stamp = line.split('$Sp')[0]
                if len(stamp.split('.')[1])!=2:
                   stamp = stamp + '0'
                #print stamp
                stamp_list.append(stamp)
           
    stamp_list_string=''.join(stamp_list) + '\n' 
    #print "stamp_list_string is %s" % stamp_list_string
    #Add stamp to first line in new data file
    open(new_file_name,'w').write(stamp_list_string)
        
    #Now do for all chans
    for chan in range(0,112):
       chan_list=[]
       for obsfile in obs_files:
           for line in open(obsfile):
               if len(line.split('$Sp'))!=2:
                 continue
               else:
                 chan_data=line.split('$Sp')[1][chan]
                 #print chan_data
                 chan_list.append(chan_data)
       #print chan_list
       chan_list_string=''.join(chan_list) + '\n' 
       open(new_file_name,'a').write(chan_list_string)

       ##Follow with 112 lines of intensity data
       #for i in range(112):
       #    open(newfile,'a').write('\n')
       #    for line in open(tmpfile):
       #        open(newfile,'a').write(line.split('$Sp')[1][i])


#Insert files for conversion here if running directly from script
if __name__ == "__main__":
    distribute(['RFEdata-2017-08-04-0-0.txt',
                'RFEdata-2017-08-04-1-1.txt'])
