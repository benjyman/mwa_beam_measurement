### CONVERSION OF RFEXPLORER RAW DATA TO F-ROWS FORMAT ###

"""To use, either import the function into a separate script or
program such as IDLE, or just edit the filenames using their
relative location into the function at the bottom of this script"""

import os, sys
from subprocess import Popen, PIPE

def distribute(obs_files, target='Converted/'):
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

    #Cycle through each file performing conversion
    for obsfile in obs_files:

        #Play around with filename here if desired (using header)
        header = open(obsfile).readline().strip()
        newfile = obsfile.split('/')[-1]
        newfile = target + newfile.strip('txt') + 'converted.txt'

        #Create tmp file clearing raw data file of any errors
        tmpfile = obsfile.strip('.txt') + '-tmp.txt'
        for line in open(obsfile):
            if len(line.split('$Sp'))!=2:
                continue
            else:
                open(tmpfile,'a').write(line)

        #Make first run extracting timestamps
        for line in open(tmpfile):
            #Grab stamp ensuring two digits after decimal
            stamp = line.split('$Sp')[0]
            if len(stamp.split('.')[1])!=2:
                stamp = stamp + '0'
            #Add stamp to first line in new data file
            open(newfile,'a').write(stamp)

        #Follow with 112 lines of intensity data
        for i in range(112):
            open(newfile,'a').write('\n')
            for line in open(tmpfile):
                open(newfile,'a').write(line.split('$Sp')[1][i])

        #Remove tmp data file and record success
        os.remove(tmpfile)
        print 'Complete: %s'%newfile

#Insert files for conversion here if running directly from script
if __name__ == "__main__":
    distribute(['RFEdata-2017-08-04-0-0.txt',
                'RFEdata-2017-08-04-1-1.txt'])
