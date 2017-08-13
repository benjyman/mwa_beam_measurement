from subprocess import call
import argparse

parser = argparse.ArgumentParser(description='Schedule some ORBCOMM observations')
parser.add_argument('--time_obs', type=int, default=10,
                    help='Lendgth of observation in seconds')
parser.add_argument('--num_tiles', type=int,
                    help='Number of tiles to schedule and run - default is 8')
parser.add_argument('--date', type=str,
                    help='Date/time of the obs in the following format: hh:mm-dd-MM-YY e.g 16:15-11-Aug-2017')

args = parser.parse_args()

time_obs = args.time_obs
num_tiles = args.num_tiles
date = args.date

time,day,month,year = date.split('-')

at_script = open('run_at_jobs_%s-%s-%s-%s.sh' %(time,day,month,year),'w+')

for tile in xrange(num_tiles):
    cmd = 'echo "python RFE_record.py --time_obs=%d --tile_index=%d --date=%s-%s-%s-%s" > tile_%02d_%s.sh' %(time_obs,tile,time,day,month,year,tile,date)
    call(cmd,shell=True)
    at_script.write('at %s %s %s %s < tile_%02d_%s.sh\n' %(time,day,month,year,tile,date))
    
at_script.close()
    
cmd = "chmod +x run_at_jobs_%s-%s-%s-%s.sh" %(time,day,month,year)
call(cmd,shell=True)
    
cmd = "./run_at_jobs_%s-%s-%s-%s.sh" %(time,day,month,year)
call(cmd,shell=True)
    
cmd = "rm tile_*_%s.sh\n" %(date)
call(cmd,shell=True)
