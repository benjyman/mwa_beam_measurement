import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from BeamPlot.Satpass    import Sateph
from BeamPlot.Settings  import lb, obs_range, obs_range_test, coords#, obs_range_test_rf0xx_01, obs_range_test_rf0xx_02, obs_range_test_rf0xx_03, obs_range_test_rf1xx,obs_range_all_rf1xx,obs_range_056xx_test
from BeamPlot.TLE import TLE_config
from BeamPlot.TimeMethods     import tm
from BeamPlot.Obstile   import Obstile
from BeamPlot.convert_to_frows_one_file import distribute
from numpy            import rad2deg
import numpy as np
import healpy as hp
np.seterr(divide='ignore')
import matplotlib.pyplot as plt
import time, datetime
import math
from reproject import reproject_from_healpix, reproject_to_healpix
from astropy.io import fits
from glob import glob

#Get info on each tile during the appropriate time period using the MWA metadata service on the web 
#http://mwa-metadata01.pawsey.org.au/metadata/obs?obs_id=1132688872
#or telescope config:
# http://mwa-metadata01.pawsey.org.au/metadata/con?obs_id=1132688872
#or get a metafits file: http://mwa-metadata01.pawsey.org.au/metadata/fits?obs_id=1132688872 

#looks like tile 51 YY has a dead dipole (lots of dead dipoles in tiles 55-58 but our data for those
#tile is no good due to off zenith pointings

#frequency for dipole and tile model
freq=137000000.0
vel_light=299792458.0
#dipole height m
height=0.3

bad_chan=86
#Define the directory of the raw data
#raw_data_dir='/data/beam/Aug2017_test/data'
#converted_frow_data_dir='./Converted'

#Set to false if data already converted
#convert_data_to_frows=True

#function to find the nearest value in an array and its index
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx,array[idx]

def nanrms(x, axis=None):
    return np.sqrt(np.nanmean(x**2, axis=axis))
    
    
def reproject_beam_to_healpix(beam_fits_filename,output_healpix_beam_filename):
    print "reprojecting %s to healpix" % beam_fits_filename
    #open fits file and get image data
    image_data=fits.open(beam_fits_filename)[0].data
    print image_data.shape
    image_data_2D=image_data[0,0,:,:]
    print image_data_2D.shape
    image_header=fits.open(beam_fits_filename)[0].header
    del image_header['NAXIS3']
    del image_header['NAXIS4']
    del image_header['CTYPE4']
    del image_header['CRVAL4']
    del image_header['CUNIT4']
    del image_header['CRPIX4']
    del image_header['CDELT4']
    del image_header['CTYPE3']
    del image_header['CRVAL3']
    del image_header['CUNIT3']
    del image_header['CRPIX3']
    del image_header['CDELT3']
    image_header['NAXIS']=2
    #Make the centre of the beam
    image_header['CRVAL1']=266.40508920
    image_header['CRVAL2']=-28.93617470
    file_name_2D=beam_fits_filename.split('.')[0] + '_2D.fits'
    fits.writeto(file_name_2D, image_data_2D, image_header,clobber=True)
    print "Using %s" % file_name_2D
    array, footprint = reproject_to_healpix(file_name_2D, 'G', hdu_in=0, order='bilinear', nested=False, nside=32)
    hp.write_map(output_healpix_beam_filename, array, overwrite=True)
    print "Made healpix image %s" % output_healpix_beam_filename
    
    
    
#define the tile name
#AUT_tile_list=["051YY","052YY","053YY","054YY","055YY","056YY","057YY","058YY","rf0YY","rf1YY"]
#AUT_tile_list=["051YY","052YY","053YY","054YY"]
#AUT_tile_list=["rf0YY"]
#AUT_tile_list=["rf0XX"]
#AUT_tile_list=["rf0YY","rf1YY"]
#AUT_tile_list=["051XX","052XX","053XX","054XX","055XX","056XX","057XX","058XX","rf0XX","rf1XX"]
#AUT_tile_list=["051XX","052XX","053XX","054XX","055XX","056XX","057XX","058XX"]

#ref_tile_list=["rf0YY","rf1YY"]
#ref_tile_list=["rf1YY","rf0YY"]
#ref_tile_list=["rf0XX"]
#ref_tile_list=["rf0XX","rf1XX"]

#set threshold for signal detection in db (noise level is around -100 to -109
#AUT_signal_threshold=-85
#AUT_signal_threshold=-50
#ref_signal_threshold=-85
#AUT_signal_threshold=-120
#ref_signal_threshold=-120
#set the how high a sat must be to be used (in deg above horizon)
#alt_threshold=30
#healpix size
nside=32

time_index_offset=0   #offset in times between ref and AUT
#ref1xx/rf0xx
#time_index_offset=-19   #offset in times between ref and AUT

#tile56xx/rf0xx (second time chunk)
#time_index_offset=-121   #offset in times between ref and AUT

plot_start_index=0
plot_end_index=4000
      
#background_level_dBm=-106
#background_level_mW=10.0**(background_level_dBm/20.0)

#fits_name="tile_map_AUT_%s_ref_%s.fits" % (AUT_tile_name,ref_tile_name)
#corrected_fits_name="tile_map_ref_corected_%s.fits" % AUT_tile_name
#rf0XX_reference_map_fits_name='tile_map_rf0XX.fits'
#fig_name="Tile_%s_map_corrected_%s.png" % (AUT_tile_name,ref_tile_name)
#corrected_fig_name="corrected_%s_tile_map.png" % (AUT_tile_name)
#sat_dictionaries_dict_filename='sat_dictionaries_dict.npy'

#sat_list=['OC-G2', 'OC-A1', 'OC-A2', 'OC-A3', 'OC-A4', 'OC-A5', 'OC-A6', 'OC-A7', 'OC-A8', 'OC-B1', 'OC-B2', 'OC-B3', 'OC-B4', 'OC-B6', 'OC-B7', 'OC-B8', 'OC-C1', 'OC-C3', 'OC-C7', 'OC-D2', 'OC-D3', 'OC-D4', 'OC-D6', 'OC-D7', 'OC-D8', 'OC-3K3', 'OC-4K4', 'OC-6K6', 'OC-7K7', 'OC-9K9', 'NOAA-15', 'NOAA-18', 'NOAA-19', 'METOP-A', 'METOP-B', 'METEOR']
#sat_list=['OC-G2', 'OC-A1', 'OC-A2', 'OC-A3', 'OC-A4', 'OC-A5', 'OC-A6', 'OC-A7', 'OC-A8', 'OC-B1', 'OC-B2', 'OC-B3', 'OC-B4', 'OC-B6', 'OC-B7', 'OC-B8', 'OC-C1', 'OC-C3', 'OC-C7', 'OC-D2', 'OC-D3', 'OC-D4', 'OC-D6', 'OC-D7', 'OC-D8', 'OC-3K3', 'OC-4K4', 'OC-6K6', 'OC-7K7', 'NOAA-15', 'NOAA-18', 'NOAA-19', 'METOP-A', 'METOP-B', 'METEOR']
useable_sat_list=[]
#sat_list=['OC-G2', 'OC-A1']
#sat_list=['OC-B2', 'OC-B6','OC-D2','OC-D3']
#sat_list=['OC-G2', 'OC-A1', 'OC-A2', 'OC-A3', 'OC-A4', 'OC-A5', 'OC-A6', 'OC-A7', 'OC-A8', 'OC-B1', 'OC-B2', 'OC-B3', 'OC-B4', 'OC-B6', 'OC-B7', 'OC-B8', 'OC-C1', 'OC-C3', 'OC-C7', 'OC-D2', 'OC-D3', 'OC-D4', 'OC-D6', 'OC-D7', 'OC-D8', 'OC-3K3', 'OC-4K4', 'OC-6K6', 'OC-7K7', 'OC-9K9', 'NOAA-15', 'NOAA-18', 'NOAA-19']
#sat_list=['OC-A1','OC-A2']
#sat_list=['NOAA-15', 'NOAA-18', 'NOAA-19','OC-G2','OC-A1']
#sat_list=['OC-A1']
#ones that look like rubbish: [OC-B3,'OC-B8','OC-C1', 'OC-C3', 'OC-C7','OC-D4','OC-D6','OC-7K7', 'OC-9K9','NOAA-18','OC-4K4']
####################only the good ones: Ben's list!
#sat_list=['OC-G2', 'OC-A1', 'OC-A2', 'OC-A3', 'OC-A4', 'OC-A5', 'OC-A6', 'OC-A7', 'OC-A8', 'OC-B1', 'OC-B2', 'OC-B4', 'OC-B6', 'OC-B7', 'OC-D2', 'OC-D3', 'OC-D7', 'OC-D8', 'OC-3K3', 'OC-6K6']
#sat_list=['OC-B2']

##########Try Jarryd's list from Sat_Data.csv #### USE this!
#sat_list=['OC-G2', 'OC-A1', 'OC-A2', 'OC-A3', 'OC-A4', 'OC-A5', 'OC-A6', 'OC-A7', 'OC-A8', 'OC-B1', 'OC-B2', 'OC-B3', 'OC-B4', 'OC-B6', 'OC-B7', 'OC-B8', 'OC-C1', 'OC-C3', 'OC-C7', 'OC-D2', 'OC-D3', 'OC-D4', 'OC-D6', 'OC-D7', 'OC-D8', 'OC-3K3', 'OC-4K4', 'OC-6K6', 'OC-7K7', 'OC-9K9', 'NOAA-15', 'NOAA-18', 'NOAA-19','METEOR']
#include 11 new sats
#sat_list=['OC-5T3','OC-8R2','OC-10T2','OC-12S3','OC-13S2','OC-14T4','OC-15R3','OC-16S1','OC-17R1','OC-18T1','OC-19S4']
#all sats
sat_list=['OC-G2', 'OC-A1', 'OC-A2', 'OC-A3', 'OC-A4', 'OC-A5', 'OC-A6', 'OC-A7', 'OC-A8', 'OC-B1', 'OC-B2', 'OC-B3', 'OC-B4', 'OC-B6', 'OC-B7', 'OC-B8', 'OC-C1', 'OC-C3', 'OC-C7', 'OC-D2', 'OC-D3', 'OC-D4', 'OC-D6', 'OC-D7', 'OC-D8', 'OC-3K3', 'OC-4K4', 'OC-6K6', 'OC-7K7', 'OC-9K9', 'NOAA-15', 'NOAA-18', 'NOAA-19','METEOR','OC-5T3','OC-8R2','OC-10T2','OC-12S3','OC-13S2','OC-14T4','OC-15R3','OC-16S1','OC-17R1','OC-18T1','OC-19S4']
#All sats except those that had no channel allocations in histograms:


#sat_list=['OC-G2', 'OC-A1', 'OC-A2', 'OC-A3']
#sat_list=['OC-A4', 'OC-A5', 'OC-A6', 'OC-A7']
#sat_list=[ 'OC-A8', 'OC-B1', 'OC-B2', 'OC-B4']
#sat_list=['OC-B6', 'OC-B7', 'OC-D2', 'OC-D3']
#bad one is in here
#sat_list=['OC-D7', 'OC-D8', 'OC-3K3']
#sat_list=['OC-3K3']
#this is the bad one!
#sat_list=['OC-D6']
#sat_list=['OC-D7']
#sat_list=['OC-D8']
#sat_list=['OC-4K4', 'OC-6K6']
#bad:
#sat_list=['OC-4K4']
#sat_list=['OC-6K6']
# 'OC-9K9' causes trouble ..?
#sat_list=['OC-7K7', 'OC-9K9']
#sat_list=['OC-7K7']
#not there:
#sat_list=['NOAA-15']
#Bad:
#sat_list=['NOAA-19']

######Bens chan dict for rf1xx/rf2/xx using all rf1xx times
#chans_dict={'NOAA-15': 66, 'NOAA-18': 89, 'NOAA-19': 24,'OC-G2':71,'OC-A1':51,'OC-A2':32,'OC-A3':34,'OC-A4':80,'OC-A5':32,'OC-A7':80,'OC-A6':51,'OC-A8':32,'OC-B1':34,'OC-B2':34,'OC-B4':51,'OC-B6':34,'OC-B7':80,'OC-C3':72,'OC-C1':36,'OC-D2':51, 'OC-D3':51,'OC-D6':68,'OC-D7':51,'OC-D8':51,'OC-3K3':36, 'OC-4K4':36, 'OC-6K6':36,'NOAA-19':24,'NOAA-15':24}

#Jarryd's Chan dict:
#chans_dict={'OC-G2':71, 'OC-A1':32, 'OC-A2':32, 'OC-A3':80, 'OC-A4':80, 'OC-A5':32, 'OC-A6':80, 'OC-A7':80, 'OC-A8':32, 'OC-B1':34, 'OC-B2':34, 'OC-B3':34, 'OC-B4':34, 'OC-B6':34, 'OC-B7':34, 'OC-B8':34, 'OC-C1':75, 'OC-C3':73, 'OC-C7':39, 'OC-D2':51, 'OC-D3':51, 'OC-D4':75, 'OC-D6':69, 'OC-D7':51, 'OC-D8':51, 'OC-3K3':36, 'OC-4K4':36, 'OC-6K6':36, 'OC-7K7':36, 'OC-9K9':36, 'NOAA-15':66, 'NOAA-18':89, 'NOAA-19':24,'METEOR':24}
#Ben chan dict with new sats
chans_dict={'OC-G2':43, 'OC-A1':np.nan, 'OC-A2':6, 'OC-A3':np.nan, 'OC-A4':np.nan, 'OC-A5':np.nan, 'OC-A6':11, 'OC-A7':41, 'OC-A8':11, 'OC-B1':np.nan, 'OC-B2':51, 'OC-B3':np.nan, 'OC-B4':46, 'OC-B6':6, 'OC-B7':11, 'OC-B8':46, 'OC-C1':46, 'OC-C3':11, 'OC-C7':11, 'OC-D2':25, 'OC-D3':23, 'OC-D4':np.nan, 'OC-D6':np.nan, 'OC-D7':23, 'OC-D8':23, 'OC-3K3':8, 'OC-4K4':8, 'OC-6K6':np.nan, 'OC-7K7':8, 'OC-9K9':8, 'NOAA-15':38, 'NOAA-18':61, 'NOAA-19':23,'METEOR':52,'OC-5T3':61,'OC-8R2':np.nan,'OC-10T2':11,'OC-12S3':41,'OC-13S2':41,'OC-14T4':11,'OC-15R3':25,'OC-16S1':41,'OC-17R1':25,'OC-18T1':11,'OC-19S4':11}


#for tile 56xx / rf0xx using rf0xx second set of times:
#chans_dict={'NOAA-15': 66, 'NOAA-18': 89, 'NOAA-19': 24,'OC-G2':71,'OC-A1':51,'OC-A2':32,'OC-A3':34,'OC-A4':80,'OC-A5':32,'OC-A7':80,'OC-A6':51,'OC-A8':32,'OC-B1':34,'OC-B2':34,'OC-B4':51,'OC-B6':34,'OC-B7':80,'OC-C3':72,'OC-C1':36,'OC-D2':51, 'OC-D3':51,'OC-D6':68,'OC-D7':51,'OC-D8':51,'OC-3K3':36, 'OC-4K4':36, 'OC-6K6':36,'NOAA-19':24,'NOAA-15':24}


def generate_pb_map(AUT_tile_name_in,ref_tile_name_in,AUT_signal_threshold_in,ref_signal_threshold_in):

   AUT_tile_name=AUT_tile_name_in
   ref_tile_name=ref_tile_name_in
   
   if ("YY" in ref_tile_name):
      polarisation="YY"
   else:
      polarisation="XX"
      
   if ("rf" in AUT_tile_name):
      AUT_signal_threshold=ref_signal_threshold_in
      ref_signal_threshold=ref_signal_threshold_in 
   else:
      AUT_signal_threshold=AUT_signal_threshold_in
      ref_signal_threshold=ref_signal_threshold_in
   
   sat_dictionaries_dict_filename='sat_dictionaries_dict_AUT_%s_ref_%s.npy' % (AUT_tile_name,ref_tile_name)
   
   fits_name="tile_map_AUT_%s_ref_%s.fits" % (AUT_tile_name,ref_tile_name)
   fits_name_no_dipole_correction="tile_map_AUT_%s_ref_%s_no_dipole_correction.fits" % (AUT_tile_name,ref_tile_name)
   #initialise healpix stuff
   AUT_tile_map_W=np.zeros(hp.nside2npix(nside))
   AUT_tile_map_data_entries_counter=np.zeros(hp.nside2npix(nside))

   ref_tile_map_W=np.zeros(hp.nside2npix(nside))
   ref_tile_map_data_entries_counter=np.zeros(hp.nside2npix(nside))
   
   dipole_model_map_W=np.zeros(hp.nside2npix(nside))
   
   #create a file to write the output to
   #with open('sat_pass_info.txt','w') as outfile:
   #  outfile.write("satellite pass info:\n")

   ###PREDICTIONS
   #initialise sateph object
   sateph = Sateph()
 

   #Set up time array, 
   #t_min=np.around(obs_range[0],0)
   #t_max=np.around(obs_range[1],0)
   t_min=np.around(obs_start_time)
   t_max=np.around(obs_end_time)
   
   
   
   #print t_min,t_max
   #predicted_time_array=np.arange(t_min,t_max)


   ###DATA
   
   #Convert to frows if not done already
   if (convert_data_to_frows):
       converted_data_filename_ref="%s_ref_%s.txt" % (ref_tile_name,converted_data_filename_base)
       print 'Looking for files %s/%s*' %(raw_data_dir,ref_tile_name)
       raw_ref_data_file_list=glob('%s/%s*' %(raw_data_dir,ref_tile_name))
       raw_ref_data_file_list=sorted(raw_ref_data_file_list)
       #print raw_ref_data_file_list
       distribute(raw_ref_data_file_list,converted_file_name=converted_data_filename_ref)             
   
   #Get the data for the ref tile
   ref_obs = Obstile([ref_tile_name],data_dir=converted_frow_data_dir)
   ref_obs.getdata()
   ref_dat = ref_obs.rdata(tile=ref_tile_name,tmin=t_min,tmax=t_max)
   
   #print ref_dat.shape
   #print t_min
   #print t_max
   
   if (convert_data_to_frows):
      converted_data_filename_AUT="%s_AUT_%s.txt" % (AUT_tile_name, converted_data_filename_base)
      raw_AUT_data_file_list=glob('%s/%s*' %(raw_data_dir,AUT_tile_name))
      raw_AUT_data_file_list=sorted(raw_AUT_data_file_list)
      #print raw_AUT_data_file_list
      distribute(raw_AUT_data_file_list,converted_file_name=converted_data_filename_AUT) 
   #select the AUT and get the data for that tile
   obs = Obstile([AUT_tile_name],data_dir=converted_frow_data_dir)
   obs.getdata()
   #print "t_min"
   #print t_min
   AUT_dat = obs.rdata(tile=AUT_tile_name,tmin=t_min,tmax=t_max)


   #print AUT_dat.shape

   data_time_array=np.array(AUT_dat[0])
   ref_time_array=np.array(ref_dat[0])
   
   #only map/plot subset of time - remove to map all time (using to plot one pass)
   #data_time_array=data_time_array[plot_start_index:plot_end_index]
   #ref_time_array=ref_time_array[plot_start_index:plot_end_index]
   
   
   
   #print data_time_array[0]

   #make a list of sats
   #sat_info=sateph.sats
   #sat_info_list=[]
   #[sat_list.append(sat.desig) for sat in sat_info]
   #print sat_info_list
   
   #print sat_list
   #create a dictionary of dictionaries, one for each satellite 
   #sat_dictionaries_dict = {sat_list[0]:{'alt':np.zeros(len(data_time_array)),'az':np.zeros(len(data_time_array)),'ref_power':[0]*len(data_time_array),'AUT_power':[0]*len(data_time_array),'chan':np.zeros(len(data_time_array))} } 
   sat_dictionaries_dict = {sat_list[0]:{'AUT_time':data_time_array,'ref_time':ref_time_array,'alt':[np.nan]*len(data_time_array),'az':[np.nan]*len(data_time_array),'ref_power':[np.nan]*len(data_time_array),'AUT_power':[np.nan]*len(data_time_array),'chan':[np.nan]*len(data_time_array)} } 
   for sat_desig in sat_list:
      #each sat has a dictionary containing an array of altitudes,azimuths, powers, chans, one for each timestep
      sat_dictionaries_dict.update({sat_desig:{'AUT_time':data_time_array,'ref_time':ref_time_array,'alt':[np.nan]*len(data_time_array),'az':[np.nan]*len(data_time_array),'ref_power':[np.nan]*len(data_time_array),'AUT_power':[np.nan]*len(data_time_array),'chan':[np.nan]*len(data_time_array)}})
    
      #populate the altitude
      sat_dictionaries_dict[sat_desig]['alt']=sateph.get_sat_alt_az(sat_desig,data_time_array)[1]
   
      #populate the azimuth
      sat_dictionaries_dict[sat_desig]['az']=sateph.get_sat_alt_az(sat_desig,data_time_array)[2]
   

   #print max(sat_dictionaries_dict['OC-A1']['alt'])
   #print max(sat_dictionaries_dict['OC-A1']['az'])

   #split the time up into chunks of about 600 timesteps 
   #length_time_array=len(data_time_array)
   #n_chunks=int(length_time_array/600)
   
   #data_time_chunks=np.array_split(data_time_array)


   #Go through each time step (got time ranges from Python.Settings, imported above)
   for index, timestep in enumerate(data_time_array):   
      #check that the ref and AUT times are almost the same?
      ref_timestamp_index,ref_timestamp=find_nearest(ref_time_array, timestep)
      #ref_timestamp_index,ref_timestamp = index, timestep
      timediff=ref_timestamp-timestep
      #go through the dict of satellites and make an ordered list of sats with alt greater than 30 deg
      sats_above_30_list = []
      sat_alt_list=[]
      for sat_desig in sat_list:
         alt=sat_dictionaries_dict[sat_desig]['alt'][index]
         #print alt
         if (alt >= alt_threshold):
            sats_above_30_list.append(sat_desig)
            sat_alt_list.append(alt)
            if (sat_desig not in useable_sat_list):
               useable_sat_list.append(sat_desig)
   
      #now find which ones have a detectable signal     
      if sats_above_30_list != []:
         print sats_above_30_list
         #sort the list in order highest alt to lowest alt
         sats_above_30_list=np.array(sats_above_30_list)
         sat_alt_list=np.array(sat_alt_list)
         inds = sat_alt_list.argsort()
         sorted_sats_above_30_list=sats_above_30_list[inds]
         sats_above_30_with_signal_list=[]
         #print sats_above_30_list
         #find the power values (for each chan) for the timestep we are in
         #timestamp_index,timestamp=find_nearest(dat[0], timestep)
         #powers=[]
         AUT_powers=[AUT_dat[1+chan][index] for chan in range(0,112)]
         ref_powers=[ref_dat[1+chan][ref_timestamp_index] for chan in range(0,112)]
         #set bad chan(s) to small power
         AUT_powers[bad_chan]=-110
         ref_powers[bad_chan]=-110
         AUT_powers[bad_chan+1]=-110
         ref_powers[bad_chan+1]=-110
         AUT_powers[bad_chan-1]=-110
         ref_powers[bad_chan-1]=-110
         #doing this doesn't work for some reason....
         #ref_powers=[ref_dat[1+chan][index] for chan in range(0,112)]
         #print sats_above_30_list
         for sat_index,sat_above_30 in enumerate(sorted_sats_above_30_list):
            #print timestep
            #print index
            #print timestamp
            #if (sat_above_30=='OC-A1' or sat_above_30=='OC-A2' or sat_above_30=='OC-A3'or sat_above_30=='OC-A4'or sat_above_30=='OC-A5'or sat_above_30=='OC-A7' or sat_above_30=='OC-A6'or sat_above_30=='OC-A8' or sat_above_30=='OC-B1'or sat_above_30=='OC-B2'or sat_above_30=='OC-B4'or sat_above_30=='OC-B7'or sat_above_30=='OC-B6' or sat_above_30=='OC-C3'or sat_above_30=='OC-C1'or sat_above_30=='OC-D2'or sat_above_30=='OC-D3'or sat_above_30=='OC-D6'or sat_above_30=='OC-D7'or sat_above_30=='OC-D8' or sat_above_30=='NOAA-19'or sat_above_30=='NOAA-15'or sat_above_30=='OC-3K3'or sat_above_30=='OC-4K4'or sat_above_30=='OC-6K6' or sat_above_30=='OC-B3'or sat_above_30=='OC-G2'or sat_above_30=='NOAA-18'or sat_above_30=='OC-B8'or sat_above_30=='OC-C7'or sat_above_30=='OC-D4'or sat_above_30=='OC-7K7'or sat_above_30=='OC-9K9'or sat_above_30=='METEOR'):
            #if (1 == 2):
            if (sat_above_30 in sat_list):
               channel=chans_dict[sat_above_30]
               if np.isnan(channel):
                  continue
               AUT_max_power=AUT_powers[channel]
               #AUT_max_power_mW=10.0**(AUT_max_power/20.0)
               #AUT_max_power_dB_rel_background=10.0*np.log10(AUT_max_power_mW/background_level_mW)
               AUT_max_power_chan=channel
               ref_max_power=ref_powers[channel]
               ref_max_power_chan=channel
               print 'Allocated chan for %s is: %s' % (sat_above_30,AUT_max_power_chan)
            else:
               continue
               #AUT_max_power=np.max(AUT_powers)
               #AUT_max_power_chan=np.argmax(AUT_powers)
               #ref_max_power=np.max(ref_powers)
               #ref_max_power_chan=np.argmax(ref_powers)
               #print 'Max power chan is: %s' % AUT_max_power_chan
            if (AUT_max_power > AUT_signal_threshold and ref_max_power>ref_signal_threshold):  
               #print 'max power:%s in chan:%s ' %  (max_power,max_power_chan)
               #populate the power in the dict of dicts!
               #print sat_dictionaries_dict[sat_above_30]['power']
               #print timestamp_index
               sat_dictionaries_dict[sat_above_30]['AUT_power'][index]=AUT_max_power
               sat_dictionaries_dict[sat_above_30]['chan'][index]=int(AUT_max_power_chan)
               sat_dictionaries_dict[sat_above_30]['ref_power'][index]=ref_max_power
               #if (ref_timestamp_index+time_index_offset < len(sat_dictionaries_dict[sat_above_30]['ref_power'])):
                  #print len(sat_dictionaries_dict[sat_above_30]['ref_power'])
                  #t= ref_timestamp_index+time_index_offset
                  #print t
                  #sat_dictionaries_dict[sat_above_30]['ref_power'][ref_timestamp_index+time_index_offset]=ref_max_power
                  
               #print 'sat %s on chan %s ' % (sat_above_30,AUT_max_power_chan)
               #set that max power to a very small value, so we can find the next highest in the next iteration
               #set one channel either side as well to make sure it is not just the same signal spread across multiple chans
               AUT_powers[AUT_max_power_chan]=-10000
               if (AUT_max_power_chan>0):
                  AUT_powers[AUT_max_power_chan-1]=-10000
                  ref_powers[AUT_max_power_chan-1]=-10000
               if (AUT_max_power_chan<len(AUT_powers)-1):
                  AUT_powers[AUT_max_power_chan+1]=-10000
                  ref_powers[AUT_max_power_chan-1]=-10000
               sats_above_30_with_signal_list.append(sat_above_30)
      
         if (sats_above_30_with_signal_list != []):
            #go through the above 30 list again and populate the healpix map for this timestep for each sat
            for sat_with_signal in sats_above_30_with_signal_list:
               #print 'sat %s on chan %s with ref power %s and AUT power %s at AUT time %s and ref time %s (time diff %s)' % (sat_above_30,AUT_max_power_chan,ref_max_power,AUT_max_power,timestep,ref_timestamp,timediff)
               #get a healpy pixel number ang2pix(alt_rad,az_rad)
               data_alt=sat_dictionaries_dict[sat_with_signal]['alt'][index]
               data_alt_rad=data_alt/180.0*np.pi
               data_theta_rad = (np.pi/2)-data_alt_rad
               data_az_deg=sat_dictionaries_dict[sat_with_signal]['az'][index]
               data_az_rad=data_az_deg/180.0*np.pi
               data_spherical_az_rad=np.pi-data_az_rad
               healpix_pixnum=hp.ang2pix(nside, data_theta_rad, data_spherical_az_rad)
               #print healpix_pixnum
               AUT_data_pt=sat_dictionaries_dict[sat_with_signal]['AUT_power'][index]
               ref_data_pt=sat_dictionaries_dict[sat_with_signal]['ref_power'][index]
               ##I think this bit was wrong; just want to put the corresponding ref time in the array at 'index' so it matches with the AUT
               #if (ref_timestamp_index+time_index_offset < len(sat_dictionaries_dict[sat_above_30]['ref_power'])):
               #   ref_data_pt=sat_dictionaries_dict[sat_with_signal]['ref_power'][ref_timestamp_index+time_index_offset]
               #print 'aut %s, ref %s' % (AUT_data_pt,ref_data_pt)
               
               #check to see why ref and AUT powers are different even if the ref is set as the AUT!
               if (AUT_data_pt != ref_data_pt):
                  print 'sat %s on chan %s with ref power %s and AUT power %s at AUT time %s and ref time %s (time diff %s)' % (sat_above_30,AUT_max_power_chan,ref_max_power,AUT_max_power,timestep,ref_timestamp,timediff)
               
               data_pt=AUT_data_pt-ref_data_pt
               #data_pt=sat_dictionaries_dict[sat_with_signal]['AUT_power'][index]
               #print data_pt
               #convert data pt from db to W
               power_W=10.0**(data_pt/10.0)
               #add the data point in W to the map
               AUT_tile_map_W[healpix_pixnum]+=power_W
               #print tile_map[healpix_pixnum]
               #keep track of how many data points you have added to this pixel
               AUT_tile_map_data_entries_counter[healpix_pixnum]+=1
               #print tile_map_data_entries_counter[healpix_pixnum]

   #save the dictionary so we can plot the data later
   np.save(sat_dictionaries_dict_filename,sat_dictionaries_dict)
   #print 'useable sats are:'
   #print useable_sat_list
   
   #divide each pixel value by the number of data entries to form the average
   #print min(tile_map)
   #print max(tile_map_data_entries_counter)
   
   #make the dipole model map
   for pixel_number in range(0,len(dipole_model_map_W)):
      dipole_model_theta_rad, dipole_model_spherical_az_rad = hp.pix2ang(nside,pixel_number)
      dipole_model_az_rad=np.pi-dipole_model_spherical_az_rad
      dipole_model_zenith_angle_rad=dipole_model_theta_rad
      
      
      ##Don't forget the groundscreen. From mwapy/mwa_tile.py:
      ##    def groundScreen(self,za,freq):
      ##  """
      ##  Calculate the groundscreen effect for an ideal infinite groundscreen
      ##  given the dipole's height above the screen and the frequency (Hz)
      ##  """
      ##  l = vel_light/freq
      ## return numpy.sin(numpy.pi * (2.0*self.height/l) * numpy.cos(za))*2.0
      
      wavelength=vel_light/freq
      ground_screen=np.sin(np.pi * (2.0*height/wavelength) * np.cos(dipole_model_zenith_angle_rad))*2.0
      
      #YY
      theta_YY=np.arccos(np.sin(dipole_model_zenith_angle_rad)*np.cos(dipole_model_az_rad))
      voltage_YY=np.sin(theta_YY)*ground_screen
      power_YY=voltage_YY**2
      
      #XX
      theta_XX=np.arccos(np.sin(dipole_model_zenith_angle_rad)*np.sin(dipole_model_az_rad))
      voltage_XX=np.sin(theta_XX)*ground_screen
      power_XX=voltage_XX**2
   
      if (polarisation=='YY'):
         dipole_model_map_W[pixel_number]=power_YY
      else:
         dipole_model_map_W[pixel_number]=power_XX
         
   #convert to db
   dipole_model_map_db=10.0*np.log10(dipole_model_map_W)
   
   #print AUT_tile_map_W
   AUT_tile_map_W_av=np.divide(AUT_tile_map_W,AUT_tile_map_data_entries_counter)
   #convert back to db
   AUT_tile_map_dB_av=10.0*np.log10(AUT_tile_map_W_av)
   
   #Multiply by the dipole model (add in log space)
   AUT_tile_map_dB_av_dipole_corrected=AUT_tile_map_dB_av+dipole_model_map_db
   
   #print AUT_tile_map_dB_av
   #if 'rf0' in tile_name:
   #   rf0XX_tile_map_W_av=tile_map_dB_av
   #print rf0XX_tile_map_av
   #if 'rf1' in tile_name:
   #   rf1XX_tile_map_av=tile_map_av
   
   #hp.gnomview(tile_map_av, coord='C',reso=60,xsize=400, title=map_title, rot=(0,90,0))  
   hp.write_map(fits_name, AUT_tile_map_dB_av_dipole_corrected,dtype=np.float32, overwrite=True)
   hp.write_map(fits_name_no_dipole_correction,AUT_tile_map_dB_av,dtype=np.float32, overwrite=True)
   
def plot_pb_map(AUT_tile_name_in,ref_tile_name_in):

   AUT_tile_name=AUT_tile_name_in
   ref_tile_name=ref_tile_name_in
   if ('XX' in ref_tile_name):
      polarisation='XX'
   else:
      polarisation='YY'
   fits_name="tile_map_AUT_%s_ref_%s.fits" % (AUT_tile_name,ref_tile_name)
   fits_name_no_dipole_correction="tile_map_AUT_%s_ref_%s_no_dipole_correction.fits" % (AUT_tile_name,ref_tile_name)
   fig_name="Tile_%s_map_corrected_%s.png" % (AUT_tile_name,ref_tile_name)
   fig_name_no_dipole_correction="Tile_%s_map_%s_no_dipole_correction.png" % (AUT_tile_name,ref_tile_name)
   AUT_tile_map_av = hp.read_map(fits_name)
   AUT_tile_map_av_no_dipole_correction = hp.read_map(fits_name_no_dipole_correction)
   map_title="AUT: Tile%s Ref:%s" % (AUT_tile_name,ref_tile_name)
   map_title_no_dipole_correction="No dipole Correction AUT: Tile%s Ref:%s" % (AUT_tile_name,ref_tile_name)
   #hp.gnomview(AUT_tile_map_av, coord='C',reso=60,xsize=400, title=map_title, rot=(0,90,0))  
   hp.orthview(map=AUT_tile_map_av,coord='E',half_sky=True,xsize=400,title=map_title,rot=(0,90,0))
   figmap = plt.gcf()
   if os.path.exists(fig_name):
      cmd = "rm -f %s " % fig_name
      os.system(cmd)
   figmap.savefig(fig_name,dpi=200)
   plt.clf()
   
   #no dipole correction
   hp.orthview(map=AUT_tile_map_av_no_dipole_correction,coord='C',half_sky=True,xsize=400,title=map_title_no_dipole_correction,rot=(0,90,0))
   figmap = plt.gcf()
   if os.path.exists(fig_name_no_dipole_correction):
      cmd = "rm -f %s " % fig_name_no_dipole_correction
      os.system(cmd)
   figmap.savefig(fig_name_no_dipole_correction,dpi=200)
   plt.clf()
   
   #hp.mollview(tile_map_av, title=map_title)   
   #hp.gnomview(tile_map_av, coord='C',reso=60,xsize=400, title=map_title, rot=(0,90,0))  
   #hp.gnomview(tile_map, coord='C',reso=60,xsize=400, title=map_title, rot=(0,90,0))

   #hp.write_map(fits_name, tile_map_av, overwrite=True)
       
   #reference_corrected_maps
   #ref_corrected_tile_map_av=tile_map_av-rf0XX_tile_map_av  #minus in log space (dB) is equivalent to division of raw numbers
   #corrected_map_title="Reference-corrected Tile Map for %s" % tile_name
   #corrected_fig_name="Tile_%s_referenceRF0_corrected.png" % (tile_name)
   #hp.gnomview(ref_corrected_tile_map_av, coord='C',reso=60,xsize=400, title=corrected_map_title, rot=(0,90,0))  
   #plt.savefig(corrected_fig_name,dpi=200)

def plot_simulated_pb_map(obs_start_time,obs_end_time,t_step,observatory,polarisation):

   fits_name="simulated_AUT_and_ref_from_%s_to_%s_at_%s_%s.fits" % (obs_start_time,obs_end_time,observatory,polarisation)
   fits_name_no_dipole_correction="simulated_AUT_and_ref_from_%s_to_%s_at_%s_%s_no_dipole_correction.fits" % (obs_start_time,obs_end_time,observatory,polarisation)
   counter_fits_name="simulated_data_pt_counts_from_%s_to_%s_at_%s_%s.fits" % (obs_start_time,obs_end_time,observatory,polarisation)
   dipole_model_name="short_dipole_model_%s.fits" % polarisation
   tile_model_name="healpix_beam_map_%s.fits" % polarisation
   fig_name="Simulated_pb_from_%s_to_%s_at_%s_%s.png" % (obs_start_time,obs_end_time,observatory,polarisation)
   fig_name_no_dipole_correction="Simulated_pb_from_%s_to_%s_at_%s_%s_no_dipole_correction.png" % (obs_start_time,obs_end_time,observatory,polarisation)
   counter_fig_name="simulated_data_pt_counts_from_%s_to_%s_at_%s_%s.png" % (obs_start_time,obs_end_time,observatory,polarisation)
   dipole_model_fig_name="short_dipole_model_%s.png" % polarisation
   tile_model_fig_name="tile_model_%s.png" % polarisation
   tile_model_dB_fig_name="tile_model_dB%s.png" % polarisation
   AUT_tile_map_av = hp.read_map(fits_name)
   AUT_tile_map_av_no_dipole_correction = hp.read_map(fits_name_no_dipole_correction)
   counter_tile_map= hp.read_map(counter_fits_name)
   dipole_model_map=hp.read_map(dipole_model_name)
   tile_model_map=hp.read_map(tile_model_name)
   tile_model_dB_map=10.0*np.log10(tile_model_map)
   map_title="Simulated_pb_%s_to_%s_at_%s_%s" % (obs_start_time,obs_end_time,observatory,polarisation)
   map_title_no_dipole_correction="No dipole Correction \n Simulated_pb_%s_to_%s_at_%s_%s" % (obs_start_time,obs_end_time,observatory,polarisation)
   counter_map_title="Simulated_data_pt_counts_%s_to_%s_at_%s_%s" % (obs_start_time,obs_end_time,observatory,polarisation)
   dipole_model_map_title="Short dipole model %s Polarisation" % polarisation
   tile_model_map_title="Tile model %s Polarisation" % polarisation
   tile_model_dB_map_title="Tile model dB %s Polarisation" % polarisation
   
   #hp.gnomview(AUT_tile_map_av, coord='C',reso=60,xsize=400, title=map_title, rot=(0,90,0))  
   hp.orthview(map=AUT_tile_map_av,coord='C',half_sky=True,xsize=400,title=map_title,rot=(0,90,0))
   figmap = plt.gcf()
   figmap.savefig(fig_name,dpi=200)
   plt.clf()
   
   #counts
   hp.orthview(map=counter_tile_map,coord='C',half_sky=True,xsize=400,title=counter_map_title,rot=(0,90,0))
   figmap = plt.gcf()
   figmap.savefig(counter_fig_name,dpi=200)
   plt.clf()
   
   #dipole model
   hp.orthview(map=dipole_model_map,coord='C',half_sky=True,xsize=400,title=dipole_model_map_title,min=-20,rot=(0,90,0))
   figmap = plt.gcf()
   figmap.savefig(dipole_model_fig_name,dpi=200)
   plt.clf()
   
   #no dipole correction map
   hp.orthview(map=AUT_tile_map_av_no_dipole_correction,coord='C',half_sky=True,xsize=400,title=map_title_no_dipole_correction,rot=(0,90,0))
   figmap = plt.gcf()
   figmap.savefig(fig_name_no_dipole_correction,dpi=200)
   plt.clf()

   #tile model (rot lon,lat,psi)
   hp.orthview(map=tile_model_map,coord='C',half_sky=True,xsize=400,title=tile_model_map_title,rot=(0,0,301.401202))
   figmap = plt.gcf()
   figmap.savefig(tile_model_fig_name,dpi=200)
   plt.clf()
   
   #tile model db (rot lon,lat,psi)
   rotated_tile_model_dB_map=hp.orthview(map=tile_model_dB_map,coord='C',half_sky=True,xsize=400,title=tile_model_map_title,rot=(0,0,301.401202),min=-40,return_projected_map=True)
   figmap = plt.gcf()
   figmap.savefig(tile_model_dB_fig_name,dpi=200)
   plt.clf()
   
   ##convert the rotated 2D numpy array back to healpix
   ##get header form original beam fits file
   #beam_filename_YY='beam_model/1162046832_cotter_zenith_hacked-image_beamYY.fits'
   #image_header=fits.open(beam_filename_YY)[0].header
   #del image_header['NAXIS3']
   #del image_header['NAXIS4']
   #del image_header['CTYPE4']
   #del image_header['CRVAL4']
   #del image_header['CUNIT4']
   #del image_header['CRPIX4']
   #del image_header['CDELT4']
   #del image_header['CTYPE3']
   #del image_header['CRVAL3']
   #del image_header['CUNIT3']
   #del image_header['CRPIX3']
   #del image_header['CDELT3']
   #image_header['NAXIS']=2
   ##Make the centre of the beam
   #image_header['CRVAL1']=266.40508920
   #image_header['CRVAL2']=-28.93617470

   #rotated_array, footprint = reproject_to_healpix((rotated_tile_model_dB_map,image_header), 'G', hdu_in=0, order='bilinear', nested=False, nside=32)
   
   ##SHould be able to plot with no rotation now....
   #hp.orthview(map=rotated_array,coord='C',half_sky=True,xsize=400,title=tile_model_map_title,rot=(0,0,0),min=-80)
   #figmap = plt.gcf()
   #figmap.savefig("tile_model_dB_rotated_hpx%s.png" % polarisation,dpi=200)
   #plt.clf()
   
   #show the rotated tile model map
   #plt.imshow(rotated_tile_model_dB_map)
   #figmap = plt.gcf()
   #figmap.savefig("tile_model_dB_rotated_%s.png" % polarisation,dpi=200)
   #plt.clf()


#Simulate a pb map from predicted satellite passes
def simulate_pb_map(obs_start_time,obs_end_time,t_step,observatory,polarisation,noise):

#   AUT_tile_name=AUT_tile_name_in
#   ref_tile_name=ref_tile_name_in
#      
   #using altitude as a proxy for power for now
   #if ("rf" in AUT_tile_name):
   #   AUT_signal_threshold=alt_threshold
   #   ref_signal_threshold=alt_threshold 
   #else:
   #   AUT_signal_threshold=alt_threshold
   #   ref_signal_threshold=alt_threshold
#   

   sat_dictionaries_dict_filename="sat_dictionaries_dict_from_%s_to_%s_at_%s_%s.npy" % (obs_start_time,obs_end_time,observatory,polarisation)
#   
   fits_name="simulated_AUT_and_ref_from_%s_to_%s_at_%s_%s.fits" % (obs_start_time,obs_end_time,observatory,polarisation)
   fits_name_no_dipole_correction="simulated_AUT_and_ref_from_%s_to_%s_at_%s_%s_no_dipole_correction.fits" % (obs_start_time,obs_end_time,observatory,polarisation)
   counter_fits_name="simulated_data_pt_counts_from_%s_to_%s_at_%s_%s.fits" % (obs_start_time,obs_end_time,observatory,polarisation)
   dipole_model_name="short_dipole_model_%s.fits" % polarisation
#   
#   #initialise healpix stuff
   AUT_tile_map_W=np.zeros(hp.nside2npix(nside))
   AUT_tile_map_data_entries_counter=np.zeros(hp.nside2npix(nside))

   ref_tile_map_W=np.zeros(hp.nside2npix(nside))
   ref_tile_map_data_entries_counter=np.zeros(hp.nside2npix(nside))
   
   dipole_model_map_W=np.zeros(hp.nside2npix(nside))


#   ###PREDICTIONS
#   #initialise sateph object
   sateph = Sateph()
 

   #Set up time array, 
   #t_min=np.around(obs_range[0],0)
   #t_max=np.around(obs_range[1],0)
   t_min=np.around(obs_start_time)
   t_max=np.around(obs_end_time)
   
   #print t_min,t_max
   predicted_time_array=np.arange(t_min,t_max,t_step)
   
   if (noise>0):
      AUT_noise_array = np.random.normal(0,noise,len(predicted_time_array))
      ref_noise_array = np.random.normal(0,noise,len(predicted_time_array))
   else:
      AUT_noise_array = np.zeros(len(predicted_time_array))
      ref_noise_array = np.zeros(len(predicted_time_array))
   
   #create a dictionary of dictionaries, one for each satellite 
   #sat_dictionaries_dict = {sat_list[0]:{'alt':np.zeros(len(data_time_array)),'az':np.zeros(len(data_time_array)),'ref_power':[0]*len(data_time_array),'AUT_power':[0]*len(data_time_array),'chan':np.zeros(len(data_time_array))} } 
   sat_dictionaries_dict = {sat_list[0]:{'AUT_time':predicted_time_array,'ref_time':predicted_time_array,'alt':[np.nan]*len(predicted_time_array),'az':[np.nan]*len(predicted_time_array),'ref_power':[np.nan]*len(predicted_time_array),'AUT_power':[np.nan]*len(predicted_time_array),'chan':[np.nan]*len(predicted_time_array),'dipole_model_power':[np.nan]*len(predicted_time_array)} } 
   for sat_desig in sat_list:
      #each sat has a dictionary containing an array of altitudes,azimuths, powers, chans, one for each timestep
      sat_dictionaries_dict.update({sat_desig:{'AUT_time':predicted_time_array,'ref_time':predicted_time_array,'alt':[np.nan]*len(predicted_time_array),'az':[np.nan]*len(predicted_time_array),'ref_power':[np.nan]*len(predicted_time_array),'AUT_power':[np.nan]*len(predicted_time_array),'chan':[np.nan]*len(predicted_time_array),'dipole_model_power':[np.nan]*len(predicted_time_array)}})
    
      #populate the altitude
      altitudes_deg=sateph.get_sat_alt_az(sat_desig,predicted_time_array)[1]
      altitudes_rad=(altitudes_deg/180.0)*np.pi
      zenith_angles_rad=(np.pi/2.0)-altitudes_rad
      sat_dictionaries_dict[sat_desig]['alt']=altitudes_deg
   
      #populate the azimuth
      azimuths_deg=sateph.get_sat_alt_az(sat_desig,predicted_time_array)[2]
      azimuths_rad=(azimuths_deg/180.0)*np.pi
      sat_dictionaries_dict[sat_desig]['az']=azimuths_deg
   
      #populate the AUT power as the altitude 
      #powers=sateph.get_sat_alt_az(sat_desig,predicted_time_array)[1]
      #low_values_indices = powers < alt_threshold
      #powers[low_values_indices]=np.nan
      #sat_dictionaries_dict[sat_desig]['AUT_power']=powers
      
      #populate the AUT power, ref power and dipole model assuming a short dipole (see "New comparison of MWA tile beams" by Benjamin McKinley on Twiki)
      #YY
      thetas_YY=np.arccos(np.sin(zenith_angles_rad)*np.cos(azimuths_rad))
      voltages_YY=np.sin(thetas_YY)
      powers_YY=voltages_YY**2
      #only populate power where alt is greater than alt_threshold
      low_values_indices = altitudes_deg < alt_threshold
      powers_YY[low_values_indices]=np.nan
      
      #XX
      thetas_XX=np.arccos(np.sin(zenith_angles_rad)*np.sin(azimuths_rad))
      voltages_XX=np.sin(thetas_XX)
      powers_XX=voltages_XX**2
      #only populate power where alt is greater than alt_threshold
      low_values_indices = altitudes_deg < alt_threshold
      powers_XX[low_values_indices]=np.nan
      
      if (polarisation=='YY'):
         sat_dictionaries_dict[sat_desig]['AUT_power']=powers_YY + AUT_noise_array
      else:
         sat_dictionaries_dict[sat_desig]['AUT_power']=powers_XX + AUT_noise_array
         
      #populate the ref power as the same as the AUT power for now (AUT is a ref)
      if (polarisation=='YY'):
         sat_dictionaries_dict[sat_desig]['ref_power']=powers_YY + ref_noise_array
      else:
         sat_dictionaries_dict[sat_desig]['ref_power']=powers_XX + ref_noise_array
      
      #dipole model
      if (polarisation=='YY'):
         sat_dictionaries_dict[sat_desig]['dipole_model_power']=powers_YY
      else:
         sat_dictionaries_dict[sat_desig]['dipole_model_power']=powers_XX


   #Go through each time step
   for index, timestep in enumerate(predicted_time_array):   
      #check that the ref and AUT times are almost the same?
      #ref_timestamp_index,ref_timestamp=find_nearest(ref_time_array, timestep)
      #ref_timestamp_index,ref_timestamp = index, timestep
      #timediff=ref_timestamp-timestep
      #go through the dict of satellites and make an ordered list of sats with alt greater than 30 deg
      sats_above_30_list = []
      sat_alt_list=[]
      for sat_desig in sat_list:
         alt=sat_dictionaries_dict[sat_desig]['alt'][index]
         #print alt
         if (alt >= alt_threshold):
            sats_above_30_list.append(sat_desig)
            sat_alt_list.append(alt)
            if (sat_desig not in useable_sat_list):
               useable_sat_list.append(sat_desig)
      
      #Don't worry about finding which ones have a detectable signal for now (if they are above the alt threshold, they do
      sats_above_30_with_signal_list=sats_above_30_list
      #.....now find which ones have a detectable signal     
      #if sats_above_30_list != []:
      #   print sats_above_30_list
      #   #sort the list in order highest alt to lowest alt
      #   sats_above_30_list=np.array(sats_above_30_list)
      #   sat_alt_list=np.array(sat_alt_list)
      #   inds = sat_alt_list.argsort()
      #   sorted_sats_above_30_list=sats_above_30_list[inds]
      #   sats_above_30_with_signal_list=[]
      #   #print sats_above_30_list
      #   #find the power values (for each chan) for the timestep we are in
      #   #timestamp_index,timestamp=find_nearest(dat[0], timestep)
      #   #powers=[]
      #   #just use altitude as a proxy for power at this point
      #   AUT_powers=[AUT_dat[1+chan][index] for chan in range(0,112)]
      #   ref_powers=[ref_dat[1+chan][ref_timestamp_index] for chan in range(0,112)]
      #   #doing this doesn't work for some reason....
      #   #ref_powers=[ref_dat[1+chan][index] for chan in range(0,112)]
      #   #print sats_above_30_list
      #   for sat_index,sat_above_30 in enumerate(sorted_sats_above_30_list):
      #      #print timestep
      #      #print index
      #      #print timestamp
      #      if (sat_above_30=='OC-A1' or sat_above_30=='OC-A2' or sat_above_30=='OC-A3'or sat_above_30=='OC-A4'or sat_above_30=='OC-A5'or sat_above_30=='OC-A7' or sat_above_30=='OC-A6'or sat_above_30=='OC-A8' or sat_above_30=='OC-B1'or sat_above_30=='OC-B2'or sat_above_30=='OC-B4'or sat_above_30=='OC-B7'or sat_above_30=='OC-B6' or sat_above_30=='OC-C3'or sat_above_30=='OC-C1'or sat_above_30=='OC-D2'or sat_above_30=='OC-D3'or sat_above_30=='OC-D6'or sat_above_30=='OC-D7'or sat_above_30=='OC-D8' or sat_above_30=='NOAA-19'or sat_above_30=='NOAA-15'or sat_above_30=='OC-3K3'or sat_above_30=='OC-4K4'or sat_above_30=='OC-6K6'):
      #         #if (1 == 2):
      #         channel=chans_dict[sat_above_30]
      #         AUT_max_power=AUT_powers[channel]
      #         #AUT_max_power_mW=10.0**(AUT_max_power/20.0)
      #         #AUT_max_power_dB_rel_background=10.0*np.log10(AUT_max_power_mW/background_level_mW)
      #         AUT_max_power_chan=channel
      #         ref_max_power=ref_powers[channel]
      #         ref_max_power_chan=channel
      #      else:
      #         AUT_max_power=np.max(AUT_powers)
      #         AUT_max_power_chan=np.argmax(AUT_powers)
      #         ref_max_power=np.max(ref_powers)
      #         ref_max_power_chan=np.argmax(ref_powers)
      #      print 'Max power chan is: %s' % AUT_max_power_chan
      #      if (AUT_max_power > AUT_signal_threshold and ref_max_power>ref_signal_threshold):  
      #         #print 'max power:%s in chan:%s ' %  (max_power,max_power_chan)
      #         #populate the power in the dict of dicts!
      #         #print sat_dictionaries_dict[sat_above_30]['power']
      #         #print timestamp_index
      #         sat_dictionaries_dict[sat_above_30]['AUT_power'][index]=AUT_max_power
      #         sat_dictionaries_dict[sat_above_30]['chan'][index]=int(AUT_max_power_chan)
      #         sat_dictionaries_dict[sat_above_30]['ref_power'][index]=ref_max_power
      #         #if (ref_timestamp_index+time_index_offset < len(sat_dictionaries_dict[sat_above_30]['ref_power'])):
      #           #print len(sat_dictionaries_dict[sat_above_30]['ref_power'])
      #           #t= ref_timestamp_index+time_index_offset
      #            #print t
      #            #sat_dictionaries_dict[sat_above_30]['ref_power'][ref_timestamp_index+time_index_offset]=ref_max_power
      #            
      #         #print 'sat %s on chan %s ' % (sat_above_30,AUT_max_power_chan)
      #         #set that max power to a very small value, so we can find the next highest in the next iteration
      #         #set one channel either side as well to make sure it is not just the same signal spread across multiple chans
      #         AUT_powers[AUT_max_power_chan]=-10000
      #         if (AUT_max_power_chan>0):
      #            AUT_powers[AUT_max_power_chan-1]=-10000
      #            ref_powers[AUT_max_power_chan-1]=-10000
      #         if (AUT_max_power_chan<len(AUT_powers)-1):
      #            AUT_powers[AUT_max_power_chan+1]=-10000
      #            ref_powers[AUT_max_power_chan-1]=-10000
      #         sats_above_30_with_signal_list.append(sat_above_30)
      
      if (sats_above_30_with_signal_list != []):
            #go through the above 30 list again and populate the healpix map for this timestep for each sat
            for sat_with_signal in sats_above_30_with_signal_list:
               #print 'sat %s on chan %s with ref power %s and AUT power %s at AUT time %s and ref time %s (time diff %s)' % (sat_above_30,AUT_max_power_chan,ref_max_power,AUT_max_power,timestep,ref_timestamp,timediff)
               #get a healpy pixel number ang2pix(alt_rad,az_rad)
               data_alt=sat_dictionaries_dict[sat_with_signal]['alt'][index]
               data_alt_rad=data_alt/180.0*np.pi
               data_theta_rad = (np.pi/2.0)-data_alt_rad
               data_az_deg=sat_dictionaries_dict[sat_with_signal]['az'][index]
               data_az_rad=data_az_deg/180.0*np.pi
               data_spherical_az_rad=np.pi-data_az_rad
               healpix_pixnum=hp.ang2pix(nside, data_theta_rad, data_spherical_az_rad)
               #print healpix_pixnum
               AUT_data_pt=sat_dictionaries_dict[sat_with_signal]['AUT_power'][index]
               ref_data_pt=sat_dictionaries_dict[sat_with_signal]['ref_power'][index]
               dipole_model_data_pt=sat_dictionaries_dict[sat_with_signal]['dipole_model_power'][index]
               ##I think this bit was wrong; just want to put the corresponding ref time in the array at 'index' so it matches with the AUT
               #if (ref_timestamp_index+time_index_offset < len(sat_dictionaries_dict[sat_above_30]['ref_power'])):
               #   ref_data_pt=sat_dictionaries_dict[sat_with_signal]['ref_power'][ref_timestamp_index+time_index_offset]
               #print 'aut %s, ref %s' % (AUT_data_pt,ref_data_pt)
               data_pt=AUT_data_pt-ref_data_pt
               #data_pt=sat_dictionaries_dict[sat_with_signal]['AUT_power'][index]
               #print data_pt
               #convert data pt from db to W
               power_W=10.0**(data_pt/20.0)
               #add the data point in W to the map
               AUT_tile_map_W[healpix_pixnum]+=power_W
               #print tile_map[healpix_pixnum]
               #keep track of how many data points you have added to this pixel
               AUT_tile_map_data_entries_counter[healpix_pixnum]+=1
               #print tile_map_data_entries_counter[healpix_pixnum]
               #add the dipole


   #save the dictionary so we can plot the data later
   np.save(sat_dictionaries_dict_filename,sat_dictionaries_dict)
   #print 'useable sats are:'
   #print useable_sat_list
   
   #divide each pixel value by the number of data entries to form the average
   #print min(tile_map)
   #print max(tile_map_data_entries_counter)
   #print AUT_tile_map_W
   AUT_tile_map_W_av=np.divide(AUT_tile_map_W,AUT_tile_map_data_entries_counter)
   
   #make the dipole model map
   for pixel_number in range(0,len(dipole_model_map_W)):
      dipole_model_theta_rad, dipole_model_spherical_az_rad = hp.pix2ang(nside,pixel_number)
      dipole_model_az_rad=np.pi-dipole_model_spherical_az_rad
      dipole_model_zenith_angle_rad=dipole_model_theta_rad
      
      ##Don't forget the groundscreen. From mwapy/mwa_tile.py:
      ##    def groundScreen(self,za,freq):
      ##  """
      ##  Calculate the groundscreen effect for an ideal infinite groundscreen
      ##  given the dipole's height above the screen and the frequency (Hz)
      ##  """
      ##  l = vel_light/freq
      ## return numpy.sin(numpy.pi * (2.0*self.height/l) * numpy.cos(za))*2.0
      
      wavelength=vel_light/freq
      ground_screen=np.sin(np.pi * (2.0*height/wavelength) * np.cos(dipole_model_zenith_angle_rad))*2.0
      
      #print ground_screen
      
      #YY
      theta_YY=np.arccos(np.sin(dipole_model_zenith_angle_rad)*np.cos(dipole_model_az_rad))
      voltage_YY=np.sin(theta_YY)*ground_screen
      power_YY=voltage_YY**2
      
      #XX
      theta_XX=np.arccos(np.sin(dipole_model_zenith_angle_rad)*np.sin(dipole_model_az_rad))
      voltage_XX=np.sin(theta_XX)*ground_screen
      power_XX=voltage_XX**2
   
      if (polarisation=='YY'):
         dipole_model_map_W[pixel_number]=power_YY
      else:
         dipole_model_map_W[pixel_number]=power_XX
         
   #convert to db
   dipole_model_map_db=10.0*np.log10(dipole_model_map_W)
   
   #convert AUT map to db
   AUT_tile_map_dB_av=10.0*np.log10(AUT_tile_map_W_av)

   #Multiply by the dipole model (add in log space)
   AUT_tile_map_dB_av_dipole_corrected=AUT_tile_map_dB_av+dipole_model_map_db
   
   
   #hp.gnomview(tile_map_av, coord='C',reso=60,xsize=400, title=map_title, rot=(0,90,0))  
   hp.write_map(fits_name, AUT_tile_map_dB_av_dipole_corrected, overwrite=True)
   hp.write_map(counter_fits_name,AUT_tile_map_data_entries_counter, overwrite=True)
   hp.write_map(dipole_model_name,dipole_model_map_db, overwrite=True)
   hp.write_map(fits_name_no_dipole_correction,AUT_tile_map_dB_av, overwrite=True)
     
#def plot_ref_corrected_pb_map():
#
#   corrected_tile_map = hp.read_map(corrected_fits_name)
#   map_title="Reference Corrected Tile Map for %s" % AUT_tile_name
#   hp.gnomview(corrected_tile_map, coord='C',reso=60,xsize=400, title=map_title, rot=(0,90,0))  
#   plt.savefig(corrected_fig_name,dpi=200)
   
       
#def reference_correct_map():
#   #uses rf0XX as the ref antenna, divides (subtract in log space) the ref from the AUT 
#   ref_map=hp.read_map(rf0XX_reference_map_fits_name)
#   AUT_map=hp.read_map(fits_name)
#   corrected_pb_map=AUT_map-ref_map
#   hp.write_map(corrected_fits_name, corrected_pb_map)

#def plot_1D_by_pass():
#   read_dictionary = np.load(sat_dictionaries_dict_filename).item()
#   for sat in sat_list:
#      AUT_time=read_dictionary[sat]['AUT_time']
#      predicted_alt=read_dictionary[sat]['alt']
#      #go through the alt array until it gets above 0
#      pass_counter=0
#      pass_dict={}
#      for alt_index,alt_value in enumerate(predicted_alt):
#         if (alt_index>0 and alt_value > 0 and predicted_alt[alt_index-1]<0):
#            start_index=alt_index       
#            pass_counter+=1
#            pass_name='pass_%s' % pass_counter
#            #pass_dict.update('pass_name':)
#        if (alt_index>0 and alt_value < 0 and predicted_alt[alt_index-1]>0):
#            end_index=alt_index
            
            
            
def plot_1D(AUT_tile_name_in,ref_tile_name_in):

   AUT_tile_name=AUT_tile_name_in
   ref_tile_name=ref_tile_name_in

   fig_name_alltime="Sat_passes_alltime_AUT_%s_ref_%s.png" % (AUT_tile_name,ref_tile_name)
   fig_name_subset_time="Sat_passes_subset_time_AUT_%s_ref_%s.png" % (AUT_tile_name,ref_tile_name)
   sat_dictionaries_dict_filename='sat_dictionaries_dict_AUT_%s_ref_%s.npy' % (AUT_tile_name,ref_tile_name)
   
   read_dictionary = np.load(sat_dictionaries_dict_filename).item()
   fig1=plt.figure(1)
   fig2=plt.figure(2)
   ax1=fig1.add_subplot(411)
   ax2=fig1.add_subplot(412,sharex=ax1)
   ax3=fig1.add_subplot(413,sharex=ax1)
   ax4=fig1.add_subplot(414,sharex=ax1)
   
   ax21=fig2.add_subplot(511)
   ax22=fig2.add_subplot(512,sharex=ax21)
   ax23=fig2.add_subplot(513,sharex=ax21)
   ax24=fig2.add_subplot(514,sharex=ax21)
   ax25=fig2.add_subplot(515,sharex=ax21)
   
   total_passes=0
   
   for sat in sat_list:
      AUT_time=read_dictionary[sat]['AUT_time']
      AUT_time_UTC=[datetime.datetime.utcfromtimestamp(unix_time) for unix_time in AUT_time]
      AUT_time_mins=AUT_time/60.0
      predicted_alt=read_dictionary[sat]['alt']
      AUT_power=read_dictionary[sat]['AUT_power']
      ref_power=read_dictionary[sat]['ref_power']
      AUT_chan=read_dictionary[sat]['chan']
      power_diff=np.array(AUT_power)-np.array(ref_power)
      print "max power diff is %s " % max(power_diff)
      
      
      time_index_plot=np.arange(0,len(predicted_alt))
      
      #only show passes higher than alt threshold
      low_values_indices = predicted_alt < alt_threshold
      predicted_alt[low_values_indices]=np.nan
      
      #Find the number of local maxima (i.e. the number of passes the sat makes)
      pass_max_index = (np.diff(np.sign(np.diff(predicted_alt))) < 0).nonzero()[0] + 1 # local max
      number_of_passes_this_sat=len(pass_max_index)
      total_passes+=number_of_passes_this_sat
      
      x_axis_hour_interval=int(math.ceil((obs_end_time-obs_start_time)/60.0/60.0/10.0))
      fig1.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%y %H:%M'))
      fig1.gca().xaxis.set_major_locator(mdates.HourLocator(interval=3))
      ax1.plot(AUT_time_UTC,predicted_alt,label=sat)
      #ax1.plot(time_index_plot,predicted_alt,label=sat)
      ax1.set_ylabel('Predicted \n altitude (deg)')

      
      ax2.plot(AUT_time_UTC,AUT_power,label=sat)
      #ax2.plot(time_index_plot,AUT_power,label=sat)
      ax2.set_ylabel('AUT Power \n (dBm)')
      
      ax3.plot(AUT_time_UTC,ref_power,label=sat)   
      #ax3.plot(time_index_plot,ref_power,label=sat)
      ax3.set_ylabel('Ref Power \n (dBm)')
      
      ax4.plot(AUT_time_UTC,power_diff,label=sat)
      #ax4.plot(time_index_plot,power_diff,label=sat)
      ax4.set_ylabel('Power diff \n (dB)')
      ax4.set_xlabel('Time UTC ')

      
      #ax1.legend(loc='upper right')
      ax1.legend(bbox_to_anchor=(1.4,1), loc='upper right', ncol=1)
      
      fig2.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%y %H:%M'))
      fig2.gca().xaxis.set_major_locator(mdates.HourLocator(interval=1))
      ax21.plot(AUT_time_mins[plot_start_index:plot_end_index],predicted_alt[plot_start_index:plot_end_index],label=sat)
      #ax21.plot(time_index_plot[plot_start_index:plot_end_index],predicted_alt[plot_start_index:plot_end_index],label=sat)
      ax21.set_ylabel('Predicted \n altitude (deg)')
      
      ax22.plot(AUT_time_mins[plot_start_index:plot_end_index],AUT_power[plot_start_index:plot_end_index],label=sat)
      #ax22.plot(time_index_plot[plot_start_index:plot_end_index],AUT_power[plot_start_index:plot_end_index],label=sat)
      ax22.set_ylabel('AUT Power \n (dBm)')
      
      ax23.plot(AUT_time_mins[plot_start_index:plot_end_index],ref_power[plot_start_index:plot_end_index],label=sat)
      #ax23.plot(time_index_plot[plot_start_index:plot_end_index],ref_power[plot_start_index:plot_end_index],label=sat)
      ax23.set_ylabel('Ref Power \n (dBm)')
      
      ax24.plot(AUT_time_mins[plot_start_index:plot_end_index],power_diff[plot_start_index:plot_end_index],label=sat)
      #ax24.plot(time_index_plot[plot_start_index:plot_end_index],power_diff[plot_start_index:plot_end_index],label=sat)
      ax24.set_ylabel('Power diff \n (dB)')
      
      ax25.plot(AUT_time_mins[plot_start_index:plot_end_index],AUT_chan[plot_start_index:plot_end_index],label=sat)
      #ax25.plot(time_index_plot[plot_start_index:plot_end_index],AUT_chan[plot_start_index:plot_end_index],label=sat)
      ax25.set_ylabel('chan no')
      ax25.set_xlabel('Time UTC') 
      #ax25.legend().draggable()
      
      #ax21.legend(loc='upper right')
      
      ax21.legend(bbox_to_anchor=(1.3,1), loc='upper right', ncol=1)
      
      #plt.figure(2)
      #f, axarr = plt.subplots(4, sharex=True)
      #axarr[0].plot(AUT_time_mins[plot_start_index:plot_end_index],predicted_alt[plot_start_index:plot_end_index])
      #axarr[0].set_ylabel('Predicted altitude (deg)')
      #axarr[0].set_title('Sat %s' % sat)
      #axarr[1].scatter(AUT_time_mins[plot_start_index:plot_end_index],AUT_power[plot_start_index:plot_end_index])
      #axarr[1].set_ylabel('AUT Power (dBm)')
      #axarr[2].scatter(AUT_time_mins[plot_start_index:plot_end_index],ref_power[plot_start_index:plot_end_index])
      #axarr[2].set_ylabel('Ref Power (dBm)')
      #axarr[3].scatter(AUT_time_mins[plot_start_index:plot_end_index],power_diff[plot_start_index:plot_end_index])
      #axarr[3].set_ylabel('Power diff (dB)')
      #axarr[3].set_xlabel('Time (minutes) ')
      
      
      #plt.figure(1)
      #plt.subplot(311)
      #plt.plot(AUT_time,predicted_alt)
      #plt.ylabel('sat alt for %s (deg)' % sat)
      
      #plt.subplot(212)
      #plt.plot(AUT_time,AUT_power)
      #plt.ylabel('AUT power for %s (dB)' % sat)
      #plt.xlabel('time (s)')
      
   #plt.legend(loc='upper right');
   #plt.show()
   #plt.legend(loc='upper right');
   
   fig1.autofmt_xdate()   
   ax1.set_title('Orbcomm Satellite Passes. Total Predicted:%s' % (total_passes))
   
   fig2.autofmt_xdate()
   ax21.set_title('Orbcomm Satellite Passes Subset')
   
   fig1.savefig(fig_name_alltime,dpi=200,bbox_inches='tight',) 
   fig2.savefig(fig_name_subset_time,dpi=200,bbox_inches='tight') 
   
   fig1.clf()
   fig2.clf()

def plot_1D_simulated(obs_start_time,obs_end_time,t_step,observatory,polarisation):

   sat_dictionaries_dict_filename="sat_dictionaries_dict_from_%s_to_%s_at_%s_%s.npy" % (obs_start_time,obs_end_time,observatory,polarisation)
   
   fig_name_alltime="Sat_passes_from_%s_to_%s_at_%s_%s.png" % (obs_start_time,obs_end_time,observatory,polarisation)

   read_dictionary = np.load(sat_dictionaries_dict_filename).item()
   fig1=plt.figure(1)

   ax1=fig1.add_subplot(411)
   ax2=fig1.add_subplot(412,sharex=ax1)
   ax3=fig1.add_subplot(413,sharex=ax1)
   ax4=fig1.add_subplot(414,sharex=ax1)
   
   total_passes=0
   
   for sat in sat_list:
      AUT_time=read_dictionary[sat]['AUT_time']
      AUT_time_UTC=[datetime.datetime.utcfromtimestamp(unix_time) for unix_time in AUT_time]
      AUT_time_mins=AUT_time/60.0
      predicted_alt=read_dictionary[sat]['alt']
      AUT_power=read_dictionary[sat]['AUT_power']
      ref_power=read_dictionary[sat]['ref_power']
      AUT_chan=read_dictionary[sat]['chan']
      power_diff=np.array(AUT_power)-np.array(ref_power)
      print "max power diff is %s " % max(power_diff)
           
      time_index_plot=np.arange(0,len(predicted_alt))
      
      #Find the number of local maxima (i.e. the number of passes the sat makes)
      pass_max_index = (np.diff(np.sign(np.diff(AUT_power))) < 0).nonzero()[0] + 1 # local max
      number_of_passes_this_sat=len(pass_max_index)
      total_passes+=number_of_passes_this_sat
      
      x_axis_hour_interval=int(math.ceil((obs_end_time-obs_start_time)/60.0/60.0/10.0))
      plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%y %H:%M'))
      plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=x_axis_hour_interval))
      ax1.plot(AUT_time_UTC,predicted_alt,label=sat)
      #ax1.plot(time_index_plot,predicted_alt,label=sat)
      #ax1.plot(AUT_time_mins[pass_max_index], AUT_power[pass_max_index], "o")
      ax1.set_ylabel('Predicted \n altitude (deg)')
      plt.gcf().autofmt_xdate()

      ax2.plot(AUT_time_UTC,AUT_power,label=sat)
      #ax2.plot(time_index_plot,AUT_power,label=sat)
      ax2.set_ylabel('AUT Power \n (dBm)')

      
      ax3.plot(AUT_time_UTC,ref_power,label=sat)   
      #ax3.plot(time_index_plot,ref_power,label=sat)
      ax3.set_ylabel('Ref Power \n (dBm)')
      
      ax4.plot(AUT_time_UTC,power_diff,label=sat)
      #ax4.plot(time_index_plot,power_diff,label=sat)
      ax4.set_ylabel('Power diff \n (dB)')
      ax4.set_xlabel('Time UTC ')
      
      #ax1.legend(loc='upper right')
      ax1.legend(bbox_to_anchor=(1.4,1), loc='upper right', ncol=1)
      
   ax1.set_title('Orbcomm Satellite Passes. Total:%s' % (total_passes))
   if os.path.exists(fig_name_alltime):
      cmd="rm -f %s" % fig_name_alltime
      os.system(cmd)
   fig1.savefig(fig_name_alltime,dpi=200,bbox_inches='tight',)    
   fig1.clf()
   
#give start and stop times in UT '%Y/%m/%d %H:%M:%S'
def plot_sat_prediction(sat_list,start_time_UT,end_time_UT):

   #initialise sateph object
   sateph = Sateph()
   
   
   #UT_time=datetime.datetime.utcfromtimestamp(int("1284101485")).strftime('%Y/%m/%d %H:%M:%S')
   start_time_UT=datetime.datetime.strptime(start_time_UT, '%Y/%m/%d %H:%M:%S')
   start_time_unix=time.mktime(start_time_UT.timetuple())
   end_time_UT=datetime.datetime.strptime(end_time_UT, '%Y/%m/%d %H:%M:%S')
   end_time_unix=time.mktime(end_time_UT.timetuple())
   
   print 'start_time unix %s' % start_time_unix
   print 'end_time unix %s' % end_time_unix
   
   
   predicted_time_array=np.arange(start_time_unix,end_time_unix)
   
   for sat_desig in sat_list:
      alt=sateph.get_sat_alt_az(sat_desig,predicted_time_array)[1]
      print "alt for sat %s is %s deg" % (sat_desig,alt)
      
def plot_chan_histogram(AUT_tile_name,ref_tile_name):
   sat_dictionaries_dict_filename='sat_dictionaries_dict_AUT_%s_ref_%s.npy' % (AUT_tile_name,ref_tile_name)
   read_dictionary = np.load(sat_dictionaries_dict_filename).item()
   for sat in sat_list:
         chans=np.array(read_dictionary[sat]['chan'])
         chans=chans[~np.isnan(chans)]
         if (np.isfinite(chans).any()):
            plt.hist(chans, bins=range(int(min(chans)), int(max(chans)) + 1, 1))
            fig_name="histogram_sat_%s.png" % sat
            plt.title("Histogram of max power chans for sat %s" % sat)
            plt.xlabel("Chan Number")
            plt.ylabel("Frequency")
            #current_fig = plt.gcf()
            #current_fig.savefig(fig_name,dpi=200)
            #plt.clf()
            plt.show()

def compare_pb_maps():
    
   #get the noise level from the r0YY - rf1YY map first
   ref_AUT_tile_name="rf0YY"
   ref_tile_name="rf1YY"
   pol='YY'
   
   ref_fits_name_no_dipole_correction="tile_map_AUT_%s_ref_%s_no_dipole_correction.fits" % (ref_AUT_tile_name,ref_tile_name)
   ref_tile_map_av_no_dipole_correction = hp.read_map(ref_fits_name_no_dipole_correction)
   #print AUT_tile_map_av_no_dipole_correction
   
   #compute the rms of the un-dipole-corrected null test
   dipole_rms=nanrms(ref_tile_map_av_no_dipole_correction)
   print "RMS of null test %s / %s: %s " % (ref_AUT_tile_name,ref_tile_name, dipole_rms)
   
   ###For loop here
   #now go through for each AU, normalise it, subtract the model and look at what's left, comaring it to the rms noise of 4.5 dB
   AUT_tile_name="051YY"
   AUT_fits_name="tile_map_AUT_%s_ref_%s.fits" % (AUT_tile_name,ref_tile_name)
   AUT_map_title="Projected normalised AUT %s ref %s" % (AUT_tile_name,ref_tile_name)
   AUT_map=hp.read_map(AUT_fits_name)
   AUT_map_W=10.0**(AUT_map/10.0)
   max_AUT_map_W=np.nanmax(AUT_map_W)
   normalised_AUT_map_W=AUT_map_W/max_AUT_map_W
   projected_normalised_AUT_map_W=hp.orthview(map=normalised_AUT_map_W,coord='E',half_sky=True,xsize=400,title=AUT_map_title,rot=(0,90,0),return_projected_map=True)
   #plt.clf()

   tile_model_name="healpix_beam_map_%s.fits" % pol
   tile_model_map_title="Tile model pol %s" % pol
   
   tile_model_map_W=hp.read_map(tile_model_name)
   #Don't convert to dB
   ##tile_model_dB_map=10.0*np.log10(tile_model_map)
   max_tile_model_map_W=np.nanmax(tile_model_map_W)
   normalised_tile_model_map_W=tile_model_map_W/max_tile_model_map_W
   projected_tile_model_map_W=hp.orthview(map=normalised_tile_model_map_W,coord='C',half_sky=True,xsize=400,title=tile_model_map_title,rot=(0,0,301.401202),return_projected_map=True)
   #plt.show()
   plt.clf()
   
   model_fractional_difference_map_W=(projected_normalised_AUT_map_W-projected_tile_model_map_W)/projected_tile_model_map_W
   model_compare_title="Fractional difference between model and Tile_%s ref_%s pol:%s" % (AUT_tile_name,ref_tile_name,pol)
   plt.imshow(model_fractional_difference_map_W,origin='lower')
   plt.title(model_compare_title)
   plt.colorbar()
   figmap = plt.gcf()
   figmap.savefig("Model_compare_Tile_%s_ref_%s_%s.png" % (AUT_tile_name,ref_tile_name,pol),dpi=200)
   plt.clf()
   

   #In dB
   model_fractional_difference_map_dB=10.0*np.log10(1.0+(model_fractional_difference_map_W))
   model_compare_title_dB="Model difference Tile_%s ref_%s pol:%s dB" % (AUT_tile_name,ref_tile_name,pol)
   plt.imshow(model_fractional_difference_map_dB,origin='lower')
   plt.title(model_compare_title_dB)
   plt.colorbar()
   figmap = plt.gcf()
   figmap.savefig("Model_compare_Tile_%s_ref_%s_%s_dB.png" % (AUT_tile_name,ref_tile_name,pol),dpi=200)
   plt.clf()
   
   
   #Since the tile model is obviously wrong near the nulls, due to beamformer delay errors etc (see other Neben paper on this)
   #just look at values away from the nulls e.g. mask the projected array for the beam model and then use that
   masked_normalised_tile_model_map_W=projected_tile_model_map_W
   masked_normalised_tile_model_map_W[masked_normalised_tile_model_map_W<0.003]=np.nan
   plt.imshow(masked_normalised_tile_model_map_W,origin='lower')
   plt.colorbar()
   figmap = plt.gcf()
   figmap.savefig("Masked_tile_model_%s.png" % (pol),dpi=200)
   plt.clf()
   
   masked_model_compare_map=(projected_normalised_AUT_map_W-masked_normalised_tile_model_map_W)/masked_normalised_tile_model_map_W
   masked_model_compare_title="Masked Model compare Tile_%s ref_%s pol:%s" % (AUT_tile_name,ref_tile_name,pol)
   plt.imshow(masked_model_compare_map,origin='lower')
   plt.title(masked_model_compare_title)
   plt.colorbar()
   figmap = plt.gcf()
   figmap.savefig("Masked_Model_compare_Tile_%s_ref_%s_%s.png" % (AUT_tile_name,ref_tile_name,pol),dpi=200)
   plt.clf()
   
   #In dB
   masked_model_fractional_difference_map_dB=10.0*np.log10(1.0+(masked_model_compare_map))
   masked_model_compare_title_dB="Masked model difference Tile_%s ref_%s pol:%s dB" % (AUT_tile_name,ref_tile_name,pol)
   plt.imshow(masked_model_fractional_difference_map_dB,origin='lower')
   plt.title(masked_model_compare_title_dB)
   plt.colorbar()
   figmap = plt.gcf()
   figmap.savefig("Masked_model_compare_Tile_%s_ref_%s_%s_dB.png" % (AUT_tile_name,ref_tile_name,pol),dpi=200)
   plt.clf()
   
   #ref_tile_name=ref_tile_name_in
   
   
   #AUT_tile_name=AUT_tile_name_in
   #ref_tile_name=ref_tile_name_in
   #if ('XX' in ref_tile_name):
   #   polarisation='XX'
   #else:
   #   polarisation='YY'
   #fits_name="tile_map_AUT_%s_ref_%s.fits" % (AUT_tile_name,ref_tile_name)
   #fits_name_no_dipole_correction="tile_map_AUT_%s_ref_%s_no_dipole_correction.fits" % (AUT_tile_name,ref_tile_name)
   #fig_name="Tile_%s_map_corrected_%s.png" % (AUT_tile_name,ref_tile_name)
   #fig_name_no_dipole_correction="Tile_%s_map_%s_no_dipole_correction.png" % (AUT_tile_name,ref_tile_name)
   #AUT_tile_map_av = hp.read_map(fits_name)
   #AUT_tile_map_av_no_dipole_correction = hp.read_map(fits_name_no_dipole_correction)
   #map_title="AUT: Tile%s Ref:%s" % (AUT_tile_name,ref_tile_name)
   #map_title_no_dipole_correction="No dipole Correction AUT: Tile%s Ref:%s" % (AUT_tile_name,ref_tile_name)
   ##hp.gnomview(AUT_tile_map_av, coord='C',reso=60,xsize=400, title=map_title, rot=(0,90,0))  
   #hp.orthview(map=AUT_tile_map_av,coord='E',half_sky=True,xsize=400,title=map_title,rot=(0,90,0))
   #figmap = plt.gcf()
   #figmap.savefig(fig_name,dpi=200)
   #plt.clf()
   
   #no dipole correction
   #hp.orthview(map=AUT_tile_map_av_no_dipole_correction,coord='C',half_sky=True,xsize=400,title=map_title_no_dipole_correction,rot=(0,90,0))
   #figmap = plt.gcf()
   #figmap.savefig(fig_name_no_dipole_correction,dpi=200)
   #plt.clf()
      

#AUT_tile_name='rf0XX'
#AUT_tile_name='052YY'
#ref_tile_name='rf1XX'
#AUT_tile_name='rf1XX'
#ref_tile_name='rf0YY'

#######beam reprojecting stuff:#########
beam_filename_XX='beam_model/1162046832_cotter_zenith_hacked-image_beamXX.fits'
beam_filename_YY='beam_model/1162046832_cotter_zenith_hacked-image_beamYY.fits'
output_healpix_beam_filename_XX='healpix_beam_map_XX.fits'
output_healpix_beam_filename_YY='healpix_beam_map_YY.fits'

#reproject_beam_to_healpix(beam_filename_YY,output_healpix_beam_filename_YY)

###########data stuff:####################
##Jarryd 2015 Nov data:
#obs_start_time=1448429025
#obs_end_time=1448447953
##Jack test data 4 Aug 2017 (start 1501833532?)

import sys,os
from optparse import OptionParser,OptionGroup

usage = 'Usage: plot_beams.py [options]'

parser = OptionParser(usage=usage)

parser.add_option('--generate_pb_map',action='store_true',dest='generate_pb_map',default=False,help='Generate the primary beam map from the RFExplorerdata and save map data [default=%default]')
parser.add_option('--plot_pb_map',action='store_true',dest='plot_pb_map',default=False,help='Make and save the plots of the beam maps [default=%default]')
parser.add_option('--plot_1D',action='store_true',dest='plot_1D',default=False,help='Make and save 1D plots of alt/power vs time [default=%default]')
parser.add_option('--obs_start_time',type='string', dest='obs_start_time',default=None,help='Start time in gps seconds of the observations e.g. --obs_start_time="1501833532" [default=%default]')
parser.add_option('--obs_end_time',type='string', dest='obs_end_time',default=None,help='Start time in gps seconds of the observations e.g. --obs_end_time="1501833542" [default=%default]')
parser.add_option('--convert_to_frows',action='store_true',dest='convert_to_frows',default=False,help='Convert data in the raw_data_dir to required f-row format [default=%default]')
parser.add_option('--raw_data_dir',type='string', dest='raw_data_dir',default='./mwa_beam_measurement/test_data',help='Directory where the raw data from the RFExplorers is stored e.g. --raw_data_dir="/data/code/git/mwa_beam_measurement/test_data" [default=%default]')
parser.add_option('--converted_data_dir',type='string', dest='converted_data_dir',default='./Converted',help='Directory where the data converted to f-row format is stored e.g. --converted_data_dir="./Converted" [default=%default]')
parser.add_option('--ref_tile_list',type='string', dest='ref_tile_list',default=None,help='list of reference antenna names e.g. --ref_ant_tile_list="rf0XX,rf1XX,rf0YY,rf1YY" [default=%default]')
parser.add_option('--AUT_tile_list',type='string', dest='AUT_tile_list',default=None,help='list of Antenna Under Test (AUT) names e.g. --AUT_tile_list="050XX,051XX,052YY,053YY" [default=%default]')
parser.add_option('--ref_signal_threshold',type='string', dest='ref_signal_threshold',default="-85",help='Threshold reference antenna power level for including data in plots in dBm e.g. --ref_signal_threshold="-85" [default=%default]')
parser.add_option('--AUT_signal_threshold',type='string', dest='AUT_signal_threshold',default="-50",help='Threshold AUT  power level for including data in plots in dBm e.g. --ref_signal_threshold="-50" [default=%default]')
parser.add_option('--alt_threshold',type='string', dest='alt_threshold',default="30",help='Threshold altitiude for how high a sat must be for inclusion e.g. --alt_threshold="30" [default=%default]')
parser.add_option('--converted_data_filename_base',type='string', dest='converted_data_filename_base',default="converted_all",help=' e.g. --converted_data_filename_base="Aug_2017_wed_all_good" [default=%default]')
parser.add_option('--use_old_TLEs',action='store_true',dest='use_old_TLEs',default=False,help='**Not yet implemented** Do not check for more recent TLEs [default=%default]')
parser.add_option('--plot_chan_histogram',action='store_true',dest='plot_chan_histogram',default=False,help='Plot the histogram of which channel has been allocated to each satellite [default=%default]')



(options, args) = parser.parse_args()

   
   
#obs_start_time=1501833532
#obs_end_time=1501833542

if options.convert_to_frows:
   convert_data_to_frows=True
   converted_data_filename_base=options.converted_data_filename_base
else:
   convert_data_to_frows=False

raw_data_dir=options.raw_data_dir
converted_frow_data_dir=options.converted_data_dir

obs_start_time=float(options.obs_start_time)
obs_end_time=float(options.obs_end_time)

ref_tile_list=options.ref_tile_list.split(',')
AUT_tile_list=options.AUT_tile_list.split(',')
ref_signal_threshold=float(options.ref_signal_threshold)
AUT_signal_threshold=float(options.AUT_signal_threshold)

alt_threshold=float(options.alt_threshold)

for ref_ant in ref_tile_list:
   for AUT in AUT_tile_list:
      if options.generate_pb_map:
         generate_pb_map(AUT,ref_ant,AUT_signal_threshold,ref_signal_threshold)
      if options.plot_pb_map:
         plot_pb_map(AUT,ref_ant)
      if (options.plot_1D):
         plot_1D(AUT,ref_ant)
      if (options.plot_chan_histogram):
         plot_chan_histogram(AUT,ref_ant)

##obsrange:(1448335711, 1448851982)
##obs_range_test =  (1448429025 , 1448447953)


#Compare the maps to the model:
#compare_pb_maps()

################simulatey stuff###############
#observatory='Greenbank'
#observatory='MWA'
#polarisation='YY'
#t_step=0.5
#noise = 1 #noise added to AUT and ref in dBm (due to amplitude stability of 3dBm for RfExplorers)

#simulate_pb_map(obs_start_time,obs_end_time,t_step,observatory,polarisation,noise)
#plot_simulated_pb_map(obs_start_time,obs_end_time,t_step,observatory,polarisation)
#plot_1D_simulated(obs_start_time,obs_end_time,t_step,observatory,polarisation)

#plot_sat_prediction(sat_list,'2016/09/11 00:00:00','2016/09/14 23:59:59')
#reference_correct_map()
#plot_ref_corrected_pb_map()




