import os, re, subprocess
import numpy as np
import matplotlib.pyplot as plt
import exoplanet as xo
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider, Button

from helper_functions import find_nearest
from helper_functions import determine_cadence

import math
import lightkurve as lk
from astropy.io import fits

def tic_id_from_simbad(other_id):
	#takes other_id (string) and queries Simbad to obtain the TIC ID
	
	from astroquery.simbad import Simbad
	import astropy
	ID_table = Simbad.query_objectids(other_id)
	
	if type(ID_table) is not astropy.table.table.Table:
		return(None)
	
	ID_table['ID'] = ID_table['ID'].astype(str)
	

	ID_pandas = ID_table.to_pandas()
	tic_id = ID_pandas[ID_pandas['ID'].str.contains("TIC")]
	
	
	return tic_id['ID'].values[0]


def tic_id_from_exoplanet_archive(other_id):
	# if SIMBAD can't get TIC ID, looks for it in exoplanet archive
	# most of this is from transit_info_from_exoplanet_archive tbh; hesitant to try to merge functions
	
	import pandas as pd


	# SIMBAD should be able to grab data for already confirmed systems
	# this is primarily for TOI candidates
	
	a = "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=toi&select=toipfx,tid,pl_tranmid,pl_orbper,pl_trandurh&format=csv"

	exoplanets = pd.read_csv(a)

	#rename columns
	column_dict = {
	'toipfx':'toi_host',    
	'tid':'tic_id',
	'pl_tranmid':'t0 [BJD]',
	'pl_orbper':'period [days]',
	'pl_trandurh':'duration [hours]',
	}

	exoplanets.rename(columns=column_dict, inplace=True)
	
	exoplanets['toi_host'] = 'TOI-' + exoplanets['toi_host'].astype(str)
	exoplanets['tic_id'] = 'TIC ' + exoplanets['tic_id'].astype(str)
	
	# replacing any chars that don't match with ones that do; standardizes user input
	
	other_id = other_id.lstrip('toi- TOI')
	other_id = 'TOI-' + other_id # there's probably a smarter way to do this
	
	
	tic_id = exoplanets['tic_id'].loc[exoplanets['toi_host'] == other_id]
	
	return tic_id.values[0]



def transit_info_from_exoplanet_archive(tic_id):
	#takes TIC ID and queries exoplanet archive to return t0, period, and duration
	
	import pyvo as vo
	import pandas as pd


	service = vo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP")
	
	a = service.search("SELECT \
					   tic_id, pl_tranmid, pl_orbper, pl_trandur\
					   FROM pscomppars")
	
   
	
	exoplanets = a.to_table()
	
	exoplanets = exoplanets.to_pandas()
	
	
	#rename columns
	column_dict = {
	'pl_tranmid':'t0 [BJD]',
	'pl_orbper':'period [days]',
	'pl_trandur':'duration [hours]',
	}
	
	exoplanets.rename(columns=column_dict, inplace=True)
	
	result = exoplanets[exoplanets['tic_id'] == tic_id]

	#if there's no row in the planetary comparison table, check TOI table
	if result.empty:
		print("Exoplanet Archive: TOI Table")
		print("----------------------------")
		a = "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=toi&select=tid,pl_tranmid,pl_orbper,pl_trandurh&format=csv"

		exoplanets = pd.read_csv(a)
		
		#rename columns
		column_dict = {
		'tid':'tic_id',
		'pl_tranmid':'t0 [BJD]',
		'pl_orbper':'period [days]',
		'pl_trandurh':'duration [hours]',
		}
		
		exoplanets.rename(columns=column_dict, inplace=True)
		exoplanets['tic_id'] = 'TIC ' + exoplanets['tic_id'].astype(str)
		
		result = exoplanets[exoplanets['tic_id'] == tic_id]
		
		
	else:
		print("Exoplanet Archive: Planet Comparison Table")
		print("------------------------------------------")
		
	

	return result
	
	




def get_transit_info(planet_id):
	#takes a id, queries Simbad to get matching TIC ID
	#then queries exoplanet archive to extract t0, period, and duration
	#if no Simbad match found, then returns None and prints error message
	#if no exoplanet archive match found, then returns None and prints error message

	tic_id = tic_id_from_simbad(planet_id)
	
	if tic_id == None:
		# search on exoplanet archive
		print("No TIC ID match found on Simbad")
		print("Checking on exoplanet archive")
		
		tic_id = tic_id_from_exoplanet_archive(planet_id)
		
		if tic_id == None:
			print("No TIC ID match found on exoplanet archive either.")
			return None
		
		else:
			print('Found ' + tic_id)
			transit_info = transit_info_from_exoplanet_archive(tic_id)
			
			if transit_info.empty:
				print("No transit info found (???? Archive TIC ID found.)")
				return None
			
			else:
				return transit_info
	
	else: # Simbad does find TIC ID
		transit_info = transit_info_from_exoplanet_archive(tic_id)
		
		if transit_info.empty:
			print("No TIC ID match found on exoplanet archive (???? Simbad TIC ID found.)")
			return None
		
		else:
			return transit_info




			
def get_light_curve(planet_id, flux_type, TESS = False, Kepler = False, 
					user_periods = None, user_t0s = None, user_durations = None,
					planet_number = 1, mask_width = 1.3):
	

	transit_info = get_transit_info(planet_id)
	if type(transit_info) == None:
		return None
	
	print('planet #,[    tic_id      ,    t0 [BJD]    ,  P [days] , tdur [hrs]')
	transit_info_list = transit_info.astype(str).values.tolist()
	for ii in range(0, len(transit_info_list)):
		print("planet " + str(ii+1) + ", " + str(transit_info_list[ii]))

	
	tic_id = str(transit_info['tic_id'].values[0])
	
	
	if user_periods != None:
		periods = user_periods
		print("using periods = " + str(user_periods))
	else:
		periods = np.array(transit_info['period [days]'].values, dtype = float)

		
	if user_t0s != None:
		t0s = user_t0s
		print("using t0s = " + str(user_t0s))
	else:
		t0s = np.array(transit_info['t0 [BJD]'].values, dtype = float)

	if user_durations != None:
		durations = user_durations
		print("using durations = " + str(user_durations))
	else:
		durations = np.array(transit_info['duration [hours]'].values, dtype = float)
	
	
	print('')
	print('')
	
	

	nplanets = len(periods)
	
	if TESS:
		#switch to TESS BJD
		t0s = t0s - 2457000
		
		if flux_type == 'qlp':
			lc_files = lk.search_lightcurve(
				tic_id, mission='TESS', author = 'qlp'
			).download_all(quality_bitmask="hardest")
		
		else:
			#pull in short cadence TESS SPOC LC
			lc_files_short_cadence = lk.search_lightcurve(
				tic_id, mission='TESS', author = 'SPOC', cadence = 'short'
			).download_all(quality_bitmask="hardest", flux_column=flux_type)

			#pull in long cadence TESS SPOC LC
			lc_files_long_cadence = lk.search_lightcurve(
				tic_id, mission='TESS', author = 'SPOC', cadence = 'long'
			).download_all(quality_bitmask="hardest", flux_column=flux_type)


			#use short cadence TESS data if if exists, else use long cadence
			if lc_files_short_cadence == []:
				lc_files = lc_files_long_cadence
			else:
				lc_files = lc_files_short_cadence
		
	if Kepler:
		#switch to Kepler BJD
		t0s = t0s - 2454833
		
		#pull in Kepler LC
		lc_files = lk.search_lightcurve(
			tic_id, mission='Kepler'
		).download_all(quality_bitmask="hardest", flux_column=flux_type)
		
	
	quarters = []
	crowding = []
	flux_fraction = []
	for file in lc_files:
		quarters.append([np.min(file.time.value),
						np.max(file.time.value)])
		
		if flux_type != 'qlp':
			crowding.append(file.CROWDSAP)
			flux_fraction.append(file.FLFRCSAP)
		
		
	
		
	lc = lc_files.stitch().remove_nans()
	
	xs = lc.time.value
	ys = lc.flux
	ys_err = lc.flux_err
	

   

	mask = np.zeros(np.shape(xs), dtype=bool)
	for ii in range(0, nplanets):
		masks = lc.create_transit_mask(period=periods[ii], 
									   duration=durations[ii]/24.*mask_width, 
									   transit_time=t0s[ii])
		mask += masks
		
	mask_fitted_planet = lc.create_transit_mask(period=periods[planet_number-1], 
												duration=durations[planet_number-1]/24.*mask_width, 
												transit_time=t0s[planet_number-1])

	
	
	
	
	
	#save the period, duration, and t0 for the planet we are fitting for...
	period = np.array([periods[planet_number-1]])
	t0 = t0s[planet_number-1]
	duration = np.array([durations[planet_number-1]])
	
	
	print('using the following params for the planet we are fitting')
	print("--------------------------------------------------------")
	
	if TESS:
		print('[  t0 [TESS BJD]  , P [days], tdur [hrs]')
	if Kepler:
		print('[ t0 [Kepler BJD] , P [days], tdur [hrs]')
	print('[' + str(t0) + ', ' + str(period[0]) + ', ' + str(duration[0]) + ']')
	
	
	
	min_time = xs.min()
	max_time = xs.max()
	t0s_all = []
	while t0 > min_time:
		t0 -= period[0]

	while t0 < max_time:
		t0s_all.append(t0)
		t0 += period[0]
	
	
	t0s_in_data = []
	for t0 in t0s_all:
		nearest_lc_time = find_nearest(xs, t0)
		if np.abs(t0 - nearest_lc_time) < 0.1: # 2.5 hours ~ 0.1 days
			t0s_in_data.append(t0)
			
			
	mu = np.median( ys )
	ys = ( ys / mu - 1 )
	ys_err = ( ys_err / mu )
	
	

	
	

	
	return \
		np.array(xs), np.array(ys), np.array(ys_err), mask, mask_fitted_planet, \
		np.array(t0s_in_data), period, duration, quarters, crowding, flux_fraction


def get_light_curve_alt(planet_id, flux_type, save_to_directory, TESS = False, Kepler = False, 
					user_periods = None, user_t0s = None, user_durations = None,
					planet_number = 1, mask_width = 1.3):

	
	# searches MAST bulk download scripts in folder "curl_scripts" to see if there is a TIC match in any sector
	

	transit_info = get_transit_info(planet_id)
	if type(transit_info) == None:
		return None
	
	print('planet #,[    tic_id      ,    t0 [BJD]    ,  P [days] , tdur [hrs]')
	transit_info_list = transit_info.astype(str).values.tolist()
	for ii in range(0, len(transit_info_list)):
		print("planet " + str(ii+1) + ", " + str(transit_info_list[ii]))

	
	tic_id = str(transit_info['tic_id'].values[0])
	
	
	if user_periods != None:
		periods = user_periods
		print("using periods = " + str(user_periods))
	else:
		periods = np.array(transit_info['period [days]'].values, dtype = float)

		
	if user_t0s != None:
		t0s = user_t0s
		print("using t0s = " + str(user_t0s))
	else:
		t0s = np.array(transit_info['t0 [BJD]'].values, dtype = float)

	if user_durations != None:
		durations = user_durations
		print("using durations = " + str(user_durations))
	else:
		durations = np.array(transit_info['duration [hours]'].values, dtype = float)
	
	
	print('')
	print('')
	
	

	nplanets = len(periods)
	
	if TESS:
		#switch to TESS BJD
		t0s = t0s - 2457000
		
		if flux_type == 'qlp':
			lc_files = lk.search_lightcurve(
				tic_id, mission='TESS', author = 'qlp'
			).download_all(quality_bitmask="hardest")
		
		else:

			# assuming we want to look at sap if not qlp
			# search bulk download scripts for matching TIC and download data
			# assuming directory with scripts is in /detrend_package/curl_scripts/
			curl_directory = './curl_scripts/'

			# working with each sector script iteratively

			curl_files = [f for f in os.listdir(curl_directory) if re.search('tesscurl_sector_\d+_lc\.sh', f)]
			#print(curl_files)
			curl_commands = [] # list of curl commands to run after all scripts have been searched through

			for curl_file in curl_files:
				#print(curl_file)
				#print(tic_id.split(' ')[1])
				# search TIC match in curl_file
				# if match, add command to curl_commands
				curl_file = open(os.path.join(curl_directory, curl_file))

				lines = curl_file.read().splitlines()
				for line in lines:
					if re.search('tess\d+-s\d+-\d+' + tic_id.split(' ')[1], line):
						# if we find a match for TIC in the line, we run the curl script in that line
						curl_commands.append(line)
						#print(line)

				curl_file.close()

			# check if there are curl commands; if not, just don't

			if curl_commands == []:
				print('No SAP flux files available for ' + planet_id + ' in sectors 1-50.')

			else: # we have curl commands to run

				for curl_command in curl_commands:
					subprocess.run(args=curl_command, cwd=save_to_directory, shell=True)
					# that should download the fits files

				# so now let's look for .fits files here
				fits_files = [f for f in os.listdir(save_to_directory) if re.search('\.fits', f)]

				lc_files = [] # getting ready to stitch lc files

				for fits_file in fits_files:
					#print(fits_file)
					lc_file = lk.read(save_to_directory+fits_file)
					lc_files.append(lc_file)

				lc_files = lk.LightCurveCollection(lc_files)

				# we should have our nice lc files now!

		
	if Kepler:
		#switch to Kepler BJD
		t0s = t0s - 2454833
		
		#pull in Kepler LC
		lc_files = lk.search_lightcurve(
			tic_id, mission='Kepler'
		).download_all(quality_bitmask="hardest", flux_column=flux_type)
		

	quarters = []
	crowding = []
	flux_fraction = []
	for file in lc_files:
		quarters.append([np.min(file.time.value),
						np.max(file.time.value)])
		
		if flux_type != 'qlp':
			crowding.append(file.CROWDSAP)
			flux_fraction.append(file.FLFRCSAP)
		
		
	
	# wow i have no idea if this will work!
	
	lc = lc_files.stitch().remove_nans()
	
	xs = lc.time.value
	ys = lc.flux
	ys_err = lc.flux_err
	

   

	mask = np.zeros(np.shape(xs), dtype=bool)
	for ii in range(0, nplanets):
		masks = lc.create_transit_mask(period=periods[ii], 
									   duration=durations[ii]/24.*mask_width, 
									   transit_time=t0s[ii])
		mask += masks
		
	mask_fitted_planet = lc.create_transit_mask(period=periods[planet_number-1], 
												duration=durations[planet_number-1]/24.*mask_width, 
												transit_time=t0s[planet_number-1])

	
	
	
	
	
	#save the period, duration, and t0 for the planet we are fitting for...
	period = np.array([periods[planet_number-1]])
	t0 = t0s[planet_number-1]
	duration = np.array([durations[planet_number-1]])
	
	
	print('using the following params for the planet we are fitting')
	print("--------------------------------------------------------")
	
	if TESS:
		print('[  t0 [TESS BJD]  , P [days], tdur [hrs]')
	if Kepler:
		print('[ t0 [Kepler BJD] , P [days], tdur [hrs]')
	print('[' + str(t0) + ', ' + str(period[0]) + ', ' + str(duration[0]) + ']')
	
	
	
	min_time = xs.min()
	max_time = xs.max()
	t0s_all = []
	while t0 > min_time:
		t0 -= period[0]

	while t0 < max_time:
		t0s_all.append(t0)
		t0 += period[0]
	
	
	t0s_in_data = []
	for t0 in t0s_all:
		nearest_lc_time = find_nearest(xs, t0)
		if np.abs(t0 - nearest_lc_time) < 0.1: # 2.5 hours ~ 0.1 days
			t0s_in_data.append(t0)
			
			
	mu = np.median( ys )
	ys = ( ys / mu - 1 )
	ys_err = ( ys_err / mu )
	
	

	
	

	
	return \
		np.array(xs), np.array(ys), np.array(ys_err), mask, mask_fitted_planet, \
		np.array(t0s_in_data), period, duration, quarters, crowding, flux_fraction
	
	

