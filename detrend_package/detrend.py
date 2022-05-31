import numpy as np
import time

from cofi_AM import *
from poly_AM import *
from poly_local import * 
from gp import *
from plot import plot_individual_outliers
from manipulate_data import *
from outlier_rejection import *

def trim_jump_times(x, y, yerr, mask, mask_fitted_planet, t0s, period, jump_times):
    #x = time 
    #y = flux 
    #yerr = flux error
    #mask = mask
    #t0s = midtransits
    #period = planet period to define plotting limit
    #jump_times = jump times to trim data around (one before and one after each transit)
    
   
    if jump_times != []:
        x_epochs = []
        y_epochs = []
        yerr_epochs = []
        mask_epochs = []
        mask_fitted_planet_epochs = []
    
        for ii in range(0, len(jump_times)-1):
            jump_start = jump_times[ii]
            jump_end = jump_times[ii+1]
            

            if ii % 2 == 0:
                jump_start = find_nearest(x, jump_start)
                jump_end = find_nearest(x, jump_end)
                

                epoch_split = [jump_start, jump_end]
                start_index = int(np.where(x == epoch_split[0])[0])
                end_index = int(np.where(x == epoch_split[1])[0])

                x_epochs.append(x[start_index:end_index])
                y_epochs.append(y[start_index:end_index])
                yerr_epochs.append(yerr[start_index:end_index])
                mask_epochs.append(mask[start_index:end_index])
                mask_fitted_planet_epochs.append(mask_fitted_planet[start_index:end_index])
                

        x_epochs = np.array(x_epochs, dtype=object)
        y_epochs = np.array(y_epochs, dtype=object)
        yerr_epochs = np.array(yerr_epochs, dtype=object)
        mask_epochs = np.array(mask_epochs, dtype=object)
        mask_fitted_planet_epochs = np.array(mask_fitted_planet_epochs, dtype=object)

    else:
        x_epochs, y_epochs, yerr_epochs, mask_epochs, mask_fitted_planet_epochs = \
        split_around_transits(x, y, yerr, mask, mask_fitted_planet, t0s, 1./2., period)
    


    
    return x_epochs, y_epochs, yerr_epochs, mask_epochs, mask_fitted_planet_epochs

def get_detrended_lc(y, detrending_model):
    '''
    input:
    y = light curve
    detrending model = stellar detrending model evaluated at same time as fluxes
    
    returns:
    detrended_lc = detrended light curve evaluated at same time as fluxes
    
    '''
    detrended_lc = (((y + 1) / (detrending_model + 1)) - 1)
    
    return np.array(detrended_lc)

def remove_trim_times(x, y, yerr, trim_times):
    
    for trim in trim_times:
        trim_min = find_nearest(x, trim[0]) 
        trim_max = find_nearest(x, trim[1]) 
        
        index_min = int(np.where(x == trim_min)[0])
        index_max = int(np.where(x == trim_max)[0])
        
        indices = np.arange(index_min, index_max)

        
        x = np.delete(x, indices)
        y = np.delete(y, indices)
        yerr = np.delete(yerr, indices)
        
    return x, y, yerr



def add_nans_for_missing_data(sap_x_local, sap_detrended_lcs, sap_yerr_local, sap_mask_local, sap_mask_fitted_planet_local, 
                              pdc_x_local, pdc_detrended_lcs, pdc_yerr_local, pdc_mask_local, pdc_mask_fitted_planet_local):
    
    
    print('pdc length in: ', str(len(pdc_x_local)))
    print('sap length in: ', str(len(sap_x_local)))
    print('---')
    
    
    for ii in range(0, len(pdc_x_local)):
        time = pdc_x_local[ii]
        yerr = pdc_yerr_local[ii]
        mask = pdc_mask_local[ii]
        mask_fitted_planet = pdc_mask_fitted_planet_local[ii]
        if time not in sap_x_local:
            for kk in range(0, len(sap_x_local)):
                if sap_x_local[kk] > time:
                    sap_x_local = np.insert(sap_x_local, kk, time)
                    sap_yerr_local = np.insert(sap_yerr_local, kk, yerr)
                    sap_mask_local = np.insert(sap_mask_local, kk, mask)
                    sap_mask_fitted_planet_local = np.insert(sap_mask_fitted_planet_local, kk, mask_fitted_planet)
                    for jj in range(0, len(sap_detrended_lcs)):
                        sap_detrended_lcs[jj] = np.insert(sap_detrended_lcs[jj], kk, np.nan)
                    
                    break
                
                elif kk + 1  == len(sap_x_local):
                    sap_x_local = np.insert(sap_x_local, kk+1, time)
                    sap_yerr_local = np.insert(sap_yerr_local, kk+1, yerr)
                    sap_mask_local = np.insert(sap_mask_local, kk+1, mask)
                    sap_mask_fitted_planet_local = np.insert(sap_mask_fitted_planet_local, kk+1, mask_fitted_planet)
                    for jj in range(0, len(sap_detrended_lcs)):
                        sap_detrended_lcs[jj] = np.insert(sap_detrended_lcs[jj], kk+1, np.nan)
                    
                    break
            
            
    
    
    for ii in range(0, len(sap_x_local)):
        time = sap_x_local[ii]
        yerr = sap_yerr_local[ii]
        mask = sap_mask_local[ii]
        mask_fitted_planet = sap_mask_fitted_planet_local[ii]
        if time not in pdc_x_local:
            for kk in range(0, len(pdc_x_local)):
                if pdc_x_local[kk] > time:
                    pdc_x_local = np.insert(pdc_x_local, kk, time)
                    pdc_yerr_local = np.insert(pdc_yerr_local, kk, yerr)
                    pdc_mask_local = np.insert(pdc_mask_local, kk, mask)
                    pdc_mask_fitted_planet_local = np.insert(pdc_mask_fitted_planet_local, kk, mask_fitted_planet)
                    for jj in range(0, len(pdc_detrended_lcs)):
                        pdc_detrended_lcs[jj] = np.insert(pdc_detrended_lcs[jj], kk, np.nan)
                    
                    break
                
                
                elif kk + 1 == len(pdc_x_local):
                    pdc_x_local = np.insert(pdc_x_local, kk+1, time)
                    pdc_yerr_local = np.insert(pdc_yerr_local, kk+1, yerr)
                    pdc_mask_local = np.insert(pdc_mask_local, kk+1, mask)
                    pdc_mask_fitted_planet_local = np.insert(pdc_mask_fitted_planet_local, kk+1, mask_fitted_planet)
                    for jj in range(0, len(pdc_detrended_lcs)):
                        pdc_detrended_lcs[jj] = np.insert(pdc_detrended_lcs[jj], kk+1, np.nan)
                    
                    break

    
                    
    
    print('pdc length out: ', str(len(pdc_x_local)))
    print('sap length out: ', str(len(sap_x_local)))
    

            
    print('')
    print('')
    print('')
    if (pdc_x_local == sap_x_local).all():
        x_detrended = pdc_x_local
        
    else:
        print("ERROR, pdc and sap x arrays aren't the same")
        

    
    yerr_detrended = np.nanmean([pdc_yerr_local, sap_yerr_local], axis=0)
    
    if (pdc_mask_local == sap_mask_local).all():
        mask_detrended = pdc_mask_local
        
    else:
        for ii in range(0, len(pdc_mask_local)):
            if pdc_mask_local[ii] != sap_mask_local[ii]:
                print(pdc_x_local[ii], pdc_mask_local[ii])
                print(sap_x_local[ii], sap_mask_local[ii])
                    
                print('')
        print("ERROR, pdc and sap mask arrays aren't the same")
        
        
    
    if (pdc_mask_fitted_planet_local == sap_mask_fitted_planet_local).all():
        mask_fitted_planet_detrended = pdc_mask_fitted_planet_local
        
    else:
        print("ERROR, pdc and sap mask for fitted planet arrays aren't the same")
        
    
                
       
    return(x_detrended, sap_detrended_lcs, pdc_detrended_lcs, yerr_detrended, mask_detrended, mask_fitted_planet_detrended)





def detrend_all_methods(x_epochs, y_epochs, yerr_epochs, mask_epochs, mask_fitted_planet_epochs, problem_times,
	t0s, period, duration, cadence, save_to_directory, show_plots):


	x_trimmed, y_trimmed, yerr_trimmed, mask_trimmed, mask_fitted_planet_trimmed = \
	trim_jump_times(x_epochs, y_epochs, yerr_epochs, mask_epochs, mask_fitted_planet_epochs, t0s, period, problem_times)

	#### polyam, gp, cofiam friendly mask arrays ####

	friendly_mask_trimmed = []
	for boolean in range(len(mask_trimmed)):
	    friendly_boolean = mask_trimmed[boolean].astype(bool)
	    friendly_mask_trimmed.append(friendly_boolean)
	
	friendly_mask_fitted_planet_trimmed = []
	for boolean in range(len(mask_fitted_planet_trimmed)):
	    friendly_boolean = mask_fitted_planet_trimmed[boolean].astype(bool)
	    friendly_mask_fitted_planet_trimmed.append(friendly_boolean)

	friendly_x_trimmed = []
	for time_array in range(len(x_trimmed)):
	    friendly_time_array = x_trimmed[time_array].astype(float)
	    friendly_x_trimmed.append(friendly_time_array)

	friendly_y_trimmed = []
	for flux_array in range(len(y_trimmed)):
	    friendly_flux_array = y_trimmed[flux_array].astype(float)
	    friendly_y_trimmed.append(friendly_flux_array)

	friendly_yerr_trimmed = []
	for flux_err_array in range(len(yerr_trimmed)):
	    friendly_flux_err_array = yerr_trimmed[flux_err_array].astype(float)
	    friendly_yerr_trimmed.append(friendly_flux_err_array)

	#################################################

	# determine local window values for later use
	# zoom in around local window
	local_x_epochs, local_y_epochs, local_yerr_epochs, \
	local_mask_epochs, local_mask_fitted_planet_epochs = \
	split_around_transits(np.concatenate(x_trimmed, axis=0, dtype=object), 
	                      np.concatenate(y_trimmed, axis=0, dtype=object),
	                      np.concatenate(yerr_trimmed, axis=0, dtype=object),
	                      np.concatenate(mask_trimmed, axis=0, dtype=object),
	                      np.concatenate(mask_fitted_planet_trimmed, axis=0, dtype=object),
	                      t0s, float(6*duration/(24.))/period, period)

	local_x = np.concatenate(local_x_epochs, axis=0, dtype=object)
	local_y = np.concatenate(local_y_epochs, axis=0, dtype=object)
	local_yerr = np.concatenate(local_yerr_epochs, axis=0, dtype=object)
	local_mask = np.concatenate(local_mask_epochs, axis=0, dtype=object)
	local_mask_fitted_planet = np.concatenate(local_mask_fitted_planet_epochs, axis=0, dtype=object)

	####################
	####################
	####################
	# local detrending
	start = time.time()
	print('')
	print('detrending via the local method')
	local_detrended = \
	local_method(x_trimmed, y_trimmed, yerr_trimmed, 
	             mask_trimmed, mask_fitted_planet_trimmed,
	             t0s, duration, period)


	# remove outliers in unmasked local detrended lc
	local_x_no_outliers, local_detrended_no_outliers = \
	reject_outliers_everywhere(local_x, local_detrended, local_yerr, 5*cadence, 5, 10)

	plot_individual_outliers(local_x, local_detrended, local_x_no_outliers, local_detrended_no_outliers,
	                         t0s, period, float(6*duration/(24.))/period, 0.009, save_to_directory)

	end = time.time()
	print('local detrending took ' + str(np.round(end - start, 2)) + ' seconds')
	




	####################
	####################
	####################
	# polyAM detrending
	start = time.time()
	print('')
	print('detrending via the polyAM method')
	poly_detrended, poly_DWs = \
	polynomial_method(friendly_x_trimmed, friendly_y_trimmed, friendly_yerr_trimmed, 
	                  friendly_mask_trimmed, friendly_mask_fitted_planet_trimmed,
	                  t0s, duration, period, local_x_epochs)


	# remove outliers in unmasked poly detrended lc
	poly_x_no_outliers, poly_detrended_no_outliers = \
	reject_outliers_everywhere(local_x, poly_detrended, local_yerr, 5*cadence, 5, 10)

	plot_individual_outliers(local_x, poly_detrended, poly_x_no_outliers, poly_detrended_no_outliers,
	                         t0s, period, float(6*duration/(24.))/period, 0.009, save_to_directory)

	end = time.time()
	print('polyAM detrending took ' + str(np.round(end - start, 2)) + ' seconds')






	####################
	####################
	####################
	# gp detrending
	start = time.time()
	print('')
	print('detrending via the GP method')
	gp_detrended = \
	gp_method(friendly_x_trimmed, friendly_y_trimmed, friendly_yerr_trimmed, 
	          friendly_mask_trimmed, friendly_mask_fitted_planet_trimmed,
	          t0s, duration, period)

	# remove outliers in unmasked gp detrended lc
	gp_x_no_outliers, gp_detrended_no_outliers = \
	reject_outliers_everywhere(local_x, gp_detrended, local_yerr, 5*cadence, 5, 10)

	plot_individual_outliers(local_x, gp_detrended, gp_x_no_outliers, gp_detrended_no_outliers,
	                         t0s, period, float(6*duration/(24.))/period, 0.009, save_to_directory)

	end = time.time()
	print('GP detrending took ' + str(np.round(end - start, 2)) + ' seconds')





	####################
	####################
	####################
	# CoFiAM detrending
	start = time.time()
	print('')
	print('detrending via the CoFiAM method')
	cofiam_detrended, cofiam_DWs = \
	cofiam_method(friendly_x_trimmed, friendly_y_trimmed, friendly_yerr_trimmed, 
	              friendly_mask_trimmed, friendly_mask_fitted_planet_trimmed,
	              t0s, duration, period, local_x_epochs)

	# remove outliers in unmasked CoFiAM detrended lc
	cofiam_x_no_outliers, cofiam_detrended_no_outliers = \
	reject_outliers_everywhere(local_x, cofiam_detrended, local_yerr, 5*cadence, 5, 10)

	plot_individual_outliers(local_x, cofiam_detrended, cofiam_x_no_outliers, cofiam_detrended_no_outliers,
	                         t0s, period, float(6*duration/(24.))/period, 0.009, save_to_directory)

	end = time.time()
	print('CoFiAM detrending took ' + str(np.round(end - start, 2)) + ' seconds')







	return local_x, local_y, local_yerr, local_mask, local_mask_fitted_planet, \
	local_detrended, local_x_no_outliers, local_detrended_no_outliers, \
	poly_detrended, poly_DWs, poly_x_no_outliers, poly_detrended_no_outliers, \
	gp_detrended, gp_x_no_outliers, gp_detrended_no_outliers, \
	cofiam_detrended, cofiam_DWs, cofiam_x_no_outliers, cofiam_detrended_no_outliers



def detrend_sap_and_pdc(sap_values, pdc_values, save_dir, pop_out_plots):

	# assumes order of sap, pdc arrays are as follows:
	# [pdc_x_epochs, pdc_y_epochs, pdc_yerr_epochs, pdc_mask_epochs, \
	# pdc_mask_fitted_planet_epochs, pdc_problem_times, pdc_t0s, pdc_period, \
	# pdc_duration, pdc_cadence]

	detrended_sap_vals = detrend_all_methods(x_epochs = sap_values[0], y_epochs = sap_values[1], yerr_epochs = sap_values[2],
		mask_epochs = sap_values[3], mask_fitted_planet_epochs = sap_values[4], problem_times = sap_values[5],
		t0s = sap_values[6], period = sap_values[7], duration = sap_values[8], cadence = sap_values[9],
		save_to_directory = save_dir, show_plots = pop_out_plots)

	detrended_pdc_vals = detrend_all_methods(x_epochs = pdc_values[0], y_epochs = pdc_values[1], yerr_epochs = pdc_values[2],
		mask_epochs = pdc_values[3], mask_fitted_planet_epochs = pdc_values[4], problem_times = pdc_values[5],
		t0s = pdc_values[6], period = pdc_values[7], duration = pdc_values[8], cadence = pdc_values[9],
		save_to_directory = save_dir, show_plots = pop_out_plots)

	return detrended_sap_vals, detrended_pdc_vals




def detrend_one_lc(lc_values, save_dir, pop_out_plots):

	# assumes order of sap, pdc arrays are as follows:
	# [pdc_x_epochs, pdc_y_epochs, pdc_yerr_epochs, pdc_mask_epochs, \
	# pdc_mask_fitted_planet_epochs, pdc_problem_times, pdc_t0s, pdc_period, \
	# pdc_duration, pdc_cadence]


	detrended_lc_vals = detrend_all_methods(x_epochs = lc_values[0], y_epochs = lc_values[1], yerr_epochs = lc_values[2],
		mask_epochs = lc_values[3], mask_fitted_planet_epochs = lc_values[4], problem_times = lc_values[5],
		t0s = lc_values[6], period = lc_values[7], duration = lc_values[8], cadence = lc_values[9],
		save_to_directory = save_dir, show_plots = pop_out_plots)

	return detrended_lc_vals




