
from datetime import date
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import exoplanet as xo
from scipy.interpolate import interp1d
from scipy.stats import median_abs_deviation
from matplotlib.widgets import Slider, Button
import sys, argparse
import os
import warnings
import ast

warnings.simplefilter('ignore', np.RankWarning)

parser = argparse.ArgumentParser(description="Looks up light curves for TESS and Kepler objects and enables labeling of jump times. \
	Detrends light curves, outputs plots, and saves detrended data as csvs.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('object', type=str, help='TESS or Kepler identifier. ex: "toi-2088"')
parser.add_argument('flux_type', type=str, help='Flux type as a string. Options: pdc, sap, or both')
parser.add_argument('mission', type=str, help='Mission data select. ex: "TESS"')
parser.add_argument('planet_num', type=int, help='Which planet to look at in system. ex: 1')
parser.add_argument('save_to_dir', type=str, help='Directory path to save csvs and figures to.')
parser.add_argument('-d', '--depth', default=0.02, help='Sets depth of detrended plots. Default is 0.02.')
parser.add_argument('-p', '--period', default=None, help='Optionally input period. Otherwise defaults to what \
	Exoplanet Archive can find.')
parser.add_argument('-t', '--t0', default=None, help='Optionally input t0. Otherwise defaults to what Exoplanet Archive can find.')
parser.add_argument('-du', '--duration', default=None, help='Optionally input duration. Otherwise defaults to what \
	Exoplanet Archive can find.')
parser.add_argument('-mw', '--mask_width', default=1.3, help='Sets mask width. Default is 1.3.')
parser.add_argument('-s', '--show_plots', default='True', help='Set whether to show non-problem-time plots.')
parser.add_argument('--dont_bin', default='False', help='if True, then jump time plots data wont be binned.')


args = vars(parser.parse_args())

# user input parameters

input_id = args['object']
flux_type = args['flux_type']
mission_select = args['mission']
input_planet_number = args['planet_num']
input_dir = args['save_to_dir']
input_depth = args['depth']
input_period = args['period']
input_t0 = args['t0']
input_duration = args['duration']
input_mask_width = args['mask_width']
input_show_plots = ast.literal_eval(args['show_plots'])
input_dont_bin = ast.literal_eval(args['dont_bin'])




# # # ---------------------------------------------- now the fun begins ! ! ! ------------------------------------------------ # # #




print(f"exoplanet.__version__ = '{xo.__version__}'")

from find_flux_jumps import *
from get_lc import *
from helper_functions import *
from outlier_rejection import *
from manipulate_data import *
from plot import plot_detrended_lc
from detrend import *


# checks user arguments
if mission_select == 'TESS': 
	tess_bool = True
	kepler_bool = False
else:
	tess_bool = False # assuming Kepler is selected
	kepler_bool = True


# determining figname
# determining today's date
today = date.today()
current_day = today.strftime('%b%d_%Y')

foldername = input_id + '.0' + str(input_planet_number)+ '_' + current_day
path = os.path.join(input_dir, foldername)

os.makedirs(path, exist_ok=True)


# check if we should run both pdc and sap flux
if flux_type == 'both':

    #pulls in light curve
    [[sap_x_epochs, sap_y_epochs, sap_yerr_epochs, sap_mask_epochs, \
    sap_mask_fitted_planet_epochs, sap_problem_times, sap_t0s, sap_period, \
    sap_duration, sap_cadence], \
    [pdc_x_epochs, pdc_y_epochs, pdc_yerr_epochs, pdc_mask_epochs, \
    pdc_mask_fitted_planet_epochs, pdc_problem_times, pdc_t0s, pdc_period, \
    pdc_duration, pdc_cadence]]  = \
    find_sap_and_pdc_flux_jumps(input_id, path + '/', show_plots = input_show_plots, TESS = tess_bool, Kepler = kepler_bool, planet_number = input_planet_number,
    user_periods = input_period, user_t0s = input_t0, user_durations = input_duration, mask_width=input_mask_width, dont_bin=input_dont_bin) 


    # now for detrending!
    [[sap_local_x, sap_local_y, sap_local_yerr, sap_local_mask, sap_local_mask_fitted_planet, \
    sap_local_detrended, sap_local_x_no_outliers, sap_local_detrended_no_outliers, \
    sap_poly_detrended, sap_poly_DWs, sap_poly_x_no_outliers, sap_poly_detrended_no_outliers, \
    sap_gp_detrended, sap_gp_x_no_outliers, sap_gp_detrended_no_outliers, \
    sap_cofiam_detrended, sap_cofiam_DWs, sap_cofiam_x_no_outliers, sap_cofiam_detrended_no_outliers], \
    [pdc_local_x, pdc_local_y, pdc_local_yerr, pdc_local_mask, pdc_local_mask_fitted_planet, \
    pdc_local_detrended, pdc_local_x_no_outliers, pdc_local_detrended_no_outliers, \
    pdc_poly_detrended, pdc_poly_DWs, pdc_poly_x_no_outliers, pdc_poly_detrended_no_outliers, \
    pdc_gp_detrended, pdc_gp_x_no_outliers, pdc_gp_detrended_no_outliers, \
    pdc_cofiam_detrended, pdc_cofiam_DWs, pdc_cofiam_x_no_outliers, pdc_cofiam_detrended_no_outliers]] = \
    detrend_sap_and_pdc(sap_values = [sap_x_epochs, sap_y_epochs, sap_yerr_epochs, sap_mask_epochs, \
    sap_mask_fitted_planet_epochs, sap_problem_times, sap_t0s, sap_period, sap_duration, sap_cadence],
    pdc_values = [pdc_x_epochs, pdc_y_epochs, pdc_yerr_epochs, pdc_mask_epochs, \
    pdc_mask_fitted_planet_epochs, pdc_problem_times, pdc_t0s, pdc_period, \
    pdc_duration, pdc_cadence],
    save_dir = path + '/', pop_out_plots = input_show_plots)

    # now to add nans, in order to make sure pdc and sap arrays are the same length
    x_detrended,\
    [sap_local_detrended2, sap_poly_detrended2, sap_cofiam_detrended2, sap_gp_detrended2],\
    [pdc_local_detrended2, pdc_poly_detrended2, pdc_cofiam_detrended2, pdc_gp_detrended2],\
    yerr_detrended, mask_detrended, mask_fitted_planet_detrended = \
    add_nans_for_missing_data(sap_local_x, [sap_local_detrended, sap_poly_detrended, sap_cofiam_detrended, sap_gp_detrended], sap_local_yerr, sap_local_mask, sap_local_mask_fitted_planet,
        pdc_local_x, [pdc_local_detrended, pdc_poly_detrended, pdc_cofiam_detrended, pdc_gp_detrended], pdc_local_yerr, pdc_local_mask, pdc_local_mask_fitted_planet)


    # # now to plot!
    green2, green1 = '#355E3B', '#18A558'
    blue2, blue1 = '#000080', '#4682B4'
    purple2, purple1 = '#2E0854','#9370DB'
    red2, red1 = '#770737', '#EC8B80'


    colors = [red1, red2,
              blue1, blue2,
              green1, green2,
              purple1, purple2]

    y_detrended = [sap_local_detrended2, pdc_local_detrended2,
                   sap_poly_detrended2, pdc_poly_detrended2,
                   sap_gp_detrended2, pdc_gp_detrended2,
                   sap_cofiam_detrended2, pdc_cofiam_detrended2]

    detrend_label = ['local SAP', 'local PDCSAP',
                     'polyAM SAP', 'polyAM PDCSAP',
                     'GP SAP', 'GP PDCSAP', 
                     'CoFiAM SAP', 'CoFiAM PDCSAP']



    plot_detrended_lc(x_detrended, y_detrended, yerr_detrended, detrend_label,
                      pdc_t0s, float(6*pdc_duration/(24.))/pdc_period, pdc_period,
                      colors, pdc_duration, depth=input_depth,
                      figname = path + '/' + 'individual_detrended.pdf')

    # plotting method_marginalized_detrended data
    y_detrended = np.array(y_detrended)
    y_detrended_transpose = y_detrended.T

    method_marg_detrended = np.nanmedian(y_detrended_transpose, axis=1)
    MAD = median_abs_deviation(y_detrended_transpose, axis=1, scale=1/1.4826, nan_policy = 'omit')

    yerr_detrended = np.sqrt(yerr_detrended.astype(float)**2 + MAD**2)


    plot_detrended_lc(x_detrended, [method_marg_detrended], yerr_detrended, ['method marg'],
                      pdc_t0s, float(6*pdc_duration/(24.))/pdc_period, pdc_period,
                      ['k'], pdc_duration, depth = input_depth,
                      figname = path + '/' + 'method_marg_detrended.pdf')


    # saving detrend data as csv

    detrend_dict = {}

    detrend_dict["time"] = x_detrended
    detrend_dict["yerr"] = yerr_detrended
    detrend_dict["mask"] = mask_detrended
    detrend_dict["method marginalized"] = method_marg_detrended



    for ii in range(0, len(y_detrended)):
        detrend = y_detrended[ii]
        label = detrend_label[ii]
        detrend_dict[label] = detrend
        
        
    detrend_df = pd.DataFrame(detrend_dict)

    detrend_df.to_csv(path + '/' + 'detrended.csv')








# check if we should run just pdc
elif flux_type == 'pdc':


    #pulls in light curve
    [pdc_x_epochs, pdc_y_epochs, pdc_yerr_epochs, pdc_mask_epochs, \
    pdc_mask_fitted_planet_epochs, pdc_problem_times, pdc_t0s, pdc_period, \
    pdc_duration, pdc_cadence]  = \
    find_flux_jumps(input_id, 'pdcsap_flux', path + '/', 
        show_plots = input_show_plots, TESS = tess_bool, Kepler = kepler_bool, 
        planet_number = input_planet_number,user_periods = input_period, 
        user_t0s = input_t0, user_durations = input_duration, 
        mask_width=input_mask_width, no_jump_times=True, dont_bin=input_dont_bin) 

    # now for detrending!


    [pdc_local_x, pdc_local_y, pdc_local_yerr, pdc_local_mask, pdc_local_mask_fitted_planet, \
    pdc_local_detrended, pdc_local_x_no_outliers, pdc_local_detrended_no_outliers, \
    pdc_poly_detrended, pdc_poly_DWs, pdc_poly_x_no_outliers, pdc_poly_detrended_no_outliers, \
    pdc_gp_detrended, pdc_gp_x_no_outliers, pdc_gp_detrended_no_outliers, \
    pdc_cofiam_detrended, pdc_cofiam_DWs, pdc_cofiam_x_no_outliers, pdc_cofiam_detrended_no_outliers] = \
    detrend_one_lc(lc_values = [pdc_x_epochs, pdc_y_epochs, pdc_yerr_epochs, pdc_mask_epochs, \
    pdc_mask_fitted_planet_epochs, pdc_problem_times, pdc_t0s, pdc_period, \
    pdc_duration, pdc_cadence],
    save_dir = path + '/', pop_out_plots = input_show_plots)


    x_detrended = pdc_local_x

    # # now to plot!

    green2, green1 = '#355E3B', '#18A558'
    blue2, blue1 = '#000080', '#4682B4'
    purple2, purple1 = '#2E0854','#9370DB'
    red2, red1 = '#770737', '#EC8B80'


    colors = [red2,
              blue2,
              green2,
              purple2]

    y_detrended = [pdc_local_detrended,
                   pdc_poly_detrended,
                   pdc_gp_detrended,
                   pdc_cofiam_detrended]

    yerr_detrended = pdc_local_yerr
    mask_detrended = pdc_local_mask

    detrend_label = ['local PDCSAP',
                     'polyAM PDCSAP',
                     'GP PDCSAP', 
                     'CoFiAM PDCSAP']



    plot_detrended_lc(x_detrended, y_detrended, yerr_detrended, detrend_label,
                      pdc_t0s, float(6*pdc_duration/(24.))/pdc_period, pdc_period,
                      colors, pdc_duration, depth=input_depth,
                      figname = path + '/' + 'individual_detrended_PDC.pdf')

    # plotting method_marginalized_detrended data
    y_detrended = np.array(y_detrended)
    y_detrended_transpose = y_detrended.T

    method_marg_detrended = np.nanmedian(y_detrended_transpose, axis=1)
    MAD = median_abs_deviation(y_detrended_transpose, axis=1, scale=1/1.4826, nan_policy = 'omit')

    yerr_detrended = np.sqrt(yerr_detrended.astype(float)**2 + MAD**2)


    plot_detrended_lc(x_detrended, [method_marg_detrended], yerr_detrended, ['method marg'],
                      pdc_t0s, float(6*pdc_duration/(24.))/pdc_period, pdc_period,
                      ['k'], pdc_duration, depth = input_depth,
                      figname = path + '/' + 'method_marg_detrended_PDC.pdf')


    # saving detrend data as csv

    detrend_dict = {}

    detrend_dict["time"] = x_detrended
    detrend_dict["yerr"] = yerr_detrended
    detrend_dict["mask"] = mask_detrended
    detrend_dict["method marginalized"] = method_marg_detrended



    for ii in range(0, len(y_detrended)):
        detrend = y_detrended[ii]
        label = detrend_label[ii]
        detrend_dict[label] = detrend
        
        
    detrend_df = pd.DataFrame(detrend_dict)

    detrend_df.to_csv(path + '/' + 'detrended_PDC.csv')







# check if we should run just sap
elif flux_type == 'sap':



    #pulls in light curve
    [sap_x_epochs, sap_y_epochs, sap_yerr_epochs, sap_mask_epochs, \
    sap_mask_fitted_planet_epochs, sap_problem_times, sap_t0s, sap_period, \
    sap_duration, sap_cadence]  = \
    find_flux_jumps(input_id, 'sap_flux', path + '/', 
        show_plots = input_show_plots, TESS = tess_bool, Kepler = kepler_bool, 
        planet_number = input_planet_number,user_periods = input_period, 
        user_t0s = input_t0, user_durations = input_duration, 
        mask_width=input_mask_width, dont_bin=input_dont_bin) 

    # now for detrending!


    [sap_local_x, sap_local_y, sap_local_yerr, sap_local_mask, sap_local_mask_fitted_planet, \
    sap_local_detrended, sap_local_x_no_outliers, sap_local_detrended_no_outliers, \
    sap_poly_detrended, sap_poly_DWs, sap_poly_x_no_outliers, sap_poly_detrended_no_outliers, \
    sap_gp_detrended, sap_gp_x_no_outliers, sap_gp_detrended_no_outliers, \
    sap_cofiam_detrended, sap_cofiam_DWs, sap_cofiam_x_no_outliers, sap_cofiam_detrended_no_outliers] = \
    detrend_one_lc(lc_values = [sap_x_epochs, sap_y_epochs, sap_yerr_epochs, sap_mask_epochs, \
    sap_mask_fitted_planet_epochs, sap_problem_times, sap_t0s, sap_period, sap_duration, sap_cadence],
    save_dir = path + '/', pop_out_plots = input_show_plots)

    x_detrended = sap_local_x

    # # now to plot!

    green2, green1 = '#355E3B', '#18A558'
    blue2, blue1 = '#000080', '#4682B4'
    purple2, purple1 = '#2E0854','#9370DB'
    red2, red1 = '#770737', '#EC8B80'


    colors = [red2,
              blue2,
              green2,
              purple2]

    y_detrended = [sap_local_detrended,
                   sap_poly_detrended,
                   sap_gp_detrended,
                   sap_cofiam_detrended]

    yerr_detrended = sap_local_yerr
    mask_detrended = sap_local_mask

    detrend_label = ['local SAP',
                     'polyAM SAP',
                     'GP SAP', 
                     'CoFiAM SAP']



    plot_detrended_lc(x_detrended, y_detrended, yerr_detrended, detrend_label,
                      sap_t0s, float(6*sap_duration/(24.))/sap_period, sap_period,
                      colors, sap_duration, depth=input_depth,
                      figname = path + '/' + 'individual_detrended_SAP.pdf')

    # plotting method_marginalized_detrended data
    y_detrended = np.array(y_detrended)
    y_detrended_transpose = y_detrended.T

    method_marg_detrended = np.nanmedian(y_detrended_transpose, axis=1)
    MAD = median_abs_deviation(y_detrended_transpose, axis=1, scale=1/1.4826, nan_policy = 'omit')

    yerr_detrended = np.sqrt(yerr_detrended.astype(float)**2 + MAD**2)


    plot_detrended_lc(x_detrended, [method_marg_detrended], yerr_detrended, ['method marg'],
                      sap_t0s, float(6*sap_duration/(24.))/sap_period, sap_period,
                      ['k'], sap_duration, depth = input_depth,
                      figname = path + '/' + 'method_marg_detrended_SAP.pdf')


    # saving detrend data as csv

    detrend_dict = {}

    detrend_dict["time"] = x_detrended
    detrend_dict["yerr"] = yerr_detrended
    detrend_dict["mask"] = mask_detrended
    detrend_dict["method marginalized"] = method_marg_detrended



    for ii in range(0, len(y_detrended)):
        detrend = y_detrended[ii]
        label = detrend_label[ii]
        detrend_dict[label] = detrend
        
        
    detrend_df = pd.DataFrame(detrend_dict)

    detrend_df.to_csv(path + '/' + 'detrended_SAP.csv')


else:
    print('ERROR!')
    print('invalid flux_type value entered...options are: pdc, sap, or both')




