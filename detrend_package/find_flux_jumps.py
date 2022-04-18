import numpy as np
import matplotlib.pyplot as plt
import exoplanet as xo
from scipy.interpolate import interp1d
from matplotlib.widgets import Slider, Button


print(f"exoplanet.__version__ = '{xo.__version__}'")

from get_lc import *
from helper_functions import *
from outlier_rejection import *
from manipulate_data import *
from plot import *



def find_flux_jumps(star_id, flux_type, TESS = False, Kepler = False, 
                    user_periods = None, user_t0s = None, user_durations = None,
                    planet_number = 1, mask_width = 1.3):



    #pulls in light curve
    time, lc, lc_err, mask, mask_fitted_planet, \
    t0s, period, duration, quarters, crowding, flux_fraction  = \
    get_light_curve(star_id, flux_type, TESS, Kepler, 
                    user_periods, user_t0s, user_durations,
                    planet_number, mask_width)



    #determine cadence of observation
    cadence = determine_cadence(time)

    #find end time of quarters
    quarters_end = [el[1] for el in quarters]


    time_out, flux_out, flux_err_out, mask_out, mask_fitted_planet_out, moving_median = \
    reject_outliers_out_of_transit(time, lc, lc_err, mask, mask_fitted_planet, 30*cadence, 4)

    plot_outliers(time, lc, time_out, flux_out, 
                  moving_median, quarters_end)

    plt.show()

    x_quarters, y_quarters, yerr_quarters, mask_quarters, mask_fitted_planet_quarters = \
    split_around_problems(time_out, flux_out, flux_err_out, 
                          mask_out, mask_fitted_planet_out, quarters_end)

    plot_split_data(x_quarters, y_quarters, t0s)
    plt.show()

    x_quarters_w_transits, y_quarters_w_transits, yerr_quarters_w_transits, \
    mask_quarters_w_transits, mask_fitted_planet_quarters_w_transits = \
    find_quarters_with_transits(x_quarters, y_quarters, yerr_quarters, 
                                mask_quarters, mask_fitted_planet_quarters, t0s)



    x_quarters_w_transits = np.concatenate(x_quarters_w_transits, axis=0, dtype='float64')
    y_quarters_w_transits = np.concatenate(y_quarters_w_transits, axis=0, dtype='float64')
    yerr_quarters_w_transits = np.concatenate(yerr_quarters_w_transits, axis=0, dtype='float64')
    mask_quarters_w_transits = np.concatenate(mask_quarters_w_transits, axis=0, dtype=bool)
    mask_fitted_planet_quarters_w_transits = np.concatenate(mask_fitted_planet_quarters_w_transits, axis=0, dtype=bool)


    x_transits, y_transits, yerr_transits, mask_transits, mask_fitted_planet_transits = split_around_transits(x_quarters_w_transits, 
                                                                                                              y_quarters_w_transits, 
                                                                                                              yerr_quarters_w_transits, 
                                                                                                              mask_quarters_w_transits, 
                                                                                                              mask_fitted_planet_quarters_w_transits, 
                                                                                                              t0s, 1./2., period)

    x_epochs = np.concatenate(x_transits, axis=0, dtype='float64')
    y_epochs = np.concatenate(y_transits, axis=0, dtype='float64')
    yerr_epochs = np.concatenate(yerr_transits, axis=0, dtype='float64')
    mask_epochs = np.concatenate(mask_transits, axis=0, dtype=bool)
    mask_fitted_planet_epochs = np.concatenate(mask_fitted_planet_transits, axis=0, dtype=bool)




    _, _, problem_times = plot_transits(x_transits, y_transits, mask_transits, t0s, period, cadence*5)



    return x_epochs, y_epochs, yerr_epochs, mask_epochs, mask_fitted_planet_epochs, problem_times, t0s, period, duration, cadence





def find_sap_and_pdc_flux_jumps(star_id, TESS = False, Kepler = False, 
                                user_periods = None, user_t0s = None, user_durations = None,
                                planet_number = 1, mask_width = 1.3):


    sap_vals = find_flux_jumps(star_id, 'sap_flux', TESS, Kepler, 
                               user_periods, user_t0s, user_durations,
                               planet_number, mask_width)    

    pdc_vals = find_flux_jumps(star_id, 'pdcsap_flux', TESS, Kepler, 
                               user_periods, user_t0s, user_durations,
                               planet_number, mask_width)


    return sap_vals, pdc_vals







#pulls in light curve
[[sap_x_epochs, sap_y_epochs, sap_yerr_epochs, sap_mask_epochs, \
sap_mask_fitted_planet_epochs, sap_problem_times, sap_t0s, sap_period, \
sap_duration, sap_cadence], \
[pdc_x_epochs, pdc_y_epochs, pdc_yerr_epochs, pdc_mask_epochs, \
pdc_mask_fitted_planet_epochs, pdc_problem_times, pdc_t0s, pdc_period, \
pdc_duration, pdc_cadence]]  = \
find_sap_and_pdc_flux_jumps('toi-1130', TESS=True, planet_number = 2, mask_width = 1.8)






