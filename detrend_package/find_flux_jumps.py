#pulls in light curve
sap_time, sap_lc, sap_lc_err, sap_mask, sap_mask_fitted_planet, \
t0s_sap, period, duration, quarters, crowding, flux_fraction  = \
get_light_curve('toi-1130', 'qlp', TESS=True, planet_number = 1, mask_width = 1.8)
                #user_periods = np.array([period]),
                #user_t0s = np.array([t0]),
                #user_durations = np.array([duration]))



#determine cadence of observation
cadence = determine_cadence(sap_time)

#find end time of quarters
quarters_end = [el[1] for el in quarters]


sap_time_out, sap_flux_out, sap_flux_err_out, sap_mask_out, sap_mask_fitted_planet_out, sap_moving_median = \
reject_outliers_out_of_transit(sap_time, sap_lc, sap_lc_err, sap_mask, sap_mask_fitted_planet, 30*cadence, 4)

plot_outliers(sap_time, sap_lc, sap_time_out, sap_flux_out, 
              sap_moving_median, quarters_end)


x_quarters, y_quarters, yerr_quarters, mask_quarters, mask_fitted_planet_quarters = \
split_around_problems(sap_time_out, sap_flux_out, sap_flux_err_out, 
                      sap_mask_out, sap_mask_fitted_planet_out, quarters_end)

plot_split_data(x_quarters, y_quarters, t0s_sap)


x_quarters_w_transits, y_quarters_w_transits, yerr_quarters_w_transits, \
mask_quarters_w_transits, mask_fitted_planet_quarters_w_transits = \
find_quarters_with_transits(x_quarters, y_quarters, yerr_quarters, 
                            mask_quarters, mask_fitted_planet_quarters, t0s_sap)



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
                                                                                                          t0s_sap, 1./2., period)

x_epochs = np.concatenate(x_transits, axis=0, dtype='float64')
y_epochs = np.concatenate(y_transits, axis=0, dtype='float64')
yerr_epochs = np.concatenate(yerr_transits, axis=0, dtype='float64')
mask_epochs = np.concatenate(mask_transits, axis=0, dtype=bool)
mask_fitted_planet_epochs = np.concatenate(mask_fitted_planet_transits, axis=0, dtype=bool)

