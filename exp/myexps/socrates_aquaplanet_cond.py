import os

import numpy as np

from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress

NCORES = 16
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = SocratesCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp = Experiment('soc_water_conservation_test_aquaplanet_5_finite_max_no_diffuse_lscale_no_evap_no_conv_current_extra_out_no_reevap_3', codebase=cb)
exp.clear_rundir()

inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')
diag.add_file('atmos_daily', 1, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'sphum_vert_int', time_avg=True)
# diag.add_field('dynamics', 'ucomp', time_avg=True)
# diag.add_field('dynamics', 'vcomp', time_avg=True)
# diag.add_field('dynamics', 'temp', time_avg=True)
# diag.add_field('dynamics', 'vor', time_avg=True)
# diag.add_field('dynamics', 'div', time_avg=True)

diag.add_field('atmosphere', 'bucket_depth', time_avg=True)
diag.add_field('atmosphere','bucket_depth_conv', time_avg=True)
diag.add_field('atmosphere','bucket_depth_cond', time_avg=True)
diag.add_field('atmosphere','bucket_depth_lh', time_avg=True)
diag.add_field('atmosphere','bucket_depth_neg_buc', time_avg=True)
diag.add_field('atmosphere','bucket_depth_neg_buc_p', time_avg=True)
diag.add_field('atmosphere','bucket_depth_neg_buc_c', time_avg=True)
diag.add_field('atmosphere','bucket_depth_neg_buc_f', time_avg=True)
diag.add_field('atmosphere','empty_bucket', time_avg=True)
diag.add_field('atmosphere','atm_water_change_cond', time_avg=True)

diag.add_field('atmosphere','mean_dt_bucket', time_avg=True)
diag.add_field('atmosphere','mean_bucket_previous', time_avg=True)
diag.add_field('atmosphere','mean_bucket_current', time_avg=True)
diag.add_field('atmosphere','mean_bucket_future', time_avg=True)
diag.add_field('atmosphere','mean_bucket_previous_post_filter', time_avg=True)
diag.add_field('atmosphere','mean_bucket_current_post_filter', time_avg=True)
diag.add_field('atmosphere','mean_bucket_future_post_filter', time_avg=True)
diag.add_field('atmosphere','mean_bucket_future_post_filter_2', time_avg=True)
diag.add_field('atmosphere','mean_bucket_future_post_filter_3', time_avg=True)
diag.add_field('atmosphere','dt_bucket_actual', time_avg=True)

#mmm add sphum outputs
diag.add_field('atmosphere', 'sphum1_current', time_avg=True)
diag.add_field('atmosphere', 'sphum1_future', time_avg=True)
diag.add_field('atmosphere', 'sphum2_current', time_avg=True)
diag.add_field('atmosphere', 'sphum2_future', time_avg=True)
diag.add_field('atmosphere', 'sphum3_current', time_avg=True)
diag.add_field('atmosphere', 'sphum3_future', time_avg=True)
diag.add_field('atmosphere', 'sphum4_current', time_avg=True)
diag.add_field('atmosphere', 'sphum4_future', time_avg=True)

#radiative tendencies
# diag.add_field('socrates', 'soc_tdt_lw', time_avg=True)
# diag.add_field('socrates', 'soc_tdt_sw', time_avg=True)
# diag.add_field('socrates', 'soc_tdt_rad', time_avg=True)

#net (up) and down surface fluxes
# diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
# diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
# diag.add_field('socrates', 'soc_surf_flux_lw_down', time_avg=True)
# diag.add_field('socrates', 'soc_surf_flux_sw_down', time_avg=True)
#net (up) TOA and downard fluxes
diag.add_field('socrates', 'soc_olr', time_avg=True)
# diag.add_field('socrates', 'soc_toa_sw', time_avg=True) 
# diag.add_field('socrates', 'soc_toa_sw_down', time_avg=True)

# additional output options commented out 
#diag.add_field('socrates', 'soc_flux_lw', time_avg=True)
#diag.add_field('socrates', 'soc_flux_sw', time_avg=True)
#diag.add_field('socrates', 'soc_co2', time_avg=True)
#diag.add_field('socrates', 'soc_ozone', time_avg=True) 
#diag.add_field('socrates', 'soc_coszen', time_avg=True) 
#diag.add_field('socrates', 'soc_spectral_olr', time_avg=True)

exp.diag_table = diag
exp.inputfiles = inputfiles

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    'socrates_rad_nml': {
        'stellar_constant':1370.,
        'lw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'),
        'sw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'),
        'do_read_ozone': True,
        'ozone_file_name':'ozone_1990',
        'ozone_field_name':'ozone_1990',
        'dt_rad':3600,
        'store_intermediate_rad':True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
        'tidally_locked':False,
        'solday': 90
    }, 
    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':2.e-4, #Ocean roughness lengths
        'roughness_heat':2.e-4,
        'roughness_moist':2.e-4,             
        'two_stream_gray': False,     #Use the grey radiation scheme
        'do_socrates_radiation': True,
        'convection_scheme': 'None', #Use simple Betts miller convection            
        'bucket':True, #Run with the bucket model
        'init_bucket_depth_land':  0.,
        'init_bucket_depth':  0.,
        'max_bucket_depth_land': 200.,
        'land_option': 'all_land',
        'finite_bucket_depth_over_land':False,
        'do_lscale_cond':True,
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':False,  
        'depth': 4,                          #Depth of mixed layer used
        'albedo_value': 0.3,                  #Albedo value used      
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':False
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  150., #Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,      
    },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 0.0,
        'reference_sea_level_press':1.0e5,
        'num_levels':40,      #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[9.8e-5],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    },

})

#Lets do a run!
if __name__=="__main__":

        cb.compile()
        #Set up the experiment object, with the first argument being the experiment name.
        #This will be the name of the folder that the data will appear in.
        exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
        for i in range(2,61):
            exp.run(i, num_cores=NCORES,overwrite_data=True)
