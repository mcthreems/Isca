#Runscript Initialization and Parameters for Isca

import os, sys
from isca import DiagTable, Experiment, Namelist, GFDL_BASE, GFDL_DATA
from isca.codebase import SocratesCodeBase
from func_sub import func_sub
from get_start import get_start

if 'sub' in sys.argv:
 jsub = True
else: 
 jsub = False

fname = os.path.basename(__file__) #retrieves name of this file
exp_name = fname[:-3] #names experiment based on filename; removes the trailing ".py"
cbase = 'soc_isca'

cb = SocratesCodeBase.from_directory(GFDL_BASE)

exp = Experiment(exp_name, codebase=cb)
exp.clear_rundir()

#Basic Experiment Parameters
NCORES = 16                             # Number of cores to use
NRUNS = 25                              # Number of runs to run the model for
RESOLUTION = 'T21', 25                  # T21 (orig T42) horizontal resolution, 15 (orig 25) levels in pressure
startrun = 1                           # Specify the starting run
#startrun = get_start(exp_name)          # get start run if restarting from previous place 
if startrun == 1:
 RESTART = False                        # flag to use restart file if not starting at run 1
else:
 RESTART = True
OVERWRITE = True                        # Whether to overwrite data files from previous runs
uselastyear = True                      # Whether to specify year cap on experiment length
lastyear = 25                           # Specify year to end on

funit = 'yearly'                        # time unit of runs (and resultant output file) (uses the -ly ending)
outunit = 'days'                        # time unit for output steps
outstep = 1                             # number of outunits between each output
rununit = 'days'                        # unit by which to measure length of one run
runlen = 30*12                          # number of rununits in one model run 
inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]
tcmax = 500.                            # Max temperature in runs (degC)
tcmin = -250.                           # Min temperature (-273 is absolute 0)
wmax = 500.                             # max valid wind speed, m/s
vsurfprop = True                        # variable surface properties
atm_bucket = 0.01                        # initial depth of water in atmosphere (m)
bucket_depth = 100 - atm_bucket        # Bucket depth, in meters
isphum = atm_bucket*9.8/100.            # initial specific humidity (kg/kg)
beta = 0.7                              # beta parameter from Fan et al
es0 = 2                               # es0 as a factor of current Earth
omega_scale = 2                         # Scaling of rotation rate relative to Earth's (rotation in days)
omega = 7.2921150e-5 / omega_scale 

dt_day = 86400.                         # Time interval of day, in sec
dt_atmos = 120.                         # Time interval of atmosphere calculations (sec)
dt_rad = 2400.                          # Time interval of radiation scheme (sec)

max_bucket = 64.*32.*bucket_depth       # Max bucket depth, set to be global starting capacity

if uselastyear:
	lastrun = lastyear
	nruns = lastrun - startrun + 1
	if nruns < 1:
		startrun = lastrun
		NRUNS = 1
	elif NRUNS > nruns:
		NRUNS = nruns

print('Starting at run #'+str(startrun))
print('Finishing at run #'+str(startrun+NRUNS-1))

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_'+funit, outstep, outunit, time_units=outunit)

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')

diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere','dt_qg_condensation', time_avg=True)
diag.add_field('atmosphere','dt_qg_convection', time_avg=True)
diag.add_field('atmosphere','convection_rain', time_avg=True)
diag.add_field('atmosphere','condensation_rain', time_avg=True)
diag.add_field('atmosphere','virga', time_avg=True)
diag.add_field('atmosphere','cape', time_avg=True)
diag.add_field('atmosphere','cin', time_avg=True)
diag.add_field('atmosphere','rh', time_avg=True)
diag.add_field('atmosphere','dt_qg_diffusion', time_avg=True)
diag.add_field('atmosphere', 'cwv', time_avg=True)
diag.add_field('atmosphere', 'dtend', time_avg=True)
diag.add_field('atmosphere', 'cvtend', time_avg=True)
diag.add_field('atmosphere', 'cdtend', time_avg=True)

diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'zsurf')

diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True)

diag.add_field('socrates', 'soc_olr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw', time_avg=True) 
diag.add_field('socrates', 'soc_toa_sw_down', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw_down', time_avg=True)
diag.add_field('socrates', 'soc_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)

#Define values for the 'core' namelist
namelist = Namelist({
    'main_nml': {
        rununit   : runlen,
        'dt_atmos':dt_atmos,
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
        'dt_rad':dt_rad,
        'store_intermediate_rad':True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
        'tidally_locked':False,
        'account_for_effect_of_water':True,
        'solday': 0
    }, 
    
    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': False,
        'roughness_mom':2.e-4, #Ocean roughness lengths
        'roughness_heat':2.e-4,
        'roughness_moist':2.e-4,      
        'two_stream_gray':False, #Don't use grey radiation
        'do_rrtm_radiation':False, #Don't use RRTM radiation
        'do_socrates_radiation': True,
        'convection_scheme':'FULL_BETTS_MILLER', #Use the specified convection scheme (DRY,SIMPLE_BETTS_MILLER,FULL_BETTS_MILLER,RAS)
        'hydro_scheme':'NONE', #Run as mixed-layer aquaplanet
        'var_surf_prop':vsurfprop, #whether to adjust albedo based on water in grid cell
        'vis_albedo':0.3, #albedo default, applies to all cells if var_surf_prop is False, to dry cells if True
        'liq_albedo':0.3,
        'evap_thresh':1.e-3
    },

    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,  
        'depth': bucket_depth,                          #Depth of mixed layer used
        'albedo_value': 0.38,                  #Albedo value used      
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': False,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':True,
        'do_simple':False,
    },

    'surface_flux_nml': {
        'use_virtual_temp':True,
        'do_simple':False,
        'old_dtaudv':False,
        'land_humidity_prefactor':beta
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':tcmin+273.15,
        'Tmax':tcmax+273.15   
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_not_calculate':False, #whether to calculate esat
        'do_simple':True,
        'tcmax_simple':tcmax, #set max temp for esat lookup table
        'tcmin_simple':tcmin, #set min
        'show_bad_value_count_by_slice':True,
        'show_all_bad_values':True #want to see if temps are too high or too low
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
        'damping_coeff':1.15740741e-4, # default value is 1.15740741e-4
        'damping_coeff_vor':-1.,       # -1 is a flag making it equal to damping_coeff   
        'damping_coeff_div':-1.,       # -1 is a flag making it equal to damping_coeff
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':40,
        'valid_range_t':[tcmin+273.15,tcmax+273.15],
        'vrange':[-1.*wmax,wmax], #mmm wind speed valid range, in m/s
        'initial_sphum':[isphum], #testing an initial sphum here to see if it affects the start
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    },

    'constants_nml': {
        'omega':omega,
        'es0':es0,
        'bparam':beta,  # beta param from Fan et al
        'use_q_max':False
    }
})

exp.diag_table = diag
exp.inputfiles = inputfiles
exp.namelist = namelist
exp.set_resolution(*RESOLUTION)

if __name__=="__main__" and jsub:
 func_sub(exp_name)
 os.system('qsub -pe dc* '+str(NCORES)+' -N out_'+exp_name+' sub_'+exp_name)

if __name__=="__main__" and not jsub:

        cb.compile()
        #Set up the experiment object, with the first argument being the experiment name.
        #This will be the name of the folder that the data will appear in.
        exp.run(startrun, use_restart=RESTART, num_cores=NCORES, overwrite_data=OVERWRITE)
        for i in range(startrun+1,startrun+NRUNS):
            exp.run(i, num_cores=NCORES,overwrite_data=OVERWRITE)

