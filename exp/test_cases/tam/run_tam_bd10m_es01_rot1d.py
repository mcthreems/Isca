#Runscript Initialization and Parameters for Isca

import os
from isca import DiagTable, Namelist, GFDL_BASE, GFDL_DATA
from func_run_sub import run_sub
from get_startrun import get_start
from holdvars import jid, hold

fname = os.path.basename(__file__) #retrieves name of this file
exp_name = fname[4:-3] #names experiment based on filename; removes the leading "run_" and trailing ".py"
cbase = 'socrates'

#Basic Experiment Parameters
NCORES = 8                              # Number of cores to use
NRUNS = 25                             # Number of runs to run the model for
RESOLUTION = 'T21', 15                  # T21 (orig T42) horizontal resolution, 15 (orig 25) levels in pressure
startrun = 1                            # Specify the starting run. If RESTART True, it is based on the present restart files for the experiment
RESTART = False                          # Whether to use the latest restart file
OVERWRITE = True                        # Whether to overwrite data files from previous runs
uselastyear = True                      # Whether to specify year cap on experiment length
lastyear = 110                          # Specify year to end on
seq = False				# Whether to submit sequential jobs
if seq:
 jid = ''
 hold = True

funit = 'yearly'                        # time unit of runs (and resultant output file) (uses the -ly ending)
outunit = 'days'                        # time unit for output steps
outstep = 1                             # number of outunits between each output
rununit = 'days'                        # unit by which to measure length of one run
runlen = 30*12                            # number of rununits in one model run 
INFILES = [os.path.join(GFDL_BASE,'exp/experiments/my_exps/input/land_all_lowres.nc'),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]
tcmax = 500                             # Max temperature in esat lookup table (degC)
tcmin = -273                            # Min temperature (absolute zero)
bucket_depth = 10.0                   # Bucket depth, in meters
es0 = 1.0				# es0 as a factor of current Earth
omega_scale = 1.0                       # Scaling of rotation rate relative to Earth's (rotation in days)
omega = 7.2921150e-5 / omega_scale 

dt_day = 86400.                         # Time interval of day, in sec
dt_atmos = 120.                         # Time interval of atmosphere calculations (sec)
dt_rad = 2400.                          # Time interval of radiation scheme (sec)

if RESTART:
	startrun = get_start(exp_name)

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

if __name__ == '__main__':
        run_sub(exp_name,NCORES,jid=jid,seq=seq,hold=hold)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_'+funit, outstep, outunit, time_units=outunit)

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'run', time_avg=True)
diag.add_field('atmosphere','sub', time_avg=True)
diag.add_field('atmosphere','infil', time_avg=True)
diag.add_field('atmosphere','gle', time_avg=True)
diag.add_field('atmosphere','res', time_avg=True)
diag.add_field('atmosphere','dis', time_avg=True)
diag.add_field('atmosphere','recharge', time_avg=True)
diag.add_field('atmosphere','table', time_avg=True)
diag.add_field('atmosphere','sens', time_avg=True)
diag.add_field('atmosphere','evap', time_avg=True)
diag.add_field('atmosphere','cape', time_avg=True)
diag.add_field('atmosphere','cin', time_avg=True)
diag.add_field('atmosphere','rh', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('socrates', 'soc_olr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw', time_avg=True) 
diag.add_field('socrates', 'soc_toa_sw_down', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw_down', time_avg=True)
diag.add_field('socrates', 'soc_flux_sw', time_avg=True)

#Define New Constants for the constants.F90 file
#constants = {
	#'ES0': es0		# Saturation vapor pressure
#}

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
        #'solday': 90
    }, 
    
    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':False,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':2.e-4, #Ocean roughness lengths
        'roughness_heat':2.e-4,
        'roughness_moist':2.e-4,      
        'two_stream_gray':False, #Don't use grey radiation
        'do_rrtm_radiation':False, #Don't use RRTM radiation
        'do_socrates_radiation': True,
        'convection_scheme':'FULL_BETTS_MILLER', #Use the specified convection scheme (SIMPLE_BETTS_MILLER,FULL_BETTS_MILLER,RAS_CONV)
        'land_option':'input', #Use land mask from input file
        'land_file_name': 'INPUT/land_all_lowres.nc', #Tell model where to find input file
        'hydro_scheme':'TAM_HYDRO', #Use the specified hydrology scheme
        'init_surf_liq':bucket_depth, #Set initial liquid depth (meters)
        'init_liq_table':bucket_depth, #Set initial liquid table (meters)
	'var_surf_prop':True #vary surface properties if liquid
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

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True,
        'tcmax_simple':tcmax, #set max temp for esat lookup table
        'tcmin_simple':tcmin, #set min
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
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    },
    
    'constants_nml': {
        'omega':omega,
        'es0':es0
    }
})
