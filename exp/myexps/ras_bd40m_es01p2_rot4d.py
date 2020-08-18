#Runscript Initialization and Parameters for Isca

import os, sys
from isca import DiagTable, Experiment, Namelist, GFDL_BASE, GFDL_DATA
from isca.codebase import SocratesCodeBase
from func_sub import func_sub

if '-jsub' in sys.argv:
 jsub = sys.argv[sys.argv.index('-jsub')+1] == 'True'
else: 
 jsub = False

fname = os.path.basename(__file__) #retrieves name of this file
exp_name = fname[:-3] #names experiment based on filename; removes the trailing ".py"
cbase = 'soc_isca'

cb = SocratesCodeBase.from_directory(GFDL_BASE)

exp = Experiment(exp_name, codebase=cb)
exp.clear_rundir()

#Basic Experiment Parameters
NCORES = 8                              # Number of cores to use
NRUNS = 25                             # Number of runs to run the model for
RESOLUTION = 'T21', 15                  # T21 (orig T42) horizontal resolution, 15 (orig 25) levels in pressure
startrun = 1                            # Specify the starting run. If RESTART True, it is based on the present restart files for the experiment
RESTART = False                          # Whether to use the latest restart file
OVERWRITE = True                        # Whether to overwrite data files from previous runs
uselastyear = True                      # Whether to specify year cap on experiment length
lastyear = 25                          # Specify year to end on

funit = 'yearly'                        # time unit of runs (and resultant output file) (uses the -ly ending)
outunit = 'days'                        # time unit for output steps
outstep = 1                             # number of outunits between each output
rununit = 'days'                        # unit by which to measure length of one run
runlen = 30*12                          # number of rununits in one model run 
inputfiles = [os.path.join(GFDL_BASE,'exp/my_inputs/land_all_lowres.nc'),os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]
tcmax = 500.                            # Max temperature in runs (degC)
tcmin = -250.                           # Min temperature (-273 is absolute 0)
wmax = 500.                             # max valid wind speed, m/s
valb = False                            # variable albedo flag
bucket_depth = 40                      # Bucket depth, in meters
es0 = 1.2                               # es0 as a factor of current Earth
omega_scale = 4                       # Scaling of rotation rate relative to Earth's (rotation in days)
omega = 7.2921150e-5 / omega_scale 

dt_day = 86400.                         # Time interval of day, in sec
dt_atmos = 120.                         # Time interval of atmosphere calculations (sec)
dt_rad = 2400.                          # Time interval of radiation scheme (sec)

max_bucket = 64.*32.*bucket_depth       # Max bucket depth, set to be global starting capacity

if RESTART:
	startrun = 1
else:
	restartpoint = 1

if not RESTART and startrun > restartpoint:
	startrun = restartpoint

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
diag.add_field('atmosphere', 'bucket_depth', time_avg=True)
diag.add_field('atmosphere','bucket_depth_conv', time_avg=True)
diag.add_field('atmosphere','bucket_depth_cond', time_avg=True)
diag.add_field('atmosphere','bucket_depth_lh', time_avg=True)
diag.add_field('atmosphere','dt_qg_condensation', time_avg=True)
diag.add_field('atmosphere','convection_rain', time_avg=True)
diag.add_field('atmosphere','condensation_rain', time_avg=True)
diag.add_field('atmosphere','virga', time_avg=True)
#diag.add_field('atmosphere','cape', time_avg=True)
#diag.add_field('atmosphere','cin', time_avg=True)
diag.add_field('atmosphere','rh', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
if valb:
 diag.add_field('mixed_layer', 'albedo', time_avg=True)
diag.add_field('mixed_layer', 'ml_heat_cap', time_avg=True)
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
#diag.add_field('ras', 'pcldb', time_avg=True)
#diag.add_field('ras', 'prec_conv', time_avg=True)
#diag.add_field('ras', 'snow_conv', time_avg=True)
#diag.add_field('ras', 'mc', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)

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
        'account_for_effect_of_water':True
        #'solday': 90
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
        'two_stream_gray':False, #Don't use grey radiation
        'do_rrtm_radiation':False, #Don't use RRTM radiation
        'do_socrates_radiation': True,
        'convection_scheme':'RAS', #Use the specified convection scheme (DRY,SIMPLE_BETTS_MILLER,FULL_BETTS_MILLER,RAS)
        'land_option':'input', #Use land mask from input file
        'land_file_name': 'INPUT/land_all_lowres.nc', #Tell model where to find input file
        'hydro_scheme':'BUCKET_HYDRO', #Run with the bucket model
        'init_bucket_depth_land':bucket_depth, #Set initial bucket depth over land (assume meters)
        'max_bucket_depth_land':max_bucket, #Set max bucket depth over land
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
        'do_simple': False,
    },

    'surface_flux_nml': {
        'use_virtual_temp': True,
        'do_simple': False,
        'old_dtaudv': False,
        'exposed_buckets': True   #flag to treat land buckets like ocean buckets    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,    
        'depth':40., #Mixed layer depth (m)
        'land_option':'input',    #Tell mixed layer to get land mask from input file
        'land_h_capacity_prefactor': 0.1, #What factor to multiply mixed-layer depth by over land. 
        'albedo_value': 0.3, #Ocean albedo value
        'land_albedo_prefactor' : 1.0, #What factor to multiply ocean albedo by over land or for static albedo
        'do_qflux' : False, #Do not use prescribed qflux formula
        'var_hcap' : True, #variable heat capacity flag
	'var_alb'  : valb  #variable albedo flag
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
        'valid_range_t':[tcmin+273.15,tcmax+273.15],
        'vrange':[-1.*wmax,wmax], #mmm wind speed valid range, in m/s
        'initial_sphum':[0.0],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.2, #Parameter that sets the vertical distribution of sigma levels
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03#,
        #'lon_max':64,
        #'lat_max':32,
        #'num_fourier':21,
        #'num_spherical':22,
        #'num_levels':15
    },
    
    'constants_nml': {
        'omega':omega,
        'es0':es0
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

