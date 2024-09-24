"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np
from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw.data import ForceDry
#from clawpack.geoclaw import fgmax_tools  #new version in same directory
import fgmax_tools
from clawpack.geoclaw import kmltools
from clawpack.geoclaw import fgout_tools

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')


#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')


    #------------------------------------------------------------------
    # GeoClaw specific parameters are set later
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    # x values should be integer multipes of 1/9"
    # y values should be integer multipes of 1/9"
    # Note: always satisfied if limits are multiples of 0.01 degree

    arcsec16 =1./(6*3600.)    # choose domain and offset edges by half a 1/3" cell so
    # cell centers are exactly at DEM grid points for both 1/9" and 1/3":

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = -135.16 - arcsec16      # west longitude
    clawdata.upper[0] = -122.16 - arcsec16      # east longitude

    clawdata.lower[1] = 38.5 - arcsec16       # south latitude
    clawdata.upper[1] = 53.5 - arcsec16         # north latitude

    # Number of grid cells: Coarsest grid = 1 degree (60 min); resolution 0.00833333 for 30" ; 1' =0.0166667
    clawdata.num_cells[0] = 13 #mx=13, 13/13 =1 ; 13/1560.00062 = 0.00833333 ; 13/780 = 0.0166667
    clawdata.num_cells[1] = 15 #my=15, 15/15 = 1 ; 15/1800.00072 = 0.00833333 ; 15/900 = 0.0166667
    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chkaaaaa'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 48
        clawdata.tfinal = 12*3600.0
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [30., 60., 300., 600.]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True
        

    clawdata.output_format = 'binary32'      # 'ascii', 'binary32', 'binary64'

    clawdata.output_q_components = 'all'   # need all
    clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1


    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.2

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75 #0.8
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'



    # --------------
    # Checkpointing:Flag
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = -2

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif np.abs(clawdata.checkpt_style) == 1:
        # Checkpoint only at tfinal.
        pass

    elif np.abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.  
        #clawdata.checkpt_times = [0.1,0.15]
        clawdata.checkpt_times = 3600.*np.arange(1,16,1)

    elif np.abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # maximum size of patches in each direction (matters in parallel):
    amrdata.max1d = 60     # default

    # max number of refinement levels:
    amrdata.amr_levels_max = 8

    # List of refinement ratios at each level (length at least mxnest-1)
   ###dx = dy = 0.00833333DEG = 30", 6", 2", 1"
    ##dx = dy = 0.016667DEG = 1', 30", 6", 2", 1"
     #dx = dy = 1DEG = 60', 6', 2', 30", 6", 2", 1/3", 1/9"
    amrdata.refinement_ratios_x = [10,3,4,5,3,6,3]
    amrdata.refinement_ratios_y = [10,3,4,5,3,6,3]
    amrdata.refinement_ratios_t = [10,3,4,5,3,6,3]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0 


    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    #rundata.regiondata.regions.append([1, 2, 0., 1e9, -360,360,-90,90])

    flagregions = rundata.flagregiondata.flagregions  # initialized to []
  
    # Levels 1 to 5, dx = dy = 1', 30", 6", 2", 1"

    # The entire domain restricted to level 1 for illustration:
    # (other regions below will force/allow more refinement)
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_domain'
    flagregion.minlevel = 1 #60'
    flagregion.maxlevel = 3 #2' ocean domain
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    # domain plus a bit so kml files look nicer:
    flagregion.spatial_region = [clawdata.lower[0] - 0.1,
                                 clawdata.upper[0] + 0.1,
                                 clawdata.lower[1] - 0.1,
                                 clawdata.upper[1] + 0.1]
    flagregions.append(flagregion)

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_source'
    flagregion.minlevel = 3 # 2'
    flagregion.maxlevel = 4 # 30" deformation region
    flagregion.t1 = 0.
    flagregion.t2 = 1*60 # ends after first minute, for dtopo initiation
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-135., -122., 39., 53.] # CSZ deform extent
    flagregions.append(flagregion)

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_coast_46_51'
    flagregion.minlevel = 3 # 2'
    flagregion.maxlevel = 4 # 30"
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle
    flagregion.spatial_region_file = os.path.abspath(\
                'input_files/RuledRectangle_Coast_46_51.data') # along coast at SJdF entrance 
    flagregions.append(flagregion)  # received from Randy LeVeque

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_northern_water 1'
    flagregion.minlevel = 3 # 2'
    flagregion.maxlevel = 4 # 30"
    flagregion.t1 = 0*60.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-124.01,-122.482,48.79,49.585] 
    flagregions.append(flagregion)
 
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_northern_water 2'
    flagregion.minlevel = 3 # 2'
    flagregion.maxlevel = 4 # 30"
    flagregion.t1 = 0*60.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-125.24,-124.01,49.23,50.125] 
    flagregions.append(flagregion)

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_SJdF'
    flagregion.minlevel = 3 # 2'
    flagregion.maxlevel = 5 # 6"
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-125.12, -123.62, 48.04, 48.72] # SJdF-DEM 2" coverage
    flagregions.append(flagregion)

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_PT_1'
    flagregion.minlevel = 3 # 2'
    flagregion.maxlevel = 5 # 6"
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-123.62,-122.91,47.91,48.79] # PT-DEM 2" coverage
    flagregions.append(flagregion)
 
    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_PT_2'# east of Region_PT_1, closer to study
    flagregion.minlevel = 4 # 30"
    flagregion.maxlevel = 6 # 2"
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-122.91,-122.2,47.91,48.79] # PT-DEM 2" coverage
    flagregions.append(flagregion)

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_PS'
    flagregion.minlevel = 4 # 30"
    flagregion.maxlevel = 6 # 2"
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-123.18,-122.16,47.01,48.19] # PS-DEM 2" coverage
    flagregions.append(flagregion)

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Bainbridge' #1/3 coverage for 1/9" topo extent
    flagregion.minlevel = 6 # 2"
    flagregion.maxlevel = 7 # 1/3"
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 1  # Rectangle
    flagregion.spatial_region = [-122.60,-122.31,47.48,47.75] # PS-DEM 2" coverage
    flagregions.append(flagregion)

    flagregion = FlagRegion(num_dim=2)
    flagregion.name = 'Region_Eagle_Harbor'
    flagregion.minlevel = 8 # fixed at 1/9"
    flagregion.maxlevel = 8 # fixed at 1/9"
    flagregion.t1 = 0.
    flagregion.t2 = 1e9
    flagregion.spatial_region_type = 2  # Ruled Rectangle 
    flagregion.spatial_region_file = os.path.abspath(\
                'input_files/RuledRectangle_fgmax_19.data')
    flagregions.append(flagregion)


    # ---------------
    # Gauges:
    # ---------------
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    #rundata.gaugedata.gauges.append([32412, -86.392, -17.975, 0., 1.e10])

    t_start_gauges = 0.

    rundata.gaugedata.gauges.append([1, -122.5386636, 47.6269394,\
                                       t_start_gauges,  1.e10]) # NE Eagle Harbor
    rundata.gaugedata.gauges.append([2, -122.5347751, 47.6236179,\
                                       t_start_gauges,  1.e10]) # Center channel 1
    rundata.gaugedata.gauges.append([3, -122.5298304, 47.6235448,\
                                       t_start_gauges,  1.e10]) # Center channel 2
    rundata.gaugedata.gauges.append([4, -122.5257788, 47.6193188,\
                                       t_start_gauges,  1.e10]) # Center channel 3
    rundata.gaugedata.gauges.append([5, -122.5193796, 47.6195408,\
                                       t_start_gauges,  1.e10]) # Center channel 4
    rundata.gaugedata.gauges.append([6, -122.5126762, 47.6191303,\
                                       t_start_gauges,  1.e10]) # Center channel 5
    rundata.gaugedata.gauges.append([7, -122.5074046, 47.6200883,\
                                       t_start_gauges,  1.e10]) # Center channel 6
    rundata.gaugedata.gauges.append([8, -122.5039281, 47.6217451,\
                                       t_start_gauges + 1*60,  1.e10]) # Center channel 7
    rundata.gaugedata.gauges.append([9, -122.5000050, 47.6207328,\
                                       t_start_gauges,  1.e10]) # Center channel 8 (DEM merge boundary)
    rundata.gaugedata.gauges.append([10, -122.4971937, 47.6195639,\
                                       t_start_gauges,  1.e10]) # Channel entrance
    rundata.gaugedata.gauges.append([11, -122.5013354, 47.6181532,\
                                       t_start_gauges,  1.e10]) # Offshore Pritchard Park
    rundata.gaugedata.gauges.append([12, -122.4958303, 47.6216316,\
                                       t_start_gauges,  1.e10]) # W of Wing Point Rd NE
    rundata.gaugedata.gauges.append([13, -122.5041835, 47.6260679,\
                                       t_start_gauges,  1.e10]) # S of Hawley Cove Park
    rundata.gaugedata.gauges.append([14, -122.5098302, 47.6223055,\
                                       t_start_gauges,  1.e10]) # BI FT 1
    rundata.gaugedata.gauges.append([15, -122.5103546, 47.6228278,\
                                       t_start_gauges,  1.e10]) # BI FT 2
    rundata.gaugedata.gauges.append([16, -122.5108047, 47.6226046,\
                                       t_start_gauges,  1.e10]) # BI FT 3
    rundata.gaugedata.gauges.append([17, -122.5140212, 47.6211288,\
                                       t_start_gauges,  1.e10]) # Offshore WSF maintenance facility
    rundata.gaugedata.gauges.append([18, -122.5153542, 47.6233859,\
                                       t_start_gauges,  1.e10]) # Offshore Waterfront Trail bridge
    rundata.gaugedata.gauges.append([19, -122.5151386, 47.6226021,\
                                       t_start_gauges,  1.e10]) # W of WSF parking lot
    rundata.gaugedata.gauges.append([20, -122.5176936, 47.6219176,\
                                       t_start_gauges,  1.e10]) # Exotic Aquatic Kayaks dock
    rundata.gaugedata.gauges.append([21, -122.5211719, 47.6211716,\
                                       t_start_gauges,  1.e10]) # Winslow Wharf Marina center
    rundata.gaugedata.gauges.append([22, -122.5214312, 47.62183,\
                                       t_start_gauges,  1.e10]) # Winslow Wharf Marina nearshore
    rundata.gaugedata.gauges.append([23, -122.5237021, 47.6206623,\
                                       t_start_gauges,  1.e10]) # Offshore Williamson Landing Marina
    rundata.gaugedata.gauges.append([24, -122.5108386, 47.6167265,\
                                       t_start_gauges,  1.e10]) # Bainbridge Island Marina center
    rundata.gaugedata.gauges.append([25, -122.5107418, 47.6157705,\
                                       t_start_gauges,  1.e10]) # Bainbridge Island Marina nearshore
    rundata.gaugedata.gauges.append([26, -122.5127642, 47.6168901,\
                                       t_start_gauges,  1.e10]) # Eagle Harbor Marina center
    rundata.gaugedata.gauges.append([27, -122.5120236, 47.6171336,\
                                       t_start_gauges,  1.e10]) # N or Eagle Harbor houseboats
    rundata.gaugedata.gauges.append([28, -122.5132347, 47.6163068,\
                                       t_start_gauges,  1.e10]) # Eagle Harbor Marina nearshore
    rundata.gaugedata.gauges.append([29, -122.5138773, 47.6168746,\
                                       t_start_gauges,  1.e10]) # W of Eagle Harbor
    rundata.gaugedata.gauges.append([30, -122.5190526, 47.6217714,\
                                       t_start_gauges,  1.e10]) # Police boat location
    rundata.gaugedata.gauges.append([31, -122.5172453, 47.6205316,\
                                       t_start_gauges,  1.e10]) # BI City dock
    rundata.gaugedata.gauges.append([32, -122.5198249, 47.6170583,\
                                       t_start_gauges,  1.e10]) # private dock center
    rundata.gaugedata.gauges.append([33, -122.5228421, 47.61776,\
                                       t_start_gauges,  1.e10]) # private dock west
    rundata.gaugedata.gauges.append([34, -122.516201, 47.6164393,\
                                       t_start_gauges,  1.e10]) # private dock east
    rundata.gaugedata.gauges.append([35, -122.4930758, 47.6222302,\
                                       t_start_gauges,  1.e10]) # E of Wing Point


    # --------------------
    # GeoClaw parameters:

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")
       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0 	# XXXX for MLW
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient =.025
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.01

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, fname]
    #topo_data.topofiles.append([2, topo_path])

    topo_path = 'input_files' 
    topofiles = topo_data.topofiles

    topofiles.append([3, topo_path + '/etopo_30sec.tt3']) #1' coverage for model domain
    topofiles.append([3, topo_path + '/topo_2secSJdF.tt3']) #2" Strait Juan de Fuca
    topofiles.append([3, topo_path + '/topo_2secPT.tt3']) # 2" Port Townsend
    topofiles.append([3, topo_path + '/topo_2secPS.tt3']) # 2" Puget Sound
    topofiles.append([3, topo_path + '/topo_19sec.tt3']) # 1/9" covering Bainbridge +

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, fname]
    #dtopo_path = os.path.join(scratch_dir, 'dtopo_usgs100227.tt3')
    #dtopo_data.dtopofiles.append([3,dtopo_path])
    
    #dtopodir = 'data1/dtopo/CSZ/'
    dtopo_data.dtopofiles.append([3, os.path.abspath\
	     ('/data1/dtopo/CSZ/CSZ_L1-extended-pmel.tt3')])
    dtopo_data.dt_max_dtopo = 1.0 # This matters for dynamic dtopos


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    #adjust sea level by dtopo:
    rundata.qinit_data.variable_eta_init = True	# for subsidence
    
    #feature to force dry land some locations below sea level:
    force_dry = ForceDry()
    force_dry.tend = 20*3600.
    force_dry.fname = 'input_files/force_dry_init_19.data' ## can only have one for now.
    rundata.qinit_data.force_dry_list.append(force_dry)

    # == setfixedgrids.data values ==  
    # rundata.fixed_grid_data removed in v5.9.0; instead use fgmax and/or fgout


    # == fgmax_grids.data values ==
    # NEW STYLE STARTING IN v5.7.0

    # set num_fgmax_val = 1 to save only max depth (also arrival info based on h)
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 5  # Save depth, arrival time 

    fgmax_grids = rundata.fgmax_data.fgmax_grids # define a shorthand name pointer to the list we need to set:

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.

    fg = fgmax_tools.FGmaxGrid() # set attributes of first fgmax grid
    fg.point_style = 4       # scattered points
    fg.npts = 0
    fg.xy_fname = os.path.abspath('input_files/fgmax_pts_eagle_harbor_19.data')
    fg.min_level_check = amrdata.amr_levels_max # which levels to monitor max on 
    fg.tstart_max = 0              # when to start monitoring max values
    fg.tend_max = 1.e10               # when to stop monitoring max values

    fg.dt_check = 3.                 # how often to update max values
    fg.arrival_tol = 0.025     # tolerance for arrival (eta) # OLD way
    fg.arrival_tols_h = [0.1016, -0.1016, -0.3048, 0.3048, 0.9144] # NEW h tols tfirstpos, neg, draw, advis, warn

    fgmax_grids.append(fg)  # append to the list we need to set


    # == fgout_grids.data values ==
    # NEW IN v5.9.0
    # Set rundata.fgout_data.fgout_grids to be a list of
    # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    #Grid around Eagle Harbor with 1/9" resolution
    fgout_dx = 1./(9*3600) #target resolution
    #fgout_dx = 1./1800 # 2" resolution test
    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1 #fgout area 1, can copy and paste to have multiple regions at different resolutions
    fgout.point_style = 2       # will specify a 2d grid of points
    fgout.output_format = 'binary32'  # 4-byte, float32
    fgout.x1 = -122.545  # specify edges (fgout pts will be cell centers) edges of grid resolution patch
    fgout.x2 = -122.4870
    fgout.y1 = 47.614
    fgout.y2 = 47.631
    fgout.nx = int(round((fgout.x2 - fgout.x1) / fgout_dx))
    fgout.ny = int(round((fgout.y2 - fgout.y1) / fgout_dx))
    fgout.tstart = 0.5*3600
    fgout.tend = 10.5*3600
    fgout.nout = 3601 # number of outputs in start and end time. 
    fgout_grids.append(fgout)    # written to fgout_grids.data


    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = True       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py
    
    return rundata
    
    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()

    kmltools.make_input_data_kmls(rundata)

