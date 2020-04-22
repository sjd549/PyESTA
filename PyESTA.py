#!/usr/bin/env python

#################################
#		Point of Contact		#
#								#
#	   Dr. Scott J. Doyle		#
#	  Scott.Doyle@Physics.org	#
#	  University of Seville		#
#  Plasma Tech & Fusion Science #
#  National Acceleration Centre #
#	     Seville, Spain		    #
#								#
#################################
#            'PyESTA'           #
# Python FIESTA Wrap & Analysis #
#################################


#====================================================================#
				 #PROGRAM FLAGS AND MODULE IMPORTS#
#====================================================================#

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--first", action="store_true", dest="install", default=False, help="Install prompt for required python modules")
(options, args) = parser.parse_args()

if 'True' in str(options): 
	import os, sys
	import os.path

	print ''
	print 'First time use requires installation of additional python modules'
	print 'Please type your password when prompted to allow installation:'
	print ''
	try:
		os.system('sudo apt-get install python-pip')
		os.system('sudo apt-get install python-matplotlib')
		os.system('sudo apt-get install python-numpy')
		os.system('sudo apt-get install python-scipy')
		os.system('sudo apt-get install ffmpeg')
		os.system('pip install tqdm')
	except:
		print ''
		print 'Error installing required packages'
		print 'Please attempt manual installation'
		print ''
	#endtry
	print ''
	print ''
#endif

#==============#

#Import core modules
import matplotlib.cm as cm
import numpy as np
import scipy as sp
import math as m
import subprocess
import os, sys
import os.path
import time

#Enforce matplotlib to avoid instancing undisplayed windows
#matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
import matplotlib
matplotlib.use('Agg')

#Import additional modules
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import savgol_filter
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from matplotlib import ticker
from scipy import ndimage
from tqdm import tqdm
from pylab import *

#====================================================================#
#====================================================================#





#====================================================================#
				  		 #LOW LEVEL INPUTS#
#====================================================================#

#Warning suppressions
np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
#Fix "can't invoke "event" command: application has been destroyed" error with PROES images
#Fix "Exception KeyError: KeyError(<weakref at 0x7fc8723ca940; to 'tqdm' at 0x7fc85cd23910>,)" error

#Define FIESTA source trunk location
FIESTATrunk = "~/Postdoc Seville/FIESTA/Source Code/FIESTA_V8.8"

#Various debug and streamlining options.
IDebugMode = False			#Produces debug outputs for relevent diagnostics.

#Set default parallelisation options
IParallel = False			#Allow multiple simultanious FIESTA simulations
NumThreads = 2				#Number of threads avaliable to each FIESTA simulation
MaxNumThreads = 4			#Maximum number of threads avaliable for all FIESTA simulations
ConvDelay = 10	#[s]		#Delay between convergence checks for race condition

#====================================================================#
#====================================================================#





#====================================================================#
					  	#PyESTA NAMELIST FILE#
#====================================================================#

#Commonly Used Settings:
#### SMART_V3Phase1 ####
#betaT    = 0.1000		#[%]
#betaP    = 0.4578		#[%]
#li2      = 1			#[]
#Ip       = 30e3		#[kA]
#Irod     = 4.7750e5	#[kA]
#B0       = 0.10		#[T]  (0.09 [T] in vacuo)
#TauPulse = 0.020		#[s]  (Timestep for pulse)

#### SMART_V3Phase2 ####
#betaT    = 0.1000		#[%]
#betaP    = 0.4578		#[%]
#li2      = 1			#[]
#Ip       = 100e3		#[kA]
#Irod     = 4.7750e5	#[kA]
#B0       = 0.20		#[T]
#TauPulse = 0.100		#[s]  (Timestep for pulse)

#### SMART_V3Phase3 ####
#betaT    = 0.1000		#[%]
#betaP    = 0.4578		#[%]
#li2      = 1			#[]
#Ip       = 500e3		#[kA]
#Irod     = 4.7750e5	#[kA]
#B0       = 1.00		#[T]  (0.09 [T] in vacuo)
#TauPulse = 0.500		#[s]  (Timestep for pulse)

#######################  DEFINE VESSEL GEOMETRY  #######################

#Define global scaling factors
RScaleVessel=1.00    #Scale all radial vessel dimensions (Relative to 2.0m)
ZScaleVessel=0.80    #Scale all axial vessel dimensions (Relative to 2.0m)
RScaleCoil=1.00      #Scale all radial coil positions (Relative to 2.0m)
ZScaleCoil=0.80      #Scale all axial coil positions (Relative to 2.0m)

#Define Vessel Outer Geometry
VesselRInnerPoint=0.15*RScaleVessel	# R min position [m]
VesselROuterPoint=0.8*RScaleVessel	# R max position [m]
VesselZMinPoint=-1.0*ZScaleVessel	# Z min position [m]
VesselZMaxPoint=1.0*ZScaleVessel	# Z max position [m]

#Define Solenoid Geometry and Parameters
nSol=210;     # Solenoid Central Windings  [-]	   	#nSol=210,800
RSol=0.13;    # R position of the solenoid [m]      #RSol=0.13,0.09
ZMinSol=-1.0*ZScaleVessel # Min Z position
ZMaxSol=1.0*ZScaleVessel  # Max Z position

#Number of Radial (R) and axial (Z) coil windings
nZDiv1=6
nRDiv1=4
nZDiv2=6
nRDiv2=4
nZPF2=6
nRPF2=4
nZPF3=6
nRPF3=4

#Define coil turn dimensions to enable cross-section calculation
width_PF=0.042  # Width of a turn (m)
height_PF=0.035 # Height of a turn (m)

#Define central location of coil sets
#R_PF1=0.9*RScaleCoil  #R position of PF1 (m)
#Z_PF1=0.3*ZScaleCoil  #Z position of PF1 (m)
R_PF2=0.9*RScaleCoil   #R position of PF2 (m)
Z_PF2=0.5*ZScaleCoil   #Z Position of PF2 (m)
R_PF3=0.9*RScaleCoil   #R Position of PF3 (m)
Z_PF3=0.8*ZScaleCoil   #Z Position of PF3 (m)
R_Div1=0.25*RScaleCoil #R Position of Div1 (m)
Z_Div1=1.05*ZScaleCoil #Z Position of Div1 (m)
R_Div2=0.55*RScaleCoil #R Position of Div2 (m)
Z_Div2=1.05*ZScaleCoil #Z Position of Div2 (m)


#######################  DEFINE INITIAL PARAMETERS  #######################

#Define any required constants
mu0 = 1.2566e-06 # Magnetic Moment      [I/m^2]

#Define Operating Conditions
Te = 250         # Electron Temperature [eV]
Ti = Te*0.1      # Ion Temperature      [eV]
BT = 0.1         # Toroidal B-Field     [T] (Defined at Rgeo)
Ip = 30e3        # Plasma current       [A]
RGeo = 0.450     # Geometrical Radius   [m]
ZGeo = 0.000     # Geometrical Axis     [m]
RSep = 0.700     # Separatrix Radius    [m]
a = RSep-RGeo    # Minor Radius         [m] (~0.25)
A = RGeo/a       # Aspect ratio         [-] (~1.8)
Kappa = 1.8      # Elongation           [-]
delta = 0.20     # Triangularity        [-] (~0.2)
li2 = 1          # Standard Value?      [-]
#q_cyl = 2.821   # Safety Factor?       [-]
#betaN = 3.529   # Normalised Beta      [%] (Obtained via VEST Excel - (2X TOO HIGH)

#Define efit Equilibrium Operating Conditions
RGeo_efit = 0.440					# Geometrical Radius	[m] (Default 0.44)
ZGeo_efit = 0.000					# Geometrical Axis		[m] (Default 0.00)
rGeo_efit = 0.238					# Minor Radius	        [m] (Default 0.44/1.85)
Aspect_efit = RGeo_efit/rGeo_efit;  # Aspect Ratio          [-] (Default 1.85)
Kappa_efit = 1.80					# Elongation			[-] (Default 1.80)
delta_efit = 0.20					# Triangularity			[-] (Default 0.20)
efit_Geometry_Init = [RGeo_efit, ZGeo_efit, rGeo_efit, Kappa_efit, delta_efit]

#Compute Further Operating Conditions
Gr_Limit = 1e20*(Ip*1e-6/(pi*a**2*Kappa))  # Greenwald Limit          [m-3]
Gr_Frac = 0.15                             # Greenwald Fraction       [-]
ne = Gr_Limit*Gr_Frac                      # Electron Density         [m-3]  ~3E19
Irod = BT*2*pi*RGeo/mu0                    # Central Rod Current      [A]
S = sqrt( (1.0+Kappa**2)/2.0 )             # Shaping factor           [-]
#deltaUp = (RGe-Rup)/a                     # Upper-Triangularity      [-]
#deltaLo = (RGe-Rlo)/a                     # Lower-Triangularity      [-]
#delta = (deltaUp+deltaLo)/2.0             # Triangularity            [-]
#betaN = (betaT*BT*a)/(Ip*1e-6*mu0)        # Normalised Beta          [%] 
#betaT = (betaN/a*(Ip*1e-6))/BT            # Beta toroidal            [%]
betaP = 3/2*ne*(Te+Ti)/(mu0*Ip/(2*pi*a))**2*2*mu0*1.6e-19*Kappa  # Beta Poloidal  [%]
BZ = -mu0*Ip/(4*pi*RGeo)*(log(8*A)+betaP+0.5*li2-3/2)      		 # Vertical field [T]

#Coil density, temperature and resistivity
coil_density = 1						# Relative Coil Density      [Arb]
coil_temp = 293.0						# Initial Coil Temperature   [K]

#Gas species analouge - H=1, He=2, Ar=11.85 (for Te < 280eV)
#https://www.webelements.com/argon/atoms.html
Z_eff=1.0								# Effective Nuclear Charge   [e-]

#Null field region radius, specifies Sensor_btheta radius
a_eff=0.10;								# Null field region radius	 [m]


###################  DEFINE SOL RAMP & COIL CURRENTS  ###################

#Notes:
#Negative coil currents attract the plasma, positive repel the plasma
#Symmetric Solenoid PrePulse and Equil currents aid power supply stability

#Solenoid coil currents [kA]		#Phase1		#Phase1-Old		#Phase2
I_Sol_Null=+775						#+0775;		#+900;			#+2200
I_Sol_MidRamp='Linear'				#Dynamic    #Dynamic		#Dynamic
I_Sol_EndRamp=-I_Sol_Null			#Dynamic    #Dynamic		#Dynamic
I_Sol_Equil=-900;					#Default -I_Sol_Null		#Efit_Equil

#PF coil currents (At Equilibrium, time(4,5,6))
I_PF1_Equil=-500;					#-500;		#-390;			#-1100
I_PF2_Equil=-500;					#-500;		#-385;			#-1700
I_Div1_Equil=+000;					#+000;						#+0000
I_Div2_Equil=+900;					#+900;						#+3300

#Define number of time-steps (vertices) in the current waveforms
nTime = 7			# Coil Waveform Timesteps	[-]
TauB = 0.020		# Buffer Timescale     		[s] Determines tstep for Ip plot
TauR = 0.050		# Ramp Timescale       		[s]
TauP = 0.020		# Pulse Timescale      		[s]
#Time   [Init      PrePulse  InitRampDown  MidRampDown  EndRampDown  MidEquil     Terminate         ]
time =  [-4*TauB,  -2*TauB,  0.0,          TauR/2.0,    TauR,        TauR+TauP,   TauR+TauP+(2*TauB)]

#Definition of time intervals:
#time(1)--> All coils and Sol initiate at zero current          Init
#time(2)--> All coils initiate null-field configuration         PrePulse
#time(3)--> All coils maintain null-field configuration         InitRampDown
#time(4)--> Sol ramps down, PF/Div coils init equilibrium       MidRampDown - InitEquil
#time(5)--> Sol completes ramp down, maintain PF/Div coils      EndRampDown - MidEquil
#time(6)--> All coils maintain equilibrium configuration        EndEquil
#time(7)--> All coils and Sol terminate at zero current         Terminate
#######
#time(3)-->time(5) lasts timescale TauR (Solenoid Ramp-Down TimeScale)
#time(5)-->time(6) lasts timescale TauP (Pulse/Discharge Timescale)
#######

#Define number of time-steps (vertices) in the current waveforms
nTime = 7;      #[Steps]
#Time   [Init      PrePulse  InitRampDown  MidRampDown  EndRampDown  MidEquil     Terminate         ]
time =  [-4*TauR,  -2*TauR,  0.0,          TauR/2.0,    TauR,        TauR+TauP,   TauR+TauP+(2*TauR)]

#Construct Sol, PF/Div coil current waveforms vertices
#Time   	     [1, 2,           3,          4,             5,             6,             7];
ISol_Waveform =  [0,  I_Sol_Null, I_Sol_Null, I_Sol_MidRamp, I_Sol_EndRamp, I_Sol_EndRamp, 0]
IPF1_Waveform =  [0,  NaN,        NaN,        NaN,           I_PF1_Equil,   I_PF1_Equil,   0]
IPF2_Waveform =  [0,  NaN,        NaN,        NaN,           I_PF2_Equil,   I_PF2_Equil,   0]
IDiv1_Waveform = [0,  I_Sol_Null, I_Sol_Null, I_Sol_MidRamp, I_Sol_Equil,   I_Sol_Equil,   0]
IDiv2_Waveform = [0,  NaN,        NaN,        NaN,           I_Div2_Equil,  I_Div2_Equil,  0]
#####
CoilWaveforms = [ISol_Waveform, IPF1_Waveform, IPF2_Waveform, IDiv1_Waveform, IDiv2_Waveform]


####################  DEFINE DIAGNOSTIC PARAMETERS  #######################

#Stability diagnostic perturbations (must be smaller than initial variable!)
deltaRGeo = 0.00;	# Small radial perturbation         [m]
deltaZGeo = 0.00;	# Small axial perturbation          [m]
deltaAspect = 0.00;	# Small aspect ratio perturbation   [-]
deltaKappa = 0.00;	# Small elongation perturbation     [-]
deltadelta = 0.00;	# Small triangiularity perturbation [-]

#====================================================================#
#====================================================================#





#====================================================================#
					  #SWITCHBOARD AND SETTINGS#
#====================================================================#

#Vessel and Coil Geometry Ranges
#'R_Div1' [x/100.0 for x in range(12,31,2)]
#'R_Div2' [x/100.0 for x in range(30,61,2)]
#'Z_PF1' [x/100.0 for x in range(00,81,2)]
#'Z_PF2' [x/100.0 for x in range(00,81,2)]

#Common General Parameter Ranges
#'Ip' [x for x in range(20000,40000,2000)]
#'Gr_Frac' [x/100.0 for x in range(9,26,3)]
#'Te' [x for x in range(50,251,50)]
#'Z_eff' [1.0,2.0,11.85]	#H, H2, Ar8+  (Old Resistivity ~34)
#'BPolEarth' [x/100000.0 for x in range(1,11,2.5)]

#Common Coil Current Ranges
#'I_Sol_Null' [x for x in range(500,1001,100)]
#'I_Sol_EndRamp' [-x for x in range(700,1001,25)]
#'I_Sol_Equil' [-x for x in range(200,851,50)]
#'I_PF1' [-x for x in range(450,651,25)]
#'I_PF2' [-x for x in range(450,651,25)]
#'I_Div1' [x for x in range(000,1501,100)]
#'I_Div2' [x for x in range(800,1151,50)]

#Coil Pulse Ranges
#'TauR' [x/1000.0 for x in range(20,61,5)]
#'TauP' [x/1000.0 for x in range(20,41,5)]

#Equilibrium Stability Ranges
#'deltaZGeo' #[x/100.0 for x in range(0,21,2)]
#'deltaRGeo' #[x/100.0 for x in range(0,21,2)]

#Common Efit Geometry Ranges
#'RGeo_efit' [x/100.0 for x in range(40,50,2)]
#'ZGeo_efit' [x/100.0 for x in range(0,0,1)]
#'a_efit' [x/10.0 for x in range(10,20,1)]
#'Kappa_efit' [x/100.0 for x in range(170,271,10)]
#'delta_efit' [-x/10.0 for x in range(0,31,2)], [x/10.0 for x in range(0,16,1)]

########################################

#Define FIESTA namelist and project directory names
FIESTAName = 'SMART_SJD.m'			#Define name of FIESTA script
ProjectName = 'S1-000002'			#Define Global Project Name (Baseline Equilibrium)
SeriesName = 'auto'					#Parameter scan series name ('auto' for automatic)

#Define simulation name structure
SimNameList = ['delta_efit','Kappa_efit','I_Sol_Null','I_PF1_Equil','I_PF2_Equil', 'I_Div1_Equil','I_Div2_Equil']

#Define if simulations are to be run
IAutorun = False			#Run requested simulation series
IParallel = False		#Enable mutli-simulations in parallel
IVerbose = True			#Verbose terminal output - not compatable with IParallel

#Define equilibrium calculation method
IEquilMethod = 'efit'					#Define equil method: 'standard','efit','feedback'
IefitCoils = ['PF1','PF2']				#Define coils for which efit, feedback is applied

#Define paramters to be varied and ranges to be varied over
ParameterVaried = 'I_Sol_EndRamp'		#Define parameter to vary - Required for diagnostics
ParameterRange = [-x for x in range(975,1001,25)]	#Define paramter range to vary over

#Define which diagnostics are to be performed
savefig_EquilStability = True		#Plots current trends in response to perturbed equilibria
savefig_EfitEquilTrends = True		#Plots efit equilibrium geometry trends from Param(equil)
savefig_UserEquilTrends = False		#Plots user defined equilibrium trends from Param(equil)
#savefig_EquilSeperatrix = False	#Plots seperatrix extrema [Rmin,Rmax,Zmin,ZMax] trends
#savefig_EquilMidplane = False		#Plots 2D Radial slice at Z=0 trends
#savefig_EquilXpoint = False		#Plots X-point location (R,Z) trends

savefig_CoilCurrentTraces = True	#Plots PF coil current timetraces for each simulation
savefig_CoilCurrentTrends = True	#Plots trends in PF coil currents over all simulations

savefig_ConnectionLength = True		#Plots trends in average connection length over all simulations
savefig_PaschenCurves = True		#Plots Paschen curves for each simulation using Lc

savefig_PlasmaCurrent = True		#Plots plasma current trends over all simulations
savefig_EddyCurrent = True			#Plots total vessel eddy current trends over all simulations


#Image overrides and tweaks
Image_TrendAxisOverride=''			#Force trend figures to use different variable

#====================================================================#
#====================================================================#





#====================================================================#
#====================================================================#

#TO DO
#IMMEDIATE FIXES
#Enable auto-detection of output folders so user doesn't have to change parameter variable

#CORE FUNCTIONALITY
#Unify PyESTA namelist inputs with .m namelist inputs 		- Ideally in external namelist file
#Rename all FIESTA output text files in unified format      - Enable CSV or Row-Wise data storing
#Save all PyESTA output data in seperate output folder      - Enable CSV or Row-Wise data storing
#Add capability to iterate towards fixed equilibrium conditions - i.e. iterate on single variable
#Add ability to change multiple variables per run
#Add ability to use multiple cores, including safety		- NEEDS TESTING!!!

#DIAGNOSTICS
#Add equilibrium Rmin,Rmax, Zmin,ZMax diagnostic showing extrema of seperatrix
#Add equilibrium 'PROES-like' diagnostics - plot Radial slice at Z=0 with parameter range
#Add equilibrium X point diagnostic, showing X-point location (R,Z) trends

#ERROR HANDLING
#Add ability to safely-eject from matlab if convergence fails 				- IMPORTANT!!!
#Add ability to produce 'nan' data files if one-or-more simulations fail	- IMPORTANT!!!
#Add general error messages to aid in debugging as the program grows larger
#DEAL WITH CORRUPTED DOUBLE LINKED LISTS - NOTE THAT | OR \ APPEARS TO SKIP TO NEXT SIMULATION!!!

#====================================================================#
#====================================================================#





#====================================================================#
				   #DEFINE COMMONLY USED FUNCTIONS#
#====================================================================#

#Takes global inputs from switchboard, returns nothing
#Alters global image options, run before any diagnostics
#Attempts to revert matplotlib changes made in 2.0 onwards.
#See: https://matplotlib.org/users/dflt_style_changes.html
def Matplotlib_GlobalOptions():

#	mpl.style.use('classic')								#Resets to classic 1.x.x format
	
	#Image options			
	mpl.rcParams['figure.figsize'] = [10.0,10.0]			#Sets default figure size
	mpl.rcParams['figure.dpi'] = 200						#Sets viewing dpi
	mpl.rcParams['savefig.dpi'] = 100						#Sets saved dpi
	mpl.rcParams['image.interpolation'] = 'bilinear'		#Applies bilinear image 'smoothing'
	mpl.rcParams['image.resample'] = True					#Resamples data before colourmapping
	mpl.rcParams['image.cmap'] = 'plasma'					#Select global colourmap 
	#'jet','plasma','gnuplot'

	#Axis options
	mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'	#View limits coencide with axis ticks
	mpl.rcParams['axes.xmargin'] = 0						#Set default x-axis padding
	mpl.rcParams['axes.ymargin'] = 0						#Set default y-axis padding
	mpl.rcParams['errorbar.capsize'] = 3					#Set error bar end cap width
	mpl.rcParams['font.size'] = 12							#Set global fontsize
	mpl.rcParams['legend.fontsize'] = 'large'				#Set legend fontsize
	mpl.rcParams['figure.titlesize'] = 'medium'				#Set title fontsize

	#Line and Colour options
#	from cycler import cycler								#See below
#	mpl.rcParams['axes.prop_cycle']=cycler(color='bgrcmyk')	#Set default colour names
	mpl.rcParams['lines.linewidth'] = 1.0					#Set Default linewidth

	#Maths and Font options
	mpl.rcParams['mathtext.fontset'] = 'cm'					#Sets 'Latex-like' maths font
	mpl.rcParams['mathtext.rm'] = 'serif'					#Sets default string font

	return()
#enddef
Matplotlib_GlobalOptions()	#MUST BE RUN BEFORE ANY DIAGNOSTICS!!!!

#=========================#

#Checks for any running process that contain given name ProcessName.
#Takes process name string input (same as process name in top or htop)
#Returns boolian, true if process is found, false if not
#By default returns false to allow for softcrash.
#Example: Bool = CheckIfProcessRunning('MATLAB')
def CheckIfProcessRunning(QueriedProcess,Bool=False):

	#Initiate required lists
	ProcessIDList,ProcessNameList = list(),list()

	#Call for all processes in 'top, htop' format and check if exists
	ProcessCall = ['ps','-A']
	Processes = subprocess.check_output(ProcessCall)
	Processes = Processes.split('\n')

	#Split processes into ID and name, and compile lists for later use
	for i in range(0,len(Processes)):
		#Check if ProcessID is first or 2nd entry and save
		try: ProcessIDList.append( float(Processes[i].split(' ')[0]) )
		except:
			try: ProcessIDList.append( Processes[i].split(' ')[1] )
			except: ProcessIDList.append( np.nan )
		#endtry

		#ProcessName is final entry in 'top' format
		ProcessNameList.append( Processes[i].split(' ')[-1] )
	#endfor

	#Check if process in ProcessNamelist
	for i in range(0,len(ProcessNameList)):
		if ProcessNameList[i] == QueriedProcess:
			Bool = True
			break
		#endif
	#endfor

	return(Bool)
#enddef

#=========================#


#Constructs and executes matlab command to run FIESTA
#Takes FIESTA .m file name and returns nothing
#Example: RunFIESTA('FIESTA.m')
def RunFIESTA(FIESTAName,Verbose=False,Parallel=False):

	#Print parallel verbosity warning
	if Verbose == True and Parallel == True:
		print ''
		print 'Warning! Parallel operation not compatable with verbose output'
		print '                 Setting Parallel = False                     '
		print ''
		Parallel = False
	#endif

	#Construct terminal command to run requested version of FIESTA
	#Example: matlab -nodisplay -nosplash -nodesktop -r "run('/path/to/FIESTA_Script');exit;"
	FIESTA_RootDir = os.getcwd()+'/'+FIESTAName
	FIESTA_Splash = '-nodisplay -nosplash -nodesktop -r '
	FIESTA_RunCMD = '\"run(\''+FIESTA_RootDir+'\');exit;\"'
	if Verbose == True or IDebugMode == True:	 FIESTA_Output = ''
	elif Verbose == False and Parallel == False: FIESTA_Output = ' > Conv.txt'
	elif Verbose == False and Parallel == True:  FIESTA_Output = ' > Conv.txt &'
	#####
	ExecuteFIESTA = 'matlab '+FIESTA_Splash+FIESTA_RunCMD+FIESTA_Output

	#Execute FIESTA script in terminal
	os.system( ExecuteFIESTA )

	return()
#enddef

#=========================#

#Takes namelist directory and namelist variable
#Locates namelist entry for variable and returns value
#Example: Value,Entry = AlterNamelistVariable(FIESTAName,ParameterVaried)
def FindNamelistVariable(Namelist_Dir,ParameterVaried):

	#Open namelist file and identify the requested variable name, line index and init value
	Namelist = open(Namelist_Dir).readlines()
	NamelistEntry = filter(lambda x:ParameterVaried in x, Namelist)[0]
	NamelistIndex = Namelist.index(NamelistEntry)
	try: NamelistValue = float(NamelistEntry.partition(';')[0].strip(ParameterVaried+' \t\n\r,='))
	except: NamelistValue = NamelistEntry.partition(';')[0].strip(ParameterVaried+' \t\n\r,=')

	return(NamelistValue,NamelistEntry)
#enddef

#=========================#

#Takes namelist directory and namelist variable and value
#Locates namelist entry for variable and alters to new value
#Returns modified namelist entry for sanity checking purposes
#WARNING -- DOESN'T APPEAR TO WORK FOR STATEMENTS INSIDE INDENTED LOOPS!!!
#Example: Init,Entry = AlterNamelistVariable(FIESTAName,ParameterVaried,VariableValue)
def AlterNamelistVariable(Namelist_Dir,ParameterVaried,VariableValue):

	#Open namelist file and identify the requested variable name, line index and init value
	Namelist = open(Namelist_Dir).readlines()
	NamelistEntry = filter(lambda x:ParameterVaried in x, Namelist)[0]
	NamelistIndex = Namelist.index(NamelistEntry)
	try: NamelistValue = float(NamelistEntry.partition(';')[0].strip(ParameterVaried+' \t\n\r,='))
	except: NamelistValue = NamelistEntry.partition(';')[0].strip(ParameterVaried+' \t\n\r,=')
	
	#Seperate variable string into 5 sub-strings of order: 'Variable,Sep2,Value,Sep,Comment'
	#Assumes single value input terminated by a semi-colon - allows for and retains comments
	Input, Sep, Comment = NamelistEntry.partition(';')
	Variable, Sep2, InitValue = Input.partition('=')
	AlteredNamelistValue = str(VariableValue)
	#Reconstruct the altered namelist entry with updated init value
	AlteredNamelistEntry = Variable+Sep2+AlteredNamelistValue+Sep+Comment

	#Replace namelist entry with altered namelist entry
	Namelist[NamelistIndex] = AlteredNamelistEntry

	#Write reconfigured namelist back into Namelist_Dir
	with open(Namelist_Dir, 'w') as file:
		file.writelines( Namelist )
	#endwith

	return(AlteredNamelistEntry)
#enddef

#=========================#

#Returns sub-folder directories within supplied simulation series folder
#Takes simulation series local directory name string
#Returns sub-folder names within series directory
#Can supply directories relative to cwd() or relative to root ('/home/...')
#Example: SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)[1]
def ExtractSubDirs(SeriesDirString,Root=True):

	#Obtain simulation series folder directories and create list for contents
	try:
		#List local directory and search for any folders associated with project
		Directorylist = os.listdir( os.getcwd() )
		for i in range(0,len(Directorylist)):
			if ProjectName in Directorylist[i]:
				#ONLY TAKES FIRST PROJECT FOLDER - NEED TO UPDATE FOR MULTI-FOLDER
				SeriesDir = Directorylist[i]
				break
			#endif
		#endif
		SimulationDirsRaw = os.listdir(SeriesDir)
		SimulationDirsCleaned = list()
	except:
		#Inform user if no simulation folders found
		print '---------------------------------'
		print 'No FIESTA output folders detected'
		print '---------------------------------'
		exit()
	#endtry
	
	#Remove any non-folder directories in SeriesDirsRaw and correct bash 'grammar'
	for i in range(0,len(SimulationDirsRaw)):

		#Define simulation series directories from root or relative to local directory
		AlwaysRoot = os.getcwd()+'/'+SeriesDir
		if Root == True: RootDir = os.getcwd()+'/'+SeriesDir
		else: RootDir = SeriesDir

		#Remove any non-folder entries - assume all folders are simulation directories
		if os.path.isdir(AlwaysRoot+'/'+SimulationDirsRaw[i]) == False:
			Directory_Is_Not_A_Folder=1.0
		elif os.path.isdir(AlwaysRoot+'/'+SimulationDirsRaw[i]) == True:
			SimulationDirsCleaned.append( ''+RootDir+'/'+SimulationDirsRaw[i]+'' )
		#endif
	#endfor

	#Maintain alphanumerical foldername ordering
	SimulationDirsRaw = sorted(SimulationDirsRaw)
	SimulationDirsCleaned = sorted(SimulationDirsCleaned)

	#Return all folder directories in requested simulation series
	return(SimulationDirsCleaned)
#enddef

#=========================#

def ReadDataFromFile(Filename,Dimension='2D',Orientation='Vertical'):
	OutputData = list()

	#If data is saved 'Row-wise', use default readin routine.
	if Orientation == 'Horizontal':
		#Determine dimensionality of profile.
		if Dimension == '2D':
			#Read in 2D data from ASCII formatted file.
			datafile = open(Filename)
			RawData = datafile.readlines()
			for m in range(0,len(RawData)):
				Row = RawData[m].split()
				for n in range(0,len(Row)):
					try: Row[n] = float(Row[n])
					except: Row[n] = str(Row[n])
				#endfor
				OutputData.append(Row)
			#endfor
		#endif

	#=====#

	#If data is saved 'column-wise', transpose the arrays to correct.
	elif Orientation == 'Vertical':
		#Determine dimensionality of profile.
		if Dimension == '2D':
			#Read in 2D data from ASCII formatted file.
			datafile = open(Filename)
			RawData = datafile.readlines()
			for m in range(0,len(RawData)):
				Row = RawData[m].split()

				#Determine how many rows of data exist.
				if len(OutputData) == 0:
					for i in range(0,len(Row)): OutputData.append(list())
				#endif
				for j in range(0,len(Row)):
					try: Row[j] = float(Row[j])
					except: Row[j] = str(Row[j])
				#endfor
				for k in range(0,len(OutputData)):
					OutputData[k].append(Row[k])
				#endfor
			#endfor
		#endif
	#endif

	#Orientation doesn't exist if 0D (scalar).
	elif Dimension == '0D':
		#Read in 0D data from ASCII formatted file.
		datafile = open(Filename)
		Row = datafile.readline().split()
		for m in range(0,len(Row)):
			OutputData.append(float(Row[m]))
		#endfor
	#endif

	return(OutputData)
#enddef

#=========================#

def ExtractFIESTAData(SeriesSubDirs,DataFileName,Dimension='2D',Orientation='Vertical',Reorder=True):

	#Create any required arrays for data storage and record HomeDir for navigation
	GlobalDataArrays,ReorderedDataArrays = list(),list()
	HomeDir = os.getcwd()

	#For all simulation directories in the requested simulation series
	for i in range(0,len(SeriesSubDirs)):
		#cd into the relevent directory and extract the data
		os.chdir(SeriesSubDirs[i]+'/RawData/')
		#GlobalDataArray organized as [Folder][Variable][Value]
		GlobalDataArrays.append(ReadDataFromFile(DataFileName,Dimension,Orientation))
	#endfor
	#cd back into PyESTA directory for continuity
	os.chdir(HomeDir)

	#Reformat GlobalDataArray to enable easy splitting of column-wise variables
	if Reorder == True:
		#Append a list for each seperate variable in the original data format
		for i in range(0,len(GlobalDataArrays[0])): ReorderedDataArrays.append(list())
		#GlobalDataArray organized as [Folder][Variable][Value]
		#ReorderedDataArrays organised as [Variable][Folder][Value]
		for i in range(0,len(ReorderedDataArrays)):
			for j in range(0,len(GlobalDataArrays)): 
				#If data array contains only one float, append as float:
				if len(GlobalDataArrays[j][i]) == 1:
					ReorderedDataArrays[i].append( GlobalDataArrays[j][i][0] )
				#if data is list, append as list:
				elif len(GlobalDataArrays[j][i]) > 1:
					ReorderedDataArrays[i].append( GlobalDataArrays[j][i] )
				#endif
			#endfor
		#endfor
		OutputDataArrays = ReorderedDataArrays
	else:
		OutputDataArrays = GlobalDataArrays
	#endif

	#Return ordered data arrays
	return(OutputDataArrays)
#enddef

#=========================#

#Creates a new folder if one does not already exist.
#Takes destination dir and namestring, returns new directory.
def CreateNewFolder(Dir,DirString):
	try:
		NewFolderDir = Dir+DirString+'/'
		os.mkdir(NewFolderDir, 0755);
	except:
		a = 1
	#endtry
	return(NewFolderDir)
#enddef

#=========================#

#Converts variable name strings into concatenated string with variable values
#Takes 1D array of variable name strings - must exist in namelist file!
#Returns 0D string of concatenated values of form: 'Var#Value '
#Example: SimulationName = CreateSimName(SimNameList)
def CreateSimName(SimNameList,VariedParameter='NaN',ParameterValue='NaN'):
	#Define empty name string
	SimulationNameString = ''

	#For each variable in SimNameList, convert to string and append value
	for i in range(0,len(SimNameList)):

		#Shorten Variable strings as much as possible by removing excess underscores
		if len(SimNameList[i].split('_')) >= 2:
			TrimmedSimName = SimNameList[i].split('_')[0]+SimNameList[i].split('_')[1]
		else:
			TrimmedSimName = SimNameList[i]
		#endif

		#Check if named parameter has been varied and use appropriate value
		if SimNameList[i] == VariedParameter:
			ParameterString = TrimmedSimName+'#'+str(ParameterValue)+' '
		else:
			ParameterString = TrimmedSimName+'#'+str(eval(SimNameList[i]))+' '
		#endif
		SimulationNameString += ParameterString
	#endfor

	#Remove final whitespace in simulation name
	SimulationNameString = SimulationNameString[:-1]

	return(SimulationNameString)
#enddef

#=========================#

#Creates a trendaxis from simulation folder names
#Takes directories of all folders in the simulation series folder
#Returns a 1D array of floating values based on the varied parameter
#Example: TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried)
def CreateTrendAxis(SimulationNames,VariableString,Image_TrendAxisOverride=''):

	#Create required list to store output
	TrendAxis = list()

	#For all simulation names
	for i in range(0,len(SimulationNames)): 
		#Split each directory folder name into substrings and identify varied parameter
		SimulationNames[i] = SimulationNames[i].split('/')[-1]	#Remove any directories
		SplitSimName = SimulationNames[i].split(' ')			#Split simulation parameters

		#Find trend variable and extract value - check override variable first
		if len(Image_TrendAxisOverride) > 0:
			TrendString = filter(lambda x: Image_TrendAxisOverride in x, SplitSimName)[0]
		else: 
			try: TrendString = filter(lambda x: VariableString in x, SplitSimName)[0]
			except: TrendString = SplitSimName[0]
		#endif

		#Attempt to convert trend value to float, if not use string.
		try: TrendValue = float( TrendString.partition('#')[2] )
		except: TrendValue = TrendString
		#endtry
		TrendAxis.append(TrendValue)
	#endfor

	return(TrendAxis)
#enddef

#=========================#


#=====================================================================#
#=====================================================================#
















#====================================================================#
						  #SOFTWARE SPLASH#
#====================================================================#

print ''
print '---------------------------------------------------------------------'
print '.______   ____    ____  _______     _______.___________.    ___      '   
print '|   _  \  \   \  /   / |   ____|   /       |           |   /   \     '
print '|  |_)  |  \   \/   /  |  |__     |   (----`---|  |----`  /  ^  \    '
print '|   ___/    \_    _/   |   __|     \   \       |  |      /  /_\  \   '
print '|  |          |  |     |  |____.----)   |      |  |     /  _____  \  '
print '| _|          |__|     |_______|_______/       |__|    /__/     \__\ '
print '                                                               V0.2.0'
print '---------------------------------------------------------------------'
print ''
print 'The following diagnostics were requested:'
print '-----------------------------------------'
if IAutorun == True:
	print'# Simulation Series Autorun'
	print''
if True in [savefig_EfitEquilTrends,savefig_UserEquilTrends]:
	#[savefig_EquilSeperatrix,savefig_EquilMidplane,savefig_EquilXpoint]
	print'# 2D Equilibrium Analysis'
if True in [savefig_PlasmaCurrent]:
	print'# 1D Plasma Current Analysis'
if True in [savefig_CoilCurrentTraces]:
	print'# 1D Coil Current Analysis'
print '-----------------------------------------'
print ''

#=====================================================================#
#=====================================================================#



#====================================================================#
					  #FIESTA AUTORUN ROUTINE#
#====================================================================#

#Auto generate series folder name if requested
if SeriesName == 'auto': SeriesName = 'Vary '+ParameterVaried
#Ensure varied parameter appears first in SimulationName
SimNameList = [SimNameList[i] for i in range(0,len(SimNameList)) if SimNameList[i]!=ParameterVaried]
SimNameList = [ParameterVaried]+SimNameList


#Autorun simulations over defined paramter range if requested
if IAutorun == True:

	#Create simulation series folder and obtain folder directories
	HomeDir = os.getcwd()
	SeriesDirString = '/'+SeriesName+' '+ProjectName+'/'
	SeriesDir = CreateNewFolder(HomeDir,SeriesDirString)

	#For all requested input parameters
	for i in range(0,len(ParameterRange)):

		#Create simulation folder for input parameter[i]
		SimulationString = CreateSimName(SimNameList,ParameterVaried,ParameterRange[i])
		SimulationDir = CreateNewFolder(SeriesDir,SimulationString)
		if IVerbose == False: print SimulationString

		#Copy FIESTA.m into simulation folder and cd into directory
		os.system('cp '+FIESTAName+' '+'\''+SimulationDir+'\'')
		os.chdir(SimulationDir)

		#Update new FIESTA.m with fixed namelist parameters
		MatlabProjectString = '\''+ProjectName+'\''
		MatlabSimulationString = '\''+SimulationString+'\''
		MatlabIEquilMethod = '\''+IEquilMethod+'\''
		AlteredEntry = AlterNamelistVariable(FIESTAName,'ProjectName',MatlabProjectString)
		AlteredEntry = AlterNamelistVariable(FIESTAName,'SimName',MatlabSimulationString)
		AlteredEntry = AlterNamelistVariable(FIESTAName,'IEquilMethod',MatlabIEquilMethod)
		AlteredEntry = AlterNamelistVariable(FIESTAName,'NumThreads',NumThreads)

		#####

		#Update new FIESTA.m with variable namelist parameters for Parameter[i]
		AlteredEntry = AlterNamelistVariable(FIESTAName,ParameterVaried,ParameterRange[i])

		#Run modified FIESTA - Verbosity determines terminal output.
		#TO IMPLIMENT::: MAXIMUM CONCURRENT RUNS = MaxNumThreads/NumThreads
		RunFIESTA(FIESTAName,Verbose=IVerbose,Parallel=IParallel)

		#Return to home directory to enable diagnostic processing
		os.chdir(HomeDir)
	#endfor

	#=================#

	#If parallel simulations have been requested:
	if IParallel == True:
		#Initial delay to allow MATLAB processes to start
		time.sleep(ConvDelay/2.0)
		TimeConv = ConvDelay/2.0

		#Parallel Race Condition Checker - Waits for all simulations to finish before analysis
		while CheckIfProcessRunning('MATLAB') == True:
			#If MATLAB process is detected, wait ConvDelay seconds and check again
			time.sleep(ConvDelay)
			TimeConv += ConvDelay
			print 'Awaiting Series Convergence:',str(TimeConv)+'[s]'
		#endwhile

		#Update user of simulation convergence
		print '------------------------------------------------'
		print 'Simulation Series Converged:',str(TimeConv)+'[s]'
		print '------------------------------------------------'
	#endfor
#endif

#=====================================================================#
#=====================================================================#



















#====================================================================#
					   #ANALYSIS AND DIAGNOSTICS#
#====================================================================#

#====================================================================#
				       #Efit EQUILIBRIUM TRENDS#
#====================================================================#

#Plot general equilibrium trends from Param(equil)
if savefig_EfitEquilTrends == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract equilibrium data from series directories
	Filename = 'Equil_Data/EquilParam.txt'
	ParamEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	ValueEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#Quick and dirty removal of most useful trends
	RGeo,ZGeo,Kappa,A,delta = list(),list(),list(),list(),list()	#efit params
	for l in range(0,len(SimulationDirs)):
		RGeo.append( ValueEquil[l][43] )		#Geomoetric Radial Length 	[m]
#		ZGeo.append( ValueEquil[l][44] )		#Geometric Axial Length 	[m]
		Kappa.append( ValueEquil[l][21] )		#Elongation 				[-]
		A.append( ValueEquil[l][16] )			#Aspect Ratio 				[-]
		delta.append( ValueEquil[l][25] )		#Triangularity (average) 	[-]
	#endfor

	#RGeo_efit = 0.44					# Geometrical Radius	[m] (0.44)
	#ZGeo_efit = 0.0					# Geometrical Axis		[m] (0.00)
	#A_efit = 1.85						# Aspect Ratio			[-] (1.85)
	#a_efit = RGeo_efit/A_efit			# Minor Radius			[m] (0.44/1.85)
	#Kappa_efit = 1.8					# Elongation			[-] (1.8)
	#delta_efit = 0.2					# Triangularity			[-] (0.2)

	#===================##===================#
	#===================##===================#

#	#Create output folder for all coil trend figures
#	EquilTrendsDir = CreateNewFolder(SeriesDirString,'/Equil_Trends/')
	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(2,2, figsize=(13,11))

	fig.suptitle('Equilibrium efit_Geometry() Trends', y=0.99, fontsize=24)

	#Plot requested equil parameter trends over full simulation series
	ax[0,0].plot(TrendAxis,RGeo,'ko-', ms=12, lw=2)
#	ax[0,0].plot((min(TrendAxis),max(TrendAxis)),(RGeo_efit,RGeo_efit), 'k--', lw=1.5)
	ax[0,0].set_ylabel('Geometric Radius R$_{\mathrm{Geo}}$ [m]', fontsize=25)
#	ax[0,0].set_xlabel('Varied Parameter: '+Parameter, fontsize=25)
	ax[0,0].tick_params(axis='x', labelsize=20)
	ax[0,0].tick_params(axis='y', labelsize=20)

	ax[1,0].plot(TrendAxis,A,'bs-', ms=12, lw=2)
#	ax[1,0].plot((min(TrendAxis),max(TrendAxis)),(A_efit,A_efit), 'b--', lw=1.5)
	ax[1,0].set_ylabel('Aspect Ratio $A$ [-]', fontsize=25)
	ax[1,0].set_xlabel('Varied Parameter: '+Parameter, fontsize=25)
	ax[1,0].tick_params(axis='x', labelsize=20)
	ax[1,0].tick_params(axis='y', labelsize=20)

	ax[0,1].plot(TrendAxis,Kappa,'r^-', ms=12, lw=2)
#	ax[0,1].plot((min(TrendAxis),max(TrendAxis)),(Kappa_efit,Kappa_efit), 'r--', lw=1.5)
	ax[0,1].set_ylabel('Elongation $\kappa$ [-]', fontsize=25)
#	ax[0,1].set_xlabel('Varied Parameter: '+Parameter, fontsize=25)
	ax[0,1].tick_params(axis='x', labelsize=20)
	ax[0,1].tick_params(axis='y', labelsize=20)

	ax[1,1].plot(TrendAxis,delta,'mh-', ms=12, lw=2)
#	ax[1,1].plot((min(TrendAxis),max(TrendAxis)),(delta_efit,delta_efit), 'm--', lw=1.5)
	ax[1,1].set_ylabel('Triangularity $\delta$ [-]', fontsize=25)
	ax[1,1].set_xlabel('Varied Parameter: '+Parameter, fontsize=25)
	ax[1,1].tick_params(axis='x', labelsize=20)
	ax[1,1].tick_params(axis='y', labelsize=20)

	plt.tight_layout(pad=2.0)
	plt.savefig(SeriesDirString+'/Equil_Trends.png')
#	plt.show()
	plt.close('all')

	print'----------------------------------'
	print'# Efit Equilibrium Trends Complete'
	print'----------------------------------'
#endif

#=====================================================================#
#=====================================================================#



#====================================================================#
				       #USER EQUILIBRIUM TRENDS#
#====================================================================#

#Plot general equilibrium trends from Param(equil)
if savefig_UserEquilTrends == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract equilibrium data from series directories
	Filename = 'Equil_Data/EquilParam.txt'
	ParamEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	ValueEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#List equilibrium parameters by index - !!! CONVERT INTO INDEX LIBRARY !!!
	#!!! MAKE LAMBDA FUNCTION TO EXTRACT REQUESTED PARAMETERS FOR PLOTTING !!!
	#!!! USER SUPPLIES WHICH TRENDS TO PLOT ON THIS FIGURE !!!
#	for i in range(len(ParamEquil[l])):
#		print i, ParamEquil[l][i]
#		requested_trend_parameter = lambda( [USER ARRAY OF TREND NAMES], ParamEquil[l])
#		requested_trend_index = ParamEquil[l].index( requested_trend_parameter )
#		requested_trend_value = ValueEquil[l][requested_trend_index]
	#endfor

	#USEFUL TRENDS TO TRACK
	#Vol(m3) #q95 #betaT #betap #betaN

	#===================##===================#
	#===================##===================#

#	#Create output folder for all coil trend figures
#	EquilTrendsDir = CreateNewFolder(SeriesDirString,'/Equil_Trends/')
	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(1, figsize=(12,10))

	#Plot requested equil parameter trends over full simulation series
	ax.plot(TrendAxis,UserEquilParameter,'ko-', ms=12, lw=2)
	####
	Legend = TrendAxis
	ax.legend(Legend, fontsize=22, frameon=False)
	ax.set_ylabel('Equilibrum Parameter [-]', fontsize=25)
	ax.set_xlabel('Varied Parameter: '+Parameter, fontsize=25)
	ax.tick_params(axis='x', labelsize=20)
	ax.tick_params(axis='y', labelsize=20)
#	ax.set_xlim(0.00,1.00)
#	ax.set_ylim(0.00,1.00)

	plt.tight_layout()
	plt.savefig(SeriesDirString+'/Equil_Trends.png')
#	plt.show()
	plt.close('all')

	print'----------------------------------'
	print'# User Equilibrium Trends Complete'
	print'----------------------------------'
#endif

#=====================================================================#
#=====================================================================#



#====================================================================#
				         #EQUILIBRIUM STABILITY#
#====================================================================#

if savefig_EquilStability == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract equilibrium data from series directories
	Filename = 'Equil_Data/EquilParam.txt'
	ParamEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	ValueEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Extract coil current data from series directories
	Filename = 'icoil_Data/efit_icoil.txt'
	ISol_efit = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	IPF1_efit = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	IPF2_efit = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	IDiv1_efit = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	IDiv2_efit = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]
	Filename = 'icoil_Data/Perturbed_icoil.txt'
	ISol_Pert = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	IPF1_Pert = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	IPF2_Pert = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	IDiv1_Pert = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	IDiv2_Pert = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#Quick and dirty removal of relevent trends
	RGeo,ZGeo,Kappa,A,delta = list(),list(),list(),list(),list()	#efit params
	for l in range(0,len(SimulationDirs)):
		RGeo.append( ValueEquil[l][43] )		#Geomoetric Radial Length 	[m]
#		ZGeo.append( ValueEquil[l][44] )		#Geometric Axial Length 	[m]
		Kappa.append( ValueEquil[l][21] )		#Elongation 				[-]
		A.append( ValueEquil[l][16] )			#Aspect Ratio 				[-]
		delta.append( ValueEquil[l][25] )		#Triangularity (average) 	[-]
	#endfor

	#Calculate ratio of coil currents between efit and perturbed equilibria
	PertFracISol,PertFracIPF1,PertFracIPF2 = list(),list(),list()
	PertFracIDiv1,PertFracIDiv2 = list(),list()
	for i in range(0,len(ISol_efit)):
		try: PertFracISol.append( ISol_Pert[i]/ISol_efit[i] )
		except: PertFracISol.append( np.nan )
		try: PertFracIPF1.append( IPF1_Pert[i]/IPF1_efit[i] )
		except: PertFracIPF1.append( np.nan )
		try: PertFracIPF2.append( IPF2_Pert[i]/IPF2_efit[i] )
		except: PertFracIPF2.append( np.nan )
		try: PertFracIDiv1.append( IDiv1_Pert[i]/IDiv1_efit[i] )
		except: PertFracIDiv1.append( np.nan )
		try: PertFracIDiv2.append( IDiv2_Pert[i]/IDiv2_efit[i] )
		except: PertFracIDiv2.append( np.nan )
	#endfor

	#min(filter(lambda v: v==v, EMinArrays[i]))

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(2, figsize=(12,14))

#	ax[0].plot(TrendAxis,ISol_Pert, 'ko-', ms=10, lw=2)
	ax[0].plot(TrendAxis,IPF1_efit, 'r--', ms=10, lw=2)
	ax[0].plot(TrendAxis,IPF1_Pert, 'r^-', ms=10, lw=2)
	ax[0].plot(TrendAxis,IPF2_efit, 'b--', ms=10, lw=2)
	ax[0].plot(TrendAxis,IPF2_Pert, 'bs-', ms=10, lw=2)
#	ax[0].plot(TrendAxis,IDiv1_Pert, 'c*-', ms=10, lw=2)
#	ax[0].plot(TrendAxis,IDiv2_Pert, 'mh-', ms=10, lw=2)

	Title = 'Variation in coil current required \n to restore equilibrium for varying '+Parameter
	ax[0].set_title(Title, fontsize=20, y=1.03)
	Legend = ['PF1 efit','PF1 Pert','PF2 efit','PF2 Pert']
	ax[0].legend(Legend, fontsize=22, frameon=False)
	ax[0].set_ylabel('Coil Current [kA]', fontsize=25)
#	ax[0].set_xlabel('Varied Parameter: '+Parameter, fontsize=25)
	ax[0].tick_params(axis='x', labelsize=20)
	ax[0].tick_params(axis='y', labelsize=20)
#	ax[0].set_xlim(0,1)		
#	ax[0].set_ylim(2,32)

	##########

	#Plot plasma current with respect to adaptive_time
#	ax[1].plot(TrendAxis,PertFracISol, 'ko-', ms=10, lw=2)
	ax[1].plot(TrendAxis,PertFracIPF1, 'r^-', ms=10, lw=2)
	ax[1].plot(TrendAxis,PertFracIPF2, 'bs-', ms=10, lw=2)
#	ax[1].plot(TrendAxis,PertFracIDiv1, 'c*-', ms=10, lw=2)
#	ax[1].plot(TrendAxis,PertFracIDiv2, 'mh-', ms=10, lw=2)

	Legend = ['PF1','PF2']	
	ax[1].legend(Legend, fontsize=22, frameon=False)
	ax[1].set_ylabel('Fractional Change In \n Coil Current [-]', fontsize=25)
	ax[1].set_xlabel('Varied Parameter: '+Parameter, fontsize=25)
	ax[1].tick_params(axis='x', labelsize=20)
	ax[1].tick_params(axis='y', labelsize=20)
#	ax[1].set_xlim(0,1)		
#	ax[1].set_ylim(0,1)

	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(SeriesDirString+'/VerticalStability_Trends.png')
#	plt.show()
	plt.close('all')

	print'--------------------------------'
	print'# Equilibrium stability Complete'
	print'--------------------------------'
#endif

#=====================================================================#
#=====================================================================#




















#====================================================================#
				  #COIL CURRENT WAVEFORM DIAGNOSTIC#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_CoilCurrentTraces == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract coil currents and time axis from series directories
	Filename = 'icoil_Data/CoilCurrents.txt'
	Time_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	ISol_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	IPF1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	IPF2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	IDiv1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]
	IDiv2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[5]

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#Rescale data for plotting: [s] to [ms]
	for i in range(0,len(Time_Arrays)):
		for j in range(0,len(Time_Arrays[i])):
			Time_Arrays[i][j] = Time_Arrays[i][j]*1000.0
		#endfor
	#endfor

	#Rescale data for plotting: [A] to [kA]
	for i in range(0,len(ISol_Arrays)):
		for j in range(0,len(ISol_Arrays[i])):
		 	#Coil currents are saved scaled by the number of windings (For reasons...)
			ISol_Arrays[i][j] = ISol_Arrays[i][j]/(1000.0*nSol)
			IPF1_Arrays[i][j] = IPF1_Arrays[i][j]/(1000.0*24)
			IPF2_Arrays[i][j] = IPF2_Arrays[i][j]/(1000.0*24)
			IDiv1_Arrays[i][j] = IDiv1_Arrays[i][j]/(1000.0*24)
			IDiv2_Arrays[i][j] = IDiv2_Arrays[i][j]/(1000.0*24)
		#endfor
	#endfor

	#Calculate dI/dt for each coil set
	DeltaIPF1,DeltaIPF2 = list(),list()
	DeltaIDiv1,DeltaIDiv2 = list(),list()
	DeltaISol = list()
	for i in range(0,len(ISol_Arrays)):
		DeltaISol.append(list())
		DeltaIPF1.append(list())
		DeltaIPF2.append(list())
		DeltaIDiv1.append(list())
		DeltaIDiv2.append(list())
		for j in range(1,len(ISol_Arrays[i])):
			Delta_t = Time_Arrays[i][j]-Time_Arrays[i][j-1]
			#
			DeltaISol[i].append( (ISol_Arrays[i][j]-ISol_Arrays[i][j-1])/Delta_t )
			DeltaIPF1[i].append( (IPF1_Arrays[i][j]-IPF1_Arrays[i][j-1])/Delta_t )
			DeltaIPF2[i].append( (IPF2_Arrays[i][j]-IPF2_Arrays[i][j-1])/Delta_t )
			DeltaIDiv1[i].append( (IDiv1_Arrays[i][j]-IDiv1_Arrays[i][j-1])/Delta_t )
			DeltaIDiv2[i].append( (IDiv2_Arrays[i][j]-IDiv2_Arrays[i][j-1])/Delta_t )
		#endfor
	#endfor

	#===================##===================#
	#===================##===================#

#	#Create output folder for all coil trend figures
#	ICoilTimeTracesDir = CreateNewFolder(SeriesDirString,'/ICoil_Trends/')
	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#For every simulation folder in the current series:
	for l in range(0,len(ISol_Arrays)):

		#Create figure for each Coil Ramp Time Trace diagnostic
		fig,ax = plt.subplots(2, figsize=(12,14), sharex=True)

		#Plot each coil current with respect to time
		ax[0].plot(Time_Arrays[l],ISol_Arrays[l], 'k-', lw=2)
		ax[0].plot(Time_Arrays[l],IPF1_Arrays[l], 'r-', lw=2)
		ax[0].plot(Time_Arrays[l],IPF2_Arrays[l], 'b-', lw=2)
		ax[0].plot(Time_Arrays[l],IDiv1_Arrays[l], 'c-', lw=2)
		ax[0].plot(Time_Arrays[l],IDiv2_Arrays[l], 'm-', lw=2)

		Range = '['+str(min(TrendAxis))+' - '+str(max(TrendAxis))+']'
		ax[0].set_title('Time-Traces of Coil Currents for '+Parameter+' in '+Range, fontsize=20, y=1.03)
		Legend = ['Sol','PF1','PF2','Div1','Div2']
		ax[0].legend(Legend, fontsize=22, ncol=2, frameon=False)
		ax[0].set_ylabel('Coil Current $I$ [kA]', fontsize=25)
#		ax[0].set_xlabel('Time $\\tau$ [ms]', fontsize=25)
#		ax[0].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#		ax[0].yaxis.set_major_locator(ticker.MultipleLocator(240))
		ax[0].tick_params(axis='x', labelsize=20)
		ax[0].tick_params(axis='y', labelsize=20)
#		ax[0].set_xlim(-50,100)		
#		ax[0].set_ylim(2,32)

		#Plot derivitive of each coil current with respect to time
		ax[1].plot(Time_Arrays[l][1::],DeltaISol[l], 'k-', lw=2)
		ax[1].plot(Time_Arrays[l][1::],DeltaIPF1[l], 'r-', lw=2)
		ax[1].plot(Time_Arrays[l][1::],DeltaIPF2[l], 'b-', lw=2)
		ax[1].plot(Time_Arrays[l][1::],DeltaIDiv1[l], 'c-', lw=2)
		ax[1].plot(Time_Arrays[l][1::],DeltaIDiv2[l], 'm-', lw=2)

		Range = '['+str(min(TrendAxis))+' - '+str(max(TrendAxis))+']'
		ax[1].set_title('Time-Traces of Delta Coil Currents for '+Parameter+' in '+Range, fontsize=20, y=1.03)
		Legend = ['Sol','PF1','PF2','Div1','Div2']
		ax[1].legend(Legend, fontsize=22, ncol=2, frameon=False)
		ax[1].set_ylabel('Change in Coil Current \n $\Delta I$ [kA ms$^{-1}$]', fontsize=25)
		ax[1].set_xlabel('Time $\\tau$ [ms]', fontsize=25)
#		ax[1].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#		ax[1].yaxis.set_major_locator(ticker.MultipleLocator(240))
		ax[1].tick_params(axis='x', labelsize=20)
		ax[1].tick_params(axis='y', labelsize=20)
#		ax[1].set_xlim(-50,100)		
#		ax[1].set_ylim(2,32)

		plt.tight_layout(pad=3.0,h_pad=1.0)
#		plt.savefig(ICoilTimeTracesDir+'/CoilRamp_'+str(TrendAxis[l])+'TimeTrace.png')
		plt.savefig(SeriesDirString+'/CoilRamp_'+str(TrendAxis[l])+'TimeTrace.png')

		plt.show()
		plt.close('all')
	#endfor
	print'----------------------------------'
	print'# Coil Current Timetraces Complete'
	print'----------------------------------'

#=====================================================================#
#=====================================================================#



#====================================================================#
				  #COIL CURRENT TRENDS DIAGNOSTIC#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_CoilCurrentTrends == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract coil currents and time axis from series directories
	Filename = 'icoil_Data/CoilCurrents.txt'
	NumCoils = len(ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical'))-1
	Time_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	ISol_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	IPF1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	IPF2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	IDiv1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]
	IDiv2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[5]


	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#Rescale data for plotting: [s] to [ms]
	for i in range(0,len(Time_Arrays)):
		for j in range(0,len(Time_Arrays[i])):
			Time_Arrays[i][j] = Time_Arrays[i][j]*1000.0
		#endfor
	#endfor

	#Rescale data for plotting: [A] to [kA]
	for i in range(0,len(ISol_Arrays)):
		for j in range(0,len(ISol_Arrays[i])):
		 	#Coil currents are saved scaled by the number of windings (For reasons...)
			ISol_Arrays[i][j] = ISol_Arrays[i][j]/(1000.0*nSol)
			IPF1_Arrays[i][j] = IPF1_Arrays[i][j]/(1000.0*24)
			IPF2_Arrays[i][j] = IPF2_Arrays[i][j]/(1000.0*24)
			IDiv1_Arrays[i][j] = IDiv1_Arrays[i][j]/(1000.0*24)
			IDiv2_Arrays[i][j] = IDiv2_Arrays[i][j]/(1000.0*24)
		#endfor
	#endfor

	#Calculate maximum coil current for each coil
	MaxIPF1,MaxIPF2 = list(),list()
	MaxIDiv1,MaxIDiv2 = list(),list()
	MaxISol = list()
	MaxIAvg = list()
	for i in range(0,len(ISol_Arrays)):
		MaxISol.append( max(ISol_Arrays[i], key=abs) )
		MaxIPF1.append( max(IPF1_Arrays[i], key=abs) )
		MaxIPF2.append( max(IPF2_Arrays[i], key=abs) )
		MaxIDiv1.append( max(IDiv1_Arrays[i], key=abs) )
		MaxIDiv2.append( max(IDiv2_Arrays[i], key=abs) )
		MaxIAvg.append( (abs(MaxISol[i])+abs(MaxIPF1[i])+abs(MaxIPF2[i])+abs(MaxIDiv1[i])+abs(MaxIDiv2[i]))/NumCoils )
	#endfor

	#Calculate dI/dt for each coil set
	DeltaIPF1,DeltaIPF2 = list(),list()
	DeltaIDiv1,DeltaIDiv2 = list(),list()
	DeltaISol = list()
	for i in range(0,len(ISol_Arrays)):
		DeltaISol.append(list())
		DeltaIPF1.append(list())
		DeltaIPF2.append(list())
		DeltaIDiv1.append(list())
		DeltaIDiv2.append(list())
		for j in range(1,len(ISol_Arrays[i])):
			Delta_t = Time_Arrays[i][j]-Time_Arrays[i][j-1]
			#
			DeltaISol[i].append( (ISol_Arrays[i][j]-ISol_Arrays[i][j-1])/Delta_t )
			DeltaIPF1[i].append( (IPF1_Arrays[i][j]-IPF1_Arrays[i][j-1])/Delta_t )
			DeltaIPF2[i].append( (IPF2_Arrays[i][j]-IPF2_Arrays[i][j-1])/Delta_t )
			DeltaIDiv1[i].append( (IDiv1_Arrays[i][j]-IDiv1_Arrays[i][j-1])/Delta_t )
			DeltaIDiv2[i].append( (IDiv2_Arrays[i][j]-IDiv2_Arrays[i][j-1])/Delta_t )
		#endfor
	#endfor

	#Calculate maximum change in current experienced for each coil set
	MaxDeltaIPF1,MaxDeltaIPF2 = list(),list()
	MaxDeltaIDiv1,MaxDeltaIDiv2 = list(),list()
	MaxDeltaISol = list()
	for i in range(0,len(DeltaISol)):
		MaxDeltaISol.append( max(DeltaISol[i], key=abs) )
		MaxDeltaIPF1.append( max(DeltaIPF1[i], key=abs) )
		MaxDeltaIPF2.append( max(DeltaIPF2[i], key=abs) )
		MaxDeltaIDiv1.append( max(DeltaIDiv1[i], key=abs) )
		MaxDeltaIDiv2.append( max(DeltaIDiv2[i], key=abs) )
	#endfor

	#===================##===================#
	#===================##===================#

#	#Create output folder for all coil trend figures
#	CurrentTrendsDir = CreateNewFolder(SeriesDirString,'/ICoil_Trends/')
	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#For every simulation folder in the current series:
	for i in range(0,NumCoils):

		#Set which coil current timetrace to compare
		if i == 0: 	
			Coil = 'ISol'
			Current_Arrays = ISol_Arrays
		if i == 1:
			Coil = 'IPF1'
			Current_Arrays = IPF1_Arrays
		if i == 2: 
			Coil = 'IPF2'	
			Current_Arrays = IPF2_Arrays
		if i == 3: 	
			Coil = 'IDiv1'
			Current_Arrays = IDiv1_Arrays
		if i == 4: 	
			Coil = 'IDiv2'
			Current_Arrays = IDiv2_Arrays
		#endif

		#Create figure to compare each coil current time-trace
		fig,ax = plt.subplots(1, figsize=(12,8))

		#Plot each coil current with respect to time
		for j in range(0,len(Current_Arrays)):
			ax.plot(Time_Arrays[j],Current_Arrays[j], lw=2)
		#endfor

		#Set title, legend and savestrings
		Range = '['+str(min(TrendAxis))+' - '+str(max(TrendAxis))+']'
		Title = 'Time-Traces of '+Coil+' Coil Currents for '+Parameter+' in '+Range
		Legend = TrendAxis
		SaveString = SeriesDirString+'/'+Coil+'_Current_Trends.png'

		ax.set_title(Title, fontsize=20, y=1.03)
		ax.legend(Legend, fontsize=22, ncol=2, frameon=False)
		ax.set_ylabel('Coil Current $I$ [kA]', fontsize=25)
		ax.set_xlabel('Time $\\tau$ [ms]', fontsize=25)
#		ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#		ax.yaxis.set_major_locator(ticker.MultipleLocator(240))
		ax.tick_params(axis='x', labelsize=20)
		ax.tick_params(axis='y', labelsize=20)
#		ax.set_xlim(-50,100)		
#		ax.set_ylim(2,32)

		plt.tight_layout(pad=3.0,h_pad=1.0)
		plt.savefig(SaveString)
		plt.show()
		plt.close('all')
	#endfor

	print'------------------------------'
	print'# Coil Current Trends Complete'
	print'------------------------------'

	#=================#

	#Create image limits
	GlobalMaxDelta = MaxDeltaISol+MaxDeltaIPF1+MaxDeltaIPF2+MaxDeltaIDiv1+MaxDeltaIDiv2
	Ylims = [min(GlobalMaxDelta),max(GlobalMaxDelta)]

	#Create figure for Coil Maximum Ramp Diagnostic
	fig,ax = plt.subplots(2, figsize=(12,14))

	#Plot derivitive of each coil current with respect to time
	ax[0].plot(TrendAxis,MaxISol, 'ko-', ms=10, lw=2)
	ax[0].plot(TrendAxis,MaxIPF1, 'r^-', ms=10, lw=2)
	ax[0].plot(TrendAxis,MaxIPF2, 'bs-', ms=10, lw=2)
	ax[0].plot(TrendAxis,MaxIDiv1, 'c*-', ms=10, lw=2)
	ax[0].plot(TrendAxis,MaxIDiv2, 'mh-', ms=10, lw=2)
	ax[0].plot(TrendAxis,MaxIAvg, 'kv:', markerfacecolor='none', ms=10, lw=2)			 #Avg ICoil
	ax[0].plot(TrendAxis[MaxIAvg.index(min(MaxIAvg))],min(MaxIAvg), 'kv', ms=14, lw=2.0) #Min Avg

	ax[0].set_title('Maximum Coil Current for Varying '+Parameter, fontsize=20, y=1.03)
	Legend = ['Sol','PF1','PF2','Div1','Div2','Avg']
	ax[0].legend(Legend, fontsize=22, ncol=2, frameon=False)
	ax[0].set_ylabel('Maximum Coil \n Current $I$ [kA]', fontsize=25)
#	ax[0].set_xlabel(Parameter, fontsize=25)
#	ax[0].xaxis.set_major_locator(ticker.MultipleLocator( (max(TrendAxis)-min(TrendAxis))/5 ))
#	ax[0].yaxis.set_major_locator(ticker.MultipleLocator(50))
	ax[0].tick_params(axis='x', labelsize=20)
	ax[0].tick_params(axis='y', labelsize=20)
	ax[0].set_xlim(TrendAxis[0],TrendAxis[-1])		
#	ax[0].set_ylim(-1.5,1.5)

	##########

	#Plot derivitive of each coil current with respect to time
	ax[1].plot(TrendAxis,MaxDeltaISol, 'ko-', ms=10, lw=2)
	ax[1].plot(TrendAxis,MaxDeltaIPF1, 'r^-', ms=10, lw=2)
	ax[1].plot(TrendAxis,MaxDeltaIPF2, 'bs-', ms=10, lw=2)
	ax[1].plot(TrendAxis,MaxDeltaIDiv1, 'c*-', ms=10, lw=2)
	ax[1].plot(TrendAxis,MaxDeltaIDiv2, 'mh-', ms=10, lw=2)

	ax[1].set_title('Maximum Delta Coil Current for Varying '+Parameter, fontsize=20, y=1.03)
	Legend = ['Sol','PF1','PF2','Div1','Div2']
	ax[1].legend(Legend, fontsize=22, ncol=2, frameon=False)
	ax[1].set_ylabel('Maximum Change in \n Current $\Delta I$ [kA ms$^{-1}$]', fontsize=25)
	ax[1].set_xlabel(Parameter, fontsize=25)
#	ax[1].xaxis.set_major_locator(ticker.MultipleLocator( (max(TrendAxis)-min(TrendAxis))/5 ))
#	ax[1].yaxis.set_major_locator(ticker.MultipleLocator(50))
	ax[1].tick_params(axis='x', labelsize=20)
	ax[1].tick_params(axis='y', labelsize=20)
	ax[1].set_xlim(TrendAxis[0],TrendAxis[-1])		
	ax[1].set_ylim(Ylims[0]*1.25,Ylims[1]*1.25)

	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(SeriesDirString+'/CoilRamp_Trends.png')
	plt.show()
	plt.close('all')

	print'----------------------------'
	print'# Coil Current Ramp Complete'
	print'----------------------------'
#endfor

#=====================================================================#
#=====================================================================#


















#====================================================================#
					  #PLASMA BREAKDOWN DIAGNOSTICS#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_ConnectionLength == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract relevent data from series directories
	Filename = 'LCon.txt'
	Lc_Array = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	Filename = 'Eta.txt'
	Eta_Perp_Array = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	Eta_Para_Array = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(1, figsize=(12,10))

	#Plot plasma current with respect to adaptive_time
	ax.plot(TrendAxis,Lc_Array, 'ko-', ms=10, lw=2)

	ax.set_title('Connection Length for Varying '+Parameter, fontsize=20, y=1.03)
	Legend = ['Lc']
	ax.legend(Legend, fontsize=22, frameon=False)
	ax.set_ylabel('Connection Length [m]', fontsize=25)
	ax.set_xlabel('Varying: '+Parameter, fontsize=25)
	ax.ticklabel_format(style='sci', axis='x', scilimits=(-2,2))
#	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#	ax.yaxis.set_major_locator(ticker.MultipleLocator(240))
	ax.tick_params(axis='x', labelsize=20)
	ax.tick_params(axis='y', labelsize=20)
#	ax.set_xlim(0,1)		
#	ax.set_ylim(2,32)

	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(SeriesDirString+'/Lc_Trends.png')
#	plt.show()
	plt.close('all')

	print'------------------------'
	print'# Lc Diagnostic Complete'
	print'------------------------'
#endif

#=====================================================================#
#=====================================================================#



#=====================================================================#
#=====================================================================#

if savefig_PaschenCurves == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract relevent data from series directories
	Filename = 'LCon.txt'
	Lc_Array = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Create arbitary pressure array over 1E-5 --> 1E-3 Torr
	Limits,Resolution = [1E-5,1E-2],10000
	PressureArray = np.linspace(Limits[0],Limits[1],Resolution).tolist()

	#Initialise required arrays EMinArrays 
	EMinArrays,MinEMinArray = list(),list()
	for i in range(0,len(Lc_Array)): EMinArrays.append( list() )

	#Construct EMinArrays for each Lcon by varying background pressure
	#N.B. !!! Lcon will vary with background pressure so this is only approximate !!!
	for i in range(0,len(EMinArrays)):
		for j in range(0,len(PressureArray)):

			#Compute minimum E-field for breakdown as cited from Song2017
			EMin = (PressureArray[j]*1.25E4)/np.log(510.0*PressureArray[j]*Lc_Array[i])
			if EMin < 0: EMin = np.nan

			EMinArrays[i].append( EMin )
		#endfor
		MinEMinArray.append( min(filter(lambda v: v==v, EMinArrays[i])) )  
	#endfor

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Round connective length for plotting
	for i in range(0,len(Lc_Array)): Lc_Array[i] = round(Lc_Array[i],1)
	#Scale pressure array for plotting
	for i in range(0,len(PressureArray)): PressureArray[i] = PressureArray[i]*1000	#[mTorr]

	#===================##===================#
	#===================##===================#

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(1, figsize=(12,10), sharex=True)

	#Plot minimum electric field for breakdown for each connection length
	for i in range(0,len(Lc_Array)):
		ax.plot(PressureArray,EMinArrays[i], '--', lw=2, ms=12)
	#####
	Legend = [Parameter+'='+str(TrendAxis[i]) for i in range(0,len(TrendAxis))]
	ax.legend(Legend, fontsize=20, loc=4, frameon=False)
	ax.set_xlabel('Prefill Pressure $P$ [mTorr]', fontsize=26)
	ax.set_ylabel('Toroidal E-Field $\\bf{E}$ [V m$^{-1}$]', fontsize=26)
	#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5E21))
	ax.tick_params(axis='x', labelsize=20)
	ax.tick_params(axis='y', labelsize=20)
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlim(1e-2,10)		#Pressure [mTorr]
	ax.set_ylim(0.25,25)		#E-Field [Vm-1]

	#Plot trend in minimum breakdown E-field with respect to varied parameter
	from mpl_toolkits.axes_grid.inset_locator import inset_axes
	left, bottom, width, height = [0.19,0.19,0.25,0.20]			#[0.62,0.27,0.23,0.23]
	ax2 = fig.add_axes([left, bottom, width, height])
	###
	ax2.plot(TrendAxis,MinEMinArray,'ko--', markerfacecolor='none', ms=8, lw=1.5)
	ax2.plot(TrendAxis[MinEMinArray.index(min(MinEMinArray))],min(MinEMinArray),'ko', ms=10)
	ax2.set_ylabel('Minimum $\\bf{E}$ [V m$^{-1}$]', labelpad=0, fontsize=14.5)
	ax2.set_xlabel('Varying: '+Parameter, fontsize=15)
	ax2.tick_params(axis='x', labelsize=14)
	ax2.tick_params(axis='y', labelsize=14)
#	ax2.set_xlim(0,1)
#	ax2.set_ylim(0.79,1.01)

	plt.tight_layout(pad=2.0)
	plt.savefig(SeriesDirString+'/Paschen_Breakdown.png')
#	plt.show()
	plt.close('all')

	print'-----------------------------'
	print'# Paschen Diagnostic Complete'
	print'-----------------------------'
#endif

#=====================================================================#
#=====================================================================#


















#====================================================================#
					  #PLASMA CURRENT DIAGNOSTICS#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_PlasmaCurrent == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract plasma current data from series directories
	Filename = 'Ip.txt'
	Time_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	Ip_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Remove any header string from the data
	for i in range(0,len(Time_Arrays)):
		Time_Arrays[i] = Time_Arrays[i][1::]
		Ip_Arrays[i] = Ip_Arrays[i][1::]
	#endfor

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#Rescale data for plotting: [s] to [ms]
	for i in range(0,len(Time_Arrays)):
		for j in range(0,len(Time_Arrays[i])):
			Time_Arrays[i][j] = Time_Arrays[i][j]*1000.0
		#endfor
	#endfor

	#Rescale data for plotting: [A] to [kA]
	for i in range(0,len(Ip_Arrays)):
		for j in range(0,len(Ip_Arrays[i])):
			Ip_Arrays[i][j] = Ip_Arrays[i][j]/1000.0
		#endfor
	#endfor

	#Calculate maximum Ip for each simulation over the full series
	Ip_MaxTrend,Ip_MinTrend = list(),list()
	for i in range(0,len(Ip_Arrays)):
		Ip_MaxTrend.append(max(Ip_Arrays[i]))
		Ip_MinTrend.append(min(Ip_Arrays[i]))
	#endfor

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(1, figsize=(12,10))

	#Plot plasma current with respect to adaptive_time
	for i in range(0,len(Ip_Arrays)): ax.plot(Time_Arrays[i],Ip_Arrays[i], lw=2)
	ax.set_title('Plasma Current Time-Trace for Varying '+Parameter, fontsize=20, y=1.03)
	ax.legend(TrendAxis, fontsize=22, loc=1, frameon=False)
	ax.set_ylabel('Plasma Current $I_{p}$ [kA]', fontsize=25)
	ax.set_xlabel('Time $\\tau$ [ms]', fontsize=25)
#	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#	ax.yaxis.set_major_locator(ticker.MultipleLocator(240))
	ax.tick_params(axis='x', labelsize=20)
	ax.tick_params(axis='y', labelsize=20)
	ax.set_xlim( min(Time_Arrays[0])*1.20,max(Time_Arrays[0])*1.50 )		
#	ax.set_ylim(2,32)

	#Plot trend in plasma current with respect to varied parameter
	from mpl_toolkits.axes_grid.inset_locator import inset_axes
	left, bottom, width, height = [0.23,0.63,0.25,0.25]			#[0.62,0.27,0.23,0.23]
	ax2 = fig.add_axes([left, bottom, width, height])
	###
	ax2.plot(TrendAxis,Ip_MaxTrend,'ko--', ms=8, lw=1.5)
	ax2.legend(['Max $I_{p}$'], fontsize=14, frameon=False)
	ax2.set_ylabel('Maximum Plasma \n Current $I_{p,max}$ [kA]', labelpad=0, fontsize=14.5)
	ax2.set_xlabel('Varied Parameter: '+Parameter, fontsize=15)
#	ax2.xaxis.set_major_locator(ticker.MultipleLocator(90))
#	ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax2.tick_params(axis='x', labelsize=14)
	ax2.tick_params(axis='y', labelsize=14)
#	ax2.set_xlim( min(TrendAxis),max(TrendAxis)*1.10 )
#	ax2.set_ylim(0.79,1.01)

	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(SeriesDirString+'/Ip_Trends.png')
#	plt.show()
	plt.close('all')

	print'-------------------------'
	print'# Ip Diagnostics Complete'
	print'-------------------------'
#endif

#=====================================================================#
#=====================================================================#






#====================================================================#
					  #EDDY CURRENT DIAGNOSTICS#
#====================================================================#

#Compare vessel eddy current profiles
if savefig_EddyCurrent == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+' '+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract plasma current data from series directories
	Filename = 'IPass.txt'
	Time_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	IPass_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Remove any header string from the data
	for i in range(0,len(Time_Arrays)):
		Time_Arrays[i] = Time_Arrays[i][1::]
		IPass_Arrays[i] = IPass_Arrays[i][1::]
	#endfor

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,Image_TrendAxisOverride)

	#Rescale data for plotting: [s] to [ms]
	for i in range(0,len(Time_Arrays)):
		for j in range(0,len(Time_Arrays[i])):
			Time_Arrays[i][j] = Time_Arrays[i][j]*1000.0
		#endfor
	#endfor

	#Rescale data for plotting: [A] to [kA]
	for i in range(0,len(IPass_Arrays)):
		for j in range(0,len(IPass_Arrays[i])):
			IPass_Arrays[i][j] = IPass_Arrays[i][j]/1000.0
		#endfor
	#endfor

	#Calculate maximum Ip for each simulation over the full series
	IPass_MaxTrend,IPass_MinTrend = list(),list()
	for i in range(0,len(IPass_Arrays)):
		IPass_MaxTrend.append(max(IPass_Arrays[i]))
		IPass_MinTrend.append(min(IPass_Arrays[i]))
	#endfor

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(1, figsize=(12,10))

	#Plot plasma current with respect to adaptive_time
	for i in range(0,len(IPass_Arrays)): ax.plot(Time_Arrays[i],IPass_Arrays[i], lw=2)
	ax.set_title('Eddy Current Time-Trace for Varying '+Parameter, fontsize=20, y=1.03)
	ax.legend(TrendAxis, fontsize=22, loc=1, frameon=False)
	ax.set_ylabel('Eddy Current $I_{p}$ [kA]', fontsize=25)
	ax.set_xlabel('Time $\\tau$ [ms]', fontsize=25)
#	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#	ax.yaxis.set_major_locator(ticker.MultipleLocator(240))
	ax.tick_params(axis='x', labelsize=20)
	ax.tick_params(axis='y', labelsize=20)
	ax.set_xlim( min(Time_Arrays[0])*1.20,max(Time_Arrays[0])*1.50 )		
#	ax.set_ylim(2,32)

	#Plot trend in plasma current with respect to varied parameter
	from mpl_toolkits.axes_grid.inset_locator import inset_axes
	left, bottom, width, height = [0.23,0.63,0.25,0.25]			#[0.62,0.27,0.23,0.23]
	ax2 = fig.add_axes([left, bottom, width, height])
	###
	ax2.plot(TrendAxis,IPass_MaxTrend,'ko--', ms=8, lw=1.5)
	ax2.legend(['Max $I_{p}$'], fontsize=14, frameon=False)
	ax2.set_ylabel('Maximum Eddy \n Current $I_{p,max}$ [kA]', labelpad=0, fontsize=14.5)
	ax2.set_xlabel('Varied Parameter: '+Parameter, fontsize=15)
#	ax2.xaxis.set_major_locator(ticker.MultipleLocator(90))
#	ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax2.tick_params(axis='x', labelsize=14)
	ax2.tick_params(axis='y', labelsize=14)
#	ax2.set_xlim( min(TrendAxis),max(TrendAxis)*1.10 )
#	ax2.set_ylim(0.79,1.01)

	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(SeriesDirString+'/PassiveCurrent_Trends.png')
#	plt.show()
	plt.close('all')

	print'-------------------------'
	print'# Ip Diagnostics Complete'
	print'-------------------------'
#endif

#=====================================================================#
#=====================================================================#






