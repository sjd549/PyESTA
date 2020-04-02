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
mu0 = 1.2566e-06; # Magnetic Moment      [I/m^2]

#Define Operating Conditions
Te = 250;         # Electron Temperature [eV]
Ti = Te*0.1;      # Ion Temperature      [eV]
BT = 0.1;         # Toroidal B-Field     [T] (Defined at Rgeo)
Ip = 30e3;        # Plasma current       [A]
TauP = 0.020;     # Pulse length         [s] (Also determines tstep for Ip plot)
RGeo = 0.450;     # Geometrical Radius   [m]
ZGeo = 0.000;     # Geometrical Axis     [m]
RSep = 0.700;     # Separatrix Radius    [m]
a = RSep-RGeo     # Minor Radius         [m]
Epsilon = RGeo/a  # Aspect ratio         [-]
Kappa = 1.8;      # Elongation           [-]
li2 = 1;          # Standard Value?      [-]
#q_cyl = 2.821;   # Safety Factor?       [-]
betaN = 3.529;    #                      [-] (Obtained via VEST Excel - LIKELY 2X TOO HIGH)

#Compute Further Operating Conditions
Gr_Limit = 1e20*(Ip*1e-6/(pi*(a**2)*Kappa))	# Greenwald Limit         [m-3]
Gr_Frac = 0.15                            	# Greenwald Fraction       [-]
ne = Gr_Limit*Gr_Frac                     	# Electron Density         [m-3]  ~3E19
Irod = BT*2*pi*RGeo/mu0                   	# Central Rod Current      [A]
S = sqrt( (1.0+Kappa**2)/2.0 )             	# Shaping factor           [-]
betaT = betaN/a*(Ip/1e6)/BT             	# Beta toroidal            [%]
betaP = 3/2*ne*(Te+Ti)/(mu0*Ip/(2*pi*a))**2*2*mu0*1.6e-19*Kappa # Beta Poloidal  [%]
BZ = -mu0*Ip/(4*pi*RGeo)*(log(8*Epsilon)+betaP+0.5*li2-3/2)    # Vertical field [T]
nT = 2.66*1e23*betaT*BT**2               	# Density Temperature Product [eV m-3]
#
coil_density = 1                        	# Relative Coil Density    [Arb]
coil_temp = 293.0                       	# Initial Coil Temperature [K]



###################  DEFINE SOL RAMP & COIL CURRENTS  ###################

#Define number of time-steps (vertices) in the current waveforms
nTime = 6      #[Steps]
tstep = TauP   #[s]
#time = [-0.05, -0.03, 0, tstep, tstep+0.03, tstep+0.05]	#Phase1_JuanJo
time =  [-0.10, -0.05, 0, tstep, tstep+0.03, tstep+0.05]	#Phase1_Daniel
#time = [-0.11, -0.05, 0, tstep, tstep+0.10, tstep+0.11]	#Phase2_JuanJo
#time = [-0.xx, -0.xx, 0, tstep, tstep+0.xx, tstep+0.xx]   	#Phase3

#Phase1 coil currents [kA] 
#ISol_Waveform = [+1300,-500,-1300]
#IPF1_Waveform = [-370]
#IPF2_Waveform = [-400]
#IDiv1_Waveform = [+000]
#IDiv2_Waveform = [+900]

#Phase1 coil currents [kA]               #SJD210
I_Sol_Start=+725       #+1250 -> +xxxx;  %+725;
I_Sol_Equil=-300       #-0000 -> -0400;  %-0000;
I_Sol_End=-725         #-1250 -> -xxxx;  %-725;
#Symmetric ISol is better for power supply
#
I_PF1=-370             #-0350 -> -xxxx;  %-370;
I_PF2=-400             #-0400 -> -xxxx;  %-400;
I_Div1=+000            #-0000 -> -xxxx;  %+000;
I_Div2=+900            #+0900 -> +xxxx;  %+900;

#Phase2 coil currents [kA]
#I_Sol_Start=+4700;    #+4700 -> +4700;     #+4700;
#I_Sol_Equil=+500;     #+500  -> +500       #+500;
#I_Sol_End=-4700;      #-4700 -> -4700;     #-4700;
#
#I_PF1=-3000;          #-3000 -> -3000      #-3000;
#I_PF2=-940;           #-940  -> -940       #-940;
#I_Div1=-000;          #-000  -> -000       #-000;
#I_Div2=+9090;         #+9090 -> +9090      #+9090;

#Phase3 coil currents [kA]
#I_Sol_Start=4700;     #4200;
#I_Sol_Equil=500;      #
#I_Sol_End=-4700;      #-5200;
#
#I_PF1=-6000;          #
#I_PF2=-3100;          #
#I_Div1=-000;          #
#I_Div2=15880;         #

#====================================================================#
#====================================================================#





#====================================================================#
					  #SWITCHBOARD AND SETTINGS#
#====================================================================#

#ParameterVaried = 'I_Sol_Start' 
#ParameterRange = [x for x in range(500,1501,100)]

#ParameterVaried = 'I_Sol_Equil' 
#ParameterRange = [x for x in range(000,451,50)]

#ParameterVaried = 'I_PF1'
#ParameterRange = [x for x in range(-700,-339,20)]

#ParameterVaried = 'I_PF2'
#ParameterRange = [x for x in range(-700,-399,20)]

#ParameterVaried = 'I_Div1'
#ParameterRange = [x for x in range(000,1501,100)]

#ParameterVaried = 'I_Div2'
#ParameterRange = [x for x in range(600,951,50)]


#ParameterVaried = 'TauP'
#ParameterRange = [x/1000.0 for x in range(10,31,2)]

#ParameterVaried = 'Ip'
#ParameterRange = [x for x in range(20000,40000,2000)]

#ParameterVaried = 'Gr_Frac'
#ParameterRange = [x/100.0 for x in range(10,27,2)]

#ParameterVaried = 'Z_eff'
#ParameterRange = [1.0,2.0,11.85]	#H, H2, Ar8+  (Old Resistivity ~34)


#ParameterVaried = 'Rgeo'
#ParameterRange = [x/10.0 for x in range(40,50,1)]

#ParameterVaried = 'Kappa'
#ParameterRange = [x/10.0 for x in range(10,20,1)]

#ParameterVaried = 'delta'
#ParameterRange = [x/100.0 for x in range(10,31,2)]

#Define feedback plasma geometry
#efit_Geometry = [RGeo, ZGeo, a, Kappa, delta]
#efit_Geometry = [0.44, 0.0, 0.44/1.85, 1.8, 0.2]

####################

#Define FIESTA namelist and project directory names
FIESTAName = 'SMART_SJD.m'			#Define name of FIESTA script
ProjectName = 'SMARTxs-P1'			#Defult global project name
SeriesName = 'auto'					#Parameter scan series name ('auto' for automatic)

#Define simulation name structure
SimNameList = ['BT','TauP','I_Sol_Start','I_PF1','I_PF2','I_Div1','I_Div2']

#Define if simulations are to be run
IAutorun = True				#Run requested simulation series
IParallel = False			#Enable mutli-simulations in parallel
IVerbose = True				#Verbose terminal output - not compatable with IParallel

#Define equilibrium calculation method
IEquilMethod = 'efit'					#Define equil method: 'standard','efit','feedback'
IefitCoils = ['PF1','PF2']				#Define coils for which efit, feedback is applied

#Define paramters to be varied and ranges to be varied over
ParameterVaried = 'I_Sol_Start'		 	#Define parameter to vary - Required for diagnostics
ParameterRange = [700,725,750,775,800]	 #Define range to vary over

#Define which diagnostics are to be performed
savefig_PlasmaCurrent = True		#Plots plasma current trends
savefig_CoilCurrents = True		#Plots maximum dI/dt in each coil

savefig_EquilTrends = True			#Plots general equilibrium trends from Param(equil)
savefig_EquilSeperatrix = False		#Plots seperatrix extrema [Rmin,Rmax,Zmin,ZMax] trends
savefig_IquilMidplane = False		#Plots 2D Radial slice at Z=0 trends
savefig_EquilXpoint = False			#Plots X-point location (R,Z) trends

#Image overrides and tweaks
Image_TrendAxisOverride=''			#Force trend figures to use different variable

#====================================================================#
#====================================================================#





#====================================================================#
#====================================================================#

#TO DO
#IMMEDIATE FIXES
#Enable concurrent diagnostic use after running simulations - requires running twice atm
# ^^^ ?This is likely due to being in the wrong directory after simulations finish? ^^^
#Enable auto-detection of output folders so user doesn't have to change parameter variable


#CORE FUNCTIONALITY
#Unify PyESTA namelist inputs with .m namelist inputs 		- Ideally in external namelist file
#Add capability to iterate towards fixed equilibrium conditions - i.e. iterate on single variable
#Rename all FIESTA output text files in unified format      - Enable CSV or Row-Wise data storing
#Save all PyESTA output data in seperate output folder      - Enable CSV or Row-Wise data storing
#Add ability to change multiple variables per run
#Add ability to use multiple cores, including safety		- NEEDS TESTING!!!

#DIAGNOSTICS
#Add breakdown diagnostic calculating path length at each point in the current waveform
#Add breakdown diagnostic calculating breakdown boolian for each simuatlion (true/false)
#Add equilibrium Rmin,Rmax, Zmin,ZMax diagnostic showing extrema of seperatrix
#Add equilibrium 'PROES-like' diagnostics - plot Radial slice at Z=0 with parameter range
#Add equilibrium X point diagnostic, showing X-point location (R,Z) trends
#Add trend diagnostics for all of the default param(equil) outputs

#ERROR HANDLING
#Add ability to safely-eject from matlab if convergence fails 				- IMPORTANT!!!
#Add ability to produce 'nan' data files if one-or-more simulations fail	- IMPORTANT!!!
#Add general error messages to aid in debugging as the program grows larger

#EXTRA IDEAS
#Add meta-convergence function, enabling self iteration towards a fixed input equilibrium
#This will require feedback between PyESTA and the output files: [Param(equil), ???, ???]


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
	elif Verbose == False and Parallel == False: FIESTA_Output = ' > Conv.out'
	elif Verbose == False and Parallel == True:  FIESTA_Output = ' > Conv.out &'
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

def ExtractFIESTAData(SeriesSubDirs,DataFileName,Dimension='2D',Orientation='Vertical'):

	#Create any required arrays for data storage and record HomeDir for navigation
	GlobalDataArray,OrderedDataArrays = list(),list()
	HomeDir = os.getcwd()

	#For all simulation directories in the requested simulation series
	for i in range(0,len(SeriesSubDirs)):
		#cd into the relevent directory and extract the data
		os.chdir(SeriesSubDirs[i]+'/RawData/')
		GlobalDataArray.append(ReadDataFromFile(DataFileName,Dimension='2D',Orientation='Vertical'))
	#endfor
	#cd back into PyESTA directory for continuity
	os.chdir(HomeDir)

	#Reformat GlobalDataArray to enable easy splitting of column-wise variables
	for i in range(0,len(GlobalDataArray[0])): OrderedDataArrays.append(list())
	#GlobalDataArray organized as [Folder][Variable][Value]
	#OrderedDataArrays organised as [Variable][Folder][Value]
	for i in range(0,len(OrderedDataArrays)):
		for j in range(0,len(GlobalDataArray)): 
			OrderedDataArrays[i].append(GlobalDataArray[j][i])
		#endfor
	#endfor

	#Return ordered data arrays
	return(OrderedDataArrays)
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
print '                                                               V0.1.0'
print '---------------------------------------------------------------------'
print ''
print 'The following diagnostics were requested:'
print '-----------------------------------------'
if IAutorun == True:
	print'# Simulation Series Autorun'
	print''
if True in [savefig_EquilTrends,savefig_EquilSeperatrix,savefig_IquilMidplane,savefig_EquilXpoint]:
	print'# 2D Equilibrium Analysis'
if True in [savefig_PlasmaCurrent]:
	print'# 1D Plasma Current Analysis'
if True in [savefig_CoilCurrents]:
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
	SeriesDirString = '/'+SeriesName+'_'+ProjectName+'/'
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
		print 'Series Convergence Achieved:',str(TimeConv)+'[s]'
		print '------------------------------------------------'
	#endfor
#endif

#=====================================================================#
#=====================================================================#



















#====================================================================#
					   #ANALYSIS AND DIAGNOSTICS#
#====================================================================#

#====================================================================#
				       #EQUILIBRIUM DIAGNOSTICS#
#====================================================================#

#Plot general equilibrium trends from Param(equil)
if savefig_EquilTrends == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+'_'+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract plasma current data from series directories
	ParamEquil = ExtractFIESTAData(SimulationDirs,'EquilParam.txt','2D','Vertical')[0]
	ValueEquil = ExtractFIESTAData(SimulationDirs,'EquilParam.txt','2D','Vertical')[1]

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

	#Quick and dirty removal of most useful trends
	RGeo,ZGeo,Kappa,Epsilon,delta = list(),list(),list(),list(),list()	#efit params
	for l in range(0,len(SimulationDirs)):
		RGeo.append( ValueEquil[l][43] )		#Geomoetric Radial Length 	[m]
#		ZGeo.append( ValueEquil[l][44] )		#Geometric Axial Length 	[m]
		Kappa.append( ValueEquil[l][21] )		#Elongation 				[-]
		Epsilon.append( ValueEquil[l][16] )		#Aspect Ratio 				[-]
		delta.append( ValueEquil[l][25] )		#Triangularity (average) 	[-]
	#endfor

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(1, figsize=(10,8))

	#Plot requested equil parameter trends over full simulation series
	ax.plot(TrendAxis,RGeo,'ko-', ms=12, lw=2)
	ax.plot(TrendAxis,Kappa,'r^-', ms=12, lw=2)
	ax.plot(TrendAxis,Epsilon,'bs-', ms=12, lw=2)
	ax.plot(TrendAxis,delta,'mh-', ms=12, lw=2)

	Legend = ['RGeo', 'Elongation', 'Aspect Ratio', 'Triangularity']
	ax.legend(Legend, fontsize=16, frameon=False)
	ax.set_ylabel('Equilibrium Parameter Trends [-]', fontsize=25)
	ax.set_xlabel('Varied Parameter: '+Parameter, fontsize=25)
	ax.tick_params(axis='x', labelsize=20)
	ax.tick_params(axis='y', labelsize=20)
#	ax.set_xlim(0.00,1.00)
#	ax.set_ylim(0.00,1.00)

	plt.tight_layout()
	plt.savefig(SeriesDirString+'/Equil_Trends.png')
#	plt.show()
	plt.close('all')

	print'-----------------------------'
	print'# Equilibrium Trends Complete'
	print'-----------------------------'
#endif

#=====================================================================#
#=====================================================================#



























#====================================================================#
					  #PLASMA CURRENT DIAGNOSTIC#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_PlasmaCurrent == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+'_'+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract plasma current data from series directories
	Time_Arrays = ExtractFIESTAData(SimulationDirs,'Ip.txt','2D','Vertical')[0]
	Ip_Arrays = ExtractFIESTAData(SimulationDirs,'Ip.txt','2D','Vertical')[1]

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
				  #COIL CURRENT WAVEFORM DIAGNOSTIC#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_CoilCurrents == True:

	#Obtain simulation folder directories for project and requested series
	SeriesDirString = SeriesName+'_'+ProjectName
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)

	#Extract coil currents and time axis from series directories
	Filename = 'CoilCurrents.txt'
	ISol_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	IPF2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	IPF3_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	IDiv1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	IDiv2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]
	Filename = 't.txt'
	Time_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]

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
			IPF2_Arrays[i][j] = IPF2_Arrays[i][j]/(1000.0*24)
			IPF3_Arrays[i][j] = IPF3_Arrays[i][j]/(1000.0*24)
			IDiv1_Arrays[i][j] = IDiv1_Arrays[i][j]/(1000.0*24)
			IDiv2_Arrays[i][j] = IDiv2_Arrays[i][j]/(1000.0*24)
		#endfor
	#endfor

	#Calculate dI/dt for each coil set
	DeltaIPF2,DeltaIPF3 = list(),list()
	DeltaIDiv1,DeltaIDiv2 = list(),list()
	DeltaISol = list()
	for i in range(0,len(ISol_Arrays)):
		DeltaISol.append(list())
		DeltaIPF2.append(list())
		DeltaIPF3.append(list())
		DeltaIDiv1.append(list())
		DeltaIDiv2.append(list())
		for j in range(1,len(ISol_Arrays[i])):
			Delta_t = Time_Arrays[i][j]-Time_Arrays[i][j-1]
			#
			DeltaISol[i].append( (ISol_Arrays[i][j]-ISol_Arrays[i][j-1])/Delta_t )
			DeltaIPF2[i].append( (IPF2_Arrays[i][j]-IPF2_Arrays[i][j-1])/Delta_t )
			DeltaIPF3[i].append( (IPF3_Arrays[i][j]-IPF3_Arrays[i][j-1])/Delta_t )
			DeltaIDiv1[i].append( (IDiv1_Arrays[i][j]-IDiv1_Arrays[i][j-1])/Delta_t )
			DeltaIDiv2[i].append( (IDiv2_Arrays[i][j]-IDiv2_Arrays[i][j-1])/Delta_t )
		#endfor
	#endfor

	#Calculate maximum change in current experienced for each coil set
	MaxDeltaIPF2,MaxDeltaIPF3 = list(),list()
	MaxDeltaIDiv1,MaxDeltaIDiv2 = list(),list()
	MaxDeltaISol = list()
	for i in range(0,len(DeltaISol)):
		MaxDeltaISol.append( max(DeltaISol[i], key=abs) )
		MaxDeltaIPF2.append( max(DeltaIPF2[i], key=abs) )
		MaxDeltaIPF3.append( max(DeltaIPF3[i], key=abs) )
		MaxDeltaIDiv1.append( max(DeltaIDiv1[i], key=abs) )
		MaxDeltaIDiv2.append( max(DeltaIDiv2[i], key=abs) )
	#endfor

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(Image_TrendAxisOverride) > 0: Parameter = Image_TrendAxisOverride
	else: Parameter = ParameterVaried
	#endif

	#Create output folder for all coil current timetraces
#	TimeTracesDir = CreateNewFolder(SeriesDirString,'/ICoil_TimeTraces/')
	#For every simulation folder in the current series:
	for l in range(0,len(ISol_Arrays)):

		#Create figure for each Coil Ramp Time Trace diagnostic
		fig,ax = plt.subplots(2, figsize=(12,14), sharex=True)

		#Plot each coil current with respect to time
		ax[0].plot(Time_Arrays[l],ISol_Arrays[l], 'k-', lw=2)
		ax[0].plot(Time_Arrays[l],IPF2_Arrays[l], 'r-', lw=2)
		ax[0].plot(Time_Arrays[l],IPF3_Arrays[l], 'b-', lw=2)
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
		ax[1].plot(Time_Arrays[l][1::],DeltaIPF2[l], 'r-', lw=2)
		ax[1].plot(Time_Arrays[l][1::],DeltaIPF3[l], 'b-', lw=2)
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
#		plt.savefig(TimeTracesDir+'/CoilRamp_'+str(TrendAxis[l])+'TimeTrace.png')
		plt.savefig(SeriesDirString+'/CoilRamp_'+str(TrendAxis[l])+'TimeTrace.png')

		plt.show()
		plt.close('all')
	#endfor
	print'----------------------------------'
	print'# Coil Current Timetraces Complete'
	print'----------------------------------'

	#=================#

	#Create image limits
	GlobalMaxDelta = MaxDeltaISol+MaxDeltaIPF2+MaxDeltaIPF3+MaxDeltaIDiv1+MaxDeltaIDiv2
	Ylims = [min(GlobalMaxDelta),max(GlobalMaxDelta)]

	#Create figure for Coil Maximum Ramp Diagnostic
	fig,ax = plt.subplots(1, figsize=(12,8))

	#Plot derivitive of each coil current with respect to time
	ax.plot(TrendAxis,MaxDeltaISol, 'ko-', ms=10, lw=2)
	ax.plot(TrendAxis,MaxDeltaIPF2, 'r^-', ms=10, lw=2)
	ax.plot(TrendAxis,MaxDeltaIPF3, 'bs-', ms=10, lw=2)
	ax.plot(TrendAxis,MaxDeltaIDiv1, 'c*-', ms=10, lw=2)
	ax.plot(TrendAxis,MaxDeltaIDiv2, 'mh-', ms=10, lw=2)

	ax.set_title('Maximum Delta Coil Current for Varying '+Parameter, fontsize=20, y=1.03)
	Legend = ['Sol','PF1','PF2','Div1','Div2']
	ax.legend(Legend, fontsize=22, ncol=2, frameon=False)
	ax.set_ylabel('Maximum Change in \n Current $\Delta I$ [kA ms$^{-1}$]', fontsize=25)
	ax.set_xlabel(Parameter, fontsize=25)
#	ax.xaxis.set_major_locator(ticker.MultipleLocator( (max(TrendAxis)-min(TrendAxis))/5 ))
#	ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
	ax.tick_params(axis='x', labelsize=20)
	ax.tick_params(axis='y', labelsize=20)
	ax.set_xlim(TrendAxis[0],TrendAxis[-1])		
	ax.set_ylim(Ylims[0],Ylims[1]*1.25)

	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(SeriesDirString+'/CoilRamp_Trends.png')
	plt.show()
	plt.close('all')

	print'--------------------------------'
	print'# Coil Ramp Diagnostics Complete'
	print'--------------------------------'
#endif

#=====================================================================#
#=====================================================================#




































