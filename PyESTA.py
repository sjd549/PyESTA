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
NumThreads = 2				#Maximum number of threads avaliable to FIESTA
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

####################


#######################  DEFINE VESSEL GEOMETRY  #######################

#Define Vessel Outer Geometry
VesselRInnerPoint=0.15; # R min position [m]
VesselROuterPoint=0.8;  # R max position [m]
VesselZMinPoint=-0.8;   # Z min position [m]
VesselZMaxPoint=0.8;    # Z max position [m]

#Define Solenoid Geometry and Parameters
nSol=800;     # number of turns of the solenoid
RSol=0.09;    # R position of the solenoid [m] (Inner Solenoid)
ZMinSol=-0.8; # Min Z position
ZMaxSol=0.8;  # Max Z position

#Define Div and PF coils
nZDiv1=6;
nRDiv1=4;
nZDiv2=6;
nRDiv2=4;
#nZPF1=6;
#nRPF1=4;
nZPF2=6;
nRPF2=4;
nZPF3=6;
nRPF3=4;

#Calculate total number of turns in each coil
nDiv1=nZDiv1*nRDiv1;
nDiv2=nZDiv2*nRDiv2;
#nPF1=nZPF1*nRPF1;
nPF2=nZPF2*nRPF2;
nPF3=nZPF3*nRPF3; 

#Define coil size to enable cross-section calculation
width_PF=0.042;  # Width of a turn (m)
height_PF=0.035; # Height of a turn (m)

#Define central location of coil sets
RScale=1.00;        #Scale all radial coil positions (Relative to 2.0m)
ZScale=0.80;        #Scale all axial coil positions (Relative to 2.0m)
#####
#R_PF1=0.9*RScale;  #R position of PF1 (m)
#Z_PF1=0.3*ZScale;  #Z position of PF1 (m)
R_PF2=0.9*RScale;   #R position of PF2 (m)
Z_PF2=0.5*ZScale;   #Z Position of PF2 (m)
R_PF3=0.9*RScale;   #R Position of PF3 (m)
Z_PF3=0.8*ZScale;   #Z Position of PF3 (m)
R_Div1=0.25*RScale; #R Position of Div1 (m)
Z_Div1=1.05*ZScale; #Z Position of Div1 (m)
R_Div2=0.55*RScale; #R Position of Div2 (m)
Z_Div2=1.05*ZScale; #Z Position of Div2 (m)


#######################  DEFINE INITIAL PARAMETERS  #######################

#Define any required constants
mu0 = 1.2566e-06; #Magnetic Moment [I/m^2]

#Define Operating Conditions
RGeo = 0.4763		# Geometrical Radius [m]
epsilon = 1.985		# Aspect ratio       [-]
a = RGeo/epsilon	# Minor radius       [m]
kappa = 1.7			# Elongation         [-]
Ip = 30e3			# Plasma current     [A]
li2 = 1				#                    [-]
#q_cyl = 2.821		# Safety Factor      [-]
betaN = 3.529		#                    [-] (Obtained via VEST Excel)
TauP = 0.020		# Pulse length       [s] (Also determines tstep for Ip plot)\r\n')

#Compute Further Operating Conditions
BT=0.1								# Toroidal B-Field            [T]
Irod=BT*2*pi*RGeo/mu0				# Central Rod Current         [kA]
S=np.sqrt( (1.0+kappa**2)/2.0 )		# Shaping factor              [-]
betaT = betaN/a*(Ip/1e6)/BT			# Beta toroidal               [%]
betaP = 25.0/(betaT/100)*(S**2)*((betaN/100)**2) # Beta Poloidal  [%]
nT = 2.66*1e20*1e3*betaT*BT**2		# Density Temperature Product [eV m-3]
n = 3e19							# Central Plasma Density      [m^-3
T = nT/n							# Central Plasma Temperature  [eV]
nGW = Ip*1e14/(pi*a**2)				# ??Central Current Density?? [kA m^2]
BZ = -mu0*Ip/(4*pi*RGeo)*( log(8*epsilon)+betaP+0.5*li2-3/2 ) #Vertical field [T]
#
coil_density = 1;					# Coil Density                [Arb]
coil_temp = 293.0;					# Coil Init Temperature       [K]


###################  DEFINE SOL RAMP & COIL CURRENTS  ###################

#Phase1 coil currents [kA]                  %SJDoyle
I_Sol_start=1300;      #+1300 -> +1500;     #+1300;
I_Sol_ramp=-900;       #-0900 -> -1100;     #-1100;
I_Sol_equil=-1300;     #-1300 -> -1500;     #-1300;
#
I_PF1=0;               #-0    -> -0         #N/A
I_PF2=-900;            #-0900 -> -1200;     #-1000;
I_PF3=-800;            #-0800 -> -0900;     #-0800;
I_Div1=-000;           #-0000 -> -0000;     #+0000;
I_Div2=+3200;          #+3200 -> +4440;     #+3200;

#Phase2 coil currents [kA]
#I_Sol_start=+4700;    #+4700 -> +4700;     #+4700;
#I_Sol_ramp=+500;      #+500  -> +500       #+500;
#I_Sol_equil=-4700;    #-4700 -> -4700;     #-4700;
#
#I_PF1=0;              #-0    -> -0         #N/A
#I_PF2=-3000;          #-3000 -> -3000      #-3000;
#I_PF3=-940;           #-940  -> -940       #-940;
#I_Div1=-000;          #-000  -> -000       #-000;
#I_Div2=+9090;         #+9090 -> +9090      #+9090;

#Phase3 coil currents [kA]
#I_Sol_start=4700;     #4200;
#I_Sol_ramp=500;       #
#I_Sol_equil=-4700;    #-5200;
#
#I_PF1=0;              #-0    -> -0         #N/A
#I_PF2=-6000;          #
#I_PF3=-3100;          #
#I_Div1=-000;          #
#I_Div2=15880;         #

#====================================================================#
#====================================================================#





#====================================================================#
					  #SWITCHBOARD AND SETTINGS#
#====================================================================#

#ParameterList = 'I_PF2' 
#ParameterRanges = [x/1.0 for x in range(-1200,-899,50)] 	

#ParameterList = 'I_Div2'
#ParameterRanges = [x/1.0 for x in range(3200,4401,100)]

#ParameterList = 'TauP'
#ParameterRanges = [x/1000.0 for x in range(10,31,2)]

#ParameterList = 'ZScale'
#ParameterRanges = [0.8,0.85,0.9,0.95,1.0]

####################

#Define FIESTA namelist and project directory names
FIESTAName = 'SMART_SJD.m'			#Define name of FIESTA script
ProjectName = 'SMARTxs-P1'			#Define global project name
SeriesName = 'VaryIPF2'				#Define parameter scan series name

#Define simulation name structure:
SimNameList = ['BT','TauP','I_Sol_start','I_PF2','I_Div2']
#SimNameList = ['ZScale','BT','TauP','I_Sol_start','I_PF2','I_Div2']
#SimNameList = ['Z_PF2','Z_PF3','Z_Div1','Z_Div2']

#Define if simulations are to be run
IAutorun = False					#Run requested simulation series
IParallel = False					#Enable mutli-simulations in parallel
IVerbose = True					#Verbose terminal output - not compatable with IParallel

#Define paramters to be varied and ranges to be varied over
ParameterList = 'I_Div2' #'TauP'
ParameterRanges = [x/1.0 for x in range(3200,4401,100)]#[x/1000.0 for x in range(10,31,2)]

#Define which diagnostics are to be performed
savefig_PlasmaCurrent = True		#Plots plasma current trends
savefig_CoilCurrents = True		#Plots maximum dI/dt in each coil

savefig_EquilSeperatrix = False		#Plots seperatrix extrema [Rmin,Rmax,Zmin,ZMax] trends
savefig_IquilMidplane = False		#Plots 2D Radial slice at Z=0 trends
savefig_EquilXpoint = False			#Plots X-point location (R,Z) trends



#Override diagnostic functions
TrendAxisOverride=''				#Force trend naming routine to recognise trend variable

#====================================================================#
#====================================================================#










#====================================================================#
#====================================================================#

#TO DO

#CORE FUNCTIONALITY
#Rename all output text files in a HPEM style format and make a standardised ASCII format
#Add capability to search for specific data files using a lambda function for data read-in
#Add ability to change multiple variables per run
#Add ability to use multiple cores 											- NEEDS TESTING!!!

#DIAGNOSTICS
#Add equilibrium Rmin,Rmax, Zmin,ZMax diagnostic showing extrema of seperatrix
#Add equilibrium 'PROES-like' diagnostics - plot Radial slice at Z=0 with parameter range
#Add equilibrium X point diagnostic, showing X-point location (R,Z) trends
#Add breakdown diagnostics showing criterion and trends with breakdown time if possible
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
#Example: Value,Entry = AlterNamelistVariable(FIESTAName,ParameterList)
def FindNamelistVariable(Namelist_Dir,ParameterList):

	#Open namelist file and identify the requested variable name, line index and init value
	Namelist = open(Namelist_Dir).readlines()
	NamelistEntry = filter(lambda x:ParameterList in x, Namelist)[0]
	NamelistIndex = Namelist.index(NamelistEntry)
	try: NamelistValue = float(NamelistEntry.partition(';')[0].strip(ParameterList+' \t\n\r,='))
	except: NamelistValue = NamelistEntry.partition(';')[0].strip(ParameterList+' \t\n\r,=')

	return(NamelistValue,NamelistEntry)
#enddef

#=========================#

#Takes namelist directory and namelist variable and value
#Locates namelist entry for variable and alters to new value
#Returns modified namelist entry for sanity checking purposes
#Example: Init,Entry = AlterNamelistVariable(FIESTAName,ParameterList,VariableValue)
def AlterNamelistVariable(Namelist_Dir,ParameterList,VariableValue):

	#Open namelist file and identify the requested variable name, line index and init value
	Namelist = open(Namelist_Dir).readlines()
	NamelistEntry = filter(lambda x:ParameterList in x, Namelist)[0]
	NamelistIndex = Namelist.index(NamelistEntry)
	try: NamelistValue = float(NamelistEntry.partition(';')[0].strip(ParameterList+' \t\n\r,='))
	except: NamelistValue = NamelistEntry.partition(';')[0].strip(ParameterList+' \t\n\r,=')
	
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
	SimulationDirsRaw = os.listdir( os.path.abspath(SeriesDirString) )
	SimulationDirsCleaned = list()
	
	#Remove any non-folder directories in SeriesDirsRaw and correct bash 'grammar'
	for i in range(0,len(SimulationDirsRaw)):

		#Define simulation series directories from root or relative to local directory
		AlwaysRoot = os.getcwd()+'/'+SeriesDirString
		if Root == True: RootDir = os.getcwd()+'/'+SeriesDirString
		else: RootDir = SeriesDirString

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
	SimulationNameString = ''
	
	#For each variable in SimNameList, convert to string and append value
	for i in range(0,len(SimNameList)):

		#Check if named parameter has been varied and use appropriate value
		if SimNameList[i] == VariedParameter:
			ParameterString = SimNameList[i]+'#'+str(ParameterValue)+' '
		else:
			ParameterString = SimNameList[i]+'#'+str(eval(SimNameList[i]))+' '
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
#Example: TrendAxis = CreateTrendAxis(SimulationNames,ParameterList)
def CreateTrendAxis(SimulationNames,VariableString,TrendAxisOverride=''):

	#Create required list to store output
	TrendAxis = list()

	#For all simulation names
	for i in range(0,len(SimulationNames)): 
		#Split each directory folder name into substrings and identify varied parameter
		SimulationNames[i] = SimulationNames[i].split('/')[-1]	#Remove any directories
		SplitSimName = SimulationNames[i].split(' ')			#Split simulation parameters
		for j in range(0,len(SplitSimName)):
			if VariableString in SplitSimName[j]:
				TrendString = SplitSimName[j]
				break
			elif len(TrendAxisOverride) > 0 and TrendAxisOverride in SplitSimName[j]:
				TrendString = SplitSimName[j]
				break
			#endif
		#endfor

		#If identifed; use TrendString. Else; default to first named parameter
		try: 
			TrendValue = TrendString.partition('#')[2]
		except:
			TrendString = SplitSimName[i].partition(',')[0]
			TrendValue = TrendString.partition('#')[2]
		#endtry
		TrendAxis.append(float(TrendValue))
	#endfor

	return(TrendAxis)
#enddef

#=========================#


#=====================================================================#
#=====================================================================#













#====================================================================#
					  #FIESTA AUTORUN ROUTINE#
#====================================================================#

#Autorun simulations over defined paramter range if requested
if IAutorun == True:

	#Create simulation series folder and obtain folder directories
	HomeDir = os.getcwd()
	SeriesDirString = '/'+SeriesName+'_'+ProjectName+'/'
	SeriesDir = CreateNewFolder(HomeDir,SeriesDirString)

	#For all requested input parameters
	for i in range(0,len(ParameterRanges)):

		#Create simulation folder for input parameter[i]
		SimulationString = CreateSimName(SimNameList,ParameterList,ParameterRanges[i])
		SimulationDir = CreateNewFolder(SeriesDir,SimulationString)

		#Copy FIESTA.m into simulation folder and cd into directory
		os.system('cp '+FIESTAName+' '+'\''+SimulationDir+'\'')
		os.chdir(SimulationDir)

		#Update new FIESTA.m with fixed parameters from PyESTA namelist
		MatlabProjectString = '\''+ProjectName+'\''
		MatlabSimulationString = '\''+SimulationString+'\''
		AlteredEntry = AlterNamelistVariable(FIESTAName,'ProjectName',MatlabProjectString)
		AlteredEntry = AlterNamelistVariable(FIESTAName,'SimName',MatlabSimulationString)

		#####

		#Update new FIESTA.m with variable parameter scan for parameter [i]
		AlteredEntry = AlterNamelistVariable(FIESTAName,ParameterList,ParameterRanges[i])

		#Run modified FIESTA - Verbosity determines terminal output.
		RunFIESTA(FIESTAName,Verbose=IVerbose,Parallel=IParallel)
	#endfor
#endif

#=====================================================================#
#=====================================================================#



















#====================================================================#
					   #ANALYSIS AND DIAGNOSTICS#
#====================================================================#

#If parallel simulations have been requested:
if IAutorun == True and IParallel == True:
	#Initial delay to allow MATLAB processes to start
	time.sleep(ConvDelay/2.0)
	TimeConv = ConvDelay/2.0

	#Parallel Race Condition Checker - Waits for all simulations to finish before analysis
	while CheckIfProcessRunning('MATLAB') == True:
		#If MATLAB process is detected, wait ConvDelay seconds and check again
		time.sleep(ConvDelay)
		TimeConv += ConvDelay
		print 'Awaiting Simulation Convergence:',TimeConv+'[s]'
	#endwhile

	print 'Simulation Convergence:',TimeConv+'[s]'
#endfor

#=========================#


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
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterList,TrendAxisOverride)

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

	####################

	#Create figure for plasma current diagnostic
	fig,ax = plt.subplots(1, figsize=(12,10))

	#Plot plasma current with respect to adaptive_time
	for i in range(0,len(Ip_Arrays)): ax.plot(Time_Arrays[i],Ip_Arrays[i], lw=2)
	ax.set_title('Plasma Current Time-Trace for Varying '+ParameterList, fontsize=20, y=1.03)
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
	ax2.set_xlabel('Varied Parameter: '+ParameterList, fontsize=15)
#	ax2.xaxis.set_major_locator(ticker.MultipleLocator(90))
#	ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
	ax2.tick_params(axis='x', labelsize=14)
	ax2.tick_params(axis='y', labelsize=14)
#	ax2.set_xlim( min(TrendAxis),max(TrendAxis)*1.10 )
#	ax2.set_ylim(0.79,1.01)

	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(SeriesDirString+'/Ip_Trends.png')
	plt.show()
	plt.close('all')
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
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterList,TrendAxisOverride)

	#Rescale data for plotting: [s] to [ms]
	for i in range(0,len(Time_Arrays)):
		for j in range(0,len(Time_Arrays[i])):
			Time_Arrays[i][j] = Time_Arrays[i][j]*1000.0
		#endfor
	#endfor

	#Rescale data for plotting: [A] to [kA]
	for i in range(0,len(ISol_Arrays)):
		for j in range(0,len(ISol_Arrays[i])):
			ISol_Arrays[i][j] = ISol_Arrays[i][j]/(1000.0*800)		#ODD SCALING!?
			IPF2_Arrays[i][j] = IPF2_Arrays[i][j]/(1000.0*24)		#ODD SCALING!?
			IPF3_Arrays[i][j] = IPF3_Arrays[i][j]/(1000.0*24)		#ODD SCALING!?
			IDiv1_Arrays[i][j] = IDiv1_Arrays[i][j]/(1000.0*24)		#ODD SCALING!?
			IDiv2_Arrays[i][j] = IDiv2_Arrays[i][j]/(1000.0*24)		#ODD SCALING!?
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

	####################	####################
	####################	####################

	#Create output folder for all coil current timetraces
#	TimeTracesDir = CreateNewFolder(SeriesDirString,'/ICoil_TimeTraces/')
	#For every simulation folder in the current series:
	for l in range(0,len(ISol_Arrays)):

		#Create figure for each Coil Ramp Time Trace diagnostic
		fig,ax = plt.subplots(2, figsize=(12,14), sharex=True)

		#Plot each coil current with respect to time
		ax[0].plot(ISol_Arrays[l], 'k-', lw=2)
		ax[0].plot(IPF2_Arrays[l], 'r-', lw=2)
		ax[0].plot(IPF3_Arrays[l], 'b-', lw=2)
		ax[0].plot(IDiv1_Arrays[l], 'c-', lw=2)
		ax[0].plot(IDiv2_Arrays[l], 'm-', lw=2)

		Range = '['+str(min(ParameterRanges))+' - '+str(max(ParameterRanges))+']'
		ax[0].set_title('Time-Traces of Coil Currents for '+ParameterList+' in '+Range, fontsize=20, y=1.03)
		Legend = ['Sol','PF2','PF3','Div1','Div2']
		ax[0].legend(Legend, fontsize=22, loc=2, ncol=2, frameon=False)
		ax[0].set_ylabel('Coil Current $I$ [kA]', fontsize=25)
#		ax[0].set_xlabel('Time $\\tau$ [ms]', fontsize=25)
#		ax[0].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#		ax[0].yaxis.set_major_locator(ticker.MultipleLocator(240))
		ax[0].tick_params(axis='x', labelsize=20)
		ax[0].tick_params(axis='y', labelsize=20)
#		ax[0].set_xlim(-50,100)		
#		ax[0].set_ylim(2,32)

		#Plot derivitive of each coil current with respect to time
		ax[1].plot(DeltaISol[l], 'k-', lw=2)
		ax[1].plot(DeltaIPF2[l], 'r-', lw=2)
		ax[1].plot(DeltaIPF3[l], 'b-', lw=2)
		ax[1].plot(DeltaIDiv1[l], 'c-', lw=2)
		ax[1].plot(DeltaIDiv2[l], 'm-', lw=2)

		Range = '['+str(min(ParameterRanges))+' - '+str(max(ParameterRanges))+']'
		ax[1].set_title('Time-Traces of Delta Coil Currents for '+ParameterList+' in '+Range, fontsize=20, y=1.03)
		Legend = ['Sol','PF2','PF3','Div1','Div2']
		ax[1].legend(Legend, fontsize=22, loc=2, ncol=2, frameon=False)
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

	####################

	#Create figure for Coil Maximum Ramp Diagnostic
	fig,ax = plt.subplots(1, figsize=(12,8))

	#Plot derivitive of each coil current with respect to time
	ax.plot(TrendAxis,MaxDeltaISol, 'ko-', ms=10, lw=2)
	ax.plot(TrendAxis,MaxDeltaIPF2, 'r^-', ms=10, lw=2)
	ax.plot(TrendAxis,MaxDeltaIPF3, 'bs-', ms=10, lw=2)
	ax.plot(TrendAxis,MaxDeltaIDiv1, 'c*-', ms=10, lw=2)
	ax.plot(TrendAxis,MaxDeltaIDiv2, 'mh-', ms=10, lw=2)

	ax.set_title('Maximum Delta Coil Current for Varying '+ParameterList, fontsize=20, y=1.03)
	Legend = ['Sol','PF2','PF3','Div1','Div2']
	ax.legend(Legend, fontsize=22, loc=1, frameon=False)
	ax.set_ylabel('Maximum Change in \n Current $\Delta I$ [kA ms$^{-1}$]', fontsize=25)
	ax.set_xlabel('Time $\\tau$ [ms]', fontsize=25)
	ax.xaxis.set_major_locator(ticker.MultipleLocator( (max(TrendAxis)-min(TrendAxis))/5 ))
#	ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
	ax.tick_params(axis='x', labelsize=20)
	ax.tick_params(axis='y', labelsize=20)
#	ax.set_xlim(0,0)		
#	ax.set_ylim(0,0)

	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(SeriesDirString+'/CoilRamp_Trends.png')
	plt.show()
	plt.close('all')
#endif

#=====================================================================#
#=====================================================================#




































