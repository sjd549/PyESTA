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
				  		 #LOW LEVEL INPUTS#
#====================================================================#

#Various debug and streamlining options.
DebugMode = False					#Produces debug outputs for relevent diagnostics.

#Warning suppressions
np.seterr(divide='ignore', invalid='ignore')		#Suppresses divide by zero errors
#Fix "can't invoke "event" command: application has been destroyed" error with PROES images
#Fix "Exception KeyError: KeyError(<weakref at 0x7fc8723ca940; to 'tqdm' at 0x7fc85cd23910>,)" error

#Set default parallelisation options
NumThreads = 2

####################

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








#====================================================================#
					  #SWITCHBOARD AND SETTINGS#
#====================================================================#

#Define local directory of FIESTA script
FIESTA_Dir = 'SMART_SJD.m'


#Define Parameter variation
VariableName = 'ZScale'
VariableRange = [0.8,0.85,0.9,0.95,1.0]




#====================================================================#
#====================================================================#



#TO DO

#CORE FUNCTIONALITY
#Add option to skip simulations and just run the diagnostics
#Add ability to change all geometric and coil variables
#Add ability to change multiple variables per run
#Add ability to use multiple cores

#DIAGNOSTICS
#Add equilibrium Rmin,Rmax, Zmin,ZMax diagnostic showing extrema of seperatrix
#Add equilibrium 'PROES-like' diagnostics - plot Radial slice at Z=0 with parameter range
#Add equilibrium X point diagnostic, showing X-point location (R,Z) trends
#Add breakdown diagnostics showing criterion and trends with breakdown time if possible
#Add coil current diagnostics, overlaying coil current profiles for chosen coils
#Add trend diagnostics for all of the default param(equil) outputs

#ERROR HANDLING
#Add ability to safely-eject from matlab if convergence fails
#Add ability to produce 'nan' data files if one-or-more simulations fail
#Add general error messages to aid in debugging as the program grows larger

#EXTRA IDEAS
#Add meta-convergence function, enabling self iteration towards a fixed output param(equil) variable


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

def ReadDataFromFile(Filename,Dimension='2D',Orientation='Horizontal'):
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

	#Lowest dimention is scalar: ==> 1D array.
	#Orientation doesn't matter if 1D.
	elif Dimension == '1D':
		#Read in 1D data from ASCII formatted file.
		datafile = open(Filename)
		Row = datafile.readline().split()
		for m in range(0,len(Row)):
			OutputData.append(float(Row[m]))
		#endfor
	#endif

	return(OutputData)
#enddef

#=========================#

#Takes namelist directory and namelist variable and value
#Locates namelist entry for variable and alters to new value
#Returns namelist initial value and modified namelist entry for sanity checking
#Example: Init,Entry = AlterNamelistVariable(FIESTA_Dir,VariableName,VariableValue)
def AlterNamelistVariable(Namelist_Dir,VariableName,VariableValue):

	#Open namelist file and identify the requested variable name, line index and init value
	Namelist = open(Namelist_Dir).readlines()
	NamelistEntry = filter(lambda x:VariableName in x, Namelist)[0]
	NamelistIndex = Namelist.index(NamelistEntry)
	NamelistInitValue = float(NamelistEntry.partition(';')[0].strip(VariableName+' \t\n\r,='))
	
	#Seperate variable string into 5 sub-strings of order: 'Variable,Sep2,Value,Sep,Comment'
	#Assumes single value input terminated by a semi-colon - allows for and retains comments
	Input, Sep, Comment = NamelistEntry.partition(';')
	Variable, Sep2, InitValue = Input.partition('=')
	InitValue = str(VariableValue)
	#Reconstruct the altered namelist entry with updated init value
	AlteredNamelistEntry = Variable+Sep2+InitValue+Sep+Comment

	#Replace namelist entry with altered namelist entry
	Namelist[NamelistIndex] = AlteredNamelistEntry

	#Write reconfigured namelist back into Namelist_Dir
	with open(Namelist_Dir, 'w') as file:
		file.writelines( Namelist )
	#endwith

	return(NamelistInitValue,AlteredNamelistEntry)
#enddef

#=========================#


#=====================================================================#
#=====================================================================#








#====================================================================#
					   #FIESTA AUTORUN ROUTINES#
#====================================================================#

#Construct terminal command to run requested version of FIESTA
#Example: matlab -nodisplay -nosplash -nodesktop -r "run('/path/to/FIESTA_Script');exit;"
FIESTA_RootDir = os.getcwd()+"/SMART_SJD.m"
SplashHandler = '-nodisplay -nosplash -nodesktop -r '
RunFIESTA = '\"run(\''+FIESTA_RootDir+'\');exit;\"'
TerminalCommand = 'matlab '+SplashHandler+RunFIESTA

#Run FIESTA over requested parameter variation
for i in range(0,len(VariableRange)):
	#Modify FIESTA namelist for requested variable in value
	VariableInit,NewEntry = AlterNamelistVariable(FIESTA_Dir,VariableName,VariableRange[i])

	#Run FIESTA, saving data to appropriate output files
	os.system( TerminalCommand )

#	#Return FIESTA namelist to default state before continuing (for safety)
#	AlterNamelistVariable(FIESTA_Dir,VariableName,VariableInit)
#endfor

#============================#

















#====================================================================#
					  #ANALYSIS AND DIAGNOSTICS#
#====================================================================#


Axes = list()
Data = list()
Filename = os.getcwd()+'/CurrentRun_SMART-V3p1/Scale Coil Geometry/'
#
Ext = 'ZScale#0.80, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/Ip_Phase_1.txt'
Data.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
Ext = 'ZScale#0.80, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/t_Phase_1.txt'
Axes.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
#
Ext = 'ZScale#0.85, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/Ip_Phase_1.txt'
Data.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
Ext = 'ZScale#0.85, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/t_Phase_1.txt'
Axes.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
#
Ext = 'ZScale#0.90, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/Ip_Phase_1.txt'
Data.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
Ext = 'ZScale#0.90, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/t_Phase_1.txt'
Axes.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
#
Ext = 'ZScale#0.95, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/Ip_Phase_1.txt'
Data.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
Ext = 'ZScale#0.95, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/t_Phase_1.txt'
Axes.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
#
Ext = 'ZScale#1.00, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/Ip_Phase_1.txt'
Data.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )
Ext = 'ZScale#1.00, BT#0.1, Tau#0.02, Sol#1500, PF2#-1175, Div2#4440,/t_Phase_1.txt'
Axes.append( ReadDataFromFile(Filename+Ext,Dimension='2D',Orientation='Vertical')[0] )


####################

#Rescale current data from [A] to [kA] for plotting
for i in range(0,len(Data)):
	for j in range(0,len(Data[i])):
		Data[i][j] = Data[i][j]/1000.0
	#endfor
#endfor

#Calculate maximum plasma current trend for each operating condition
Data_MaxTrend,Data_MinTrend = list(),list()
for i in range(0,len(Data)):
	Data_MaxTrend.append(max(Data[i]))
#endfor
TrendAxis = [0.80,0.85,0.90,0.95,1.00]


#=====================================================================#
#=====================================================================#

fig,ax = plt.subplots(figsize=(11,7))

#Plot waveform shape
ax.plot(Axes[0],Data[0], 'k-', lw=2)
ax.plot(Axes[1],Data[1], 'r-', lw=2)
ax.plot(Axes[2],Data[2], 'b-', lw=2)
ax.plot(Axes[3],Data[3], 'c-', lw=2)
ax.plot(Axes[4],Data[4], 'm-', lw=2)
#####
ax.legend(['$Z_{Scale}$ 0.80','$Z_{Scale}$ 0.85','$Z_{Scale}$ 0.90','$Z_{Scale}$ 0.95','$Z_{Scale}$ 1.00'], fontsize=22, loc=2, frameon=False)
ax.set_ylabel('Plasma Current $I_{p}$ [kA]', fontsize=25)
ax.set_xlabel('Time $\\tau$ [s]', fontsize=25)
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#ax.yaxis.set_major_locator(ticker.MultipleLocator(240))
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_xlim(-0.02,0.08)		
ax.set_ylim(-2,32)

from mpl_toolkits.axes_grid.inset_locator import inset_axes
left, bottom, width, height = [0.50, 0.27, 0.25, 0.25]
ax2 = fig.add_axes([left, bottom, width, height])
###
ax2.plot(TrendAxis,Data_MaxTrend,'ko--', ms=5, lw=1.5)
###
#ax2.legend(['$\epsilon_{s}$','$\eta_{dc}/V_{\mathrm{pp}}$'], loc=6, fontsize=10, frameon=False)
ax2.set_ylabel('Maximum Plasma \n Current $I_{p,max}$ [kA]', labelpad=0, fontsize=14.5)
ax2.set_xlabel('Coil Scaling Parameter $Z_{Scale}$ [-]', fontsize=15)
#ax2.xaxis.set_major_locator(ticker.MultipleLocator(90))
#ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax2.set_ylim(22,27)
ax2.set_xlim(0.79,1.01)


plt.tight_layout(pad=3.0,h_pad=1.0)
#plt.savefig('PlasmaCurrentTrends.eps')
#plt.savefig('PlasmaCurrentTrends.pdf')
plt.show()
plt.close('all')

#=====================================================================#
#=====================================================================#



























