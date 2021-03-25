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
import matplotlib
import numpy as np
import scipy as sp
import math as m
import subprocess
import os, sys
import os.path
import time

#Enforce matplotlib to avoid instancing undisplayed windows
#matplotlib-tcl-asyncdelete-async-handler-deleted-by-the-wrong-thread
#matplotlib.use('Agg')			!!! CAUSES FIGURES TO NOT PLOT !!!

#Import additional modules
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import savgol_filter
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from matplotlib import ticker
from scipy import ndimage
from cycler import cycler			#Enables easy modification of mpl.rcParams.color()
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

#====================================================================#
					  #SWITCHBOARD AND SETTINGS#
#====================================================================#

#Requested simulation series:
SeriesDirString = "S3-000005"

#Requested Equil Variables
TrendAxisVariables=''				#Force trend figures to use different variable		[BREAKS!!!]

#Various Diagnostic inputs
PaschenPressure = [1E-7,1E-2]		#Operating pressure range in Torr [Min,Max]

#Requested diagnostics and plotting routines.
savefig_1DEquilProfiles = False		#Plots 1D profiles through equilibrium midplane		#TO DO!!!
savefig_2DEquilPlots = False		#Plots 2D images of the target equilibria
savefig_EquilTrends = False			#Plots efit equilibrium geometry trends from Param(equil)

savefig_CoilCurrentTrends = True	#Plots trends in PF coil currents over all simulations
savefig_CoilVoltageTrends = False	#Plots trends in PF coil voltages over all simulations
savefig_PlasmaCurrent = False		#Plots plasma current trends over all simulations
savefig_EddyCurrent = False			#Plots vessel net eddy current trends over all simulations
savefig_Breakdown = False			#Plots Paschen curves and connection length for each simulation

savefig_VerticalStability = False	#Plots vertical stability growth rates and perturbed coil currents


#Image plotting options.
image_extension = '.eps'				#Extensions ('.png', '.jpg', '.eps')
image_aspectratio = [10,10]				#[x,y] in cm [Doesn't rotate dynamically]
image_radialcrop = []					#[R1,R2] in cm
image_axialcrop = []					#[Z1,Z2] in cm
image_cbarlimit = []					#[min,max] colourbar limits	

image_contourplot = True				#Plot contourlines onto 2D images
image_normalise = False					#normalise image/profiles to local max
image_plotgrid = False					#Plot major/minor gridlines on profiles
image_logplot = False					#Plot ln(Data), against linear axis.
image_rotate = False					#Rotate 2D figures 90 degrees to the left

#Image overrides and tweaks
titleoverride = []						#TBC
legendoverride = []						#TBC
xaxisoverride = []						#TBC
xlabeloverride = []						#TBC
ylabeloverride = []						#TBC
cbaroverride = []						#TBC

#====================================================================#
#====================================================================#





#====================================================================#
#====================================================================#

#TO DO
#IMMEDIATE FIXES
#Unify PyESTA namelist inputs with .m namelist inputs 		- Ideally in external namelist file

#CORE FUNCTIONALITY
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
				#DEFINE COMMONLY USED I/O FUNCTIONS#
#====================================================================#

#Takes global inputs from switchboard, returns nothing
#Alters global image options, run before any diagnostics
#Attempts to revert matplotlib changes made in 2.0 onwards.
#See: https://matplotlib.org/users/dflt_style_changes.html
def Matplotlib_GlobalOptions():

#	mpl.style.use('classic')								#Resets to classic 1.x.x format
	
	#Image options			
	mpl.rcParams['figure.figsize'] = [10.0,10.0]				#Sets default figure size
	mpl.rcParams['figure.dpi'] = 100							#Sets viewing dpi
	mpl.rcParams['savefig.dpi'] = 100							#Sets saved dpi
	mpl.rcParams['image.interpolation'] = 'bilinear'			#Applies bilinear image 'smoothing'
	mpl.rcParams['image.resample'] = True						#Resamples data before colourmapping
	mpl.rcParams['axes.prop_cycle'] = cycler(color='krbgcmy')	#Select global line colour cycle
	mpl.rcParams['image.cmap'] = 'plasma'						#Select global colourmap 
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
			if SeriesDirString in Directorylist[i]:
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
			RawData,Header = RawData[1::],RawData[0]

			#Enlarge output data array by number of columns
			NumColumns = len(RawData[0].split())
			for m in range(0,NumColumns):
				OutputData.append(list())
			#endfor

			for i in range(0,len(RawData)):
				#For each column, split row and turn to float
				Row = RawData[i].split()
				for j in range(0,len(Row)):
					try: Row[j] = float(Row[j])
					except: Row[j] = str(Row[j])
				#endfor
				for k in range(0,NumColumns):
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

def ExtractFIESTAData(SeriesSubDirs,DataFilename,Dimension='2D',Orientation='Vertical',Reorder=True):

	#Create any required arrays for data storage and record HomeDir for navigation
	GlobalDataArrays,ReorderedDataArrays = list(),list()
	HomeDir = os.getcwd()

	#For all simulation directories in the requested simulation series
	for i in range(0,len(SeriesSubDirs)):
		#cd into the relevent directory and extract the data
		os.chdir(SeriesSubDirs[i])
		#GlobalDataArray organized as [Folder][Variable][Value]
		GlobalDataArrays.append(ReadDataFromFile(DataFilename,Dimension,Orientation))
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

def ExtractFIESTAEquil(SeriesSubDirs,DataFilename):

	#Create any required arrays for data storage and record HomeDir for navigation
	OutputData,RGrids,ZGrids = list(),list(),list()
	HomeDir = os.getcwd()

	#For all simulation directories in the requested simulation series
	for i in range(0,len(SeriesSubDirs)):
		#cd into the relevent directory and extract the data
		os.chdir(SeriesSubDirs[i]+'/RawData/')
		datafile = open(DataFilename)
		RawData = datafile.readlines()

		#Isolate header and extract grid resolution
		Header = RawData[0].split()
		Date, Unknown = str(Header[1]), int(Header[2])
		RGrid, ZGrid = int(Header[3]), int(Header[4])
		DataLim = ((RGrid*ZGrid)/5)+1

		#Extract equlibrium 2D flux surface, excluding header and footer
		Equil1DArray = list()
		for m in range(1,DataLim):
			Row = RawData[m]
			#Equil data is saved as double precision (%1.16f) in blocks of 5 values
			for n in range(0,5): 
				try: 
					Value = float(Row[16*n:16*(n+1)])
					Equil1DArray.append(Value)
				except: 
					print 'ERROR: Equilibrium data corrupted at line: '+str(m+1)
					break
				#endtry
			#endfor
		#endfor
		#Reshape array into 2D using grid resolution and 'roll' image 20 cells left 
		Equil2DArray = np.reshape(Equil1DArray,(ZGrid,RGrid))
		Equil2DArray = np.roll(Equil2DArray, -20, axis=1)		#Fixes alignment issue

		#Skim values that are beyond threshold limits - improves image contrast
		for i in range(0,5):
			for j in range(0,len(Equil2DArray[i])):
				Equil2DArray[i][j] = 0.0
			#endfor
		#endfor

		OutputData.append(Equil2DArray)
		RGrids.append(RGrid)
		ZGrids.append(ZGrid)
	#endfor

	#cd back into PyESTA directory for continuity
	os.chdir(HomeDir)
	return(OutputData,RGrids,ZGrids)
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
def CreateTrendAxis(SimulationNames,VariableString,TrendAxisVariables=''):

	#Create required list to store output
	TrendAxis = list()

	#For all simulation names
	for i in range(0,len(SimulationNames)): 
		#Split each directory folder name into substrings and identify varied parameter
		SimulationNames[i] = SimulationNames[i].split('/')[-1]	#Remove any directories
		SplitSimName = SimulationNames[i].split(' ')			#Split simulation parameters

		#Find trend variable and extract value - check override variable first
		if len(TrendAxisVariables) > 0:
			TrendString = filter(lambda x: TrendAxisVariables in x, SplitSimName)[0]
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

#=====================================================================#
#=====================================================================#






#====================================================================#
				  #COMMONLY USED PLOTTING FUNCTIONS#
#====================================================================#

#Create figure of desired size and with variable axes.
#Returns figure and axes seperately.
#fig,ax = figure(image_aspectratio,1,shareX=False)
def figure(aspectratio=[],subplots=1,shareX=False):
	if len(aspectratio) == 2:
		fig, ax = plt.subplots(subplots, figsize=(aspectratio[0],aspectratio[1]),sharex=shareX)
	else:
		fig, ax = plt.subplots(subplots, figsize=(10,10), sharex=shareX)
	#endif
	return(fig,ax)
#enddef

#=========================#

#Create figure and plot a 1D graph with associated image plotting requirements.
#Returns plotted axes and figure if new ones were created.
#Else plots to existing figure and returns the image object.
#ImagePlotter1D(Zlineout,Zaxis,image_aspectratio,fig,ax[0]):
def ImagePlotter1D(axis,profile,aspectratio,fig=111,ax=111):

	#Generate new figure if required. {kinda hacky...}
	if fig == 111 and ax == 111:
		fig,ax = figure(aspectratio)
	elif fig == 111:
		fig = figure(aspectratio)
	#endif

	#Apply any required numerical changes to the profile.
	if image_logplot == True:
		profile = np.log(profile)
	if image_normalise == True:
		profile = normalise(profile)[0]
	#endif

	#Plot profile and return.
	im = ax.plot(axis,profile, lw=2)

	try: return(fig,ax,im)
	except: return()
#enddef

#=========================#

#Create figure and plot a 2D image with associated image plotting requirements.
#Returns plotted image, axes and figure after applying basic data restructuring.
#fig,ax,im,Image = ImagePlotter2D(Image,extent,image_aspectratio,variablelist[l],fig,ax[0])
def ImagePlotter2D(Image,extent,aspectratio=image_aspectratio,fig=111,ax=111):

	#Generate new figure if required. {kinda hacky...}
	if fig == 111 and ax == 111:
		fig,ax = figure(aspectratio)
	elif fig == 111:
		fig = figure(aspectratio)
	#endif

	#Rotate image if required
	if image_rotate == True:
		Image = np.asarray(Image)
		Image = Image.transpose().tolist()
	#endif

	#Apply any required numerical changes to the image.
	if image_logplot == True:
		Image = np.log(Image)
	elif image_normalise == True:
		Image = normalise(Image)[0]
	#endif

	#Plot image with or without contour plots, (contour scale = 90% of cbar scale)
	if image_contourplot == True:
		im = ax.contour(Image,extent=extent,origin="lower")
		im.set_clim(np.min(Image)*0.90,np.max(Image)*0.90)
		im = ax.imshow(Image,extent=extent,origin="lower")
	else:
		im = ax.imshow(Image,extent=extent,origin="lower")
	#endif
	return(fig,ax,im,Image)
#enddef

#=========================#

#Applies plt.options to current figure based on user input.
#Returns nothing, open figure is required, use figure().
#For best results call immediately before saving/displaying figure.
#ImageOptions(fig,plt.gca(),Xlabel,Ylabel,Title,Legend,Crop=False,Rotate=False)
def ImageOptions(fig,ax,Xlabel='',Ylabel='',Title='',Legend=[],Crop=False,Rotate=False):

	#Apply user overrides to plots.
	if len(titleoverride) > 0:
		Title = titleoverride
	if len(legendoverride) > 0:
		Legend = legendoverride
	if len(xlabeloverride) > 0:
		Xlabel = xlabeloverride[0]
	if len(ylabeloverride) > 0:
		Ylabel = ylabeloverride[0]
	#endif

	#Set title and legend if one is supplied.
	if len(Title) > 0:
		ax.set_title(Title, fontsize=18, y=1.03)
	if len(Legend) > 0:
		ax.legend(Legend, fontsize=16, frameon=False)		#loc = automatic
	#endif

	#Set labels and ticksize.
	ax.set_xlabel(Xlabel, fontsize=24)
	ax.set_ylabel(Ylabel, fontsize=24)
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)

	#Force scientific notation for all axes, accounting for non-scalar x-axes.
	try: ax.xaxis.get_major_locator().set_params(style='sci',scilimits=(-2,3),axis='both')
	except: Axes_Contain_Strings = True
#	try: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='both')	#Old tickformat.
#	except: ax.ticklabel_format(style='sci',scilimits=(-2,3),axis='y')	#Old tickformat.
	#endtry

	#Set grid, default is off.
	if image_plotgrid == True: ax.grid(True)
	#endif

	#Crop image dimensions, use provided dimensions or default if not provided.
	if isinstance(Crop, (list, np.ndarray) ) == True:
		CropImage(ax,Extent=Crop,Rotate=Rotate)
	elif any( [len(image_radialcrop),len(image_axialcrop)] ) > 0:
		if Crop == True:
			CropImage(ax,Rotate=Rotate)
		#endif
	#endif

	#Arrange figure such that labels, legends and titles fit within frame.
	fig.tight_layout()

	return()
#enddef

#=========================#

#Creates and plots a colourbar with given label and binsize.
#Takes image axis, label string, number of ticks and limits
#Allows pre-defined colourbar limits in form [min,max].
#Returns cbar axis if further changes are required.
#cbar = Colourbar(ax[0],'Label',5,Lim=[0,1])
def Colourbar(ax,im,Label='',Ticks=5,Lim=[]):

	#Set default font and spacing options and modify if required
	Rotation,Labelpad = 270,30
	LabelFontSize,TickFontsize = 24,18
	if '\n' in Label: Labelpad += 25		#Pad label for multi-line names

	#Create and define colourbar axis
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)

	#Set number of ticks, label location and define scientific notation.
	cbar.set_label(Label, rotation=Rotation,labelpad=Labelpad,fontsize=LabelFontSize)
	cbar.formatter.set_powerlimits((-2,3))
	cbar.locator = ticker.MaxNLocator(nbins=Ticks)
	cbar.ax.yaxis.offsetText.set(size=TickFontsize)
	yticks(fontsize=TickFontsize)

	#Apply colourbar limits if specified.  (lim=[min,max])
	if len(Lim) == 2: im.set_clim(vmin=Lim[0], vmax=Lim[1])

	return(cbar)
#enddef

#=========================#

#Creates an invisible colourbar to align subplots without colourbars.
#Takes image axis, returns colourbar axis if further edits are required
#cax = InvisibleColourbar(ax[0])
def InvisibleColourbar(ax='NaN'):
	if ax == 'NaN': ax = plt.gca()

	#Create colourbar axis, ideally should 'find' values of existing cbar! 
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.1)

	#Set new cax to zero size and remove ticks.
	try: cax.set_facecolor('none')				#matplotlib v2.x.x method
	except: cax.set_axis_bgcolor('none')		#matplotlib v1.x.x method
	for axis in ['top','bottom','left','right']:
		cax.spines[axis].set_linewidth(0)
	cax.set_xticks([])
	cax.set_yticks([])

	return(cax)
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
print '                                                               V0.4.0'
print '---------------------------------------------------------------------'
print ''
print 'The following diagnostics were requested:'
print '-----------------------------------------'
if True in [savefig_2DEquilPlots,savefig_EquilTrends]:
	print'# 2D Equilibrium Analysis'
if True in [savefig_PlasmaCurrent,savefig_EddyCurrent]:
	print'# 1D Dynamic Current Analysis'
if True in [savefig_CoilCurrentTrends,savefig_CoilVoltageTrends]:
	print'# 1D Coil Waveform Analysis'
if True in [savefig_Breakdown]:
	print'# 1D Breakdown Analysis'
print '-----------------------------------------'
print ''

#=====================================================================#
#=====================================================================#





















#====================================================================#
				   	   #2D EQUILIBRIUM FIGURES#
#====================================================================#

#Plot 2D equilibria from Equil.txt files
if savefig_2DEquilPlots == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

	#Extract equilibrium data from series directories
	Filename = 'Equil_Data/Equil.txt'
	Equil_Arrays,RGrids,ZGrids = ExtractFIESTAEquil(SimulationDirs,Filename)

	#Create trendaxis from folder names
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

	#===================##===================#
	#===================##===================#


	#Plot each equilibrium		#NEEDS COMPLETING!!!
	for l in range(0,NumFolders):

		Extent = [0.01,1.1, -1.3,1.3]		#[GridSize_R, GridSize_Z]
		fig,ax,im,Image = ImagePlotter2D(Equil_Arrays[l],Extent,image_aspectratio)
#		SMARTVessel()
#		SMARTCoils()

		#Apply figure labels, title, legend and cropping
		Title, Legend = 'SMART Target Equilibrium', TrendAxis
		Xlabel, Ylabel = 'Radius $R$ [m]', 'Height $Z$ [m]'
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False,Rotate=False)
		cbar = Colourbar(ax, im, 'Flux Surface Function $\Phi(R,Z)$',5)
		
		#Hacky, need to update cropping function for 1D figures
		ax.set_xlim(0.00, 1.1)		
		ax.set_ylim(-1.0, 1.0)

		plt.tight_layout(pad=3.0,h_pad=1.0)
		plt.savefig('TargetEquilibrium'+image_extension)
#		plt.show()
		plt.close('all')
	#endfor
#endif

#=====================================================================#
#=====================================================================#



#====================================================================#
				    #EQUILIBRIUM PARAMETER TRENDS#	---- TO BE REFACTORED
#====================================================================#

#Plot general equilibrium trends from Param(equil)
if savefig_EquilTrends == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

	#Extract equilibrium data from series directories
	Filename = 'Equil_Data/EquilParam.txt'
	ParamEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	ValueEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Create trendaxis from folder names
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

	#Quick and dirty removal of most useful trends
	RGeo,ZGeo,Kappa,AspectRatio,delta = list(),list(),list(),list(),list()	#efit params
	for l in tqdm(range(0,len(SimulationDirs))):
		RGeo.append( ValueEquil[l][13] )		#Magnetic Radius		 	[m]
		ZGeo.append( ValueEquil[l][14] )		#Magnetic Axis			 	[m]
		Kappa.append( ValueEquil[l][21] )		#Elongation 				[-]
		AspectRatio.append( ValueEquil[l][16] )	#Aspect Ratio 				[-]
		delta.append( ValueEquil[l][25] )		#Triangularity (average) 	[-]
		#####
		#!!!THESE NEED TO BE EXTRACTED FOR ELI - COMPARE SHIFTS TO EFIT INPUTS!!!#
		rGeo = RGeo[-1]/AspectRatio[-1]			#Minor Radius				[m]
		PlasVolume = ValueEquil[l][35]			#3D Plasma Volume			[m3]
		SurfaceArea = ValueEquil[l][36]			#2D Plasma X-section area	[m2]
		GradShift = RGeo[-1]-0.42				#Shafranov Shift			[m]
		VertShift = ZGeo[-1]-0.00				#Vertical Shift				[m]
	#endfor

	#===================##===================#
	#===================##===================#

#	#Create output folder for all coil trend figures
#	EquilTrendsDir = CreateNewFolder(SeriesDirString,'/Equil_Trends/')
	#Organize figure labelling variables
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
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

	ax[1,0].plot(TrendAxis,AspectRatio,'bs-', ms=12, lw=2)
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
	plt.savefig('Equil_Trends'+image_extension)
#	plt.show()
	plt.close('all')
#endif

#=====================================================================#
#=====================================================================#


#	---- TO BE MERGED WITH DIAGNOSTIC ABOVE
#====================================================================#
				       #USER EQUILIBRIUM TRENDS#
#====================================================================#
#	---- TO BE MERGED WITH DIAGNOSTIC ABOVE
savefig_UserEquilTrends = False
#Plot general equilibrium trends from Param(equil)
if savefig_UserEquilTrends == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

	#Extract equilibrium data from series directories
	Filename = 'Equil_Data/EquilParam.txt'
	ParamEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	ValueEquil = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]

	#Create trendaxis from folder names
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

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
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
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
	plt.savefig('Equil_Trends'+image_extension)
#	plt.show()
	plt.close('all')

	print'----------------------------------'
	print'# Efit Equilibrium Trends Complete'
	print'----------------------------------'
#endif

#=====================================================================#
#=====================================================================#





















#====================================================================#
				  #COIL CURRENT TRENDS DIAGNOSTIC#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_CoilCurrentTrends == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

	#Extract absolute coil currents and time axis from series directories
	Filename = 'icoil_Data/CoilCurrents.txt'
	NumCoils = len(ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical'))-1
	Time_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	ISol_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	IPF1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	IPF2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	IDiv1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]
	IDiv2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[5]

	#Extract delta coil currents and time axis from series directories
	Filename = 'icoil_Data/DeltaCoilCurrents.txt'
	NumCoils = len(ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical'))-1
	DeltaTime_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	DeltaISol_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	DeltaIPF1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	DeltaIPF2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	DeltaIDiv1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]
	DeltaIDiv2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[5]

	#Create trendaxis from folder names
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

	#Rescale absolute coil currents: [A] to [kA]
	for i in range(0,len(ISol_Arrays)):
		for j in range(0,len(ISol_Arrays[i])):	
			ISol_Arrays[i][j] = ISol_Arrays[i][j]/1000.0				#[kA]
			IPF1_Arrays[i][j] = IPF1_Arrays[i][j]/1000.0				#[kA]
			IPF2_Arrays[i][j] = IPF2_Arrays[i][j]/1000.0				#[kA]
			IDiv1_Arrays[i][j] = IDiv1_Arrays[i][j]/1000.0				#[kA]
			IDiv2_Arrays[i][j] = IDiv2_Arrays[i][j]/1000.0				#[kA]
		#endfor
	#endfor

	#Rescale Delta coil currents: [A] to [kA]
	for i in range(0,len(DeltaISol_Arrays)):
		for j in range(0,len(DeltaISol_Arrays[i])):	
			DeltaISol_Arrays[i][j] = DeltaISol_Arrays[i][j]/1000.0		#[kA/ms]
			DeltaIPF1_Arrays[i][j] = DeltaIPF1_Arrays[i][j]/1000.0		#[kA/ms]
			DeltaIPF2_Arrays[i][j] = DeltaIPF2_Arrays[i][j]/1000.0		#[kA/ms]
			DeltaIDiv1_Arrays[i][j] = DeltaIDiv1_Arrays[i][j]/1000.0	#[kA/ms]
			DeltaIDiv2_Arrays[i][j] = DeltaIDiv2_Arrays[i][j]/1000.0	#[kA/ms]
		#endfor
	#endfor

	#Calculate maximum coil current for each coil
	MaxIPF1,MaxIPF2,MaxIDiv1,MaxIDiv2 = list(),list(),list(),list()
	MaxISol,MaxIAvg = list(),list()
	for i in range(0,len(ISol_Arrays)):
		MaxISol.append( max(ISol_Arrays[i], key=abs) )			#[kA]
		MaxIPF1.append( max(IPF1_Arrays[i], key=abs) )			#[kA]
		MaxIPF2.append( max(IPF2_Arrays[i], key=abs) )			#[kA]
		MaxIDiv1.append( max(IDiv1_Arrays[i], key=abs) )		#[kA]
		MaxIDiv2.append( max(IDiv2_Arrays[i], key=abs) )		#[kA]
		Tot = abs(MaxISol[i])+abs(MaxIPF1[i])+abs(MaxIPF2[i])+abs(MaxIDiv1[i])+abs(MaxIDiv2[i])
		MaxIAvg.append( Tot/NumCoils )							#[kA]
	#endfor

	#Calculate maximum change in current experienced for each coil set
	MaxDeltaIPF1,MaxDeltaIPF2 = list(),list()
	MaxDeltaIDiv1,MaxDeltaIDiv2 = list(),list()
	MaxDeltaISol = list()
	for i in range(0,len(DeltaISol_Arrays)):
		MaxDeltaISol.append( max(DeltaISol_Arrays[i], key=abs) )
		MaxDeltaIPF1.append( max(DeltaIPF1_Arrays[i], key=abs) )
		MaxDeltaIPF2.append( max(DeltaIPF2_Arrays[i], key=abs) )
		MaxDeltaIDiv1.append( max(DeltaIDiv1_Arrays[i], key=abs) )
		MaxDeltaIDiv2.append( max(DeltaIDiv2_Arrays[i], key=abs) )
	#endfor

	#===================##===================#
	#===================##===================#

#	#Create output folder for all coil trend figures
#	CurrentTrendsDir = CreateNewFolder(SeriesDirString,'/ICoil_Trends/')
	#Organize figure labelling variables
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
	else: Parameter = ParameterVaried
	#endif

	#For every simulation folder in the current series:
	for i in range(0,NumCoils):

		#Set which coil current timetrace to compare
		Coilset = ['ISol','IPF1','IPF2','IDiv1','IDiv2']
		if i == 0: Coil,ICoil_Arrays = 'ISol',ISol_Arrays
		if i == 1: Coil,ICoil_Arrays = 'IPF1',IPF1_Arrays
		if i == 2: Coil,ICoil_Arrays = 'IPF2',IPF2_Arrays
		if i == 3: Coil,ICoil_Arrays = 'IDiv1',IDiv1_Arrays
		if i == 4: Coil,ICoil_Arrays = 'IDiv2',IDiv2_Arrays

		#Plot each coil current with respect to time
		fig,ax = figure(image_aspectratio,1,shareX=False)
		for j in tqdm(range(0,len(ICoil_Arrays))):
			ImagePlotter1D(Time_Arrays[j],ICoil_Arrays[j],image_aspectratio,fig,ax)
		#endfor

		#Apply figure labels, title, legend and cropping
		Range = '['+str(min(TrendAxis))+' - '+str(max(TrendAxis))+']'
		Title = 'Time-Traces of '+Coil+' Coil Currents for '+Parameter+' in '+Range
		Xlabel, Ylabel = 'Coil Current $I$ [kA]', 'Time $\\tau$ [ms]'
		Legend = TrendAxis
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False,Rotate=False)

		plt.tight_layout(pad=3.0,h_pad=1.0)
		plt.savefig(SeriesDirString+'/'+Coil+'_Current_Trends'+image_extension)
#		plt.show()
		plt.close('all')
	#endfor

	#===================##===================#
	#===================##===================#

	#Create image limits
	GlobalMaxDelta = MaxDeltaISol+MaxDeltaIPF1+MaxDeltaIPF2+MaxDeltaIDiv1+MaxDeltaIDiv2
	Ylims = [min(GlobalMaxDelta),max(GlobalMaxDelta)]

	#Create fig and plot...
	fig,ax = figure(image_aspectratio,2,shareX=False)
	#...ax[0] coil currents for each simulation...
	ImagePlotter1D(TrendAxis,MaxISol,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxIPF1,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxIPF2,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxIDiv1,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxIDiv2,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxIAvg,image_aspectratio,fig,ax[0])
	ax[0].plot(TrendAxis[MaxIAvg.index(min(MaxIAvg))],min(MaxIAvg), 'kv', ms=14, lw=2.0) #Min Avg
	#...and [1] delta coil currents for each simulation
	ImagePlotter1D(TrendAxis,MaxDeltaISol,image_aspectratio,fig,ax[1])
	ImagePlotter1D(TrendAxis,MaxDeltaIPF1,image_aspectratio,fig,ax[1])
	ImagePlotter1D(TrendAxis,MaxDeltaIPF2,image_aspectratio,fig,ax[1])
	ImagePlotter1D(TrendAxis,MaxDeltaIDiv1,image_aspectratio,fig,ax[1])
	ImagePlotter1D(TrendAxis,MaxDeltaIDiv2,image_aspectratio,fig,ax[1])

	#Apply figure labels, title, legend and cropping
	Title = 'Maximum Coil Current for Varying '+Parameter
	Xlabel1, Ylabel1 = '', 'Maximum Coil \n Current $I$ [kA]'
	Xlabel2, Ylabel2 = 'Varying: '+Parameter, 'Maximum Change in \n Current $\Delta I$ [kA ms$^{-1}$]'
	Legend = Coilset
	ImageOptions(fig,ax[0],Xlabel1,Ylabel1,Title,Legend,Crop=False,Rotate=False)
	ImageOptions(fig,ax[1],Xlabel2,Ylabel2,'',Legend,Crop=False,Rotate=False)

	#Clean up and save figure to home directory
	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig('Coil_CurrentTrends'+image_extension)
#	plt.show()
	plt.close('all')

	print'------------------------------'
	print'# Coil Current Trends Complete'
	print'------------------------------'
#endfor

#=====================================================================#
#=====================================================================#





#====================================================================#
				  #COIL VOLTAGE TRENDS DIAGNOSTIC#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_CoilVoltageTrends == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

	#Extract absolute coil voltages and time axis from series directories
	Filename = 'icoil_Data/CoilVoltages.txt'
	NumCoils = len(ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical'))-1
	Time_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	VSol_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	VPF1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	VPF2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	VDiv1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]
	VDiv2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[5]

	#Extract delta coil voltages and time axis from series directories
	Filename = 'icoil_Data/DeltaCoilVoltages.txt'
	NumCoils = len(ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical'))-1
	DeltaTime_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[0]
	DeltaVSol_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[1]
	DeltaVPF1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[2]
	DeltaVPF2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[3]
	DeltaVDiv1_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[4]
	DeltaVDiv2_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical')[5]

	#Create trendaxis from folder names
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

	#Rescale absolute coil currents: [V] to [kV]
	for i in range(0,len(VSol_Arrays)):
		for j in range(0,len(VSol_Arrays[i])):	
			VSol_Arrays[i][j] = VSol_Arrays[i][j]/1000.0				#[kV]
			VPF1_Arrays[i][j] = VPF1_Arrays[i][j]/1000.0				#[kV]
			VPF2_Arrays[i][j] = VPF2_Arrays[i][j]/1000.0				#[kV]
			VDiv1_Arrays[i][j] = VDiv1_Arrays[i][j]/1000.0				#[kV]
			VDiv2_Arrays[i][j] = VDiv2_Arrays[i][j]/1000.0				#[kV]
		#endfor
	#endfor

	#Rescale Delta coil currents: [V] to [kV]
	for i in range(0,len(DeltaVSol_Arrays)):
		for j in range(0,len(DeltaVSol_Arrays[i])):	
			DeltaVSol_Arrays[i][j] = DeltaVSol_Arrays[i][j]/1000.0		#[kV/ms]
			DeltaVPF1_Arrays[i][j] = DeltaVPF1_Arrays[i][j]/1000.0		#[kV/ms]
			DeltaVPF2_Arrays[i][j] = DeltaVPF2_Arrays[i][j]/1000.0		#[kV/ms]
			DeltaVDiv1_Arrays[i][j] = DeltaVDiv1_Arrays[i][j]/1000.0	#[kV/ms]
			DeltaVDiv2_Arrays[i][j] = DeltaVDiv2_Arrays[i][j]/1000.0	#[kV/ms]
		#endfor
	#endfor

	#Calculate maximum coil current for each coil
	MaxVPF1,MaxVPF2,MaxVDiv1,MaxVDiv2 = list(),list(),list(),list()
	MaxVSol,MaxVAvg = list(),list()
	for i in range(0,len(VSol_Arrays)):
		MaxVSol.append( max(VSol_Arrays[i], key=abs) )			#[kV]
		MaxVPF1.append( max(VPF1_Arrays[i], key=abs) )			#[kV]
		MaxVPF2.append( max(VPF2_Arrays[i], key=abs) )			#[kV]
		MaxVDiv1.append( max(VDiv1_Arrays[i], key=abs) )		#[kV]
		MaxVDiv2.append( max(VDiv2_Arrays[i], key=abs) )		#[kV]
		Tot = abs(MaxVSol[i])+abs(MaxVPF1[i])+abs(MaxVPF2[i])+abs(MaxVDiv1[i])+abs(MaxVDiv2[i])
		MaxVAvg.append( Tot/NumCoils )							#[kV]
	#endfor

	#Calculate maximum change in current experienced for each coil set
	MaxDeltaVPF1,MaxDeltaVPF2 = list(),list()
	MaxDeltaVDiv1,MaxDeltaVDiv2 = list(),list()
	MaxDeltaVSol = list()
	for i in range(0,len(DeltaVSol_Arrays)):
		MaxDeltaVSol.append( max(DeltaVSol_Arrays[i], key=abs) )
		MaxDeltaVPF1.append( max(DeltaVPF1_Arrays[i], key=abs) )
		MaxDeltaVPF2.append( max(DeltaVPF2_Arrays[i], key=abs) )
		MaxDeltaVDiv1.append( max(DeltaVDiv1_Arrays[i], key=abs) )
		MaxDeltaVDiv2.append( max(DeltaVDiv2_Arrays[i], key=abs) )
	#endfor

	#===================##===================#
	#===================##===================#

#	#Create output folder for all coil trend figures
#	CurrentTrendsDir = CreateNewFolder(SeriesDirString,'/ICoil_Trends/')
	#Organize figure labelling variables
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
	else: Parameter = ParameterVaried
	#endif

	#For every simulation folder in the current series:
	for i in range(0,NumCoils):

		#Set which coil current timetrace to compare
		Coilset = ['VSol','VPF1','VPF2','VDiv1','VDiv2']
		if i == 0: Coil,VCoil_Arrays = 'ISol',VSol_Arrays
		if i == 1: Coil,VCoil_Arrays = 'IPF1',VPF1_Arrays
		if i == 2: Coil,VCoil_Arrays = 'IPF2',VPF2_Arrays
		if i == 3: Coil,VCoil_Arrays = 'IDiv1',VDiv1_Arrays
		if i == 4: Coil,VCoil_Arrays = 'IDiv2',VDiv2_Arrays

		#Plot each coil current with respect to time
		fig,ax = figure(image_aspectratio,1,shareX=False)
		for j in tqdm(range(0,len(VCoil_Arrays))):
			ImagePlotter1D(Time_Arrays[j],VCoil_Arrays[j],image_aspectratio,fig,ax)
		#endfor

		#Apply figure labels, title, legend and cropping
		Range = '['+str(min(TrendAxis))+' - '+str(max(TrendAxis))+']'
		Title = 'Time-Traces of '+Coil+' Coil Voltages for '+Parameter+' in '+Range
		Xlabel, Ylabel = 'Coil Voltage $V$ [kV]', 'Time $\\tau$ [ms]'
		Legend = TrendAxis
		ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False,Rotate=False)

		plt.tight_layout(pad=3.0,h_pad=1.0)
		plt.savefig(SeriesDirString+'/'+Coil+'_Voltage_Trends'+image_extension)
#		plt.show()
		plt.close('all')
	#endfor

	#===================##===================#
	#===================##===================#

	#Create image limits
	GlobalMaxDelta = MaxDeltaISol+MaxDeltaIPF1+MaxDeltaIPF2+MaxDeltaIDiv1+MaxDeltaIDiv2
	Ylims = [min(GlobalMaxDelta),max(GlobalMaxDelta)]

	#Create fig and plot...
	fig,ax = figure(image_aspectratio,2,shareX=False)
	#...ax[0] coil currents for each simulation...
	ImagePlotter1D(TrendAxis,MaxVSol,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxVPF1,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxVPF2,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxVDiv1,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxVDiv2,image_aspectratio,fig,ax[0])
	ImagePlotter1D(TrendAxis,MaxVAvg,image_aspectratio,fig,ax[0])
	ax[0].plot(TrendAxis[MaxVAvg.index(min(MaxVAvg))],min(MaxVAvg), 'kv', ms=14, lw=2.0) #Min Avg
	#...and [1] delta coil currents for each simulation
	ImagePlotter1D(TrendAxis,MaxDeltaVSol,image_aspectratio,fig,ax[1])
	ImagePlotter1D(TrendAxis,MaxDeltaVPF1,image_aspectratio,fig,ax[1])
	ImagePlotter1D(TrendAxis,MaxDeltaVPF2,image_aspectratio,fig,ax[1])
	ImagePlotter1D(TrendAxis,MaxDeltaVDiv1,image_aspectratio,fig,ax[1])
	ImagePlotter1D(TrendAxis,MaxDeltaVDiv2,image_aspectratio,fig,ax[1])

	#Apply figure labels, title, legend and cropping
	Title = 'Maximum Coil Current for Varying '+Parameter
	Xlabel1, Ylabel1 = '', 'Maximum Coil \n Voltage $V$ [kV]'
	Xlabel2, Ylabel2 = 'Varying: '+Parameter, 'Maximum Change in \n Voltage $\Delta V$ [kV ms$^{-1}$]'
	Legend = Coilset
	ImageOptions(fig,ax[0],Xlabel1,Ylabel1,Title,Legend,Crop=False,Rotate=False)
	ImageOptions(fig,ax[1],Xlabel2,Ylabel2,'',Legend,Crop=False,Rotate=False)

	#Clean up and save figure to home directory
	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig('Coil_VoltageTrends'+image_extension)
#	plt.show()
	plt.close('all')

	print'------------------------------'
	print'# Coil Voltage Trends Complete'
	print'------------------------------'

#=====================================================================#
#=====================================================================#




























#====================================================================#
					  #PLASMA BREAKDOWN DIAGNOSTICS#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_Breakdown == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

	#Extract relevent data from series directories
	Lc_Arrays = ExtractFIESTAData(SimulationDirs,'Lc.txt','2D','Vertical',False)
	EMF_Arrays = ExtractFIESTAData(SimulationDirs,'VLoop.txt','2D','Vertical',False)

	#Strip and organize breakdown parameters into seperate variables
	VLoop_Arrays,ELoop_Arrays,ELoopEff_Arrays = list(),list(),list()
	DeltaPhi_Arrays = list()
	for l in range(0,NumFolders):
		#HACKY:: Extra [0][0] removes outer array(s) for single variable readin - FIX THIS!!!
		Lc_Arrays[l] = Lc_Arrays[l][0][0]				
		VLoop_Arrays.append(EMF_Arrays[l][0][0])
		DeltaPhi_Arrays.append(EMF_Arrays[l][1][0])
		ELoop_Arrays.append(EMF_Arrays[l][2][0])
		ELoopEff_Arrays.append(EMF_Arrays[l][3][0])
	#endfor

	#Create trendaxis from folder names
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
	else: Parameter = ParameterVaried
	#endif

	#Plot connection length trends for each detected simulation folder
	fig,ax,im = ImagePlotter1D(TrendAxis,Lc_Arrays,image_aspectratio)
	#endfor
	
	#Apply figure labels, title, legend and cropping
	Xlabel, Ylabel = 'Varying: '+Parameter, 'Connection Length $L_{c}$ [m]'
	Title = 'Connection Length for Varying '+str(Parameter)
	Legend = ['Lc']
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False,Rotate=False)

	#Clean up and save figure to home directory
	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig(os.getcwd()+'/'+SeriesDirString+'/Lc_Trends'+image_extension)
#	plt.show()
	plt.close('all')

	#===================##===================##===================#
	#===================##===================##===================#

	#Create arbitary pressure array over requested pressure range
	Resolution = 25000
	PressureArray = np.linspace(PaschenPressure[0],PaschenPressure[1],Resolution).tolist()

	#Initialise Paschen minimum Efield arrays 
	EMinArrays,MinEMinArray = list(),list()
	for i in range(0,len(Lc_Arrays)): EMinArrays.append( list() )

	#Construct EMinArrays for each connection length over provided background pressure range
	for i in range(0,len(EMinArrays)):
		for j in range(0,len(PressureArray)):
			#Compute minimum E-field for breakdown employing townsend coefficients for Duterium (Song2017)
			EMin = (PressureArray[j]*1.25E4)/np.log(510.0*PressureArray[j]*Lc_Arrays[i])	#[V/m]
			if EMin < 0: EMin = np.nan
			EMinArrays[i].append( EMin )
		#endfor
		MinEMinArray.append( min(filter(lambda v: v==v, EMinArrays[i])) )  					#[V/m]
	#endfor

	#Stretch ELoop_Arrays to the same length as EMinArrays for plotting
	for l in range(0,NumFolders):
		ELoop = ELoop_Arrays[l]
		ELoop_Arrays[l] = list()
		for i in range(0,len(EMinArrays[0])):
			ELoop_Arrays[l].append(ELoop)
		#endfor
	#endfor

	#Create trendaxis from folder names and organize figure labelling variables
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
	else: Parameter = ParameterVaried
	#endif

	#Round connective length and scale pressure array for plotting
	for i in range(0,len(Lc_Arrays)): Lc_Arrays[i] = round(Lc_Arrays[i],1)			#[m]
	for i in range(0,len(PressureArray)): PressureArray[i] = PressureArray[i]*1000	#[mTorr]

	#===================##===================#
	#===================##===================#

	#Create fig and plot plasma current for each detected simulation folder
	fig,ax = figure(image_aspectratio,1,shareX=False)
	ColourCycle = ['k','r','b','g','c','m','y']
	for l in tqdm(range(0,NumFolders)): ax.plot(PressureArray,EMinArrays[l], ColourCycle[l]+'-', lw=2, ms=12)
	for l in range(0,NumFolders):	ax.plot(PressureArray,ELoop_Arrays[l], ColourCycle[l]+'--', lw=1.5, ms=12)
	#endfor
	
	#Apply figure labels, title, legend and cropping
	Xlabel, Ylabel = 'Prefill Pressure $P$ [mTorr]', 'Toroidal E-Field $\\bf{E}$ [V m$^{-1}$]'
	Title = 'Plasma Current Time-Trace for Varying '+str(Parameter)
	Legend = [Parameter+'='+str(TrendAxis[i]) for i in range(0,len(TrendAxis))]
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False,Rotate=False)
	ax.set_xscale("log")
	ax.set_yscale("log")

	#Hacky, need to update cropping function for 1D figures
	ax.set_xlim(1e-3,5)			#Pressure [mTorr]
	ax.set_ylim(0.05,25)		#E-Field [Vm-1]

	#Clean up and save figure to home directory
	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig('Paschen_Breakdown'+image_extension)
#	plt.show()
	plt.close('all')

	#Plot trend in minimum breakdown E-field with respect to varied parameter
#	from mpl_toolkits.axes_grid.inset_locator import inset_axes
#	left, bottom, width, height = [0.19,0.19,0.25,0.20]			#[0.62,0.27,0.23,0.23]
#	ax2 = fig.add_axes([left, bottom, width, height])
	###
#	ax2.plot(TrendAxis,MinEMinArray,'ko--', markerfacecolor='none', ms=8, lw=1.5)
#	ax2.plot(TrendAxis[MinEMinArray.index(min(MinEMinArray))],min(MinEMinArray),'ko', ms=10)
#	ax2.set_ylabel('Minimum $\\bf{E}$ [V m$^{-1}$]', labelpad=0, fontsize=14.5)
#	ax2.set_xlabel('Varying: '+Parameter, fontsize=15)
#	ax2.ticklabel_format(axis="x", style="sci", scilimits=(-2,3))
#	ax2.tick_params(axis='x', labelsize=14)
#	ax2.tick_params(axis='y', labelsize=14)
#	ax2.set_xlim(0,1)
#	ax2.set_ylim(0.79,1.01)

	print'--------------------------------'
	print'# Breakdown Diagnostics Complete'
	print'--------------------------------'
#endif

#=====================================================================#
#=====================================================================#



















#====================================================================#
		 			  #PLASMA CURRENT DIAGNOSTICS#
#====================================================================#

#Compare optimised plasma current profiles
if savefig_PlasmaCurrent == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

	#Extract plasma current data from series directories
	Filename = 'RZIP_Data/Ip.txt'
	Ip_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical',False)

	#Create trendaxis from folder names
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

	#Rescale data for plotting: [A] to [kA]
	for l in range(0,NumFolders):
		for i in range(0,len(Ip_Arrays[l][0])):
			Ip_Arrays[l][1][i] = Ip_Arrays[l][1][i]/1000.0
		#endfor
	#endfor

	#Calculate maximum Ip for each simulation over the full series
	Ip_MaxTrend,Ip_MinTrend = list(),list()
	for l in range(0,NumFolders):
		Ip_MaxTrend.append(max(Ip_Arrays[l][1]))
		Ip_MinTrend.append(min(Ip_Arrays[l][1]))
	#endfor

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
	else: Parameter = ParameterVaried
	#endif

	#Create fig and plot plasma current for each detected simulation folder
	fig,ax = figure(image_aspectratio,1,shareX=False)
	for l in tqdm(range(0,NumFolders)):
		ImagePlotter1D(Ip_Arrays[l][0],Ip_Arrays[l][1],image_aspectratio,fig,ax)
	#endfor
	
	#Apply figure labels, title, legend and cropping
	Xlabel, Ylabel = 'Time $\\tau$ [ms]', 'Plasma Current $I_{p}$ [kA]'
	Title = 'Plasma Current Time-Trace for Varying '+str(Parameter)
	Legend = TrendAxis
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False,Rotate=False)

	#Hacky, need to update cropping function for 1D figures
	ax.set_xlim( min(Ip_Arrays[l][0])*1.20	,max(Ip_Arrays[l][0])*1.30   )		
	ax.set_ylim( -5							,max(max(Ip_Arrays[l]))*1.05 )

	#Clean up and save figure to home directory
	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig('Ip_Trends'+image_extension)
#	plt.show()
	plt.close('all')

	#Plot trend in plasma current with respect to varied parameter
#	from mpl_toolkits.axes_grid.inset_locator import inset_axes
#	left, bottom, width, height = [0.23,0.63,0.25,0.25]			#[0.62,0.27,0.23,0.23]
#	ax2 = fig.add_axes([left, bottom, width, height])
	###
#	ax2.plot(TrendAxis,Ip_MaxTrend,'ko--', ms=8, lw=1.5)
#	ax2.legend(['Max $I_{p}$'], fontsize=14, frameon=False)
#	ax2.set_ylabel('Maximum Plasma \n Current $I_{p,max}$ [kA]', labelpad=0, fontsize=14.5)
#	ax2.set_xlabel('Varied Parameter: '+Parameter, fontsize=15)
#	ax2.xaxis.set_major_locator(ticker.MultipleLocator(90))
#	ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
#	ax2.tick_params(axis='x', labelsize=14)
#	ax2.tick_params(axis='y', labelsize=14)
#	ax2.set_xlim( min(TrendAxis),max(TrendAxis)*1.10 )
#	ax2.set_ylim(0.79,1.01)
#endif

#=====================================================================#
#=====================================================================#



#====================================================================#
					  #EDDY CURRENT DIAGNOSTICS#
#====================================================================#

#Compare vessel eddy current profiles
if savefig_EddyCurrent == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

	#Extract plasma current data from series directories
	Filename = 'RZIP_Data/IPass.txt'
	IPassive_Arrays = ExtractFIESTAData(SimulationDirs,Filename,'2D','Vertical',False)

	#Create trendaxis from folder names
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

	#Rescale data for plotting: [A] to [kA]
	for l in range(0,NumFolders):
		for i in range(0,len(IPassive_Arrays[l][0])):
			IPassive_Arrays[l][1][i] = IPassive_Arrays[l][1][i]/1000.0
		#endfor
	#endfor

	#Calculate maximum net eddy current for each simulation over the full series
	IPassive_MaxTrend,IPassive_MinTrend = list(),list()
	for l in range(0,NumFolders):
		IPassive_MaxTrend.append(max(IPassive_Arrays[l][1]))
		IPassive_MinTrend.append(min(IPassive_Arrays[l][1]))
	#endfor

	#===================##===================#
	#===================##===================#

	#Organize figure labelling variables
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
	else: Parameter = ParameterVaried
	#endif

	#Create fig and plot passive currents for each detected simulation folder
	fig,ax = figure(image_aspectratio,1,shareX=False)
	for l in tqdm(range(0,NumFolders)):
		ImagePlotter1D(IPassive_Arrays[l][0],IPassive_Arrays[l][1],image_aspectratio,fig,ax)
	#endfor
	
	#Apply figure labels, title, legend and cropping
	Xlabel, Ylabel = 'Time $\\tau$ [ms]', 'Passive Current $I_{passive}$ [kA]'
	Title = 'Passive Current Time-Trace for Varying '+str(Parameter)
	Legend = TrendAxis
	ImageOptions(fig,ax,Xlabel,Ylabel,Title,Legend,Crop=False,Rotate=False)

	#Hacky, need to update cropping function for 1D figures
	ax.set_xlim( min(IPassive_Arrays[l][0])*1.20,	max(IPassive_Arrays[l][0])*1.30 )		
	ax.set_ylim( -5,								max(max(IPassive_Arrays[l]))*1.05 )

	#Clean up and save figure to home directory
	plt.tight_layout(pad=3.0,h_pad=1.0)
	plt.savefig('PassiveCurrent_Trends'+image_extension)
#	plt.show()
	plt.close('all')

	#Plot trend in plasma current with respect to varied parameter
#	from mpl_toolkits.axes_grid.inset_locator import inset_axes
#	left, bottom, width, height = [0.23,0.63,0.25,0.25]			#[0.62,0.27,0.23,0.23]
#	ax2 = fig.add_axes([left, bottom, width, height])
	###
#	ax2.plot(TrendAxis,IPass_MaxTrend,'ko--', ms=8, lw=1.5)
#	ax2.legend(['Max $I_{p}$'], fontsize=14, frameon=False)
#	ax2.set_ylabel('Maximum Eddy \n Current $I_{p,max}$ [kA]', labelpad=0, fontsize=14.5)
#	ax2.set_xlabel('Varied Parameter: '+Parameter, fontsize=15)
#	ax2.xaxis.set_major_locator(ticker.MultipleLocator(90))
#	ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
#	ax2.tick_params(axis='x', labelsize=14)
#	ax2.tick_params(axis='y', labelsize=14)
#	ax2.set_xlim( min(TrendAxis),max(TrendAxis) )
#	ax2.set_ylim(0.79,1.01)

	print'------------------------------'
	print'# Current Diagnostics Complete'
	print'------------------------------'
#endif

#=====================================================================#
#=====================================================================#



















#====================================================================#
				         #VERTICAL STABILITY#	---- TO BE REFACTORED
#====================================================================#

if savefig_VerticalStability == True:

	#Obtain simulation folder directories for project and requested series
	SimulationNames = ExtractSubDirs(SeriesDirString,Root=False)
	SimulationDirs = ExtractSubDirs(SeriesDirString,Root=True)
	NumFolders = len(SimulationDirs)

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

	#Remove any header string from the data
	for i in range(0,len(ISol_efit)):
		ISol_efit[i] = ISol_efit[i][1::]
		IPF1_efit[i] = IPF1_efit[i][1::]
		IPF2_efit[i] = IPF2_efit[i][1::]
		IDiv1_efit[i] = IDiv1_efit[i][1::]
		IDiv2_efit[i] = IDiv2_efit[i][1::]
		###
		ISol_Pert[i] = ISol_Pert[i][1::]
		IPF1_Pert[i] = IPF1_Pert[i][1::]
		IPF2_Pert[i] = IPF2_Pert[i][1::]
		IDiv1_Pert[i] = IDiv1_Pert[i][1::]
		IDiv2_Pert[i] = IDiv2_Pert[i][1::]
	#endfor

	#Create trendaxis from folder names
	ParameterVaried = 'test'
	TrendAxis = CreateTrendAxis(SimulationNames,ParameterVaried,TrendAxisVariables)

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
	if len(TrendAxisVariables) > 0: Parameter = TrendAxisVariables
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
	plt.savefig('VerticalStability_Trends.png')
#	plt.show()
	plt.close('all')

	print'-----------------------------'
	print'# Stability Analysis Complete'
	print'-----------------------------'
#endif

#=====================================================================#
#=====================================================================#





















