%% SMall Aspect Ratio Tokamak (SMART), V3p1 init.nam

clear 
clc
close all

%%%%%%%%%%

%Add FIESTA trunk path, include path to any extra functions.
FIESTATrunk = "~/Postdoc Seville/FIESTA/Source Code/FIESTA V8.8";
FunctionsLocal = "Functions";
FunctionsRemote = "../../Functions";
addpath(genpath(FIESTATrunk),genpath(FunctionsLocal),genpath(FunctionsRemote));

%Set maximum number of concurrent CPU threads in use
NumThreads = 2;
NumThreads = maxNumCompThreads(NumThreads);

%%%%%%%%%%%%%%%%%%%  DEFINE DATA OUTPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%

%Define figure extension
FigExt = '.png'; 		%'.png','.eps','.pdf'

%Define project and series names
ProjectName = 'S1-000008';		%Define global project name
SeriesName = 'Default';         %Define parameter scan series name

%Create global output folders for saved data and figures
ASCIIDir = 'RawData/'; mkdir(ASCIIDir);
%FigDir = 'Figures/'; mkdir(FigDir);		%Not Currently Implimented

%Create simulation name based upon relevant run parameters
SimName = 'DefaultSimName';
disp([ 'SimName: ' SimName ]);
disp([ ' ' ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DEFINE REACTOR GEOMETRY                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Vessel Wall Thickness
VWall_Inboard=0.004;	%Inboard Wall Thickness		[m]
VWall_Outboard=0.008;	%Outboard Wall Thickness	[m]
VWall_Upper=0.015;		%Top Wall Thickness			[m]
VWall_Lower=0.015;		%Bottom Wall Thickness		[m]

%Define Vessel Internal Geometry (Does not include wall thickness)
VesselRMinInner=+0.15+VWall_Inboard;	% R min position [m] 	%Inboard wall 'fixed' by outer edge.
VesselRMaxInner=+0.80;					% R max position [m]
VesselZMinInner=-0.80;					% Z min position [m]
VesselZMaxInner=+0.80;					% Z max position [m]

%Define center points of vessel walls (Inner Geometry + half wall thickness)
ZMinCentre=VesselZMinInner-(VWall_Lower/2);		% Lower Wall 'grows' outwards (-Z direction)
ZMaxCentre=VesselZMaxInner+(VWall_Upper/2);		% Upper Wall 'grows' outwards (+Z direction)
RMinCentre=VesselRMinInner-(VWall_Inboard/2);	% Inboard wall 'grows' inwards (-R direction)
RMaxCentre=VesselRMaxInner+(VWall_Outboard/2);	% Outboard wall 'grows outwards (+R direction)

%%%%%%%%%%%%%%%%%%%%%%%  DEFINE COIL GEOMETRY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Solenoid Geometry and Parameters
nSol=210;					 		% Number of Solenoid Windings
RSolInner=0.115; RSolOuter=0.145;   % Inner and Outer solenoid radii    [m] 
RSol=(RSolInner+RSolOuter)/2;       % Central radius of solenoid (0.13) [m]
ZMinSol=ZMinCentre-(VWall_Upper/2); % Solenoid Min Z position           [m]
ZMaxSol=ZMaxCentre+(VWall_Lower/2); % Solenoid Max Z position           [m]

%Number of Radial (R) and axial (Z) PF coil windings
nZDiv1=6; nRDiv1=4;
nZDiv2=6; nRDiv2=4;
nZPF1=6; nRPF1=4;
nZPF2=6; nRPF2=4;

%Calculate total number of windings in each coil
nDiv1=nZDiv1*nRDiv1;
nDiv2=nZDiv2*nRDiv2;
nPF1=nZPF1*nRPF1;
nPF2=nZPF2*nRPF2;

%Define coil turn dimensions to enable cross-section calculation
width_PF=0.042;  % Width of a turn (m)
height_PF=0.035; % Height of a turn (m)

%{
%Define central location of coil sets
R_PF1=(0.90)+ww_R/2;  %R position of PF1 (m)		%0.90m
Z_PF1=(0.30)+ww_Z/2;  %Z Position of PF1 (m)		%0.30m
R_PF2=(0.90)+ww_R/2;  %R Position of PF2 (m)		%0.90m
Z_PF2=(0.60)+ww_Z/2;  %Z Position of PF2 (m)		%0.60m
R_Div1=(0.15)+ww_R/2; %R Position of Div1 (m)		%0.15m	(Originally 0.25m)
Z_Div1=(0.85)+ww_Z/2; %Z Position of Div1 (m)		%0.85m
R_Div2=(0.45)+ww_R/2; %R Position of Div2 (m)		%0.45m	(Originally 0.55m)
Z_Div2=(0.85)+ww_Z/2; %Z Position of Div2 (m)		%0.85m
%}

%Define central location of coil sets
R_PF1=0.90;  %R position of PF1 (m)		%0.90m
Z_PF1=0.30;  %Z Position of PF1 (m)		%0.30m
R_PF2=0.90;  %R Position of PF2 (m)		%0.90m
Z_PF2=0.60;  %Z Position of PF2 (m)		%0.60m
R_Div1=0.15; %R Position of Div1 (m)	%0.15m	(Originally 0.25m)
Z_Div1=0.85; %Z Position of Div1 (m)	%0.85m
R_Div2=0.45; %R Position of Div2 (m)	%0.45m	(Originally 0.55m)
Z_Div2=0.85; %Z Position of Div2 (m)	%0.85m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      DEFINE OPERATIONAL PARAMETERS                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%  DEFINE INITIAL PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%

%Define any required constants
mu0 = 1.2566e-06; % Magnetic Moment      	[I/m^2]
BEarth = 5.0E-5;  % Earth's Magnetic Field	[T]

%Define initial operating conditions (primarily used for Topeol2)
Te = 250;			% Electron Temperature [eV]
Ti = Te*0.1;		% Ion Temperature      [eV]
BT = 0.1;			% Toroidal B-Field     [T] (Defined at Rgeo)
Ip = 30e3;			% Plasma current       [A]
RGeo = 0.450;		% Geometrical Radius   [m]
ZGeo = 0.000;		% Geometrical Axis     [m]
RSep = 0.700;		% Separatrix Radius    [m] (~0.70)
rGeo = RSep-RGeo;	% Minor Radius         [m] (~0.25)
Aspect = RGeo/rGeo;	% Aspect ratio         [-] (~1.85)
Kappa = 1.80;		% Elongation           [-] (~1.8)
delta = 0.20;		% Triangularity        [-] (~0.2)
li2 = 1;			% Inductance	       [-]

%Compute further operating conditions (primarily used for Topeol2)
Gr_Limit = 1e20*(Ip*1e-6/(pi*Kappa*rGeo^2));     % Greenwald Limit       [m-3]
Gr_Frac = 0.15;                            % Greenwald Fraction       [-]
ne = Gr_Limit*Gr_Frac;                     % Electron Density         [m-3]  ~3E19
Irod = (BT*2*pi*RGeo)/mu0;                 % Central Rod Current      [A]
S = sqrt( (1.0+Kappa^2)/2.0 );             % Shaping factor           [-]
%deltaUp = (RGe-Rup)/a;                    % Upper-Triangularity      [-]
%deltaLo = (RGe-Rlo)/a;                    % Lower-Triangularity      [-]
%delta = (deltaUp+deltaLo)/2.0;            % Triangularity            [-]
%betaN = (betaT*BT*a)/(Ip*1e-6*mu0)        % Normalised Beta          [%] 
%betaT = (betaN/a*(Ip*1e-6))/BT;           % Beta toroidal            [%]
betaP = 3/2*ne*(Te+Ti)/(mu0*Ip/(2*pi*rGeo))^2*2*mu0*1.6e-19*Kappa; 	% Beta Poloidal  [%]
BZ = -mu0*Ip/(4*pi*RGeo)*(log(8*Aspect)+betaP+0.5*li2-3/2);    		% Vertical field [T]

%Define efit Equilibrium Operating Conditions
RGeo_efit = 0.440;					% Geometrical Radius	[m] (Default 0.44)
ZGeo_efit = 0.000;					% Geometrical Axis		[m] (Default 0.00)
rGeo_efit = 0.238;                  % Minor Radius	        [m] (Default 0.44/1.85)
Aspect_efit = RGeo_efit/rGeo_efit;  % Aspect Ratio          [-] (Default 1.85)
Kappa_efit = 1.80;					% Elongation			[-] (Default 1.80)
delta_efit = 0.20;					% Triangularity			[-] (Default 0.20)
efitGeometry_Init = [RGeo_efit, ZGeo_efit, rGeo_efit, Kappa_efit, delta_efit];

%Define feedback stability perturbations
deltaRGeo = 0.00;	% Small radial perturbation         [m]
deltaZGeo = 0.00;	% Small axial perturbation          [m]
deltaAspect = 0.00;	% Small aspect ratio perturbation   [-]
deltaKappa = 0.00;	% Small elongation perturbation     [-]
deltadelta = 0.00;	% Small triangiularity perturbation [-]
PertGeometry_Init = [deltaRGeo,deltaZGeo,deltaAspect,deltaKappa,deltadelta];

%Define Coil density, temperature and resistivity
coil_density = 1;                       % Relative Coil Density      [Arb]
coil_temp = 293.0;                      % Initial Coil Temperature   [K]
resistivity = copper_resistivity_at_temperature(coil_temp);

%Gas species analouge - H=1, He=2, Ar=11.85 (for Te < 280eV) https://www.webelements.com/argon/atoms.html
%May need to increase to account for impurities in plasma (Carbon Divertor)
Z_eff=1.0;                              % Effective Nuclear Charge   [e-]
%Calculate perpendicular and parallel plasma resistivity using Spitzer model
Lambda=(12*pi*((8.854E-12*1.6E-19*Te)^3/(ne*(1.6E-19)^6))^(1/2));
PlasmaResistPerp=(0.74*1.65E-9*Z_eff*log(Lambda))/((Te*1E-3)^(3/2));
PlasmaResistPara=PlasmaResistPerp/1.96;

%Define null field region radius (also specifies sensor_btheta radius)
a_eff=0.15;								% Null field region radius	 [m]


%%%%%%%%%%%%%%%%%%%  DEFINE SOL RAMP & COIL CURRENTS  %%%%%%%%%%%%%%%%%%%%%

%Notes:
%Negative coil currents attract the plasma, positive repel the plasma
%Symmetric Solenoid PrePulse and Equil currents aid power supply stability
%Rod current (Irod) sets the toroidal component of the magnetic field.

%Definition of CoilWaveform time intervals:
%time(1)--> All coils and Sol initiate at zero current          Init
%time(2)--> All coils initiate null-field configuration         PrePulse
%time(3)--> All coils maintain null-field configuration         InitRampDown
%time(4)--> Sol ramps down, PF/Div coils init equilibrium       MidRampDown - InitEquil
%time(5)--> Sol completes ramp down, maintain PF/Div coils      EndRampDown - MidEquil
%time(6)--> All coils maintain equilibrium configuration        EndEquil
%time(7)--> All coils and Sol terminate at zero current         Terminate
%%%%%%%
%time(3)-->time(5) lasts timescale TauR (Solenoid Ramp-Down TimeScale)
%time(5)-->time(6) lasts timescale TauP (Pulse/Discharge Timescale)
%%%%%%%

%Solenoid coil currents [kA]		%Phase1		%Phase1HomoWalls
I_Sol_Null=+750;					%+0750;		%+0850;
I_Sol_MidRamp='Linear';				%Dynamic    %Dynamic
I_Sol_Equil=-I_Sol_Null;			%-0750;     %-0850
I_Sol_EndEquil=-725; %-I_Sol_Null;  %-0725;     %-0775

%PF coil currents (At Equilibrium, time(4,5,6))
I_PF1_Equil=-500;					%-500;		%-500;
I_PF2_Equil=-500;					%-500;		%-500;
I_Div1_Equil=+000;					%+000;		%+ISol;
I_Div2_Equil=+900;					%+900;      %+900;

%Define number of time-steps (vertices) in the current waveforms
TauB = 0.020;			% Buffer Timescale     		[s] Determines tstep for Ip plot
TauR = 0.050;			% Ramp Timescale       		[s]
TauP = 0.020;			% Pulse Timescale      		[s]
%Time   [Init      PrePulse  InitRampDown  MidRampDown  EndRampDown  MidEquil     Terminate         ];
time =  [-4*TauB   -2*TauB   0.0           TauR/2.0     TauR         TauR+TauP    TauR+TauP+(2*TauB)];
nTime = length(time);	% Coil Waveform Timesteps	[-]

%Fit any dynamic coil currents {Pre-ramp, mid-ramp, end-ramp}
I_Sol_MidRamp = FitSolenoidRamp({I_Sol_Null,I_Sol_MidRamp,I_Sol_Equil},time);

%Construct Sol, PF/Div coil current waveforms vertices
%					  %!Null-Field!             %!Breakdown!   %!Efit Icoil!
%Time   	     [1,  2,            3,          4,             5,             6,             7];
ISol_Waveform =  [0,  I_Sol_Null,   I_Sol_Null, I_Sol_MidRamp, I_Sol_Equil,   I_Sol_EndEquil,0];
IPF1_Waveform =  [0,  NaN,          NaN,        NaN,           I_PF1_Equil,   I_PF1_Equil,   0];
IPF2_Waveform =  [0,  NaN,          NaN,        NaN,           I_PF2_Equil,   I_PF2_Equil,   0];
IDiv1_Waveform = ISol_Waveform;     %IDiv1 in Series with Solenoid
IDiv2_Waveform = [0,  NaN,          NaN,        NaN,           I_Div2_Equil,  I_Div2_Equil,  0];
%%%%%
CoilWaveforms = [ISol_Waveform; IPF1_Waveform; IPF2_Waveform; IDiv1_Waveform; IDiv2_Waveform];

%Define dynamic coils (i.e. which coil currents are fit by efit)
global efitCoils; efitCoils = {'PF1','PF2'};


%%%%%%%%%%%%%%%%%%  DISPLAY VARIABLE OUTPUT TO USER  %%%%%%%%%%%%%%%%%%%%%%

disp([ ' ' ]);
disp([ '%===== Initial Operating Parameters =====%' ]);
disp([ 'TauP = ' num2str(TauP*1000) ' [ms]' ]);
disp([ 'Ip = ' num2str(Ip/1000) ' [kA]' ]);
disp([ 'IRod = ' num2str(Irod) ' [kA]' ]);
disp([ 'BT = ' num2str(BT) ' [T]' ]);
disp([ 'BZ = ' num2str(BZ) ' [T]' ]);
%disp([ 'betaN = ' num2str(betaN) ' [%]' ]);
%disp([ 'betaT = ' num2str(betaT) ' [%]' ]);
%disp([ 'betaP = ' num2str(betaP) ' [%]' ]);
disp([ 'ne = ' num2str(ne) ' [m-3]' ]);
disp([ 'Te = ' num2str(Te) ' [eV]' ]);
disp([ 'Ti = ' num2str(Ti) ' [eV]' ]);

disp([ 'RGeo = ' num2str(RGeo) ' [m]' ]);
disp([ 'ZGeo = ' num2str(ZGeo) ' [m]' ]);
disp([ 'RSep = ' num2str(RSep) ' [m]' ]);
disp([ 'Minor Radius = ' num2str(rGeo) ' [m]' ]);
disp([ 'AspectRatio = ' num2str(Aspect) ' [-]' ]);
disp([ 'Elongation = ' num2str(Kappa) ' [-]' ]);
disp([ 'Shaping Factor = ' num2str(S) ' [-]' ]);
%disp([ 'Triangularity = ' num2str(delta) ' [-]' ]);

disp([ ' ' ]);
disp([ '%===== Initial Coil Currents =====%' ]);
disp([ 'I_Sol_PrePulse = ' num2str(I_Sol_Null/1000) ' [kA]' ]);
disp([ 'I_Sol_MidRamp = Dynamic [kA]' ]);
disp([ 'I_Sol_Equil = ' num2str(I_Sol_Equil/1000) ' [kA]' ]);
disp([ ' ' ]);
disp([ 'I_PF1_Equil = ' num2str(I_PF1_Equil/1000) ' [kA]' ]);
disp([ 'I_PF2_Equil = ' num2str(I_PF2_Equil/1000) ' [kA]' ]);
disp([ 'I_Div1_Equil = ' num2str(I_Div1_Equil/1000) ' [kA]' ]);
disp([ 'I_Div2_Equil = ' num2str(I_Div2_Equil/1000) ' [kA]' ]);
disp([ ' ' ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    INITIATE VESSEL AND COIL OBJECTS                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define FIESTA simulation grid limits and resolution
GridSize_R = [0.03, 1.0];	%[m]
GridSize_Z = [-1.3, 1.3];	%[m]
GridCells_R = 200;          %[Cells]
GridCells_Z = 251;          %[Cells]
GridRes_R = max(GridSize_R)/GridCells_R;	%[m/Cell] GridRes_R = 0.00500
GridRes_Z = max(GridSize_Z)/GridCells_Z;	%[m/Cell] GridRes_Z = 0.00518

%Generate fiesta_grid object over which equilibrum simulation will be performed
global Grid;
Grid = fiesta_grid(GridSize_R(1),GridSize_R(2),GridCells_R, GridSize_Z(1),GridSize_Z(2),GridCells_Z);

%Extract vectors of R and Z grid points for use in further diagnostics
rGrid=get(Grid,'r'); %1*200
zGrid=get(Grid,'z'); %1*251


%%%%%%%%%%%%%%%%%%  INITIATE VACUUM VESSEL FILAMENTS  %%%%%%%%%%%%%%%%%%%

%Define four vertices defined as the centre of each vessel corner
Vertice1=[RMinCentre ZMaxCentre];	%Top Left
Vertice2=[RMaxCentre ZMaxCentre];	%Top Right
Vertice3=[RMaxCentre ZMinCentre];	%Bottom Right
Vertice4=[RMinCentre ZMinCentre];	%Bottom Left
%Define wall corners going clockwise from top left (R,Z)
WallCorners=[[Vertice1]; [Vertice2]; [Vertice3]; [Vertice4]];
WallThickness=[[VWall_Upper]; [VWall_Outboard]; [VWall_Lower]; [VWall_Inboard]];

%Define constant vessel cell surface area (required for RZIp inductance matrix stability)
%NEED TO DETERMINE A BASELINE AREA WHICH GIVES GOOD RESOLUTION ON ALL WALLS (5.0e-5?)
BaselineArea = max(WallThickness)^2;     %[m^2] 
%Define normalisation factors for each wall to achieve baseline area
WallNormFactor1 = BaselineArea/(WallThickness(1)^2);
WallNormFactor2 = BaselineArea/(WallThickness(2)^2);
WallNormFactor3 = BaselineArea/(WallThickness(3)^2);
WallNormFactor4 = BaselineArea/(WallThickness(4)^2);
WallNormFactors = [[WallNormFactor1]; [WallNormFactor2]; [WallNormFactor3]; [WallNormFactor4]];
%WallNormFactors = [[1]; [1]; [1]; [1]];      % !!!FUDGE TO AVOID AREA SCALING!!!

%Construct filament arrays for each section of vessel wall:
%Radial Top Wall, Z=Zmax :: Top Left to Top Right (Vertice1 to Vertice2)
n_fil1_2 = (RMaxCentre-RMinCentre+2*WallThickness(1))/WallThickness(1);   %Set number of filaments within vessel wall
n_fil1_2 = round( n_fil1_2/WallNormFactors(1) );                          %Scale number of filaments if required
R_lin1_2 = linspace(RMinCentre,RMaxCentre,n_fil1_2);            %Create evenly spaced array of R and Z coordinates
Z_lin1_2 = ZMaxCentre*ones(1,n_fil1_2);                         %Create evenly spaced array of R and Z coordinates
Length = WallThickness(1)*WallNormFactors(1);                   %Define thickness parallel to wall direction
Width = WallThickness(1);                                       %Define thickness perpendicular to wall direction
R_VWall1_2 = linspace(Length,Length,n_fil1_2);                  %Create Radial wall thickness array of size n_fil
Z_VWall1_2 = linspace(Width,Width,n_fil1_2);                    %Create Axial wall thickness array of size n_fil

%Axial Outboard Wall, R=Rmax :: Top Right to Bottom Right (Vertice2 to Vertice3)
n_fil2_3 = (ZMaxCentre-ZMinCentre+2*WallThickness(2))/WallThickness(2);   %Set number of filaments within vessel wall
n_fil2_3 = round( n_fil2_3/WallNormFactors(2) );                          %Scale number of filaments if required
R_lin2_3 = RMaxCentre*ones(1,n_fil2_3);
Z_lin2_3 = linspace(ZMaxCentre,ZMinCentre,n_fil2_3);
R_VWall2_3 = linspace(VWall_Outboard,VWall_Outboard,n_fil2_3);
Z_VWall2_3 = linspace(VWall_Outboard,VWall_Outboard,n_fil2_3);
Length = WallThickness(2)*WallNormFactors(2);                   %Define thickness parallel to wall direction
Width = WallThickness(2);                                       %Define thickness perpendicular to wall direction
R_VWall2_3 = linspace(Width,Width,n_fil2_3);                    %Create Radial wall thickness array of size n_fil
Z_VWall2_3 = linspace(Length,Length,n_fil2_3);                  %Create Axial wall thickness array of size n_fil

%Radial Bottom Wall, Z=Zmin :: Bottom Right to Bottom Left (Vertice3 to Vertice4)
n_fil3_4 = (RMaxCentre-RMinCentre+2*WallThickness(3))/WallThickness(3);   %Set number of filaments within vessel wall
n_fil3_4 = round( n_fil3_4/WallNormFactors(3) );                          %Scale number of filaments if requiredr
Z_lin3_4 = ZMinCentre*ones(1,n_fil3_4);
R_lin3_4 = linspace(RMaxCentre,RMinCentre,n_fil3_4);
Length = WallThickness(3)*WallNormFactors(3);                    %Define thickness parallel to wall direction
Width = WallThickness(3);                                        %Define thickness perpendicular to wall direction
R_VWall3_4 = linspace(Length,Length,n_fil3_4);                   %Create Radial wall thickness array of size n_fil
Z_VWall3_4 = linspace(Width,Width,n_fil3_4);                     %Create Axial wall thickness array of size n_fil

%Axial Inboard Wall, R=Rmin :: Bottom Left to Top Left (Vertice4 to Vertice1)
n_fil4_1 = (ZMaxCentre-ZMinCentre+2*WallThickness(4))/WallThickness(4);   %Set number of filaments within vessel wall
n_fil4_1 = round( n_fil4_1/WallNormFactors(4) );                          %Scale number of filaments if required
R_lin4_1 = RMinCentre*ones(1,n_fil4_1);
Z_lin4_1 = linspace(ZMinCentre,ZMaxCentre,n_fil4_1);    
Length = WallThickness(4)*WallNormFactors(4);                    %Define thickness parallel to wall direction
Width = WallThickness(4);                                        %Define thickness perpendicular to wall direction
R_VWall4_1 = linspace(Width,Width,n_fil4_1);                     %Create Radial wall thickness array of size n_fil
Z_VWall4_1 = linspace(Length,Length,n_fil4_1);                   %Create Axial wall thickness array of size n_fil

%Assemble the vessel wall filament position arrays
R_Lin_Array=[R_lin1_2 R_lin2_3 R_lin3_4 R_lin4_1]';
Z_Lin_Array=[Z_lin1_2 Z_lin2_3 Z_lin3_4 Z_lin4_1]';
RVWall_Array=[R_VWall1_2 R_VWall2_3 R_VWall3_4 R_VWall4_1]';
ZVWall_Array=[Z_VWall1_2 Z_VWall2_3 Z_VWall3_4 Z_VWall4_1]';

%Remove duplicate cells at wall corners
%This method should ensure that each wall owns it's starting vertex
dup = (abs(diff(R_Lin_Array))+abs(diff(Z_Lin_Array))) > 0;
R_Lin_Array = R_Lin_Array(dup);         %length 173     %696 - NoNorm
Z_Lin_Array = Z_Lin_Array(dup);         %length 173     %696 - NoNorm
RVWall_Array = RVWall_Array(dup);       %length 173     %696 - NoNorm
ZVWall_Array = ZVWall_Array(dup);       %length 173     %696 - NoNorm

%Construct vessel wall FIESTA filaments using position arrays
%Inputs(R,Z,r_thick,z_thick,1,0,0) where {R=MajorRadius, Z=Height, r=MinorRadius, z=MinorHeight}
for i=1:length(R_Lin_Array)
    vessel_filament(i) = fiesta_filament(R_Lin_Array(i),Z_Lin_Array(i),RVWall_Array(i),ZVWall_Array(i),1,0,0);
end
%Enable induced currents in vessel wall filaments - used only to calculate eddy currents
%The vessel density and resistivity are set within fiesta_passive.m, may be settable here!
global passive; passive = fiesta_passive('STVesselPas',vessel_filament,'g');
global vessel; vessel = fiesta_vessel( 'STVessel',passive);

%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIATE PF COILS  %%%%%%%%%%%%%%%%%%%%%%%%%%

%Define and initiate PF coils - Arbitrary numbering of coils
global iSol; iSol = 1;       %Central Inducting Solenoid
global iPF1; iPF1 = 2;       %Upper Plasma Forming Coil
global iPF2; iPF2 = 3;       %Lower Plasma Forming Coil
global iDiv1; iDiv1 = 4;     %Inboard Divertor Coil
global iDiv2; iDiv2 = 5;     %Outboard Divertor Coil

%Create array containing number of coil windings - Used to generate coil objects
global coilturns; coilturns=[];
coilturns(iSol) = nSol; 
coilturns(iDiv1) = nDiv1;
coilturns(iDiv2) = nDiv2;
coilturns(iPF1) = nPF1;
coilturns(iPF2) = nPF2;
nPF = 5; 				%Total number of coils including solenoid

%Create coil set from parameters defined above. (Function made by Carlos Soria)
%Function createVESTPFCircuit creates two PF coils. One in (R, Z) and another in (R, -Z)
PF1  = createVestPFCircuit('PF1',R_PF1,Z_PF1,width_PF,height_PF,coilturns(iPF1),nZPF1,nRPF1,true, coil_temp, resistivity, coil_density);
PF2  = createVestPFCircuit('PF2',R_PF2,Z_PF2,width_PF,height_PF,coilturns(iPF2),nZPF2,nRPF2,true, coil_temp, resistivity, coil_density);
Div1 = createVestPFCircuit('Div1',R_Div1,Z_Div1,width_PF,height_PF,coilturns(iDiv1),nZDiv1,nRDiv1,true, coil_temp, resistivity, coil_density); 
Div2 = createVestPFCircuit('Div2',R_Div2,Z_Div2,width_PF,height_PF,coilturns(iDiv2),nZDiv2,nRDiv2,true, coil_temp, resistivity, coil_density);


%%%%%%%%%%%%%%%%%%%%%%  INITIATE CENTRAL SOLENOID  %%%%%%%%%%%%%%%%%%%%%%

%Number of filaments of the inductor (coil = number of turns)
nfil_ind_coil = coilturns(iSol); 

clear('coil_filaments');
Z_filament = linspace(ZMinSol,ZMaxSol,nfil_ind_coil);
%Construct central solenoid filaments - solenoid is treated as 'vessel wall' with nonzero current
for iFilament=1:nfil_ind_coil
	Constant=sqrt(70e-6);
    coil_filaments(iFilament) = fiesta_filament( RSol, Z_filament(iFilament), Constant, Constant ); 
end
Sol_Coil = fiesta_coil( 'psh_coil', coil_filaments, 'Blue', resistivity, coil_density );
Sol_circuit = fiesta_circuit( 'Sol', [1], [Sol_Coil] );

%Collate completed coilset and create FIESTA icoil object for equilibrium computation
global coilset; coilset = fiesta_coilset('SMARTcoilset',[Sol_circuit,PF1,PF2,Div1,Div2],false,R_Lin_Array',Z_Lin_Array');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  END SIMULATION INITIAL SET-UP  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           COMPUTE INITIAL TARGET AND NULL-FIELD EQUILIBRIUA             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TimeArray index for target equilibrium 
TimeIndex_Discharge = 5;                          %default time(5)
TimeIndex_NullField = 2;                          %default time(2)

%Create initial icoil object at requested TimeIndex
CoilWaveforms_Init = CoilWaveforms;
global icoil_init; icoil_init = fiesta_icoil(coilset);
%Assign equilibrium coil currents to icoil object [kA]
icoil_init.Sol=CoilWaveforms_Init(iSol,TimeIndex_Discharge);   %Solenoid Equilibrium Current
icoil_init.PF1=CoilWaveforms_Init(iPF1,TimeIndex_Discharge);   %PF1 Equilibrium Current
icoil_init.PF2=CoilWaveforms_Init(iPF2,TimeIndex_Discharge);   %PF2 Equilibrium Current
icoil_init.Div1=CoilWaveforms_Init(iDiv1,TimeIndex_Discharge); %Div1 Equilibrium Current
icoil_init.Div2=CoilWaveforms_Init(iDiv2,TimeIndex_Discharge); %Div2 Equilibrium Current


%%%%%%%%%%%%%%%%%%%%  COMPUTE DISCHARGE EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%

%Compute Jprofile from betaP and Ip employing Topeol2 Solver
jprofile = fiesta_jprofile_topeol2( 'Topeol2', betaP, 1, li2, Ip );

%Compute equilibrium (Psi(R,Z)) from the supplied jprofile, icoil and geometry
%Returns target equilibrium and CoilWaveforms for PF1 and PF2 at requested time_Index
[Equil,EquilParams,CoilWaveforms,efitGeometry,config] = ...                 %% THIS CONFIG SHOULD BE GLOBAL %%
    efitInverse(jprofile,Irod,CoilWaveforms_Init,efitGeometry_Init,TimeIndex_Discharge);

%Save discharge coil currents for all coils at TimeIndex_Discharge
CoilCurrentsEfit = transpose(CoilWaveforms(:,TimeIndex_Discharge));
icoil_efit = fiesta_icoil(coilset, CoilCurrentsEfit);

%Initiate virtual B-field sensors centered on Rgeo
sensor_btheta = InitiateBSensors(EquilParams,a_eff);

%RZIP computes coefficients [A B C D] using the null field sensors
%Output C is used to compute the null-field PF coil currents 
%Outputs curlyM and curlyR are used to compute the plasma and eddy currents
rzip_config = fiesta_rzip_configuration( 'RZIP', config, vessel, {sensor_btheta} );
[A, B, C, D, curlyM, curlyR, gamma, plasma_parameters, index, label_index, state] = ...
    response(rzip_config, Equil, 'rp', PlasmaResistPerp);
%Homogenious Vessel Walls
%Lp0prime = 1.1841e-06      %Condition number of curlyM: 497781.0 
%Differential Vessel Walls
%Lp0prime = 1.1923e-06      %Condition number of curlyM: 30321857.9 

%%%%%%%%%%%%%%%%%%%%%  COMPUTE OPTIMISED NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%%

%Update CoilWaveforms array with null-field values (Using NaN Mask)
RZIP_C = C;
CoilWaveforms = NullFieldWaveforms(CoilWaveforms,RZIP_C,sensor_btheta,TimeIndex_NullField);

%Save null-field coil currents for all coils at TimeIndex_NullField
CoilCurrentsNull = transpose(CoilWaveforms(:,TimeIndex_NullField)); 
icoil_null = fiesta_icoil(coilset, CoilCurrentsNull);

%Compute null-field equilibrium using null-field coil currents
equil_null = fiesta_equilibrium('SMART-Null', config, Irod, icoil_null );
EquilParams_Null = parameters(equil_null);

%%%%%%%%%%%%%%%%%%  COMPUTE BREAKDOWN FOR NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%

%Extract the null poloidal and toroidal B-field vector arrays
[BrData_Null,BzData_Null,BPhiData_Null,BpolData_Null,BtorData_Null] = ExtractBField(equil_null);

%Average null poloidal and toroidal fields over region of area a_eff^2
[BpolAvg_Null,BtorAvg_Null] = ExtractNullBMin(EquilParams,BpolData_Null,BtorData_Null,a_eff);

%Compute the average connection length within the null-field sensor region
%Lc = ConnectionLength(EquilParams,VesselRMaxInner,VesselRMinInner,VesselZMaxInner,VesselZMinInner);
Lc = 0.25*a_eff*(BtorAvg_Null/BpolAvg_Null);

%%%%%%%%%%%%%%%  COMPUTE DYNAMIC PLASMA & EDDY CURRENTS  %%%%%%%%%%%%%%%%%%

%Compute dynamic plasma and vessel eddy currents with new coil waveforms
[time_linear,time_adaptive,I_PF_output,V_PF_output,Ip_output,Vp_output,I_Passive] = ...
    DynamicCurrents(CoilWaveforms, time, curlyM, curlyR);

%Extract Vessel Eddy Currents during discharge (time='false' for absolute max)
VesselEddyCurrents = ExtractPassiveCurrents(I_Passive,time_adaptive,time(TimeIndex_Discharge));


%%%%%%%%%%%%%%%%%%%%%%%% PLOT TARGET EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%%%

%Plot target equilibrium following convergence
Title = {'SMART Target Equilibrium iter(0)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_TargetEquilibrium';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({Equil},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%  PLOT VIRTUAL SENSORS ONTO EQUILIBRIUM  %%%%%%%%%%%%%%%%%

Title = {'SMART Virtual Sensors',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '__VirtualBSensors';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({Equil,sensor_btheta},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%  PLOT NULL-FIELD PHI SURFACES  %%%%%%%%%%%%%%%%%%%%%

%Plot the optimised null-field phi
Title = {'SMART Null-field Equilibrium iter(0)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_NullPhi';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({equil_null},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%% PLOT NULL-FIELD BPOL  %%%%%%%%%%%%%%%%%%%%%

%Log poloidal and toroidal magnetic fields to show details (Sol Obscures)
logBpolData_Null = log(BpolData_Null);
logBtorData_Null = log(BtorData_Null);
%Plot the optimised null-field phi
Title = {'SMART Null-field iter(0)',' '};
CbarLabel = 'Null-field B_{\theta} ln([T])';
Filename = '__NullBpol';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({logBpolData_Null},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%% PLOT COIL CURRENT WAVEFORMS %%%%%%%%%%%%%%%%%%%%%%%% 

%Plot figure showing dynamic coil currents
figure;
plot(time_adaptive*1000, I_PF_output/1000);
title(gca,'SMART Initial Coil Current Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Coil Current (kA)');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
Filename = '_CurrentWaveforms';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%% PLOT PLASMA CURRENT %%%%%%%%%%%%%%%%%%%%%%%%%%%            

%Plot plasma current over full timescale
close all
plot(time_adaptive*1000, Ip_output/1000)
title(gca,'SMART Plasma Current iter(0)');
legend(gca,'Plasma Current'); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Plasma Current (kA)');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
Filename = '_PlasmaCurrent';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%% PLOT TOTAL EDDY CURRENT %%%%%%%%%%%%%%%%%%%%%%%%% 

%Sum all filaments (row-wise) to get total net passive current
Net_IPassive = sum(I_Passive,2);
%Plot net passive current density over full timescale
close all
plot(time_adaptive*1000, Net_IPassive/1000)
title(gca,'Net SMART Eddy Current iter(0)');
legend(gca,'Net Eddy Current'); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Net Vessel Current I_{Eddy} [kA]');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca,'YLim',[min(Net_IPassive/1000)*1.15 max(Net_IPassive/1000)*1.20]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
Filename = '_1DEddyCurrent';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%% PLOT SPATIALLY RESOLVED EDDY CURRENTS %%%%%%%%%%%%%%%%%%   

%Obtain the filament variables r and z
ptmp = get(vessel,'passives');
ftmp = get(ptmp,'filaments');
RR = get(ftmp(:),'r'); %dim 1*number of filaments
ZZ = get(ftmp(:),'z'); %dim 1*number of filaments
%Plot eddy currents within a cross-section of the vessel
close all
figure; hold on; axis equal;
plot(coilset);
scatter3(RR,ZZ,VesselEddyCurrents/1000,100,VesselEddyCurrents/1000,'filled');
title('SMART Vessel Eddy Currents iter(0)');
view(2) %2D view
colormap(plasma);
cbar = colorbar;
cbar.Label.String = 'Eddy-Current I_{eddy} [kA]';
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_EddyCurrent';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  DETERMINE FORCES UPON THE VESSEL                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Obtain the equilibrium B-field in R,Z and Phi
[BrData,BzData,BPhiData,BpolData,BtorData] = ExtractBField(Equil);

%Interpolate the B-fields onto a grid that aligns with the vessel grid 
%These are the values of the B-field at the vessel grid points
%Br_interp aligns with the vessel filament cells (RR, ZZ) from before
Br_interp = @(r,z) interpn(zGrid,rGrid,BrData,z,r);
Bz_interp = @(r,z) interpn(zGrid,rGrid,BzData,z,r);
Bphi_interp = @(r,z) interpn(zGrid,rGrid,BPhiData,z,r);

%Extract B-field at vessel walls - meshes are aligned so indexes are the same
Br_vessel=Br_interp(RR,ZZ);
Bz_vessel=Bz_interp(RR,ZZ);
Bphi_vessel=Bphi_interp(RR,ZZ);

%Combine Bfield values in [R,Phi,Z] into a single array for the vessel
%size(number of filaments*3), each row is the vector field on one filament
B_vessel=[Br_vessel' Bphi_vessel' Bz_vessel']; 
%The maximum current on each vessel filament is I_Passive_fil (size 1*number of filaments)
%Current vector is in the phi direction so only take magnitude [0, 1*I_Passive(phi), 0]
VesselEddyCurrentVec=VesselEddyCurrents'*[0 1 0]; 		%size [number of filaments*3]

%The force upon all the filament would be 2piR*Force_fil_cross. R is stores
%in RR, which contains all the R values in a vector form with number of fil components. 
%Force_fil_cross is a vector of 3 components. It would be difficult to
%multiply them, but we do not need to, right now, because to obtain the
%pressure R cancles out, since the areas are 2piR*anchura (or altura). We
%assimilate the 3D filament as a 2D filament, so that it has no width in
%the R axis, s that its surface is 2piR*altura

%Compute J X B force acting on each filament (J X B computed for all directions) 
Force_fil=cross(VesselEddyCurrentVec,B_vessel);	%[N] %size [number of filaments*3]
%Take magntude of all forces as some will be negative (only care about the maximum force)
[Force_max, index]=max(abs(Force_fil));			%[N] %Also obtain index of each force

%Pressure acting on vessel wall is force over unit filiment area
%These are absolute numbers - don't include any directionality
PressureR=abs(Force_max(1))/(height_PF);	%[Pa]
PressureZ=abs(Force_max(3))/(height_PF);	%[Pa]

%Stress acting on vessel wall is the combined force divided by the unit filiment area
%These are directional, some are negative and some are positive
StressR=(Force_fil(:,1))/(height_PF);		%[Pa]
StressZ=(Force_fil(:,3))/(height_PF);		%[Pa]
%Obtain maximum radial and axial stresses - either positive or negative
StressR_max=max(abs(StressR));				%[Pa]
StressZ_max=max(abs(StressZ));				%[Pa]


%%%%%%%%%%%%%%%%%%%%%% PLOT VESSEL EDDY STRESSES %%%%%%%%%%%%%%%%%%%%%%%%

%Scale stresses from [Pa] to [Atm]
StressR=StressR/2.0e5;
StressZ=StressZ/2.0e5;

%Plot figure showing vessel eddy stresses
close all
figure; hold on; axis equal;
plot(coilset);
%plot(vessel);
quiver(R_Lin_Array,Z_Lin_Array,StressR,StressZ,'color',[1 0 0],'AutoScale','off');
title('SMART Vessel Eddy-Stresses');
view(2) %2D view
legend(gca,'hide');
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_EddyStresses';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      COMPUTE PERTURBED EQUILIBRIA                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Apply small perturbation(s) to the efit_Geometry_init values
RGeo_Pert = efitGeometry_Init(1)+deltaRGeo;
ZGeo_Pert = efitGeometry(2)+deltaZGeo;
rGeo_Pert = efitGeometry_Init(3)+deltaAspect;
Kappa_Pert = efitGeometry_Init(4)+deltaKappa;
delta_Pert = efitGeometry_Init(5)+deltadelta;
%Construct new perturbed efit geometry
efitGeometry_Pert = [RGeo_Pert, ZGeo_Pert, rGeo_Pert, Kappa_Pert, delta_Pert];

%Recompute betaP and jprofile for the perturbed equilibrium - NOT USED DIRECTLY HERE
betaP_Pert = 3/2*ne*(Te+Ti)/(mu0*Ip/(2*pi*rGeo_Pert))^2*2*mu0*1.6e-19*Kappa_Pert;	% [%]
jprofile_Pert = fiesta_jprofile_topeol2( 'Topeol2', betaP_Pert, 1, li2, Ip );       % [-]
Irod_Pert = BT*2*pi*RGeo_Pert/mu0;                                                  % [A]

%EMPLOY FEEDBACK FOR PERTURBED EQUILIBRIUM?
%feedback = shape_controller(config, {'PF1','PF2','Div1','Div2'}, RGeo, ZGeo, rGeo, Kappa, delta);
%[efit_config, signals, weights, index] = efit_shape_controller(config, {'PF1','PF2','Div1','Div2'}, efitGeometry_Init);
%Calculate equilibrium fitting coil currents to provided jprofile
%equil = set(equil, config, 'feedback', feedback);
%EquilParams=parameters(equil);

%Compute perturbed equilibrium using the original jprofile and perturbed efit_Geometry
%PRETTY SURE THIS IS BACKWARDS, SHOULD USE PERTURBED JPROFILE WITH ORIGINAL efit_Geometry
%[efit_config, signals, weights, index] = efit_shape_controller(config, {'PF1','PF2'}, efit_Geometry_Pert);
%equil_pert = fiesta_equilibrium('ST', config, Irod, jprofile, control, efit_config, icoil_init, signals, weights);
%EquilParams_Pert = parameters(equil_pert);
equil_pert = Equil;                 % NOT CURRENTLY IN USE
CoilWaveforms = CoilWaveforms;      % NOT CURRENTLY IN USE    

%Save perturbed coil currents for all coils at time(5)
CoilCurrentsPert = transpose(CoilWaveforms(:,5));
icoil_pert = fiesta_icoil(coilset, CoilCurrentsPert);

%%%%%%%%%%%%%%%%%%%%%% PLOT PERTURBED EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%

%Plot perturbed equilibrium following convergence
%N.B. THIS IS NOT TECHNICALLY THE PERTURBED EQUILIBRIUM !!!!!
%BUT RATHER THE EFIT EQUILIBRIUM REQUIRED TO RETURN IT TO NORMAL
Title = {'SMART Perturbed Equilibrium \Psi(R,Z)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_PerturbedEquilibrium';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({equil_pert},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          DATA I/O MANAGEMENT                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create subdirectory for equilibrium related data
EquilDir = strcat(ASCIIDir,'Equil_Data/'); mkdir(EquilDir);

%Write 2D qeqdsk equilibrium file
Filename = strcat(EquilDir,'Equil.txt');
geqdsk_write_BUXTON(config, Equil, Filename);
%Write 1D equilibrium qprofile parameters file
qProfile = qprofile(Equil);
qProfileVariables = fieldnames(qProfile);
qProfileParams = struct2cell(qProfile(:,1));
Filename = strcat(EquilDir,'EquilProfiles.txt');
fileID=fopen(Filename,'w');
for i = 1:length(qProfileVariables)
    fprintf(fileID,'%s\r\n', string(qProfileVariables(i)));
    for j = 1:length(qProfileParams{i})
        try fprintf(fileID,'%1.13f\r\n', qProfileParams{i}(j));
        catch fprintf(fileID,'%1.13f\r\n', 'NaN');
        end 
    end
    fprintf(fileID,'%s\r\n', '*');
end
%Write 0D equilibrium parameters file
EquilParams=parameters(Equil);
ParamVariables = fieldnames(EquilParams);
ParamValues = struct2cell(EquilParams(:,1));
Filename = strcat(EquilDir,'EquilParam.txt');
fileID=fopen(Filename,'w');
for i = 1:length(ParamValues)
    try fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); ParamValues(i)]);
    catch fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); 'NaN']);
    end
end


%Write 2D qeqdsk perturbed equilibrium file
Filename = strcat(EquilDir,'Pert_Equil.txt');
geqdsk_write_BUXTON(config, Equil, Filename);
%Write 1D perturbed equilibrium qprofile parameters file
qProfile = qprofile(Equil);
qProfileVariables = fieldnames(qProfile);
qProfileParams = struct2cell(qProfile(:,1));
Filename = strcat(EquilDir,'Pert_EquilProfiles_Pert.txt');
fileID=fopen(Filename,'w');
for i = 1:length(qProfileVariables)
    fprintf(fileID,'%s\r\n', string(qProfileVariables(i)));
    for j = 1:length(qProfileParams{i})
        try fprintf(fileID,'%1.13f\r\n', qProfileParams{i}(j));
        catch fprintf(fileID,'%1.13f\r\n', 'NaN');
        end 
    end
    fprintf(fileID,'%s\r\n', '*');
end
%Write 0D perturbed equilibrium parameters file
EquilParams=parameters(Equil);
ParamVariables = fieldnames(EquilParams);
ParamValues = struct2cell(EquilParams(:,1));
Filename = strcat(EquilDir,'Pert_EquilParam.txt');
fileID=fopen(Filename,'w');
for i = 1:length(ParamValues)
    try fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); ParamValues(i)]);
    catch fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); 'NaN']);
    end
end


%Write 2D null-field equilibrium file
Filename = strcat(EquilDir,'Null_Equil.txt');
fileID=fopen(Filename,'w');
%Vacuum equilibria can't use geqdsk format - save as 2D array
Psi_Null = struct2cell(get(equil_null,'Psi_vac')); 
Psi_Null = Psi_Null(3); Psi_Null = Psi_Null{1,1};               %len(50200) = 200*251
Psi_Null = reshape(Psi_Null,[length(zGrid),length(rGrid)]);     %[rGrid,zGrid] = [200,251]
fprintf(fileID,'%s %s %s \n', '     Psi_Null', string(length(rGrid)), string(length(zGrid)));
for i = 1:size(Psi_Null,1)
    fprintf(fileID,'%g\t',Psi_Null(i,:));
    fprintf(fileID,'\n');
end
%Write Null Bpol as 2D array
Filename = strcat(EquilDir,'Null_Bpol.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s \n', '     Bpol_Null', string(length(rGrid)), string(length(zGrid)));
for i = 1:size(BpolData_Null,1)
    fprintf(fileID,'%g\t',BpolData_Null(i,:));
    fprintf(fileID,'\n');
end
%Write Null Btor as 2D array
Filename = strcat(EquilDir,'Null_Bpol.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s \n', '     Btor_Null', string(length(rGrid)), string(length(zGrid)));
for i = 1:size(BtorData_Null,1)
    fprintf(fileID,'%g\t',BtorData_Null(i,:));
    fprintf(fileID,'\n');
end
%{
%Write 1D null-field equilibrium qprofile parameters file
qProfile = qprofile(equil_optimised_null)
qProfileVariables = fieldnames(qProfile);
qProfileParams = struct2cell(qProfile(:,1));
Filename = strcat(EquilDir,'Null_EquilProfiles.txt');
fileID=fopen(Filename,'w');
for i = 1:length(qProfileVariables)
    fprintf(fileID,'%s\r\n', string(qProfileVariables(i)));
    for j = 1:length(qProfileParams{i})
        try fprintf(fileID,'%1.13f\r\n', qProfileParams{i}(j));
        catch fprintf(fileID,'%1.13f\r\n', 'NaN');
        end 
    end
    fprintf(fileID,'%s\r\n', '*');
end
%Write 0D null-field equilibrium parameters file
EquilParam=parameters(equil_optimised_null);
ParamVariables = fieldnames(Null_EquilParam);
ParamValues = struct2cell(EquilParam(:,1));
Filename = strcat(EquilDir,'EquilParam_Null.txt');
fileID=fopen(Filename,'w');
for i = 1:length(ParamValues)
    try fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); ParamValues(i)]);
    catch fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); 'NaN']);
    end
end
%}

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Write efit geometry and perturbed efit geometry (if applicable)
Filename = strcat(EquilDir,'efit_Geometry_Init.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'RGeo','ZGeo','a','kappa','delta');
fprintf(fileID,'%1.12f %1.12f %1.12f %1.12f %1.12f\r\n',[efitGeometry_Init(1)'; efitGeometry_Init(2)'; efitGeometry_Init(3)'; efitGeometry_Init(4)'; efitGeometry_Init(5)']);

Filename = strcat(EquilDir,'efit_Geometry_Equil.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'RGeo','ZGeo','a','kappa','delta');
fprintf(fileID,'%1.12f %1.12f %1.12f %1.12f %1.12f\r\n',[efitGeometry(1)'; efitGeometry(2)'; efitGeometry(3)'; efitGeometry(4)'; efitGeometry(5)']);
    
Filename = strcat(EquilDir,'efit_Geometry_Pert.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'RGeo','ZGeo','a','kappa','delta');
fprintf(fileID,'%1.12f %1.12f %1.12f %1.12f %1.12f\r\n',[efitGeometry_Pert(1)'; efitGeometry_Pert(2)'; efitGeometry_Pert(3)'; efitGeometry_Pert(4)'; efitGeometry_Pert(5)']);

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Create subdirectory for coil current related data
%%
icoilDir = strcat(ASCIIDir,'icoil_Data/'); mkdir(icoilDir);

Filename = strcat(icoilDir,'icoil_position.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s\r\n', 'Coil','R [m]  ','Z [m]');
fprintf(fileID,'%s %0.5f %0.5f\r\n', 'Sol ',RSol,ZMaxSol);
fprintf(fileID,'%s %0.5f %0.5f\r\n', 'PF1 ',R_PF1,Z_PF1);
fprintf(fileID,'%s %0.5f %0.5f\r\n', 'PF2 ',R_PF2,Z_PF2);
fprintf(fileID,'%s %0.5f %0.5f\r\n', 'Div1',R_Div1,Z_Div1);
fprintf(fileID,'%s %0.5f %0.5f\r\n', 'Div2',R_Div2,Z_Div2);

Filename = strcat(icoilDir,'Init_icoil.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol','PF1','PF2','Div1','Div2');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_init.Sol'; icoil_init.PF1'; icoil_init.PF2'; icoil_init.Div1'; icoil_init.Div2']);

Filename = strcat(icoilDir,'efit_icoil.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol','PF1','PF2','Div1','Div2');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_init.Sol'; icoil_efit.PF1'; icoil_efit.PF2'; icoil_efit.Div1'; icoil_efit.Div2']);

Filename = strcat(icoilDir,'Perturbed_icoil.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol','PF1','PF2','Div1','Div2');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_pert.Sol'; icoil_pert.PF1'; icoil_pert.PF2'; icoil_pert.Div1'; icoil_pert.Div2']);

Filename = strcat(icoilDir,'Null_icoil.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol','PF1','PF2','Div1','Div2');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_null.Sol'; icoil_null.PF1'; icoil_null.PF2'; icoil_null.Div1'; icoil_null.Div2']);

%Extract coil current time-traces
Sol=I_PF_output(:,1);   %Sol
PF1=I_PF_output(:,2);   %PF1
PF2=I_PF_output(:,3);   %PF2
Div1=I_PF_output(:,4);  %Div1
Div2=I_PF_output(:,5);  %Div2
Filename = strcat(icoilDir,'CoilCurrents.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s %s\r\n', 'time_adaptive','I_Sol','I_PF1','I_PF2','I_Div1','I_Div2');
fprintf(fileID,'%1.12f %0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[time_adaptive'; Sol'; PF1'; PF2'; Div1'; Div2']);

%Extract coil voltage time-traces
Sol=V_PF_output(:,1);     %Sol
PF1=V_PF_output(:,2);     %PF1
PF2=V_PF_output(:,3);     %PF2
Div1=V_PF_output(:,4);     %Div1
Div2=V_PF_output(:,5);     %Div2
Filename = strcat(icoilDir,'CoilVoltages.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s %s\r\n', 'time_adaptive','I_Sol','I_PF1','I_PF2','I_Div1','I_Div2');
fprintf(fileID,'%1.12f %0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[time_adaptive'; Sol'; PF1'; PF2'; Div1'; Div2']);

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Create subdirectory for dynamic current data (RZIP)
RZIPDir = strcat(ASCIIDir,'RZIP_Data/'); mkdir(RZIPDir);

Filename = strcat(RZIPDir,'time.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s\r\n','time_adaptive');
fprintf(fileID,'%1.12f\r\n',time_adaptive);

Filename = strcat(RZIPDir,'Ip.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'time_adaptive','Ip_output');
fprintf(fileID,'%1.12f %1.12f\r\n',[time_adaptive'; Ip_output']);

Filename = strcat(RZIPDir,'IPass.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'time_adaptive','Net_I_Passive');
fprintf(fileID,'%1.12f %1.12f\r\n',[time_adaptive'; Net_IPassive']);

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Misc outputs, save unordered in main RawData directory

Filename = strcat(ASCIIDir,'Eta.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'Eta_Perp', 'Eta_Para');
fprintf(fileID,'%1.12f %1.12f\r\n', PlasmaResistPerp', PlasmaResistPara');

Filename = strcat(ASCIIDir,'Bpol.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'Null_Bpol', 'Null_Btor');
fprintf(fileID,'%1.12f %1.12f\r\n', BpolAvg_Null', BtorAvg_Null');

Filename = strcat(ASCIIDir,'betaP.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'betaP', 'betaP_Pert');
fprintf(fileID,'%1.12f %1.12f\r\n', betaP', betaP_Pert');

Filename = strcat(ASCIIDir,'LCon.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s \r\n', 'Lc');
fprintf(fileID,'%1.12f\r\n', Lc');

Filename = strcat(ASCIIDir,'MaxStress.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'StressR_max', 'StressZ_max');
fprintf(fileID,'%1.12f %1.12f\r\n', 'Stress_Rmax', StressR_max', StressZ_max');

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Done!
disp([ 'Done!' ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% ADDING EDDIES TO DISCHARGE EQUIL %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%THIS SECTION MUST BE RUN AFTER THE INITIAL SECTION
%THIS WILL BECOME THE EDDY LOOP - PROBABLY EASIER THAN FULL LOOP

%%%%%%%%%%%%%%%%%%%  RECOMPUTE DISCHARGE EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%

%Book-keeping for the start of each loop
close all       %Ensure all open figures are closed

%Extract previously calculated efit coil currents without eddies
CoilCurrents = transpose(CoilWaveforms(:,TimeIndex_Discharge)); %n=5, coil filaments 
%NEED A CATCH-RETRY SECTION FOR UPDATING THE GUESSES (General rule - Slightly increase IDiv2 and retry)
CoilCurrents(iDiv2) = I_Div2_Equil;     %900;
%Combine efit coil currents and vessel eddy currents into new array
CoilAndVesselCurrents = [CoilCurrents, VesselEddyCurrents];     %n=313, coil + vessel filaments

%Create new fiesta coilset and configuration to re-perform efit
global coilvesselset
coilvesselset = fiesta_loadassembly(coilset, vessel);                          %Creates object of size coil and vessel currents ????    
config_passive = fiesta_configuration('Config_Passive', Grid, coilvesselset);  %Creates config including coil and vessel filaments
icoil_passive = fiesta_icoil( coilvesselset, CoilAndVesselCurrents );          %Creates icoil object with coil and vessel currents

%Compute equilibrium (Psi(R,Z)) from the supplied jprofile, icoil and geometry
%Returns target equilibrium and CoilWaveforms for PF1 and PF2 at requested time_Index
%[Equil_Passive,EquilParams_Passive,CoilWaveforms_Passive,efitGeometry_Passive,config_passive] = ... %% ITERATION 2 CONFIG %%
%    efitInverse(jprofile,Irod,CoilWaveforms,efitGeometry_Init,TimeIndex_Discharge);
%NEED TO EITHER SUPPLY THE COILSET TO THE EFIT FUNCTION OR COMPUTE IT OUTSIDE OF THE FUNCTION FOR ALL CASES
%FIRST JOB IS TO GET THE EFIT FUNCTION CAPABLE OF WORKING WITH DIFFERENT CONFIGS, COILSETS AND VESSEL FILAMENTS

%Recompute efit including eddy currents using efit_config_passive and icoil_passive (maintain original efit_Geometry_Init)
control = fiesta_control('diagnose',true, 'quiet',false, 'convergence',1e-5, 'boundary_method',2);
[efit_config_passive, signals_passive, weights_passive, index_passive] = efit_shape_controller(config_passive, {'PF1','PF2'}, efitGeometry_Init);
Equil_Passive = fiesta_equilibrium('SMART_Passive', config_passive, Irod, jprofile, control, efit_config_passive, icoil_passive, signals_passive, weights_passive);
EquilParams_Passive = parameters(Equil_Passive);

%Extract the new coil currents from the efit-equilibrium:
icoil_efit_passive = get(Equil_Passive,'icoil'); 
CoilCurrentsEfit_Passive = get(icoil_efit_passive,'currents');
%NEED TO UPDATE THE COIL WAVEFORMS HERE TOO - USUALLY INSIDE EFIT FUNCTION

%Initiate virtual B-field sensors centered on Rgeo
sensor_btheta_passive = InitiateBSensors(EquilParams_Passive,a_eff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSSIBLE FIX FOR THE RZIP ISSUE??
%
%CoilCurrents = CoilCurrentsEfit_Passive(1:nPF)
%VesselCurrents = CoilCurrentsEfit_Passive(nPF:313)
%
%MAKE ALL VESSEL FILAMENTS INTO COILS WITH THEIR RESPECTIVE CURRENTS
%icoil_Vessel =  fiesta_icoil( nPF + length(vessel)) %313 'coils'
%for i=1:313
%   icoil_vessel(i) = CoilCurrentsEfit_Passive(i)
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NEED TO RE-CALCULATE PLASMA AND EDDY WITH NEW EFIT EQUILIBRIUM
%ADDING EDDY CURRENTS VARIES THE PLASMA CURRENT - 
%rzip_config = fiesta_rzip_configuration( 'RZIP',config_passive,vessel,{sensor_btheta_passive} );   %CONFIG ACCOUNTS FOR VESSEL FILAMENTS
%A_Passive,B_Passive,C_Passive,D_Passive,curlyM_Passive,curlyR_Passive,gamma_Passive,plasma_parameters_Passive,index_Passive,label_index_Passive,state_Passive] = ...
%    response(rzip_config, Equil_Passive, 'rp', PlasmaResistPerp);

%%%%%%%%%%%%%%%%%%%%%  COMPUTE OPTIMISED NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%%

%Update CoilWaveforms array with null-field values (Using NaN Mask)
%   RZIP_C = C;
%   CoilWaveforms = NullFieldWaveforms(CoilWaveforms, RZIP_C, sensor_btheta_passive, TimeIndex_NullField);
%ISSUE :: RZIP_C is computed from RZIP - Need recomputed RZIP with eddy currents

%Extract previously calculated efit coil currents without eddys
CoilCurrentsNull = transpose(CoilWaveforms(:,TimeIndex_NullField)); %Null-field coil currents without eddys
CoilAndVesselCurrents = [CoilCurrentsNull, VesselEddyCurrents];     %n=313, coil + vessel filaments
icoil_null_passive = fiesta_icoil(coilvesselset, CoilAndVesselCurrents);

%Compute null-field equilibrium using null-field coil and vessel eddy currents
equil_null_passive = fiesta_equilibrium('SMART-Null', config_passive, Irod, icoil_null_passive);
EquilParams_Null_Passive = parameters(equil_null_passive);
%ISSUE? :: THIS GIVES VERY UNUSUAL COIL CURRENTS REQUIRED FOR NULL-FIELD

%Extract the new coil currents from the null-field equilibrium:
icoil_null_passive = get(equil_null_passive,'icoil'); 
CoilCurrentsNull_Passive = get(icoil_null_passive,'currents');

%%%%%%%%%%%%%%%%%%  COMPUTE BREAKDOWN FOR NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%

%Extract the null poloidal and toroidal B-field vector arrays
[BrData_Null_Passive,BzData_Null_Passive,BPhiData_Null_Passive,BpolData_Null_Passive,BtorData_Null_Passive] = ...
    ExtractBField(equil_null_passive);

%Minimum null poloidal and toroidal fields, averaged over region of area a_eff^2
[BpolAvg_Null_Passive,BtorAvg_Null_Passive] = ...
    ExtractNullBMin(EquilParams_Passive,BpolData_Null_Passive,BtorData_Null_Passive,a_eff);
    %NOTE - USES EQUILPARAMS_PASSIVE ONLY TO EXTRACT RGeo and ZGeo (Would be less confusing to use EQUILPARAMS_NULL)

%Compute the average connection length within the null-field sensor region
%Lc_Passive = ConnectionLength(EquilParams,VesselRMaxInner,VesselRMinInner,VesselZMaxInner,VesselZMinInner);
Lc_Passive = 0.25*a_eff*(BtorAvg_Null_Passive/BpolAvg_Null_Passive);

%%%%%%%%%%%%%%%  COMPUTE DYNAMIC PLASMA & EDDY CURRENTS  %%%%%%%%%%%%%%%%%%

%Recompute dynamic plasma and vessel eddy currents with new coil waveforms
%ISSUE :: curlyM and curlyR are computed from RZIP - Need recomputed RZIP with eddy currents

%Update dynamic plasma and vessel eddy currents
%[time_linear_passive,time_adaptive_passive,I_PF_output_Passive,V_PF_output_Passive,Ip_output_Passive,Vp_output_Passive,I_Passive_Passive] = ...
%    DynamicCurrents(CoilWaveforms_Passive, time, curlyM_Passive, curlyR_Passive);

%Extract Vessel Eddy Currents during discharge (time='false' for absolute max)
%VesselEddyCurrents_Passive = ...
%    ExtractPassiveCurrents(I_Passive_Passive,time_adaptive_Passive,time(TimeIndex_Discharge));

%%%%%%%% !!!! RESTART THE EDDY LOOP SECTION AT THIS POINT !!!! %%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% PLOT TARGET EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%%%

%Plot target equilibrium following convergence
Title = {'SMART Target Equilibrium iter(1)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_TargetEquilibrium_Passive';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({Equil_Passive},{rGrid,zGrid},Title,CbarLabel,SaveString);

CoilCurrentsEfit(1:nPF)
CoilCurrentsEfit_Passive(1:nPF)

%%%%%%%%%%%%%%%%%%%%  PLOT NULL-FIELD PHI SURFACES  %%%%%%%%%%%%%%%%%%%%%

%Plot the optimised null-field phi
Title = {'SMART Null-field Equilibrium iter(1)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_NullPhi_Passive';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({equil_null_passive},{rGrid,zGrid},Title,CbarLabel,SaveString);

CoilCurrentsNull(1:nPF)
CoilCurrentsNull_Passive(1:nPF)

%%%%%%%%%%%%%%%%%%%%%% PLOT NULL-FIELD BPOL  %%%%%%%%%%%%%%%%%%%%%

%Log poloidal and toroidal magnetic fields to show details (Sol Obscures)
logBpolData_Null_Passive = log(BpolData_Null_Passive);
logBtorData_Null_Passive = log(BtorData_Null_Passive);
%Plot the optimised null-field phi
Title = {'SMART Null-field iter(1)',' '};
CbarLabel = 'Null-field B_{\theta} ln([T])';
Filename = '__NullBpol_Passive';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({logBpolData_Null_Passive},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%      FUNCTIONS      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Equilibrium,EquilParams,OutputCoilWaveforms,efitGeometry,config]= ...
    efitInverse(Jprofile,Irod,InputCoilWaveforms,InputGeometry,TimeIndex)

    %Obtain required global variables
    global Grid;
    global coilset;
    global efitCoils;
    global iPF1; global iPF2;
    global iDiv1; global iDiv2;
    global iSol;
    
    %Initiate icoil object and set currents to desired TimeIndex
    CoilCurrents = transpose(InputCoilWaveforms(:,TimeIndex));
    icoil = fiesta_icoil(coilset, CoilCurrents);
    
    %Initiate efit configuration - default configuration fixed for now
    config = fiesta_configuration('SMART', Grid, coilset);
    control = fiesta_control('diagnose',true, 'quiet',false, 'convergence',1e-5, 'boundary_method',2);
    
    %Efit outputs coil currents resulting from the supplied jprofile, icoil and geometry
	%Returns new currents for the requested coils: {'Coil1, {...}, 'Coiln'}
	[efit_config, signals, weights, index] = efit_shape_controller(config, efitCoils, InputGeometry);

	%Inverse equilibrium, outputs coil currents resulting in the supplied jprofile and icoil config
	Equilibrium = fiesta_equilibrium('SMART', config, Irod, Jprofile, control, efit_config, icoil, signals, weights);
	EquilParams = parameters(Equilibrium);

	%Extract achieved efit equilibrium geometry values 
    RGeo = EquilParams.r0_geom; ZGeo = EquilParams.z0_geom;
    AspectRatio = EquilParams.aspectratio; Rminor = (RGeo/AspectRatio);
    Elongation = EquilParams.kappa; Triangularity = EquilParams.delta;
    %Compile efit output geometry
	efitGeometry = [RGeo,ZGeo,Rminor,Elongation,Triangularity];

	%Extract the new coil currents from the efit-equilibrium:
    OutputCoilWaveforms = InputCoilWaveforms;                %Initiate output current waveforms
	icoil_efit = get(Equilibrium,'icoil'); 	CoilCurrents_efit = get(icoil_efit,'currents');
	OutputCoilWaveforms(2,5:6) = CoilCurrents_efit(iPF1);	 %Assumes IPF1 is flat over equilibrium
	OutputCoilWaveforms(3,5:6) = CoilCurrents_efit(iPF2);    %Assumes IPF2 is flat over equilibrium
%	OutputCoilWaveforms(4,5:6) = CoilCurrents_efit(iDiv1);   %Need to auto-select which coils to update
	OutputCoilWaveforms(5,5:6) = CoilCurrents_efit(iDiv2);   %Assumes IDiv2 is flat over equilibrium
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute time_resolved plasma and vessel eddy currents for given input coil waveforms
function [Time_Linear,Time_Adaptive,I_PF_output,V_PF_output,Ip_output,Vp_output,I_Passive]= ...
    DynamicCurrents(CoilWaveforms,TimeIndices,CurlyM,CurlyR)

    %Obtain required global variables
    global coilturns;
    global iPF1; global iPF2;
    global iDiv1; global iDiv2;
    global iSol;
    
    %Name and colour coils for plotting
    coil_names{iSol} = 'Sol'; PF_colors{iSol} = 'Red';
    coil_names{iPF1} = 'PF1'; PF_colors{iPF1} = 'Magenta';
    coil_names{iPF2} = 'PF2'; PF_colors{iPF2} = 'Black';
    coil_names{iDiv1} = 'Div1'; PF_colors{iDiv1} = 'Cyan';
    coil_names{iDiv2} = 'Div2'; PF_colors{iDiv2} = 'Green';
    
    %Determine number of coil waveforms and time points within each
    SizeCoilArrays = size(CoilWaveforms);
    nCoils = SizeCoilArrays(1);             %Number of PF/Div coils
    nTime = SizeCoilArrays(2);              %Number of TimeVertics

    %Initiate RZip PF Arrays used to calculate plasma and eddy currents
    IPFinput_Discrete = transpose(CoilWaveforms);        %Discrete Coil currents from efit
    VPFinput_Discrete = NaN(nTime,nCoils);               %Discrete Coil voltages initiated to zero

    %Convert from discrete time-points to a continuous time-axis
    TemporalResolution = 1000;
    Time_Linear = linspace(min(TimeIndices),max(TimeIndices),TemporalResolution);
    %Construct input PF/Div coil waveforms in continuous time-axis
    IPFinput_Continous = NaN(length(Time_Linear),nCoils);
    for i=1:nCoils
        IPFinput_Continous(:,i) = interp1(TimeIndices,IPFinput_Discrete(:,i),Time_Linear);
    end
    VPFinput_Continous = NaN*IPFinput_Continous;

    %Initiate Plasma Currrent and Plasma Potential arrays:
    Ip_long = zeros(size(Time_Linear));         %Sets Ip_long to zero array
    Vp_long = NaN(size(Time_Linear));           %Sets Vp_long to 'NaN' array (Current Driven)

    %Compute dynamic coil currents employing Current Criven Ip
	%CurlyM and CurlyR are large inductance and resistance matrices.
    [V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, Time_Adaptive ] = ...
        state_space_including_passive_elements_v4( CurlyM, CurlyR, Time_Linear, IPFinput_Continous, VPFinput_Continous, Ip_long, Vp_long, 'adaptive_timesteping',true );

    %Set breakdown time and prepare Ip_long and Vp_long for voltage driven Ip
    %!!! NEED Time_Breakdown TO INCLUDE BREAKDOWN AND BURNTHROUGH TIME !!!
    Time_Breakdown = 0;                                          %Set time for plasma breakdown (default 0)
    Time_Plasma = Time_Adaptive > Time_Breakdown;                %Set times for which plasma exists
    Vp_output(Time_Plasma) = 0;                                  %Set voltage to zero when plasma exists
    Vp_long = interp1(Time_Adaptive, Vp_output, Time_Linear);    %Sets Vp_long = 0 when Time_Linear > 0.
    Ip_long = NaN*Vp_long;                                       %Sets Ip_long to 'NaN' array (Voltage Driven)

    %Compute dynamic coil currents employing Voltage Criven Ip
    [ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, Time_Adaptive ] = ...
        state_space_including_passive_elements_v4( CurlyM, CurlyR, Time_Linear, IPFinput_Continous, VPFinput_Continous, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names', coil_names, 'show_plot',true, 'turns',coilturns, 'currentScale',1e3, 'PF_colors',PF_colors );
    close all
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CoilWaveformsOutput]=NullFieldWaveforms(CoilWaveformsInput,RZIP_C,sensor_btheta,TimeIndex)
   
    %Obtain required global variables
    global iPF1; global iPF2;
    global iDiv1; global iDiv2;
    global iSol; 
    
    %Determine number of coil waveforms and time points within each
    SizeCoilArrays = size(CoilWaveformsInput);
    nCoils = SizeCoilArrays(1);             %Number of PF/Div coils
    nTime = SizeCoilArrays(2);              %Number of TimeVertics

    %Extract scaling factors for null-field coil currents
    %Cn is the part of the matrix C related to the sensors (see response)
    C_temp = RZIP_C(end-get(sensor_btheta,'n')+1:end,1:nCoils);
    C1 = C_temp(:,1);          %Elements of C_temp(Cn) for Sol coil
    
    D1_PF1 = C_temp(:,iPF1);   %Elements of C_temp(Cn) for PF1 coil
    D1_PF2 = C_temp(:,iPF2);   %Elements of C_temp(Cn) for PF2 coil
    D1_Div1 = C_temp(:,iDiv1); %Elements of C_temp(Cn) for Div1 coil
    D1_Div2 = C_temp(:,iDiv2); %Elements of C_temp(Cn) for Div2 coil
    
    %Determine if Div1 is in series with solenoid or not and optimise null-field accordingly
    if isnan(CoilWaveformsInput(iDiv1,TimeIndex)) == true                  %If Div1 NOT in series Sol
        %Scale ALL null-field coil currents relative to Solenoid current
        D1=[D1_PF1, D1_PF2, D1_Div1, D1_Div2];                              %Optimise for PF1,2 & Div1,2
        ISolNullField = CoilWaveformsInput(iSol,TimeIndex);                 %Extract Solenoid Current 
        IPF_null = -pinv(D1) * (C1*ISolNullField);                          %Scale null-field currents
        IPF_null = [ISolNullField,transpose(IPF_null(:))];                  %Add Sol into IPF_null
    
    elseif isnan(CoilWaveformsInput(iDiv1,TimeIndex)) == false               %If Div1 IS in series with Sol
        %Scale ALL EXCEPT Div1 relative to Solenoid current:
        D1=[D1_PF1, D1_PF2, D1_Div2];                                       %Optimise for PF1,2 & Div2
        ISolNullField = CoilWaveformsInput(iSol,TimeIndex);                 %Extract Solenoid Current 
        IPF_null = -pinv(D1) * (C1*ISolNullField + D1_Div1*ISolNullField);  %Scale null-field currents
        IPF_null = [ISolNullField,IPF_null(1),IPF_null(2),ISolNullField,IPF_null(3)]; %Add Sol and Div1 into I_PF_null
    end
    
    %Update CoilWaveforms array with null-field values (Using NaN Mask)
    for i = 1:nCoils
        for j = 1:nTime
            %Determine if coil 'i' at timestep 'j' requires null-field
            if isnan(CoilWaveformsInput(i,j))
                CoilWaveformsInput(i,j) = IPF_null(i);
            end
        end
    end
    %Update output coil waveform array
    CoilWaveformsOutput = CoilWaveformsInput;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initiate virtual sensors within null-field region
function SensorBTheta=InitiateBSensors(EquilParams,a_eff)

    %Determine number of sensors (constant for now)
    NumSensors = 10;

    %Define null field region as defined by the equilibrium RGeo and ZGeo
    %Null-field region radius taken as effective minor radius (a_eff)
    RGeo = EquilParams.r0_geom; ZGeo=EquilParams.z0_geom;
    BP_virt_R = linspace(RGeo-a_eff,RGeo+a_eff,NumSensors);
    BP_virt_Z = linspace(ZGeo-a_eff,ZGeo+a_eff,NumSensors);

    %Create null field region grid
    [BP_virt_R,BP_virt_Z] = meshgrid(BP_virt_R,BP_virt_Z);
    BP_virt_R = BP_virt_R(:)';
    BP_virt_Z = BP_virt_Z(:)';

    %Create sensors over the null field region
    BP_virt_theta = zeros(1,length(BP_virt_R));
    nSensors = length(BP_virt_theta);

    %Create array of unique sensor names (technically required...)
    BP_virt_names = {};
    for iSensor=1:nSensors
        BP_virt_names{iSensor} = ['Radial Bp Virtual Sensor #' num2str(iSensor) ];
    end

    %R and Z of Dim[1*200] and cyclical (i.e. BP_virt_R[201] == BP_virt_R [1]) 
    BP_virt_R = [BP_virt_R  BP_virt_R];
    BP_virt_Z = [BP_virt_Z  BP_virt_Z];
    %Theta of Dim[1*200] and: ???"The first 100 have 0, and the second has pi/2"???
    BP_virt_theta = [BP_virt_theta  BP_virt_theta+pi/2];

    %Name each sensor and create a FIESTA sensor object for use in RZIP
    for iSensor=nSensors+1:2*nSensors
        BP_virt_names{iSensor} = ['Vertical Bp Virtual Sensor #' num2str(iSensor) ];
    end
    
    %Initiate sensors with sensor names over defined region
    SensorBTheta = fiesta_sensor_btheta('BSensors',BP_virt_R,BP_virt_Z,BP_virt_theta,BP_virt_names);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extract Passive Vessel currents at desired time [s] during the pulse
%If discharge time is supplied as 'false' then absolute maximum currents are extracted
function VesselEddyCurrents=ExtractPassiveCurrents(I_Passive,time_adaptive,DischargeTime)

    %NOTES
    %INPUT  :: Coilwaveforms with optimised null-field and efit discharge currents
    %       :: Computed from RZIP using equil_null and efit using equil_efit, respectively
    %INPUT  :: Variables curlyM and curlyR contain vessel information
    %       :: Computed from RZIP using equil_efit
    %OUTPUT :: I_Passive contains eddy current of the nf filaments at each instant of time.
    %       :: I_Passive filaments are [WallThickness x WallThickness] square by default 
    %       :: (len(I_Passive) = 3811*nfilaments == len(time_adaptive) = 3811*1

    if isfloat(DischargeTime) == false
        %Time intervals intersected with number of filaments (time intervals*number of filaments)
        SizeIPas=size(I_Passive);
        
        %For each vessel filament extract the greatest absolute current
        for i=1:SizeIPas(2)
            %Obtain largest positive and negative for each filament
            positive=max(I_Passive(:,i));           %I_Passive(Timepoint,Filament)          
            negative=min(I_Passive(:,i));           %I_Passive(Timepoint,Filament)
            %Keep the largest absolute value
            if abs(positive)> abs(negative)
                I_Passive_Abs(i)=positive;
            else
                I_Passive_Abs(i)=negative;
            end
        end
        VesselEddyCurrents = I_Passive_Abs;         %1D array of max vessel eddy currents
        
    %%%%%     %%%%%     %%%%%     %%%%%
        
    elseif isfloat(DischargeTime) == true
        %Time intervals intersected with number of filaments (time intervals*number of filaments)
        SizeIPas=size(I_Passive);

        %Find adaptive_time index closest to the user requested time
        DischargeIndices = find(time_adaptive>DischargeTime);    %Allows for non-identical times
        DischargeIndex = DischargeIndices(1);                    %First index is close enough
        %Extract filiment eddy currents at desired time
        for i=1:SizeIPas(2)                                      %SizeIPas(2) = 696 vessel filaments 
            I_Passive_AtTime(i)=I_Passive(DischargeIndex,i);     %I_Passive(Timepoint,Filament)
        end
        VesselEddyCurrents = I_Passive_AtTime;      %1D array of vessel eddy currents at specific time
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILITY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BpolAvg,BtorAvg]=ExtractNullBMin(EquilParams,BpolData,BtorData,a_eff)

    %Obtain required global variables
    global Grid;

    %Null-field region centre is at (RGeo, ZGeo) to align with sensor_btheta
    RGeo = EquilParams.r0_geom; ZGeo=EquilParams.z0_geom;
    
    %Convert from SI to cell index notation to enable averaging
    RAxis = get(Grid,'r'); ZAxis = get(Grid,'z');
    IndexRGeo = find(RGeo < RAxis); IndexRGeo = IndexRGeo(1);
    IndexZGeo = find(ZGeo < ZAxis); IndexZGeo = IndexZGeo(1);

    %Determine null-field index range; radius of a_eff around (R, Z)
    CellRange = ceil(a_eff/0.0055)/2.0;     %Null-field region radius in cells
    CellRangeR = [IndexRGeo-CellRange, IndexRGeo+CellRange];
    CellRangeZ = [IndexZGeo-CellRange, IndexZGeo+CellRange];
    
    %Resize Bpol and Btor arrays to null-region
    BpolData_NullRegion = BpolData(CellRangeZ(1):CellRangeZ(2), CellRangeR(1):CellRangeR(2));
    BtorData_NullRegion = BtorData(CellRangeZ(1):CellRangeZ(2), CellRangeR(1):CellRangeR(2));

    %Find min, max and average values within null-field region
    BpolMin = min(BpolData_NullRegion,[],'all'); BpolMax = max(BpolData_NullRegion,[],'all');
    BtorMin = min(BtorData_NullRegion,[],'all'); BtorMax = max(BtorData_NullRegion,[],'all');
    BpolAvg = mean(BpolData_NullRegion,'all');  %[T]
    BtorAvg = mean(BtorData_NullRegion,'all');  %[T]
    
    %LEGACY CODE - MAINTAINED WITH A VERY SMALL BMin
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BMin = 1.0E-10; %[T]
    %Enforce lower limit for BpolMin (default to Earth's B-field (5.0E-5 [T])
    %Song2017 suggests Bpolmin as 0.2mT -> 1mT (2e-4 -- 1e-3)
    if BpolAvg < BMin;
        BpolAvg = BpolAvg+BMin;	
    elseif BtorAvg < BMin;
        BtorAvg = BtorAvg+BMin;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BrData,BzData,BphiData,BpolData,BtorData]=ExtractBField(Equilibrium)

    %Extract the null poloidal field
    BrData = get(get(Equilibrium,'Br'),'data');
    BzData = get(get(Equilibrium,'Bz'),'data');
    BphiData = get(get(Equilibrium,'Bphi_vac'),'data');
    rGrid = get(get(get(Equilibrium,'Br'),'grid'),'r');
    zGrid = get(get(get(Equilibrium,'Br'),'grid'),'z');

    %Reshape into 2D to match simulation mesh grid (Z,R) 
    BrData = reshape( BrData, length(zGrid), length(rGrid) );
    BzData = reshape( BzData, length(zGrid), length(rGrid) );
    BphiData = reshape( BphiData, length(zGrid), length(rGrid) );

    %Extract poloidal and toroidal magnetic vectors (1D arrays)
    BpolData = sqrt(BzData.^2+BrData.^2);       %Compute BPoloidal vector
    BtorData = sqrt(BphiData.^2);               %Compute BToroidal vector
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! NEEDS TESTING !!!
function Lc=ConnectionLength(EquilParams,VesselRMax,VesselRMin,VesselZMax,VesselZMin)

    %Define Location of Inner Vessel Walls (four corners)
    RadialCorners = [VesselRMin, VesselRMax, VesselRMax, VesselRMin];
    AxialCorners = [VesselZMin, VesselZMin, VesselZMax, VesselZMax];
    
    %Define starting location for connection length (default RGeo,ZGeo) 
    RPoint = EquilParams.r0;
    ZPoint = EquilParams.z0;
    
    %Convert StartLoc into fiesta_point and InnerWall into fiesta_line
    StartLoc = fiesta_point('Start', RPoint, ZPoint);
    InnerWall = fiesta_line('InnerWall', RadialCorners, AxialCorners);

    %Compute connection length from StartLoc to intersection at any point on InnerWall
    %!!!! Undefined function 'line_intersect' for input arguments of type 'fiesta_line' !!!!
    [length_3d, length_2d, connection, phi, path_3d, path_2d] = ...
        connection_length2(equil_optimised_null, StartLoc, InnerWall)
    %!!!! Undefined function 'line_intersect' for input arguments of type 'fiesta_line' !!!!
    
    %Update output variables
    Lc = connection;
end
% !!! NEEDS TESTING !!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Specify or extrapolate coil current for a given current ramp timescale
function CoilRampCurrent=FitSolenoidRamp(CoilRampCurrents,TimeVertices)

    %Extract relevent coil currents
    CoilStartCurrent = CoilRampCurrents{1};
    CoilRampCurrent = CoilRampCurrents{2};
    CoilEndCurrent = CoilRampCurrents{3};

    %Extrapolate a linear ramp-down:
    if strcmp(CoilRampCurrent, 'Linear');
        %Maintain a linear solenoid ramp-down from time(4), through time(5) to time (6)
        %Apply a linear fit to the solenoid ramp-down profile between PrePulse to Equil
        [coef] = polyfit([TimeVertices(3), TimeVertices(5)], [CoilStartCurrent, CoilEndCurrent], 1);
        %Extrapolate solenoid current when PF and Div coils reach equilibrium values, at time(4)
        CoilRampCurrent = (coef(1)*TimeVertices(4)) + coef(2);
   
    %Employ user specified value if requested
    elseif isfloat(CoilRampCurrent);
        CoilRampCurrent = double(CoilRampCurrent);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig=PlotEquilibrium(Arrays,Grid,Title,CbarLabel,SaveString)

    %USEFUL FIESTA CLASSES
    %class(sensor_btheta) = 'fiesta_sensor_btheta'
    %class(equil) = 'fiesta_equilibrium'
    
    %Obtain required global variables
    global vessel; Vessel = vessel;
    global coilset; Coilset = coilset;

    %Initiate a Clean Figure
    close all
    figure; hold on; axis equal;
    
    %for each supplied sub-array (Arrays{i})
    for i=1:length(Arrays);
        %If data is a FIESTA class then use in-built function
        if isa(Arrays{i},'fiesta_equilibrium') == true;
            plot(Arrays{i});
        elseif isa(Arrays{i},'fiesta_sensor_btheta') == true;
            plot(Arrays{i});
        %Else plot as a regular contour plot using the supplied grid
        else
            contourf(Grid{1},Grid{2},Arrays{i});
        end
    end
    
    %Plot Vessel and Coilset
    plot(Vessel);
    plot(Coilset);
    
    %Colourmap
    colormap(plasma);
    cbar = colorbar;
    cbar.Label.String = CbarLabel;
    
    %Title, Legend, Labels, etc...
    title(gca,Title);
    legend(gca,'hide');
    set(gca,'XLim',[0 1.1]);
    set(gca,'YLim',[-1.1 1.1]);
    set(gca, 'FontSize', 13, 'LineWidth', 0.75);
    xlabel(gca,'R (m)');
    ylabel(gca,'Z (m)');
    
    %Data output
    saveas(gcf, SaveString);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cm_data=plasma(m)
cm = [[  5.03832136e-02,   2.98028976e-02,   5.27974883e-01],
       [  6.35363639e-02,   2.84259729e-02,   5.33123681e-01],
       [  7.53531234e-02,   2.72063728e-02,   5.38007001e-01],
       [  8.62217979e-02,   2.61253206e-02,   5.42657691e-01],
       [  9.63786097e-02,   2.51650976e-02,   5.47103487e-01],
       [  1.05979704e-01,   2.43092436e-02,   5.51367851e-01],
       [  1.15123641e-01,   2.35562500e-02,   5.55467728e-01],
       [  1.23902903e-01,   2.28781011e-02,   5.59423480e-01],
       [  1.32380720e-01,   2.22583774e-02,   5.63250116e-01],
       [  1.40603076e-01,   2.16866674e-02,   5.66959485e-01],
       [  1.48606527e-01,   2.11535876e-02,   5.70561711e-01],
       [  1.56420649e-01,   2.06507174e-02,   5.74065446e-01],
       [  1.64069722e-01,   2.01705326e-02,   5.77478074e-01],
       [  1.71573925e-01,   1.97063415e-02,   5.80805890e-01],
       [  1.78950212e-01,   1.92522243e-02,   5.84054243e-01],
       [  1.86212958e-01,   1.88029767e-02,   5.87227661e-01],
       [  1.93374449e-01,   1.83540593e-02,   5.90329954e-01],
       [  2.00445260e-01,   1.79015512e-02,   5.93364304e-01],
       [  2.07434551e-01,   1.74421086e-02,   5.96333341e-01],
       [  2.14350298e-01,   1.69729276e-02,   5.99239207e-01],
       [  2.21196750e-01,   1.64970484e-02,   6.02083323e-01],
       [  2.27982971e-01,   1.60071509e-02,   6.04867403e-01],
       [  2.34714537e-01,   1.55015065e-02,   6.07592438e-01],
       [  2.41396253e-01,   1.49791041e-02,   6.10259089e-01],
       [  2.48032377e-01,   1.44393586e-02,   6.12867743e-01],
       [  2.54626690e-01,   1.38820918e-02,   6.15418537e-01],
       [  2.61182562e-01,   1.33075156e-02,   6.17911385e-01],
       [  2.67702993e-01,   1.27162163e-02,   6.20345997e-01],
       [  2.74190665e-01,   1.21091423e-02,   6.22721903e-01],
       [  2.80647969e-01,   1.14875915e-02,   6.25038468e-01],
       [  2.87076059e-01,   1.08554862e-02,   6.27294975e-01],
       [  2.93477695e-01,   1.02128849e-02,   6.29490490e-01],
       [  2.99855122e-01,   9.56079551e-03,   6.31623923e-01],
       [  3.06209825e-01,   8.90185346e-03,   6.33694102e-01],
       [  3.12543124e-01,   8.23900704e-03,   6.35699759e-01],
       [  3.18856183e-01,   7.57551051e-03,   6.37639537e-01],
       [  3.25150025e-01,   6.91491734e-03,   6.39512001e-01],
       [  3.31425547e-01,   6.26107379e-03,   6.41315649e-01],
       [  3.37683446e-01,   5.61830889e-03,   6.43048936e-01],
       [  3.43924591e-01,   4.99053080e-03,   6.44710195e-01],
       [  3.50149699e-01,   4.38202557e-03,   6.46297711e-01],
       [  3.56359209e-01,   3.79781761e-03,   6.47809772e-01],
       [  3.62553473e-01,   3.24319591e-03,   6.49244641e-01],
       [  3.68732762e-01,   2.72370721e-03,   6.50600561e-01],
       [  3.74897270e-01,   2.24514897e-03,   6.51875762e-01],
       [  3.81047116e-01,   1.81356205e-03,   6.53068467e-01],
       [  3.87182639e-01,   1.43446923e-03,   6.54176761e-01],
       [  3.93304010e-01,   1.11388259e-03,   6.55198755e-01],
       [  3.99410821e-01,   8.59420809e-04,   6.56132835e-01],
       [  4.05502914e-01,   6.78091517e-04,   6.56977276e-01],
       [  4.11580082e-01,   5.77101735e-04,   6.57730380e-01],
       [  4.17642063e-01,   5.63847476e-04,   6.58390492e-01],
       [  4.23688549e-01,   6.45902780e-04,   6.58956004e-01],
       [  4.29719186e-01,   8.31008207e-04,   6.59425363e-01],
       [  4.35733575e-01,   1.12705875e-03,   6.59797077e-01],
       [  4.41732123e-01,   1.53984779e-03,   6.60069009e-01],
       [  4.47713600e-01,   2.07954744e-03,   6.60240367e-01],
       [  4.53677394e-01,   2.75470302e-03,   6.60309966e-01],
       [  4.59622938e-01,   3.57374415e-03,   6.60276655e-01],
       [  4.65549631e-01,   4.54518084e-03,   6.60139383e-01],
       [  4.71456847e-01,   5.67758762e-03,   6.59897210e-01],
       [  4.77343929e-01,   6.97958743e-03,   6.59549311e-01],
       [  4.83210198e-01,   8.45983494e-03,   6.59094989e-01],
       [  4.89054951e-01,   1.01269996e-02,   6.58533677e-01],
       [  4.94877466e-01,   1.19897486e-02,   6.57864946e-01],
       [  5.00677687e-01,   1.40550640e-02,   6.57087561e-01],
       [  5.06454143e-01,   1.63333443e-02,   6.56202294e-01],
       [  5.12206035e-01,   1.88332232e-02,   6.55209222e-01],
       [  5.17932580e-01,   2.15631918e-02,   6.54108545e-01],
       [  5.23632990e-01,   2.45316468e-02,   6.52900629e-01],
       [  5.29306474e-01,   2.77468735e-02,   6.51586010e-01],
       [  5.34952244e-01,   3.12170300e-02,   6.50165396e-01],
       [  5.40569510e-01,   3.49501310e-02,   6.48639668e-01],
       [  5.46157494e-01,   3.89540334e-02,   6.47009884e-01],
       [  5.51715423e-01,   4.31364795e-02,   6.45277275e-01],
       [  5.57242538e-01,   4.73307585e-02,   6.43443250e-01],
       [  5.62738096e-01,   5.15448092e-02,   6.41509389e-01],
       [  5.68201372e-01,   5.57776706e-02,   6.39477440e-01],
       [  5.73631859e-01,   6.00281369e-02,   6.37348841e-01],
       [  5.79028682e-01,   6.42955547e-02,   6.35126108e-01],
       [  5.84391137e-01,   6.85790261e-02,   6.32811608e-01],
       [  5.89718606e-01,   7.28775875e-02,   6.30407727e-01],
       [  5.95010505e-01,   7.71902878e-02,   6.27916992e-01],
       [  6.00266283e-01,   8.15161895e-02,   6.25342058e-01],
       [  6.05485428e-01,   8.58543713e-02,   6.22685703e-01],
       [  6.10667469e-01,   9.02039303e-02,   6.19950811e-01],
       [  6.15811974e-01,   9.45639838e-02,   6.17140367e-01],
       [  6.20918555e-01,   9.89336721e-02,   6.14257440e-01],
       [  6.25986869e-01,   1.03312160e-01,   6.11305174e-01],
       [  6.31016615e-01,   1.07698641e-01,   6.08286774e-01],
       [  6.36007543e-01,   1.12092335e-01,   6.05205491e-01],
       [  6.40959444e-01,   1.16492495e-01,   6.02064611e-01],
       [  6.45872158e-01,   1.20898405e-01,   5.98867442e-01],
       [  6.50745571e-01,   1.25309384e-01,   5.95617300e-01],
       [  6.55579615e-01,   1.29724785e-01,   5.92317494e-01],
       [  6.60374266e-01,   1.34143997e-01,   5.88971318e-01],
       [  6.65129493e-01,   1.38566428e-01,   5.85582301e-01],
       [  6.69845385e-01,   1.42991540e-01,   5.82153572e-01],
       [  6.74522060e-01,   1.47418835e-01,   5.78688247e-01],
       [  6.79159664e-01,   1.51847851e-01,   5.75189431e-01],
       [  6.83758384e-01,   1.56278163e-01,   5.71660158e-01],
       [  6.88318440e-01,   1.60709387e-01,   5.68103380e-01],
       [  6.92840088e-01,   1.65141174e-01,   5.64521958e-01],
       [  6.97323615e-01,   1.69573215e-01,   5.60918659e-01],
       [  7.01769334e-01,   1.74005236e-01,   5.57296144e-01],
       [  7.06177590e-01,   1.78437000e-01,   5.53656970e-01],
       [  7.10548747e-01,   1.82868306e-01,   5.50003579e-01],
       [  7.14883195e-01,   1.87298986e-01,   5.46338299e-01],
       [  7.19181339e-01,   1.91728906e-01,   5.42663338e-01],
       [  7.23443604e-01,   1.96157962e-01,   5.38980786e-01],
       [  7.27670428e-01,   2.00586086e-01,   5.35292612e-01],
       [  7.31862231e-01,   2.05013174e-01,   5.31600995e-01],
       [  7.36019424e-01,   2.09439071e-01,   5.27908434e-01],
       [  7.40142557e-01,   2.13863965e-01,   5.24215533e-01],
       [  7.44232102e-01,   2.18287899e-01,   5.20523766e-01],
       [  7.48288533e-01,   2.22710942e-01,   5.16834495e-01],
       [  7.52312321e-01,   2.27133187e-01,   5.13148963e-01],
       [  7.56303937e-01,   2.31554749e-01,   5.09468305e-01],
       [  7.60263849e-01,   2.35975765e-01,   5.05793543e-01],
       [  7.64192516e-01,   2.40396394e-01,   5.02125599e-01],
       [  7.68090391e-01,   2.44816813e-01,   4.98465290e-01],
       [  7.71957916e-01,   2.49237220e-01,   4.94813338e-01],
       [  7.75795522e-01,   2.53657797e-01,   4.91170517e-01],
       [  7.79603614e-01,   2.58078397e-01,   4.87539124e-01],
       [  7.83382636e-01,   2.62499662e-01,   4.83917732e-01],
       [  7.87132978e-01,   2.66921859e-01,   4.80306702e-01],
       [  7.90855015e-01,   2.71345267e-01,   4.76706319e-01],
       [  7.94549101e-01,   2.75770179e-01,   4.73116798e-01],
       [  7.98215577e-01,   2.80196901e-01,   4.69538286e-01],
       [  8.01854758e-01,   2.84625750e-01,   4.65970871e-01],
       [  8.05466945e-01,   2.89057057e-01,   4.62414580e-01],
       [  8.09052419e-01,   2.93491117e-01,   4.58869577e-01],
       [  8.12611506e-01,   2.97927865e-01,   4.55337565e-01],
       [  8.16144382e-01,   3.02368130e-01,   4.51816385e-01],
       [  8.19651255e-01,   3.06812282e-01,   4.48305861e-01],
       [  8.23132309e-01,   3.11260703e-01,   4.44805781e-01],
       [  8.26587706e-01,   3.15713782e-01,   4.41315901e-01],
       [  8.30017584e-01,   3.20171913e-01,   4.37835947e-01],
       [  8.33422053e-01,   3.24635499e-01,   4.34365616e-01],
       [  8.36801237e-01,   3.29104836e-01,   4.30905052e-01],
       [  8.40155276e-01,   3.33580106e-01,   4.27454836e-01],
       [  8.43484103e-01,   3.38062109e-01,   4.24013059e-01],
       [  8.46787726e-01,   3.42551272e-01,   4.20579333e-01],
       [  8.50066132e-01,   3.47048028e-01,   4.17153264e-01],
       [  8.53319279e-01,   3.51552815e-01,   4.13734445e-01],
       [  8.56547103e-01,   3.56066072e-01,   4.10322469e-01],
       [  8.59749520e-01,   3.60588229e-01,   4.06916975e-01],
       [  8.62926559e-01,   3.65119408e-01,   4.03518809e-01],
       [  8.66077920e-01,   3.69660446e-01,   4.00126027e-01],
       [  8.69203436e-01,   3.74211795e-01,   3.96738211e-01],
       [  8.72302917e-01,   3.78773910e-01,   3.93354947e-01],
       [  8.75376149e-01,   3.83347243e-01,   3.89975832e-01],
       [  8.78422895e-01,   3.87932249e-01,   3.86600468e-01],
       [  8.81442916e-01,   3.92529339e-01,   3.83228622e-01],
       [  8.84435982e-01,   3.97138877e-01,   3.79860246e-01],
       [  8.87401682e-01,   4.01761511e-01,   3.76494232e-01],
       [  8.90339687e-01,   4.06397694e-01,   3.73130228e-01],
       [  8.93249647e-01,   4.11047871e-01,   3.69767893e-01],
       [  8.96131191e-01,   4.15712489e-01,   3.66406907e-01],
       [  8.98983931e-01,   4.20391986e-01,   3.63046965e-01],
       [  9.01807455e-01,   4.25086807e-01,   3.59687758e-01],
       [  9.04601295e-01,   4.29797442e-01,   3.56328796e-01],
       [  9.07364995e-01,   4.34524335e-01,   3.52969777e-01],
       [  9.10098088e-01,   4.39267908e-01,   3.49610469e-01],
       [  9.12800095e-01,   4.44028574e-01,   3.46250656e-01],
       [  9.15470518e-01,   4.48806744e-01,   3.42890148e-01],
       [  9.18108848e-01,   4.53602818e-01,   3.39528771e-01],
       [  9.20714383e-01,   4.58417420e-01,   3.36165582e-01],
       [  9.23286660e-01,   4.63250828e-01,   3.32800827e-01],
       [  9.25825146e-01,   4.68103387e-01,   3.29434512e-01],
       [  9.28329275e-01,   4.72975465e-01,   3.26066550e-01],
       [  9.30798469e-01,   4.77867420e-01,   3.22696876e-01],
       [  9.33232140e-01,   4.82779603e-01,   3.19325444e-01],
       [  9.35629684e-01,   4.87712357e-01,   3.15952211e-01],
       [  9.37990034e-01,   4.92666544e-01,   3.12575440e-01],
       [  9.40312939e-01,   4.97642038e-01,   3.09196628e-01],
       [  9.42597771e-01,   5.02639147e-01,   3.05815824e-01],
       [  9.44843893e-01,   5.07658169e-01,   3.02433101e-01],
       [  9.47050662e-01,   5.12699390e-01,   2.99048555e-01],
       [  9.49217427e-01,   5.17763087e-01,   2.95662308e-01],
       [  9.51343530e-01,   5.22849522e-01,   2.92274506e-01],
       [  9.53427725e-01,   5.27959550e-01,   2.88883445e-01],
       [  9.55469640e-01,   5.33093083e-01,   2.85490391e-01],
       [  9.57468770e-01,   5.38250172e-01,   2.82096149e-01],
       [  9.59424430e-01,   5.43431038e-01,   2.78700990e-01],
       [  9.61335930e-01,   5.48635890e-01,   2.75305214e-01],
       [  9.63202573e-01,   5.53864931e-01,   2.71909159e-01],
       [  9.65023656e-01,   5.59118349e-01,   2.68513200e-01],
       [  9.66798470e-01,   5.64396327e-01,   2.65117752e-01],
       [  9.68525639e-01,   5.69699633e-01,   2.61721488e-01],
       [  9.70204593e-01,   5.75028270e-01,   2.58325424e-01],
       [  9.71835007e-01,   5.80382015e-01,   2.54931256e-01],
       [  9.73416145e-01,   5.85761012e-01,   2.51539615e-01],
       [  9.74947262e-01,   5.91165394e-01,   2.48151200e-01],
       [  9.76427606e-01,   5.96595287e-01,   2.44766775e-01],
       [  9.77856416e-01,   6.02050811e-01,   2.41387186e-01],
       [  9.79232922e-01,   6.07532077e-01,   2.38013359e-01],
       [  9.80556344e-01,   6.13039190e-01,   2.34646316e-01],
       [  9.81825890e-01,   6.18572250e-01,   2.31287178e-01],
       [  9.83040742e-01,   6.24131362e-01,   2.27937141e-01],
       [  9.84198924e-01,   6.29717516e-01,   2.24595006e-01],
       [  9.85300760e-01,   6.35329876e-01,   2.21264889e-01],
       [  9.86345421e-01,   6.40968508e-01,   2.17948456e-01],
       [  9.87332067e-01,   6.46633475e-01,   2.14647532e-01],
       [  9.88259846e-01,   6.52324832e-01,   2.11364122e-01],
       [  9.89127893e-01,   6.58042630e-01,   2.08100426e-01],
       [  9.89935328e-01,   6.63786914e-01,   2.04858855e-01],
       [  9.90681261e-01,   6.69557720e-01,   2.01642049e-01],
       [  9.91364787e-01,   6.75355082e-01,   1.98452900e-01],
       [  9.91984990e-01,   6.81179025e-01,   1.95294567e-01],
       [  9.92540939e-01,   6.87029567e-01,   1.92170500e-01],
       [  9.93031693e-01,   6.92906719e-01,   1.89084459e-01],
       [  9.93456302e-01,   6.98810484e-01,   1.86040537e-01],
       [  9.93813802e-01,   7.04740854e-01,   1.83043180e-01],
       [  9.94103226e-01,   7.10697814e-01,   1.80097207e-01],
       [  9.94323596e-01,   7.16681336e-01,   1.77207826e-01],
       [  9.94473934e-01,   7.22691379e-01,   1.74380656e-01],
       [  9.94553260e-01,   7.28727890e-01,   1.71621733e-01],
       [  9.94560594e-01,   7.34790799e-01,   1.68937522e-01],
       [  9.94494964e-01,   7.40880020e-01,   1.66334918e-01],
       [  9.94355411e-01,   7.46995448e-01,   1.63821243e-01],
       [  9.94140989e-01,   7.53136955e-01,   1.61404226e-01],
       [  9.93850778e-01,   7.59304390e-01,   1.59091984e-01],
       [  9.93482190e-01,   7.65498551e-01,   1.56890625e-01],
       [  9.93033251e-01,   7.71719833e-01,   1.54807583e-01],
       [  9.92505214e-01,   7.77966775e-01,   1.52854862e-01],
       [  9.91897270e-01,   7.84239120e-01,   1.51041581e-01],
       [  9.91208680e-01,   7.90536569e-01,   1.49376885e-01],
       [  9.90438793e-01,   7.96858775e-01,   1.47869810e-01],
       [  9.89587065e-01,   8.03205337e-01,   1.46529128e-01],
       [  9.88647741e-01,   8.09578605e-01,   1.45357284e-01],
       [  9.87620557e-01,   8.15977942e-01,   1.44362644e-01],
       [  9.86509366e-01,   8.22400620e-01,   1.43556679e-01],
       [  9.85314198e-01,   8.28845980e-01,   1.42945116e-01],
       [  9.84031139e-01,   8.35315360e-01,   1.42528388e-01],
       [  9.82652820e-01,   8.41811730e-01,   1.42302653e-01],
       [  9.81190389e-01,   8.48328902e-01,   1.42278607e-01],
       [  9.79643637e-01,   8.54866468e-01,   1.42453425e-01],
       [  9.77994918e-01,   8.61432314e-01,   1.42808191e-01],
       [  9.76264977e-01,   8.68015998e-01,   1.43350944e-01],
       [  9.74443038e-01,   8.74622194e-01,   1.44061156e-01],
       [  9.72530009e-01,   8.81250063e-01,   1.44922913e-01],
       [  9.70532932e-01,   8.87896125e-01,   1.45918663e-01],
       [  9.68443477e-01,   8.94563989e-01,   1.47014438e-01],
       [  9.66271225e-01,   9.01249365e-01,   1.48179639e-01],
       [  9.64021057e-01,   9.07950379e-01,   1.49370428e-01],
       [  9.61681481e-01,   9.14672479e-01,   1.50520343e-01],
       [  9.59275646e-01,   9.21406537e-01,   1.51566019e-01],
       [  9.56808068e-01,   9.28152065e-01,   1.52409489e-01],
       [  9.54286813e-01,   9.34907730e-01,   1.52921158e-01],
       [  9.51726083e-01,   9.41670605e-01,   1.52925363e-01],
       [  9.49150533e-01,   9.48434900e-01,   1.52177604e-01],
       [  9.46602270e-01,   9.55189860e-01,   1.50327944e-01],
       [  9.44151742e-01,   9.61916487e-01,   1.46860789e-01],
       [  9.41896120e-01,   9.68589814e-01,   1.40955606e-01],
       [  9.40015097e-01,   9.75158357e-01,   1.31325517e-01]];
   
if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    hsv(153:end,1)=hsv(153:end,1)+1; % hardcoded
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data(cm_data(:,1)>1,1)=cm_data(cm_data(:,1)>1,1)-1;
    cm_data=hsv2rgb(cm_data);
  
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
