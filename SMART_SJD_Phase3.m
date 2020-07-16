%% SMall Aspect Ratio Tokamak (SMART)
%#################################
%#		Point of Contact		 #
%#								 #
%#	   Dr. Scott J. Doyle		 #
%#	  Scott.Doyle@Physics.org	 #
%#	  University of Seville		 #
%#  Plasma Tech & Fusion Science #
%#  National Acceleration Centre #
%#	     Seville, Spain		     #
%#								 #
%#################################
%%

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
ProjectName = 'S3-000002';		%Define global project name
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
nSol = 230                                  % Number of Axial Solenoid Windings
RSolInner = 0.115; RSolOuter = 0.145;       % Inner and Outer solenoid radii    [m]
Width_Sol = 0.011; Height_Sol = 0.011;      % Width and height of Sol Winding   [m]
RSolCentre = (RSolInner+RSolOuter)/2;       % Geometric Centre of Sol (0.13)    [m]
RSolCentreWinding = RSolOuter-Width_Sol;    % Winding Centre of Sol (0.1340)    [m] 
ZMinSol = -0.7750; %ZMinCentre;             % Solenoid Min Z position           [m]
ZMaxSol = +0.7750; %ZMaxCentre;             % Solenoid Max Z position           [m]

%Number of Radial (R) and axial (Z) PF coil windings
nZDiv1=6; nRDiv1=6;     % Square
nZDiv2=6; nRDiv2=4;     % Axial Rectangle
nZPF1=6; nRPF1=4;       % Axial Rectangle
nZPF2=6; nRPF2=4;       % Axial Rectangle

%Calculate total number of windings in each coil
nDiv1=nZDiv1*nRDiv1;    % 36 Total Windings
nDiv2=nZDiv2*nRDiv2;    % 24 Total Windings
nPF1=nZPF1*nRPF1;       % 24 Total Windings
nPF2=nZPF2*nRPF2;       % 24 Total Windings

%Define coil total cross-sectional dimensions
width_PF1 = 0.075; height_PF1 = 0.050;    %[m]
width_PF2 = 0.075; height_PF2 = 0.050;    %[m]
width_Div1 = 0.075; height_Div1 = 0.075;  %[m]
width_Div2 = 0.075; height_Div2 = 0.050;  %[m]

%Define central location of coil sets                                          %NOTES  (Outer midplane flange diameter = 176.8mm (180 mm))
R_PF1 = 0.940;  %R position of PF1 (m)	%0.940m     (MINIMUM OF 938mm)
Z_PF1 = 0.200;  %Z Position of PF1 (m)	%0.200m     (MINIMUM OF 200mm)         %Closer together is optimal for +d and -d   (35cm seperation)
R_PF2 = 0.700;  %R Position of PF2 (m)	%0.700m     (MINIMUM OF 938mm)         %Closer to wall is optimal for +d and -d
Z_PF2 = 0.575;  %Z Position of PF2 (m)	%0.575m     (MINIMUM OF 608mm)         %Lower is better for -d, but reduces volume
R_Div1 = 0.250; %R Position of Div1 (m)	%0.250m     (MINIMUM OF 236mm)
Z_Div1 = 0.700; %Z Position of Div1 (m)	%0.700m     (MINIMUM OF 890mm)         %Higher increases aspect ratio (at expense of triang)
R_Div2 = 0.500; %R Position of Div2 (m)	%0.500m     (MINIMUM OF 458mm)
Z_Div2 = 0.700; %Z Position of Div2 (m)	%0.700m     (MINIMUM OF 890mm)         %Lower is better for +d, but reduces volume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      DEFINE OPERATIONAL PARAMETERS                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%  DEFINE INITIAL PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%

%Define any required constants
global mu0; mu0 = 1.2566e-06;    % Magnetic Moment              [I/m^2]
global BEarth; BEarth = 1.0E-4;  % Earth's B-Field (Def 5e-5)	[T]       

%Define initial operating conditions (primarily used for Topeol2)
Te = 250;			% Electron Temperature [eV]
Ti = Te*0.1;		% Ion Temperature      [eV]
BT = 1.0;			% Toroidal B-Field     [T] (Defined at Rgeo)
Ip = 300e3;			% Plasma current       [A]
RGeo = 0.420;		% Geometrical Radius   [m] (~0.420)
ZGeo = 0.000;		% Geometrical Axis     [m] (~0.000)
RSep = 0.700;		% Separatrix Radius    [m] (~0.700)
rGeo = RSep-RGeo;	% Minor Radius         [m] (~0.250)
Aspect = RGeo/rGeo;	% Aspect ratio         [-] (~1.850)
Kappa = 1.80;		% Elongation           [-] (~1.800)
delta = 0.20;		% Triangularity        [-] (~0.200)
li2 = 1;			% Inductance	       [-]

%Compute further operating conditions (primarily used for Topeol2)
Gr_Limit = 1e20*(Ip*1e-6/(pi*Kappa*rGeo^2));     % Greenwald Limit    [m-3]
Gr_Frac = 0.70;                            % Greenwald Fraction       [-]
ne = Gr_Limit*Gr_Frac;                     % Electron Density         [m-3]
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
RGeo_efit = 0.420;					% Geometrical Radius	[m] (Default 0.420) ::
ZGeo_efit = 0.000;					% Geometrical Axis		[m] (Default 0.000) ::
Aspect_efit = 1.85;                 % Aspect Ratio          [-] (Default 1.850) :: RGeo/rGeo
rGeo_efit = RGeo_efit/Aspect_efit;  % Minor Radius	        [m] (Default 0.238) :: RGeo/Aspect
Kappa_efit = 1.80;					% Elongation			[-] (Default 1.800) ::
delta_efit = 0.00;					% Triangularity			[-] (-1.00-> +0.00 -> +1.00) ::
efitGeometry_Init = [RGeo_efit, ZGeo_efit, rGeo_efit, Kappa_efit, delta_efit];

%Define feedback stability perturbations
deltaRGeo = 0.00;	% Small radial perturbation         [m]
deltaZGeo = 0.01;	% Small axial perturbation          [m]
deltaAspect = 0.00;	% Small aspect ratio perturbation   [-]
deltaKappa = 0.00;	% Small elongation perturbation     [-]
deltadelta = 0.00;	% Small triangiularity perturbation [-]
PertGeometry_Init = [deltaRGeo,deltaZGeo,deltaAspect,deltaKappa,deltadelta];

%Define Coil density, temperature and resistivity
coil_density = 1;                       % Relative Coil Density      [Arb]
coil_temp = 293.0;                      % Initial Coil Temperature   [K]
resistivity = copper_resistivity_at_temperature(coil_temp);

%Gas species analouge - H=1, He=2, Ar=11.85 (for Te < 280eV) https://www.webelements.com/argon/atoms.html
%H discharge, Z_eff increased to 2 to allow for impurities in the plasma (Carbon wall tiles)
Z_eff = 2.0;                            % Effective Nuclear Charge      %[e-]
%Calculate perpendicular and parallel plasma resistivity using Spitzer model
Lambda=(12*pi*((8.854E-12*1.6E-19*Te)^3/(ne*(1.6E-19)^6))^(1/2));
PlasmaResistPerp=(0.74*1.65E-9*Z_eff*log(Lambda))/((Te*1E-3)^(3/2));    %[Ohms]
PlasmaResistPara=PlasmaResistPerp/1.96;                                 %[Ohms]

%Define null field region (sensor_btheta) radius
R_Null = 0.15;                      	% Null field region radius      %[m]


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
                                    %TauR1=30ms %TauR=30ms      %TauR=30ms
                                    %RGeo=0.42  %RGeo=0.xx      %RGeo=0.xx
%Solenoid coil currents [kA]		%Phase3     %Phase3NegTri   %Phase3PosTri
I_Sol_Null=+10000;					%+10000;    %+xx000;        %+xx000;
I_Sol_MidRamp=+0000;				%+00000;    %+00000;        %+00000;
I_Sol_Equil=-1000;                  %-01000;    %-xx000;        %-xx000;
I_Sol_EndEquil=-7000;           	%-07000;    %-xx000;        %-xx000;

%PF coil currents (At Equilibrium, time(4,5,6))
I_PF1_Equil=-4000;					%-4000;     %-xx00;         %-xx00;
I_PF2_Equil=-4000;					%-4000;     %-xx00;         %-xx00;     (NEG FOR +delta, POS FOR -delta, after efit) 
I_Div1_Equil=+4500;					%+4500;     %-xx00;         %+xx00;     (HIGH FOR +delta, LOW FOR -delta, before efit)
I_Div2_Equil=+0000;                 %+0000;     %+0000;         %+0000;

%Define number of time-steps (vertices) in the current waveforms
TauN  = 0.030;			% Null-Field Timescale      [s] Determines null-field decay timescale
TauR1 = 0.030;			% Breakdown Ramp Timescale  [s] Determines max loop voltage
TauR2 = 0.100;			% PF & Div Ramp Timescale   [s] Determines max PF/Div current ramp
TauR  = TauR1+TauR2;    % Total Ramp Timescale      [s] 
TauP  = 0.500;			% Pulse Timescale      		[s] Determines flat-top timescale
%Time   [Init      PrePulse   InitRampDown  MidRampDown  EndRampDown  MidEquil     Terminate         ];
time =  [-4*TauN   -TauN      0.0           TauR1        TauR         TauR+TauP    TauR+TauP+(4*TauN)];
nTime = length(time);	% Coil Waveform Timesteps	[-]

%Fit any dynamic coil currents, set with 'linear', {pre-ramp, mid-ramp, end-ramp}
%I_Sol_MidRamp = FitSolenoidRamp({I_Sol_Null,I_Sol_MidRamp,I_Sol_Equil},time);

%Construct Sol, PF/Div coil current waveforms vertices
%					              %!Null-Field! %!Breakdown!   %!Efit Icoil!
%Time   	     [1,  2,          3,            4,             5,             6,             7];
ISol_Waveform =  [0,  I_Sol_Null, I_Sol_Null,   I_Sol_MidRamp, I_Sol_Equil,   I_Sol_EndEquil,0];
IPF1_Waveform =  [0,  NaN,        NaN,          NaN,           I_PF1_Equil,   I_PF1_Equil,   0];
IPF2_Waveform =  [0,  NaN,        NaN,          NaN,           I_PF2_Equil,   I_PF2_Equil,   0];
IDiv1_Waveform = [0,  NaN,        NaN,          NaN,           I_Div1_Equil,  I_Div1_Equil,  0];
IDiv2_Waveform = [0,  NaN,        NaN,          NaN,           I_Div2_Equil,  I_Div2_Equil,  0];
%%%%%
CoilWaveforms = [ISol_Waveform; IPF1_Waveform; IPF2_Waveform; IDiv1_Waveform; IDiv2_Waveform];

%Define dynamic coils (i.e. which coil currents are fit by efit)
global efitCoils; efitCoils = {'PF1','PF2'};                        % Default PF1, PF2
global feedbackCoils; feedbackCoils = {'Div2'};                     % Default Div2


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
disp([ 'I_Sol_MidRamp = ' num2str(I_Sol_MidRamp/1000) ' [kA] ']);
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
GridSize_R = [0.01, 1.1];	%[m]        (0.01, 1.1)
GridSize_Z = [-1.3, 1.3];	%[m]        (-1.3, 1.3)
GridCells_R = 300;          %[Cells]    (300)   (>300 would be good, causes RAM issues)
GridCells_Z = 251;          %[Cells]    (251)   (~251 seems good, higher becomes unstable)
global GridRes_R; global GridRes_Z;
GridRes_R = max(GridSize_R)/GridCells_R;	%[m/Cell] GridRes_R = 0.00366
GridRes_Z = max(GridSize_Z)/GridCells_Z;	%[m/Cell] GridRes_Z = 0.01036

%Generate fiesta_grid object over which equilibrum simulation will be performed
global Grid;
Grid = fiesta_grid(GridSize_R(1),GridSize_R(2),GridCells_R, GridSize_Z(1),GridSize_Z(2),GridCells_Z);

%Extract vectors of R and Z grid points for use in further diagnostics
rGrid=get(Grid,'r'); %1*300
zGrid=get(Grid,'z'); %1*251


%%%%%%%%%%%%%%%%%%  INITIATE VACUUM VESSEL FILAMENTS  %%%%%%%%%%%%%%%%%%%

%Define vessel corners, thickness and filament cross-sectional area 
VesselDimensions = [RMinCentre, RMaxCentre, ZMinCentre, ZMaxCentre];       %[m]
WallThickness = [VWall_Upper, VWall_Outboard, VWall_Lower, VWall_Inboard]; %[m]
%Lower filament areas give higher passive current resolution (Note :: FilamentArea affects convergence)
FilamentArea = 1.50e-4; %(2.5e-4 > A > 1.5e-4 or RZIp M,R matrices fail)   %[m^2]

%Construct SMART vessel wall filaments ("Static"=fixed fil area, "Diff"=scaled fil area)
[vessel_filament,R_Fil_Array,Z_Fil_Array] = ... 
    CreateRectilinearVessel(VesselDimensions,WallThickness,FilamentArea,"Diff");

%Construct passive vessel components and arrange into a vessel object
%Resistivity and density are set for stainless steel
VesselResistivity = 6.9e-7;  VesselDensity = 7.8e3;	 %[Ohm]; [Kg/m3]; Stainless Steel (GRADE)
global passive; passive = fiesta_passive('STVesselPas',vessel_filament,'g',VesselResistivity,VesselDensity);
global vessel; vessel = fiesta_vessel( 'STVessel',passive);

%Compute characteristic magnetic field penetration timescale (Amoskov2005)
TauVessel = (mu0*max(WallThickness)^2)/VesselResistivity;       %[s]

%%%%%%%%%%%%%%%%%%%%%%  INITIATE SOL & PF COILS  %%%%%%%%%%%%%%%%%%%%%%%%%%

%Define and initiate PF coils - Arbitrary numbering of coils
global iSol; iSol = 1;       %Central Inducting Solenoid
global iPF1; iPF1 = 2;       %Upper Plasma Forming Coil
global iPF2; iPF2 = 3;       %Lower Plasma Forming Coil
global iDiv1; iDiv1 = 4;     %Inboard Divertor Coil
global iDiv2; iDiv2 = 5;     %Outboard Divertor Coil

%Create array containing number of coil windings - Used to generate coil objects
global coilturns; coilturns=[];
coilturns(iPF1) = nPF1; coilturns(iPF2) = nPF2;
coilturns(iDiv1) = nDiv1; coilturns(iDiv2) = nDiv2;
coilturns(iSol) = nSol; nSolR = 1;
nPF = 5;

%Create coil set from parameters defined above. (Function made by Carlos Soria)
%Function createVESTPFCircuit creates two PF coils. One in (R, Z) and another in (R, -Z)
Sol = CreateSMARTSolenoidCircuit('Sol',RSolOuter,RSolInner,ZMaxSol,ZMinSol,coilturns(iSol),nSolR,coil_temp,resistivity,coil_density);
PF1  = CreateSMARTCoilCircuit('PF1',R_PF1,Z_PF1,width_PF1,height_PF1,coilturns(iPF1),nZPF1,nRPF1,true,coil_temp,resistivity,coil_density);
PF2  = CreateSMARTCoilCircuit('PF2',R_PF2,Z_PF2,width_PF2,height_PF2,coilturns(iPF2),nZPF2,nRPF2,true,coil_temp,resistivity,coil_density);
Div1 = CreateSMARTCoilCircuit('Div1',R_Div1,Z_Div1,width_Div1,height_Div1,coilturns(iDiv1),nZDiv1,nRDiv1,true,coil_temp,resistivity,coil_density); 
Div2 = CreateSMARTCoilCircuit('Div2',R_Div2,Z_Div2,width_Div2,height_Div2,coilturns(iDiv2),nZDiv2,nRDiv2,true,coil_temp,resistivity,coil_density);

%Collate global coilset containing Solenoid, PF and Div coil circuits (expects a row aligned filament array)
R_Fil_Array = transpose(R_Fil_Array); Z_Fil_Array = transpose(Z_Fil_Array);     
global coilset; coilset = fiesta_coilset('SMARTcoilset',[Sol,PF1,PF2,Div1,Div2],false,R_Fil_Array',Z_Fil_Array');
coilset_init = coilset;

%Reduce solenoid current by number of windings
CoilWaveforms(1,:) = CoilWaveforms(1,:)/nSolR;
% ISSUE :: Need to reduce solenoid current by number of radial filaments
% ISSUE :: This doesn't appear to be required for the PF coils!?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  END SIMULATION INITIAL SET-UP  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TO DO %%%

%%%   GET ALL VESSEL, COIL AND CURRENT WAVEFORM DATA UPDATED ON SMART REPO
%%%   GET ALL FIGURES INTO FUNCTIONS - MAKE VESSEL/COIL SUB-FUNCTION
%%%   GET I/O ALL INTO FUNCTIONS AND GET NEW SAVING ROUTINES SORTED OUT
%%%   MIGRATE THE CreateSMARTCoilCircuit FUCTION TO THIS VERSION OF CODE
%%%   FINISH COMMENTS ON THE NEW CreateSMARTSolenoidCircuit FUNCTION
%%%   FINISH COMMENTS ON ALL NEW FUNCTIONS

%%%   ARCHIVE ALL NEW BIBLO AND ADD TO MENDELEY ASAP

%%%   ENABLE TOGGLEABLE UPPER AND LOWER SINGLE NULL CONFIGURATION 

%%%   GET THE FEEDBACK SYSTEM WORKING FOR VERTICAL AND HORIZONTAL STABILITY

%%%   FIX THE CORNER OF THE DIFF VESSEL WALLS (LAST FILAMENT IS LARGER)

%%%   GET RZIP ABLE TO TAKE BOTH COIL AND VESSEL FILAMENTS !!!!
%%%   findboundary.m function contains rules for LCFS boundary









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           COMPUTE INITIAL TARGET AND NULL-FIELD EQUILIBRIUA             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TimeArray index for target equilibrium 
TimeIndex_Discharge = 5;                          %default time(5)  End of Sol Ramp
TimeIndex_NullField = 3;                          %default time(3)  Prior to Sol Ramp

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
global config; [Equil,EquilParams,CoilWaveforms,efitGeometry,config] = ...
    efit(jprofile,Irod,'nullconfig',efitGeometry_Init,CoilWaveforms_Init,[],TimeIndex_Discharge);

%Save discharge coil currents for all coils at TimeIndex_Discharge
CoilCurrentsEfit = transpose(CoilWaveforms(:,TimeIndex_Discharge));
icoil_efit = fiesta_icoil(coilset, CoilCurrentsEfit);

%Initiate virtual B-field sensors centered on Rgeo
sensor_btheta = InitiateBSensors(EquilParams.r0_geom,EquilParams.z0_geom,R_Null);

%RZIP computes coefficients [A B C D] using the null field sensors
%Output C is used to compute the null-field PF coil currents 
%Outputs curlyM and curlyR are used to compute the plasma and eddy currents
rzip_config = fiesta_rzip_configuration( 'RZIP', config, vessel, {sensor_btheta} );
[A, B, C, D, curlyM, curlyR, Gamma, plasma_parameters, index, label_index, state] = ...
    response(rzip_config, Equil, 'rp', PlasmaResistPerp);            


%%%%%%%%%%%%%%%%%%%%%  COMPUTE OPTIMISED NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%%

%Update CoilWaveforms array with null-field values (Using NaN Mask)
RZIP_C = C;
CoilWaveforms = NullFieldWaveforms(CoilWaveforms,RZIP_C,sensor_btheta,TimeIndex_NullField);

%Save null-field coil currents for all coils at TimeIndex_NullField
CoilCurrentsNull = transpose(CoilWaveforms(:,TimeIndex_NullField));
icoil_Null = fiesta_icoil(coilset, CoilCurrentsNull);

%Compute null-field equilibrium using null-field coil currents
Equil_Null = fiesta_equilibrium('SMART-Null', config, Irod, icoil_Null );
EquilParams_Null = parameters(Equil_Null);


%%%%%%%%%%%%%%%%%%  COMPUTE BREAKDOWN FOR NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%

%Extract the null poloidal and toroidal B-field vector arrays
[BrData_Null,BzData_Null,BPhiData_Null,BpolData_Null,BtorData_Null] = ExtractBField(Equil_Null);

%Average null poloidal and toroidal fields over region of area R_null^2
[BpolAvg_Null,BtorAvg_Null] = ExtractNullBMin(EquilParams,BpolData_Null,BtorData_Null,R_Null);

%Compute the average connection length within the null-field sensor region
%InnerVesselDimensions=[VesselRMaxInner,VesselRMinInner,VesselZMaxInner,VesselZMinInner]
%Lc = ConnectionLength(EquilParams,InnerVesselDimensions);             % NOTE :: Issue with inbuilt connection length functions
a_eff = min([abs(EquilParams.r0_geom-VesselRMinInner),abs(EquilParams.r0_geom-VesselRMaxInner)]);
Lc = 0.25*a_eff*(BtorAvg_Null/BpolAvg_Null);

%Compute maximum loop voltage and E-field during solenoid ramp-down
[Vloop,Eloop,DeltaPhiSol] = LoopVoltage(CoilWaveforms,time,RSolCentreWinding,ZMaxSol,ZMinSol,EquilParams.rin); 
Eloop_eff = abs(Eloop)*(BtorAvg_Null/BpolAvg_Null); %[V/m]   %Rough estimate threshold Eloop condition for breakdown
%Generally Eloop_eff > 100 [V/m] for startup with ECRH      (An2015)
%Generally Eloop_eff > 1000 [V/m] for solenoid only startup (Lloyd1991)

%Compute breakdown (avalanche) timescale at given pressure (!!! NEEDS TESTING !!!)
Pressure = EquilParams.P0*(7.5e-7);    %[Torr]  (~2e-4 Torr)
[TauAvalanche,Pressure,Alpha,Vde]=AvalancheTimescale(Pressure,Eloop,Lc,ne,1.0,true);

%%%%%%%%%%%%%%%  COMPUTE DYNAMIC PLASMA & EDDY CURRENTS  %%%%%%%%%%%%%%%%%%

%Compute dynamic plasma and vessel eddy currents with new coil waveforms
[time_linear,time_adaptive,I_PF_output,V_PF_output,Ip_output,Vp_output,I_Passive] = ...
    DynamicCurrents(CoilWaveforms, time, curlyM, curlyR);

%Compute maximum change in each coil current - !!! MAKE THIS INTO A FUNCTION !!!
for j=1:length(I_PF_output(1,:))
    for i=2:length(time_adaptive)-1
        Delta_Ip_output(i) = ( (Ip_output(i)-Ip_output(i-1))/(time_adaptive(i)-time_adaptive(i-1)) )/1000;  %[A/ms]
        Delta_IPFoutput(i,j) = ( (I_PF_output(i,j)-I_PF_output(i-1,j))/(time_adaptive(i)-time_adaptive(i-1)) )/1000;  %[A/ms]
        Delta_VPFoutput(i,j) = ( (V_PF_output(i,j)-V_PF_output(i-1,j))/(time_adaptive(i)-time_adaptive(i-1)) )/1000;  %[V/ms]
    end
end
MinDelta_Ipoutput = min(Delta_Ip_output); MaxDelta_Ipoutput = max(Delta_Ip_output);    %[A/ms]
MinDelta_IPFoutput = min(Delta_IPFoutput); MaxDelta_VPFoutput = max(Delta_VPFoutput);  %[A/ms]
MinDelta_VPFoutput = min(Delta_VPFoutput); MaxDelta_IPFoutput = max(Delta_IPFoutput);  %[V/ms]

%Compute solenoid OH flux swing employing the Ejima-Wesley coefficient
Cew = 0.40*EquilParams.aspectratio;                           %[-]      (Gryaznevich2006)
DeltaPhiSolew = Cew*mu0*EquilParams.r0_geom*max(Ip_output);   %[Vs]     (Menard2016 pg36)

%Extract Vessel Eddy Currents during discharge and null-field (time='false' for absolute max)
VesselEddyCurrents = ExtractPassiveCurrents(I_Passive,time_adaptive,time(TimeIndex_Discharge));
VesselEddyCurrentsNull = ExtractPassiveCurrents(I_Passive,time_adaptive,time(TimeIndex_NullField));
% ISSUE :: EDDY CURRENTS CHANGE OVER THE AVALANCHE TIMESCALE - TAKING THEM AT TimeIndex_NullField IS SLIGHTLY WRONG


%%%%%%%%%%%%%%% PLOT VESSEL AND COIL FILAMENT OVERVIEW %%%%%%%%%%%%%%%%%% 

figure; axes;
set(gca, 'DataAspectRatio', [1,1,1], 'NextPlot', 'add')
coil = plot(coilset);
fil = plot(vessel);
set(fil, 'EdgeColor', 'k')
set(coil, 'EdgeColor', 'k')
tau=get(vessel, 'tau');
r=get(vessel, 'r'); z=get(vessel, 'z');
%set(line(r(i), z(i)), 'LineStyle', 'none', 'marker' , '*', 'markersize', 10);
%set(line(r(j), z(j)), 'LineStyle', 'none', 'marker' , 'o', 'markersize', 10);
set(gca,'XLim',[0.00 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_VesselFilaments';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

figure; axes;
set(gca, 'DataAspectRatio', [1,1,1], 'NextPlot', 'add')
coil = plot(coilset);
fil = plot(vessel);
set(fil, 'EdgeColor', 'k')
set(coil, 'EdgeColor', 'k')
tau=get(vessel, 'tau');
r=get(vessel, 'r'); z=get(vessel, 'z');
%set(line(r(i), z(i)), 'LineStyle', 'none', 'marker' , '*', 'markersize', 10);
%set(line(r(j), z(j)), 'LineStyle', 'none', 'marker' , 'o', 'markersize', 10);
set(gca,'XLim',[0.05 0.45]);
set(gca,'YLim',[0.45 1.00]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_VesselFilaments_Closeup';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

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
Filename = '_VirtualBSensors';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({Equil,sensor_btheta},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%  PLOT NULL-FIELD PHI SURFACES  %%%%%%%%%%%%%%%%%%%%%

%Plot the optimised null-field phi
Title = {'SMART Null-field Equilibrium iter(0)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_NullPhi';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({Equil_Null},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%% PLOT NULL-FIELD BPOL  %%%%%%%%%%%%%%%%%%%%%

%Log poloidal and toroidal magnetic fields to show details (Sol Obscures)
logBpolData_Null = log10(BpolData_Null);
logBtorData_Null = log10(BtorData_Null);
%Plot the optimised null-field phi
Title = {'SMART Null-field iter(0)',' '};
CbarLabel = 'Null-field B_{\theta} log_{10}([T])';
Filename = '_NullBpol';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({logBpolData_Null},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%% PLOT COIL CURRENT WAVEFORMS %%%%%%%%%%%%%%%%%%%%%%%% 

%Plot figure showing dynamic coil currents
figure('units','inch','position',[10 10 12 12]);
subplot(2,1,1)
plot(time_adaptive*1000, I_PF_output/1000, 'LineWidth',2);
title(gca,'SMART Coil Current Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Coil Current I [kA]');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
%%%%%
subplot(2,1,2)
plot(time_adaptive(1:end-1)*1000,Delta_IPFoutput/1000, 'LineWidth',2)
title(gca,'SMART Delta Coil Current Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Delta Coil Current \Delta I (kA ms^{-1})');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
%%%%%
Filename = '_CurrentWaveforms';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%% PLOT COIL VOLTAGE WAVEFORMS %%%%%%%%%%%%%%%%%%%%%%%% 

%Plot figure showing dynamic coil currents
figure('units','inch','position',[10 10 12 12]);
subplot(2,1,1)
plot(time_adaptive*1000, V_PF_output/1000, 'LineWidth',2);
title(gca,'SMART Coil Voltage Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Coil Voltage V [kV]');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
%%%%%
subplot(2,1,2)
plot(time_adaptive(1:end-1)*1000,Delta_VPFoutput/1000, 'LineWidth',2)
title(gca,'SMART Delta Coil Voltage Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Delta Coil Voltage \Delta V (kV ms^{-1})');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
%%%%%
Filename = '_VoltageWaveforms';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%% PLOT PLASMA CURRENT %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot plasma current over full timescale
figure('units','inch','position',[10 10 12 12]);
subplot(2,1,1); hold on; grid on;
plot(time_adaptive*1000, Ip_output/1000, 'LineWidth',2)
title(gca,'SMART Plasma Current iter(0)');
legend(gca,'I_{p}', 'FontSize',16); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Plasma Current I_{p} (kA)');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
%%%%%
subplot(2,1,2); hold on; grid on;
plot(time_adaptive(1:end-1)*1000, Delta_Ip_output/1000, 'LineWidth',2)
title(gca,'SMART Plasma Current iter(0)');
legend(gca,'dI_{p}/dt', 'FontSize',16); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Delta Plasma Current (kA/ms)');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%set(gca,'XLim',[time(3)*1e3, time(5)*1e3]);        %Breakdown Ramp Closeup
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
%%%%%
Filename = '_PlasmaCurrent';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%% PLOT TOTAL EDDY CURRENT %%%%%%%%%%%%%%%%%%%%%%%%% 

%Sum all filaments (row-wise) to get total net passive current
Net_IPassive = sum(I_Passive,2);
%Plot net passive current density over full timescale
close all
figure('units','inch','position',[12 12 8 8]); hold on; grid on;
plot(time_adaptive*1000, Net_IPassive/1000, 'LineWidth',2)
title(gca,'Net SMART Eddy Current iter(0)');
legend(gca,'Net Eddy Current'); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Net Vessel Current I_{Eddy} [kA]');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca,'YLim',[min(Net_IPassive/1000)*1.15 max(Net_IPassive/1000)*1.20]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
Filename = '_1DEddyCurrent';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%% PLOT 2D RESOLVED EDDY CURRENTS %%%%%%%%%%%%%%%%%%%%%%   

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

%%%%%     %%%%%     %%%%%     %%%%%     %%%%%     %%%%%

%Obtain the filament variables r and z
ptmp = get(vessel,'passives');
ftmp = get(ptmp,'filaments');
RR = get(ftmp(:),'r'); %dim 1*number of filaments
ZZ = get(ftmp(:),'z'); %dim 1*number of filaments
%Plot eddy currents within a cross-section of the vessel
close all
figure; hold on; axis equal;
plot(coilset);
scatter3(RR,ZZ,VesselEddyCurrentsNull/1000,100,VesselEddyCurrentsNull/1000,'filled');
title('SMART Vessel Null Eddy Currents iter(0)');
view(2) %2D view
colormap(plasma);
cbar = colorbar;
cbar.Label.String = 'Eddy-Current I_{eddy} [kA]';
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_EddyCurrentNull';
saveas(gcf, strcat(ProjectName,Filename,FigExt));


%%%%%%%%%%%%%%%%%%%%%%  PLOT PASCHEN CURVES  %%%%%%%%%%%%%%%%%%%%%

%Calculate Paschen Curve and determine breakdown
PressureLimits = [1e-6, 1e-3]; Resolution = 25000;
PressureArray = linspace(PressureLimits(1),PressureLimits(2),Resolution);

for i=1:length(PressureArray)
    PaschenCurve(i) = (PressureArray(i)*1.25e4)/(log(510.0*PressureArray(i)*Lc));
    if PaschenCurve(i) < 0
        PaschenCurve(i) = nan;
    end
end

for i=1:length(PressureArray)
    EloopArray(i) = abs(Eloop);
end

close all
figure('units','inch','position',[12 12 8 8]); hold on; grid on;
title(gca,'SMART Paschen Breakdown');
LegendString = {'Paschen Curve','Max E_{loop}'};
plot(PressureArray,PaschenCurve, 'k', 'LineWidth',2)
plot(PressureArray,EloopArray, 'r', 'LineWidth',2)
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Pressure P (Torr)');
ylabel(gca,'E_{min} for Breakdown (V/m)');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'XLim',[2.5e-6 1e-3]);
set(gca,'YLim',[0.1 100]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
Filename = '_PaschenBreakdown';
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
PressureR=abs(Force_max(1))/(FilamentArea);	%[Pa]
PressureZ=abs(Force_max(3))/(FilamentArea);	%[Pa]

%Stress acting on vessel wall is the combined force divided by the unit filiment area
%These are directional, some are negative and some are positive
StressR=(Force_fil(:,1))/(FilamentArea);	%[Pa]
StressZ=(Force_fil(:,3))/(FilamentArea);	%[Pa]
%Obtain maximum radial and axial stresses - either positive or negative
StressR_max=max(abs(StressR));				%[Pa]
StressZ_max=max(abs(StressZ));				%[Pa]


%%%%%%%%%%%%%%%%%%%%%% PLOT VESSEL EDDY STRESSES %%%%%%%%%%%%%%%%%%%%%%%%

%Plot figure showing vessel eddy stresses
close all
figure; hold on; axis equal;
plot(coilset);
%plot(vessel);
quiver(R_Fil_Array,Z_Fil_Array,StressR,StressZ,'color',[1 0 0],'AutoScale','on');
title('SMART Vessel Eddy-Stresses iter(1)');
view(2) %2D view
legend(gca,'hide');
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_EddyStresses';
saveas(gcf, strcat(ProjectName,Filename,FigExt));
close('all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      COMPUTE PERTURBED EQUILIBRIA                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Apply small perturbation(s) to the efitGeometry values
RGeo_Pert = efitGeometry(1)+PertGeometry_Init(1);
ZGeo_Pert = efitGeometry(2)+PertGeometry_Init(2);
rGeo_Pert = efitGeometry(3)+PertGeometry_Init(3);
Kappa_Pert = efitGeometry(4)+PertGeometry_Init(4);
delta_Pert = efitGeometry(5)+PertGeometry_Init(5);
PertGeometry = [RGeo_Pert, ZGeo_Pert, rGeo_Pert, Kappa_Pert, delta_Pert];

%Calculate feedbackCoil currents required to fix the applied perturbation
feedbackCoils = {'PF2','Div2'};
%[feedback_config, signals, weights, index] = efit_shape_controller(config, feedbackCoils, PertGeometry);
feedback = shape_controller(config, feedbackCoils, RGeo_Pert, ZGeo_Pert, rGeo_Pert, Kappa_Pert, delta_Pert);
Equil_Pert = set(Equil, config, 'feedback', feedback);                                 
EquilParams_Pert = parameters(Equil_Pert);
% ISSUE :: Equil_Pert is not shaped the same as Equil  
% ISSUE :: Coil currents are unreliable until Equil_Pert ~ Equil

%Extract coil currents required to offset perturbation
icoil_pert = get(Equil_Pert,'icoil'); 
CoilCurrents_Pert = get(icoil_pert,'currents');
%}

%Obtain the Real Vertical Growth Rate from RZIp  (Gamma = eig(-curlyM\curlyR), sort for positive values)
if length(Gamma(Gamma>0)) > 0; Gamma_Real = Gamma(Gamma>0); else Gamma_Real = 0; end     %[s-1]

%If feedback fails, overwrite with default equilibrium for now
Equil_Pert = Equil; EquilParams_Pert = EquilParams;
PertGeometry = efitGeometry;
icoil_pert = icoil_efit; CoilCurrents_Pert = CoilCurrentsEfit;
%end

%%%%%%%%%%%%%%%%%%%%%% PLOT PERTURBED EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%

%Plot perturbed equilibrium following convergence
Title = {'SMART Perturbed Equilibrium \Psi(R,Z)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_PerturbedEquilibrium';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({Equil_Pert},{rGrid,zGrid},Title,CbarLabel,SaveString);
close('all')

%%%%%%%%%%%%%%%%%%%%%% PLOT VERTICAL GROWTH RATES  %%%%%%%%%%%%%%%%%%%%%%

%Plot the vertical growth rates as calculated by RZIp and manually
figure('units','inch','position',[12 12 8 8]); hold on;
plot(Gamma,'LineWidth',2);
plot(ones(length(Gamma)),'k--','LineWidth',1.5);
title('SMART Vertical Growth Rates');
LegendString = {strcat('RZIp Gamma: ',string(Gamma_Real),' s^{-1}')};
legend(LegendString);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'Eigenvalue (Sorted) [-]');
ylabel(gca,'Growth Rate \gamma [s^{-1}]');
Filename = '_VerticalGrowthRates';
saveas(gcf, strcat(ProjectName,Filename,FigExt));
close('all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% ADDING EDDIES TO DISCHARGE EQUIL %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%  RECOMPUTE DISCHARGE EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%

%Book-keeping for the start of each loop
close all

%Apply alterations to efit CoilWaveforms to increase simulation stability
CoilWaveforms(:,TimeIndex_Discharge) = CoilWaveforms_Init(:,TimeIndex_Discharge);   %Initial discharge guesses are more stable
CoilWaveforms(iDiv1,TimeIndex_Discharge) = I_Div1_Equil;                            %Slightly increase IDiv1 and retry if required

%Compute equilibrium (Psi(R,Z)) from the supplied jprofile, icoil and geometry
%Returns target equilibrium and CoilWaveforms for PF1 and PF2 at requested time_Index
if isa(coilset,'fiesta_loadassembly') == 0; coilset = fiesta_loadassembly(coilset, vessel); end
[Equil_Passive,EquilParams_Passive,CoilWaveforms_Passive,efitGeometry_Passive,config_passive] = ...
    efit(jprofile,Irod,'nullconfig',efitGeometry_Init,CoilWaveforms,VesselEddyCurrents,TimeIndex_Discharge);
%NOTE :: PERHAPS REPLACE CoilWaveforms WITH CoilWaveforms_Init TO REPLACE NAN NULL-FIELD VALUES FOR NullFieldWaveforms FUNCTION?

%Save discharge coil currents for all coils at TimeIndex_Discharge
CoilCurrentsEfit_Passive = transpose(CoilWaveforms_Passive(:,TimeIndex_Discharge));
icoil_efit_passive = fiesta_icoil(coilset, [CoilCurrentsEfit_Passive,VesselEddyCurrents]);

%Initiate virtual B-field sensors centered on Rgeo
sensor_btheta_passive = InitiateBSensors(EquilParams_Passive.r0_geom,EquilParams_Passive.z0_geom,R_Null);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NEED TO RE-CALCULATE PLASMA AND EDDY WITH NEW EFIT EQUILIBRIUM - ADDING EDDY CURRENTS VARIES THE PLASMA CURRENT 
%rzip_config = fiesta_rzip_configuration( 'RZIP', config_passive, vessel, {sensor_btheta_passive} );   %CONFIG ACCOUNTS FOR VESSEL FILAMENTS
%[A_Passive,B_Passive,C_Passive,D_Passive,curlyM_Passive,curlyR_Passive,Gamma_Passive,plasma_parameters_Passive,index_Passive,label_index_Passive,state_Passive,condRZIp_Passive] = ...
%    response(rzip_config, Equil_Passive, 'rp', PlasmaResistPerp);
%ISSUE :: CV = greens(coilset, vessel) returns NaNs, likely because the vessel filaments are overlapping the vessel 'coils' in Equil_Passive
%      :: Equil_Passive contains vessel filament 'coils', while RZIP is being fed 'vessel' filaments, both of which have the same coordinates
%      :: SET RZIP UP WITH ONLY PF COILS FOR ALL FILAMENTS

%%%%%%%%%%%%%%%%%%%%%  COMPUTE OPTIMISED NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%%

%Update CoilWaveforms array with null-field values (Using NaN Mask)
%   RZIP_C = C;
%   CoilWaveforms = NullFieldWaveforms(CoilWaveforms_Passive, RZIP_C, sensor_btheta_passive, TimeIndex_NullField);
%ISSUE  :: RZIP_C is computed from RZIP - Need recomputed RZIP with eddy currents
%NOTE   :: REMEMBER THAT NULLFIELDWAVEFORMS EXPECTS NaN VALUES FOR COILS TO BE UPDATED

%Extract previously calculated efit coil currents without eddys
CoilCurrentsNull = transpose(CoilWaveforms_Passive(:,TimeIndex_NullField)); %Null-field coil currents without eddys
CoilAndVesselCurrents = [CoilCurrentsNull, VesselEddyCurrentsNull];         %n=5+n_fil; coil + vessel filaments
icoil_Null_Passive = fiesta_icoil(coilset, CoilAndVesselCurrents);
%NOTE :: CoilWaveforms_Passive here does NOT have the re-optimised coil currents!

%Compute null-field equilibrium using null-field coil and vessel eddy currents
Equil_Null_Passive = fiesta_equilibrium('SMART-Null', config_passive, Irod, icoil_Null_Passive);
EquilParams_Null_Passive = parameters(Equil_Null_Passive);

%Extract the new coil currents from the null-field equilibrium:
icoil_Null_Passive = get(Equil_Null_Passive,'icoil'); 
CoilCurrentsNull_Passive = get(icoil_Null_Passive,'currents');

%%%%%%%%%%%%%%%%%%  COMPUTE BREAKDOWN FOR NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%

%Extract the null poloidal and toroidal B-field vector arrays
[BrData_Null_Passive,BzData_Null_Passive,BPhiData_Null_Passive,BpolData_Null_Passive,BtorData_Null_Passive] = ...
    ExtractBField(Equil_Null_Passive);

%Minimum null poloidal and toroidal fields, averaged over region of area R_null^2
[BpolAvg_Null_Passive,BtorAvg_Null_Passive] = ...
    ExtractNullBMin(EquilParams_Passive,BpolData_Null_Passive,BtorData_Null_Passive,R_Null);
    %NOTE - USES EQUILPARAMS_PASSIVE ONLY TO EXTRACT RGeo and ZGeo (Would be less confusing to use EQUILPARAMS_NULL)

%Compute the average connection length within the null-field sensor region
%InnerVesselDimensions=[VesselRMaxInner,VesselRMinInner,VesselZMaxInner,VesselZMinInner]
%Lc = ConnectionLength(EquilParams,InnerVesselDimensions);
a_eff = min([abs(EquilParams_Passive.r0_geom-VesselRMinInner),abs(EquilParams_Passive.r0_geom-VesselRMaxInner)]);
Lc_Passive = 0.25*a_eff*(BtorAvg_Null_Passive/BpolAvg_Null_Passive);

%Compute maximum loop voltage and E-field during solenoid ramp-down
[Vloop_Passive,Eloop_Passive,DeltaPhiSol_Passive] = LoopVoltage(CoilWaveforms_Passive,time,RSolCentreWinding,ZMaxSol,ZMinSol,EquilParams_Passive.rin);
Eloop_eff_Passive = abs(Eloop_Passive)*(BtorAvg_Null_Passive/BpolAvg_Null_Passive); %[V/m]   %Rough estimate threshold Eloop condition for breakdown
%Generally Eloop_eff > 100 [V/m] for startup with ECRH      (An2015)
%Generally Eloop_eff > 1000 [V/m] for solenoid only startup (Lloyd1991)

%Compute breakdown (avalanche) timescale at given pressure (!!! NEEDS TESTING !!!)
Pressure_Passive = EquilParams_Passive.P0*(7.5e-7);    %[Torr]  (~2e-4 Torr)
[TauAvalanche_Passive,Pressure_Passive,Alpha_Passive,Vde_Passive] = ...
    AvalancheTimescale(Pressure_Passive,Eloop_Passive,Lc_Passive,ne,1.0,true);


%%%%%%%%%%%%%%%  COMPUTE DYNAMIC PLASMA & EDDY CURRENTS  %%%%%%%%%%%%%%%%%%

%Recompute dynamic plasma and vessel eddy currents with new coil waveforms
%ISSUE :: curlyM and curlyR are computed from RZIP - Need recomputed RZIP with eddy currents

%Update dynamic plasma and vessel eddy currents
%[time_linear_passive,time_adaptive_passive,I_PF_output_Passive,V_PF_output_Passive,Ip_output_Passive,Vp_output_Passive,I_Passive_Passive] = ...
%    DynamicCurrents(CoilWaveforms_Passive, time, curlyM_Passive, curlyR_Passive);

%Compute maximum change in each coil current --- !!! MAKE THIS INTO A FUNCTION !!!
% !!! ONCE MADE INTO A FUNCTIONm ENSURE OUTPUT PASSIVE VARIABLES AT THIS POINT !!!
%for j=1:length(I_PF_output(1,:))
%    for i=2:length(time_adaptive)-1
%        Delta_IPFoutput(i,j) = ( (I_PF_output(i,j)-I_PF_output(i-1,j))/(time_adaptive(i)-time_adaptive(i-1)) )/1000;  %[A/ms]
%        Delta_VPFoutput(i,j) = ( (V_PF_output(i,j)-V_PF_output(i-1,j))/(time_adaptive(i)-time_adaptive(i-1)) )/1000;  %[V/ms]
%    end
%end
%MinDelta_IPFoutput = min(Delta_IPFoutput); MaxDelta_VPFoutput = max(Delta_VPFoutput);  %[A/ms]
%MinDelta_VPFoutput = min(Delta_VPFoutput); MaxDelta_IPFoutput = max(Delta_IPFoutput);  %[V/ms]

%Compute solenoid OH flux swing via the Ejima-Wesley coefficient
%Cew = 0.40*EquilParams_Passive.aspectratio;                                           %[-]      (Gryaznevich2006)
%DeltaPhiSolew_Passive = Cew*mu0*EquilParams_Passive.r0_geom*max(Ip_output_Passive);   %[Vs]     (Menard2016 pg36)

%Extract Vessel Eddy Currents during discharge (time='false' for absolute max)
%VesselEddyCurrents_Passive = ...
%    ExtractPassiveCurrents(I_Passive_Passive,time_adaptive_Passive,time(TimeIndex_Discharge));
%VesselEddyCurrentsNull_Passive...
%    = ExtractPassiveCurrents(I_Passive_Passive,time_adaptive_Passive,time(TimeIndex_NullField));

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
PlotEquilibrium({Equil_Null_Passive},{rGrid,zGrid},Title,CbarLabel,SaveString);

CoilCurrentsNull(1:nPF)
CoilCurrentsNull_Passive(1:nPF)

%%%%%%%%%%%%%%%%%%%%%% PLOT NULL-FIELD BPOL  %%%%%%%%%%%%%%%%%%%%%

%Log poloidal and toroidal magnetic fields to show details (Sol Obscures)
logBpolData_Null_Passive = log10(BpolData_Null_Passive);
logBtorData_Null_Passive = log10(BtorData_Null_Passive);
%Plot the optimised null-field phi
Title = {'SMART Null-field iter(1)',' '};
CbarLabel = 'Null-field B_{\theta} log_{10}([T])';
Filename = '_NullBpol_Passive';
SaveString = strcat(ProjectName,Filename,FigExt);
PlotEquilibrium({logBpolData_Null_Passive},{rGrid,zGrid},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%%  PLOT PASCHEN CURVES  %%%%%%%%%%%%%%%%%%%%%

%Calculate Paschen Curve and determine breakdown
PressureLimits = [1e-6, 1e-3]; Resolution = 25000;
PressureArray = linspace(PressureLimits(1),PressureLimits(2),Resolution);

for i=1:length(PressureArray)
    PaschenCurve_Passive(i) = (PressureArray(i)*1.25e4)/(log(510.0*PressureArray(i)*Lc_Passive));
    if PaschenCurve_Passive(i) < 0
        PaschenCurve_Passive(i) = nan;
    end
end

for i=1:length(PressureArray)
    EloopArray_Passive(i) = abs(Eloop_Passive);
end

close all
figure('units','inch','position',[12 12 8 8]); hold on; grid on;
title(gca,'SMART Paschen Breakdown');
LegendString = {'Paschen Curve','Max E_{loop}'};
plot(PressureArray,PaschenCurve_Passive, 'k', 'LineWidth',2)
plot(PressureArray,EloopArray_Passive, 'r', 'LineWidth',2)
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Pressure P (Torr)');
ylabel(gca,'E_{min} for Breakdown (V/m)');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca,'XLim',[2.5e-6 1e-3]);
set(gca,'YLim',[0.1 100]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
Filename = '_PaschenBreakdown_Passive';
saveas(gcf, strcat(ProjectName,Filename,FigExt));
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          DATA I/O MANAGEMENT                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ISSUE :: Need to impliment these functions into the eddy loop once complete
% ISSUE :: Functions overwrite data with most recent (which is fine) but
%          having an optional method to save a 'iter_movie' file would be good
% ISSUE :: fileID does not mean anything at the moment, may be needed later
%          with regards to the above 'iter_movie' file.

%{
%Stuff required for I/O if eddy currents are not computed
Equil_Passive = Equil; Equil_Null_Passive = Equil_Null;
EquilParams_Passive = EquilParams; config_passive = config;
icoil_efit_passive = icoil_efit; icoil_Null_Passive = icoil_Null;
Vloop_Passive = Vloop; Vloop_Lc_Passive = Vloop_Lc;
Eloop_Passive = Eloop; Eloop_eff_Passive = Eloop_eff;
Lc_Passive = Lc;
BpolAvg_Null_Passive = BpolAvg_Null; BtorAvg_Null_Passive = BtorAvg_Null;
%}

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%
%Create subdirectory for equilibrium related data
EquilDir = strcat(ASCIIDir,'Equil_Data/'); mkdir(EquilDir);

%Write 2D, 1D and 0D equilibrium values to text files once per iteration
[fileID] = WriteEquilibrium(Equil_Passive, config_passive, EquilDir, '', false);
[fileID] = WriteEquilibrium(Equil_Null_Passive, config_passive, EquilDir, '', true);
[fileID] = WriteEquilibrium(Equil_Pert, config, EquilDir, '_Pert', false);

%Write initial target geometry, efit geometry and perturbed geometry
[fileID] = WriteGeometry(efitGeometry_Init, EquilDir, 'efit_Geometry_Init.txt');
[fileID] = WriteGeometry(efitGeometry, EquilDir, 'efit_Geometry_Equil.txt');
[fileID] = WriteGeometry(PertGeometry, EquilDir, 'efit_Geometry_Pert.txt');

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Create subdirectory for coil current related data
icoilDir = strcat(ASCIIDir,'icoil_Data/'); mkdir(icoilDir);

Filename = strcat(icoilDir,'icoil_position.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s\r\n', 'Coil','R [m]  ','Z [m]');
fprintf(fileID,'%s %0.5f %0.5f\r\n', 'Sol ',RSolCentre,ZMaxSol);
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
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_efit_passive.Sol'; icoil_efit_passive.PF1'; icoil_efit_passive.PF2'; icoil_efit_passive.Div1'; icoil_efit_passive.Div2']);

Filename = strcat(icoilDir,'Perturbed_icoil.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol','PF1','PF2','Div1','Div2');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_pert.Sol'; icoil_pert.PF1'; icoil_pert.PF2'; icoil_pert.Div1'; icoil_pert.Div2']);

Filename = strcat(icoilDir,'Null_icoil.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol','PF1','PF2','Div1','Div2');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_Null_Passive.Sol'; icoil_Null_Passive.PF1'; icoil_Null_Passive.PF2'; icoil_Null_Passive.Div1'; icoil_Null_Passive.Div2']);

%Extract coil current time-traces
ISol=I_PF_output(:,1);   %ISol
IPF1=I_PF_output(:,2);   %IPF1
IPF2=I_PF_output(:,3);   %IPF2
IDiv1=I_PF_output(:,4);  %IDiv1
IDiv2=I_PF_output(:,5);  %IDiv2
Filename = strcat(icoilDir,'CoilCurrents.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s %s\r\n', 'time_adaptive [ms]','I_Sol [A]','I_PF1 [A]','I_PF2 [A]','I_Div1 [A]','I_Div2 [A]');
fprintf(fileID,'%1.12f %0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[time_adaptive'*1000; ISol'; IPF1'; IPF2'; IDiv1'; IDiv2']);

%Extract delta coil current time-traces
dISol=Delta_IPFoutput(:,1);   %ISol
dIPF1=Delta_IPFoutput(:,2);   %IPF1
dIPF2=Delta_IPFoutput(:,3);   %IPF2
dIDiv1=Delta_IPFoutput(:,4);  %IDiv1
dIDiv2=Delta_IPFoutput(:,5);  %IDiv2
Filename = strcat(icoilDir,'DeltaCoilCurrents.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s %s\r\n', 'time_adaptive [ms]','Del_I_Sol [A/ms]','Del_I_PF1 [A/ms]','Del_I_PF2 [A/ms]','Del_I_Div1 [A/ms]','Del_I_Div2 [A/ms]');
fprintf(fileID,'%1.12f %0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[time_adaptive(1:end-1)'*1000; dISol'; dIPF1'; dIPF2'; dIDiv1'; dIDiv2']);

%Extract coil voltage time-traces
VSol=V_PF_output(:,1);     %VSol
VPF1=V_PF_output(:,2);     %VPF1
VPF2=V_PF_output(:,3);     %VPF2
VDiv1=V_PF_output(:,4);     %VDiv1
VDiv2=V_PF_output(:,5);     %VDiv2
Filename = strcat(icoilDir,'CoilVoltages.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s %s\r\n', 'time_adaptive [ms]','V_Sol [V]','V_PF1 [V]','V_PF2 [V]','V_Div1 [V]','V_Div2 [V]');
fprintf(fileID,'%1.12f %0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[time_adaptive'*1000; VSol'; VPF1'; VPF2'; VDiv1'; VDiv2']);

%Extract delta coil voltage time-traces
dVSol=Delta_VPFoutput(:,1);   %ISol
dVPF1=Delta_VPFoutput(:,2);   %IPF1
dVPF2=Delta_VPFoutput(:,3);   %IPF2
dVDiv1=Delta_VPFoutput(:,4);  %IDiv1
dVDiv2=Delta_VPFoutput(:,5);  %IDiv2
Filename = strcat(icoilDir,'DeltaCoilVoltages.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s %s\r\n', 'time_adaptive [ms]','Del_V_Sol [V/ms]','Del_V_PF1 [V/ms]','Del_V_PF2 [V/ms]','Del_V_Div1 [V/ms]','Del_V_Div2 [V/ms]');
fprintf(fileID,'%1.12f %0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[time_adaptive(1:end-1)'*1000; dVSol'; dVPF1'; dVPF2'; dVDiv1'; dVDiv2']);

%Extract max/min delta coil currents
Filename = strcat(icoilDir,'Delta_IPF.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'I_Sol [A/ms]','I_PF1 [A/ms]','I_PF2 [A/ms]','I_Div1 [A/ms]','I_Div2 [A/ms]');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[MaxDelta_IPFoutput(1)'; MaxDelta_IPFoutput(2)'; MaxDelta_IPFoutput(3)'; MaxDelta_IPFoutput(4)'; MaxDelta_IPFoutput(5)']);
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[MinDelta_IPFoutput(1)'; MinDelta_IPFoutput(2)'; MinDelta_IPFoutput(3)'; MinDelta_IPFoutput(4)'; MinDelta_IPFoutput(5)']);

%Extract max/min delta coil voltages
Filename = strcat(icoilDir,'Delta_VPF.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'V_Sol [V/ms]','V_PF1 [V/ms]','V_PF2 [V/ms]','V_Div1 [V/ms]','V_Div2 [V/ms]');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[MaxDelta_VPFoutput(1)'; MaxDelta_VPFoutput(2)'; MaxDelta_VPFoutput(3)'; MaxDelta_VPFoutput(4)'; MaxDelta_VPFoutput(5)']);
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[MinDelta_VPFoutput(1)'; MinDelta_VPFoutput(2)'; MinDelta_VPFoutput(3)'; MinDelta_VPFoutput(4)'; MinDelta_VPFoutput(5)']);

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Create subdirectory for dynamic current data (RZIP)
RZIPDir = strcat(ASCIIDir,'RZIP_Data/'); mkdir(RZIPDir);

Filename = strcat(RZIPDir,'time.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s\r\n','time_adaptive [ms]');
fprintf(fileID,'%1.12f\r\n',time_adaptive'*1000);

Filename = strcat(RZIPDir,'Ip.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'time_adaptive [ms]','Ip_output [A]');
fprintf(fileID,'%1.12f %1.12f\r\n',[time_adaptive'*1000; Ip_output']);

Filename = strcat(RZIPDir,'IPass.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'time_adaptive [ms]','Net_I_Passive [A]');
fprintf(fileID,'%1.12f %1.12f\r\n',[time_adaptive'*1000; Net_IPassive']);

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Misc outputs, save unordered in main RawData directory

Filename = strcat(ASCIIDir,'Eta.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'Eta_Perp [Ohm]', 'Eta_Para [Ohm]');
fprintf(fileID,'%1.12f %1.12f\r\n', PlasmaResistPerp', PlasmaResistPara');

Filename = strcat(ASCIIDir,'Bpol.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'Null_Bpol [T]', 'Null_Btor [T]');
fprintf(fileID,'%1.12f %1.12f\r\n', BpolAvg_Null_Passive', BtorAvg_Null_Passive');

Filename = strcat(ASCIIDir,'IRod.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s\r\n', 'IRod [A]');
fprintf(fileID,'%1.12f\r\n', EquilParams_Passive.irod');

Filename = strcat(ASCIIDir,'TauVessel.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s\r\n', 'TauVessel [ms]');
fprintf(fileID,'%1.12f\r\n', TauVessel*1000');

Filename = strcat(ASCIIDir,'betaP.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'betaP [%]', 'betaP_Pert [%]');
fprintf(fileID,'%1.12f %1.12f\r\n', EquilParams_Passive.betap', EquilParams_Pert.betap');

Filename = strcat(ASCIIDir,'TauAvalanche.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s\r\n', 'TauAvalanche [ms]');
fprintf(fileID,'%1.12f\r\n', TauAvalanche_Passive*1000');

Filename = strcat(ASCIIDir,'Lc.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s \r\n', 'Lc [m]');
fprintf(fileID,'%1.12f\r\n', Lc_Passive');

Filename = strcat(ASCIIDir,'VLoop.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s\r\n', 'Vloop [V]', '     DeltaPhi [Vs]', ' Eloop [V/m]  ', ' Eloop_eff [V/m]');
fprintf(fileID,'%1.12f %1.12f %1.12f %1.12f\r\n', Vloop_Passive', DeltaPhiSol_Passive', Eloop_Passive', Eloop_eff_Passive');

Filename = strcat(ASCIIDir,'MaxStress.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'StressR_max [Pa]', 'StressZ_max [Pa]');
fprintf(fileID,'%1.12f %1.12f\r\n', StressR_max', StressZ_max');

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Done!
disp([ ' ' ]);
disp([ 'Done!' ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%      FUNCTIONS      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%AVALANCHE FUNCTION TESTING AREA
%NOTES ::
%IF AVALANCHE APPROACHES FROM NEGATIVE -- THEN ALPHA^{-1} < Lc AND BREAKDOWN IS NOT POSSIBLE
%   IN THIS CASE, TAUAVALANCHE WILL CROSS THE X-AXIS AT ne0 = BreakdownFrac*ne
%   NOTABLY, THE GROWTH RATE IS NEGATIVE, SO TAUAVALANCHE INCREASES WITH ne0 
%   THIS DATA IS MOSTLY USELESS, SAVE TO SAY THAT THERE IS NO 'REAL' AVALANCHE TIMESCALE
%IF AVALANCHE APPROACHES FROM POSITIVE -- THEN ALPHA^{-1} > Lc AND BREAKDOWN IS ACHIEVED
%   IN THIS CASE, TAUAVALANCHE WILL APPROACH 0 WITH INCREASING ne0
%   NOTABLY, THE GROWTH RATE IS POSITIVE, SO TAUAVALANCHE DECREASES WITH ne0
%   THIS DATA CAN BE USED TO SET A LIMIT FOR TAUR1 FOR A GIVEN PRE-IONISATION DENSITY
%{
%%
ne0 = linspace(1e0,1e18,10000);                 %ne = 1.1772e18 [m-3] Phase1
DashedLine = transpose(zeros(length(ne0),1));

for i=1:length(ne0)
    %Compute breakdown (avalanche) timescale at given pressure (!!! NEEDS TESTING !!!)
    Pressure = EquilParams.P0*(7.5e-7);    %[Torr]  (~2e-4 Torr)
    Pressure = 1.0e-4;
    [TauAvalanche,Pressure,Alpha,Vde]=AvalancheTimescale(Pressure,Eloop,Lc,ne,ne0(i),false);
    TauAvalancheArray(i) = TauAvalanche;
end

for i=1:length(ne0)
    %Compute breakdown (avalanche) timescale at given pressure (!!! NEEDS TESTING !!!)
    Pressure_Passive = EquilParams_Passive.P0*(7.5e-7);    %[Torr]  (~2e-4 Torr)
    [TauAvalanche_Passive,Pressure_Passive,Alpha_Passive,Vde_Passive] = ...
        AvalancheTimescale(Pressure_Passive,Eloop_Passive,Lc_Passive,ne,ne0(i),false);
    TauAvalancheArray_Passive(i) = TauAvalanche_Passive;
end

figure; hold on;
plot(ne0,TauAvalancheArray*1000, 'k-', 'LineWidth',2)
plot(ne0,TauAvalancheArray_Passive*1000, 'r-', 'LineWidth',2)
plot(ne0,DashedLine, 'k--', 'LineWidth',1)
LegendString = {'No Eddys', 'With Eddys'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Pre-Ionisation Density (m^{-3})');
ylabel(gca,'Avalanche Timescale [ms]');
set(gca, 'xScale', 'log')
%set(gca,'YLim',[-0.5,10]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
%%
%} 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Equilibrium,EquilParams,efitCoilWaveforms,efitGeometry,config]= ...
    efit(Jprofile,Irod,InputConfig,InputGeometry,InputCoilWaveforms,InputVesselWaveforms,TimeIndex)

    %Obtain required global variables
    global Grid; global coilset;
    global efitCoils;
    global iPF1; global iPF2;
    global iDiv1; global iDiv2;
    global iSol;

    %Use FIESTA config file if supplied - default control is fixed
    if isa(InputConfig,'fiesta_configuration') == true;
        control = fiesta_control('diagnose',true, 'quiet',false, 'convergence',1e-4, 'boundary_method',2);
        config = InputConfig;
    else
        %Else, initiate efit configuration - default configuration and control are fixed
        control = fiesta_control('diagnose',true, 'quiet',false, 'convergence',1e-4, 'boundary_method',2);
        config = fiesta_configuration('SMART_config', Grid, coilset);
    end

    %Initiate icoil object and extract currents at desired time index
    if isa(InputVesselWaveforms,'double') & length(InputVesselWaveforms) == 0
        CoilCurrents = transpose(InputCoilWaveforms(:,TimeIndex));              %Extracts coil currents from desired time index
        icoil = fiesta_icoil(coilset, CoilCurrents);                            %Creates icoil object with coil currents
    else
        CoilCurrents = transpose(InputCoilWaveforms(:,TimeIndex));              %Extracts coil currents from desired time index
        icoil = fiesta_icoil(coilset, [CoilCurrents, InputVesselWaveforms]);    %Creates icoil object with coil & vessel currents
    end
    
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
    efitCoilWaveforms = InputCoilWaveforms;                   %Initiate output current waveforms
	icoil_efit = get(Equilibrium,'icoil'); 	CoilCurrents_efit = get(icoil_efit,'currents');
	efitCoilWaveforms(iPF1,5:6) = CoilCurrents_efit(iPF1);    %Assumes IPF1 is flat over equilibrium
	efitCoilWaveforms(iPF2,5:6) = CoilCurrents_efit(iPF2);    %Assumes IPF2 is flat over equilibrium
	efitCoilWaveforms(iDiv1,5:6) = CoilCurrents_efit(iDiv1);  %Assumes IDiv1 is flat over equilibrium
	efitCoilWaveforms(iDiv2,5:6) = CoilCurrents_efit(iDiv2);  %Assumes IDiv2 is flat over equilibrium
    
    %Clean up before returning to main code
    close all
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

    %Compute dynamic coil currents employing Current Driven Ip
	%CurlyM and CurlyR are large inductance and resistance matrices.
    [V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, Time_Adaptive ] = ...
        state_space_including_passive_elements_v4( CurlyM, CurlyR, Time_Linear, IPFinput_Continous, VPFinput_Continous, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names',coil_names, 'show_plot',false, 'turns',coilturns, 'currentScale',1e3, 'PF_colors',PF_colors );
    
    %Set breakdown time and prepare Ip_long and Vp_long for voltage driven Ip
    %!!! NEED Time_Breakdown TO INCLUDE BREAKDOWN AND BURNTHROUGH TIME !!!
    Time_Breakdown = 0;                                          %Set time for plasma breakdown (default 0)
    Time_Plasma = Time_Adaptive > Time_Breakdown;                %Set times for which plasma exists
    Vp_output(Time_Plasma) = 0;                                  %Set voltage to zero when plasma exists
    Vp_long = interp1(Time_Adaptive, Vp_output, Time_Linear);    %Sets Vp_long = 0 when Time_Linear > 0.
    Ip_long = NaN*Vp_long;                                       %Sets Ip_long to 'NaN' array (Voltage Driven)
    
    %Compute dynamic coil currents employing Voltage Driven Ip
    [ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, Time_Adaptive ] = ...
        state_space_including_passive_elements_v4( CurlyM, CurlyR, Time_Linear, IPFinput_Continous, VPFinput_Continous, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names',coil_names, 'show_plot',false, 'turns',coilturns, 'currentScale',1e3, 'PF_colors',PF_colors );

    %Clean up before returning to main code
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
    if isnan(CoilWaveformsInput(iDiv1,TimeIndex)) == true                   %If Div1 NOT in series Sol
        %Scale ALL null-field coil currents relative to Solenoid current
        D1=[D1_PF1, D1_PF2, D1_Div1, D1_Div2];                              %Optimise for PF1,2 & Div1,2
        ISolNullField = CoilWaveformsInput(iSol,TimeIndex);                 %Extract Solenoid Current 
        IPF_null = -pinv(D1) * (C1*ISolNullField);                          %Scale null-field currents
        IPF_null = [ISolNullField,transpose(IPF_null(:))];                  %Add Sol into IPF_null
    
    elseif isnan(CoilWaveformsInput(iDiv1,TimeIndex)) == false              %If Div1 IS in series with Sol
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
function SensorBTheta=InitiateBSensors(RGeo,ZGeo,R_null)

    %Determine number of sensors (constant for now)
    NumSensors = 10;

    %Define null field region, of size R_null^2, centred on equilibrium RGeo, ZGeo
    BP_virt_R = linspace(RGeo-R_null,RGeo+R_null,NumSensors);
    BP_virt_Z = linspace(ZGeo-R_null,ZGeo+R_null,NumSensors);

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

%Constructs a set of rectilinear vessel wall filaments of variable thickness
%Set FilamentArea to 0.0 to autocompute filament area
%Set WallNorm to 'Static' to enable variable filament area - not recommended
function [vessel_filaments,R_Fil_Array,Z_Fil_Array]=...
    CreateRectilinearVessel(VesselDimensions,WallThickness,FilamentArea,WallNorm)
    
    %Initiate any required data or control arrays
    VesselFaces = ["Horizontal", "Vertical", "Horizontal", "Vertical"];
    R_Fil_Array = [];  Z_Fil_Array = [];
    dR_Fil_Array = []; dZ_Fil_Array = [];
	
    %Unpack wall corners into single parameters
    RMinCentre = VesselDimensions(1); RMaxCentre = VesselDimensions(2); %[m]
    ZMinCentre = VesselDimensions(3); ZMaxCentre = VesselDimensions(4); %[m]
 	
    %Determine constant filament area if not otherwise specified
    if FilamentArea == 0.0  
        FilamentArea = max(WallThickness)^2;    %[m^2]
    end
	
    %Define four vessel vertices going clockwise from inboard upper (V1)
    Vertice1 = [RMinCentre ZMaxCentre];	%Top Left
    Vertice2 = [RMaxCentre ZMaxCentre];	%Top Right
    Vertice3 = [RMaxCentre ZMinCentre];	%Bottom Right
    Vertice4 = [RMinCentre ZMinCentre];	%Bottom Left
    WallCorners = [[Vertice1]; [Vertice2]; [Vertice3]; [Vertice4]];
 	
    %Define four wall edges going clockwise from inboard upper vertex (V1)
    %Top and bottom walls 'own' their vertices and are sized to match radial wall extents
    WallEdge1 = [RMinCentre+WallThickness(4)/2, RMaxCentre-WallThickness(3)/2];  %Upper Wall
    WallEdge2 = [ZMaxCentre, ZMinCentre];  %Outboard Wall
    WallEdge3 = [RMaxCentre-WallThickness(3)/2, RMinCentre+WallThickness(4)/2];  %Lower Wall
    WallEdge4 = [ZMinCentre, ZMaxCentre];  %Inboard Wall
    WallEdges = [[WallEdge1]; [WallEdge2]; [WallEdge3]; [WallEdge4]];
    
    %Define normalisation factors for each wall to achieve filament area
    WallNormFactor1 = FilamentArea/(WallThickness(1)^2);
    WallNormFactor2 = FilamentArea/(WallThickness(2)^2);
    WallNormFactor3 = FilamentArea/(WallThickness(3)^2);
    WallNormFactor4 = FilamentArea/(WallThickness(4)^2);
    WallNormFactors = [WallNormFactor1, WallNormFactor2, WallNormFactor3, WallNormFactor4];

    %If requested, set variable filament area (not recommended)
    if WallNorm == "Static"
        WallNormFactors = [1, 1, 1, 1];
    end
    
    %For each wall, create a linspace of filament coordinates in (R,Z)
    %and a corresponding linspace of scaled width and height values.
    for i=1:length(VesselFaces)
        
        %If vessel face is horizontal, scale width of filament to maintain area
        if VesselFaces(i) == "Horizontal"
            Height = WallThickness(i);                               %Define thickness parallel to wall direction
            Width = WallThickness(i)*WallNormFactors(i);             %Define thickness perpendicular to wall direction
            NumFil = (RMaxCentre-RMinCentre+2*WallThickness(i))/WallThickness(i);
            NumFil = floor( NumFil/WallNormFactors(i) );             %Scale number of filaments to maintain total width

            R_fil = linspace(WallEdges(i,1),WallEdges(i,2),NumFil);  %Create evenly spaced array of wall R coordinates
            Z_fil = WallCorners(i,2)*ones(1,NumFil);                 %Create evenly spaced array of wall Z coordinates
            dR_Wall = linspace(Height,Height,NumFil);                %Create Axial wall thickness array of size NumFil
            dZ_Wall = linspace(Width,Width,NumFil);                  %Create Radial wall thickness array of size NumFil

            %Top and bottom walls 'own' their vertices - i.e. full width of filaments
            R_Fil_Array = [R_Fil_Array, R_fil];                      %Append filament R coordinates to array
            Z_Fil_Array = [Z_Fil_Array, Z_fil];                      %Append filament Z coordinates to array
            dR_Fil_Array = [dR_Fil_Array, dR_Wall];                  %Append filament radial widths to array
            dZ_Fil_Array = [dZ_Fil_Array, dZ_Wall];                  %Append filament axial heights to array

        %If vessel face is vertical, scale height of filament to maintain area
        elseif VesselFaces(i) == "Vertical"
            Height = WallThickness(i)*WallNormFactors(i);            %Define thickness parallel to wall direction
            Width = WallThickness(i);                                %Define thickness perpendicular to wall direction
            NumFil = (ZMaxCentre-ZMinCentre+2*WallThickness(i))/WallThickness(i);
            NumFil = floor( NumFil/WallNormFactors(i) );             %Scale number of filaments to maintain total height

            Z_fil = linspace(WallEdges(i,1),WallEdges(i,2),NumFil);  %Create evenly spaced array of wall Z coordinates
            R_fil = WallCorners(i,1)*ones(1,NumFil);                 %Create evenly spaced array of wall R coordinates
            dR_Wall = linspace(Width,Width,NumFil);                  %Create Axial wall thickness array of size NumFil
            dZ_Wall = linspace(Height,Height,NumFil);                %Create Radial wall thickness array of size NumFil

            %Radial walls don't own the vertices - remove first and last filaments
            R_Fil_Array = [R_Fil_Array, R_fil(2:end-1)];             %Append filament R coordinates to array (removing corners)
            Z_Fil_Array = [Z_Fil_Array, Z_fil(2:end-1)];             %Append filament Z coordinates to array (removing corners)
            dR_Fil_Array = [dR_Fil_Array, dR_Wall(2:end-1)];         %Append filament radial widths to array (removing corners)
            dZ_Fil_Array = [dZ_Fil_Array, dZ_Wall(2:end-1)];         %Append filament axial heights to array (removing corners)
        end 
    end

    %Construct vessel wall employing passive FIESTA filaments using position arrays
    for i=1:length(R_Fil_Array)
        vessel_filaments(i) = fiesta_filament(R_Fil_Array(i),Z_Fil_Array(i),dR_Fil_Array(i),dZ_Fil_Array(i),1,0,0);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function is a little dodgy at the moment - SolTurnsR > 1 requires CoilWaveforms Modifications
function [Sol_Circuit]=...
    CreateSMARTSolenoidCircuit(SolName,RMaxSol,RMinSol,ZMaxSol,ZMinSol,SolTurnsZ,SolTurnsR,SolTemp,SolResistivity,SolDensity)

    %SolTurnsZ :: Solenoid axial windings (true solenoid turns)
    %SolTurnsR :: Solenoid radial windings (proxy for radial resolution)

    %Initiate Z-coordinates for solenoid filaments and calculate width/height
    Z_filaments = linspace(ZMinSol,ZMaxSol,SolTurnsZ); clear('coil_filaments');
    
    %Calculate solenoid filament width, height and radial spacing
    SolFilWidth = (RMaxSol-RMinSol)/SolTurnsR;    %Filament width  [m]
    SolFilHeight = (2*ZMaxSol)/SolTurnsZ;             %Filament height [m]
    SolFilRadii = linspace(RMinSol+SolFilWidth/2, RMaxSol-SolFilWidth/2, SolTurnsR );   %[m]
    
    %Construct central solenoid filaments - solenoid is treated as 'vessel wall' with nonzero current
    for i=1:length(SolFilRadii)
        %Save each vertical filament stack as a seperate coil - moving from inner radius outwards
        for iFilament=1:SolTurnsZ
            Sol_filaments(iFilament) = fiesta_filament( SolFilRadii(i), Z_filaments(iFilament), SolFilWidth, SolFilHeight ); 
        end
        Sol_Coil(i) = fiesta_coil( 'Sol_Zcoil', Sol_filaments, 'Blue', SolResistivity, SolDensity );
    end
    
    %Compile vertical coil stacks into a single circuit
    Sol_Circuit = fiesta_circuit( SolName, ones(1,length(SolFilRadii)), Sol_Coil(:) );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BpolAvg,BtorAvg]=ExtractNullBMin(EquilParams,BpolData,BtorData,R_null)

    %Obtain required global variables
    global GridRes_R; global GridRes_Z;
    global Grid;
    global BEarth;

    %Null-field region centre is at (RGeo, ZGeo) to align with sensor_btheta
    RGeo = EquilParams.r0_geom; ZGeo=EquilParams.z0_geom;
    
    %Convert from SI to cell index notation to enable averaging
    RAxis = get(Grid,'r'); ZAxis = get(Grid,'z');
    IndexRGeo = find(RGeo < RAxis); IndexRGeo = IndexRGeo(1);
    IndexZGeo = find(ZGeo < ZAxis); IndexZGeo = IndexZGeo(1);

    %Determine null-field index range; radius of R_null around (R, Z)
    RCells = ceil(R_null/GridRes_R)/2.0;     %Null-field region radius in cells
    ZCells = ceil(R_null/GridRes_Z)/2.0;     %Null-field region height in cells
    CellRangeR = [IndexRGeo-RCells, IndexRGeo+RCells];
    CellRangeZ = [IndexZGeo-ZCells, IndexZGeo+ZCells];
    
    %Resize Bpol and Btor arrays to null-region
    BpolData_NullRegion = BpolData(CellRangeZ(1):CellRangeZ(2), CellRangeR(1):CellRangeR(2));
    BtorData_NullRegion = BtorData(CellRangeZ(1):CellRangeZ(2), CellRangeR(1):CellRangeR(2));

    %Find min, max and average values within null-field region
    BpolMin = min(BpolData_NullRegion,[],'all'); BpolMax = max(BpolData_NullRegion,[],'all');
    BtorMin = min(BtorData_NullRegion,[],'all'); BtorMax = max(BtorData_NullRegion,[],'all');
    BpolAvg = mean(BpolData_NullRegion,'all');  %[T]
    BtorAvg = mean(BtorData_NullRegion,'all');  %[T]

    %Enforce lower limit for BpolMin (default to a multiple of Earth's B-field)
    BMin = BEarth; %[T]     (Default BEarth = 1e-5 T = 0.5 G)
    if BpolAvg < BMin
        BpolAvg = BpolAvg+BMin;	
    elseif BtorAvg < BMin
        BtorAvg = BtorAvg+BMin;
    end
    %Lloyd1991 suggests DIII-D Tokamak Bpolmin as 0.2 -> 1.2mT (2 -> 12G)
    %Song2017 suggests HL-2A Tokamak Bpolmin as 0.2 -> 1mT (2 -> 10G)
    %Chung2013a suggests VEST Tokamak Bpolmin as ???? 
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
function Lc=ConnectionLength(EquilParams,InnerVesselDimensions)

    %Extract inner vessel dimensions for further processing
    VesselRMax=InnerVesselDimensions(1); VesselRMin=InnerVesselDimensions(2);
    VesselZMax=InnerVesselDimensions(3); VesselZMin=InnerVesselDimensions(4);
    
    %Define Location of Inner Vessel Walls (four corners)
    RadialCorners = [VesselRMin, VesselRMax, VesselRMax, VesselRMin];
    AxialCorners = [VesselZMin, VesselZMin, VesselZMax, VesselZMax];
    
    %Define starting location for connection length (default RGeo,ZGeo) 
    RPoint = EquilParams.r0_geom;
    ZPoint = EquilParams.z0_geom;
    
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

%Induced loop voltage and E-field calculator 
%Assumes induced field arises purely from solenoid ramp-down
function [Vloop,Eloop,DeltaPhi]=LoopVoltage(CoilWaveforms,time,RSol,ZMaxSol,ZMinSol,RGeo)

    %Obtain required global variables
    global coilturns; global iSol;
    global mu0;
    
    %Compute the maximum average change in Sol current over the full discharge
    for i=2:length(time)
        dI = CoilWaveforms(iSol,i)-CoilWaveforms(iSol,i-1);
        dt = time(i)-time(i-1);
        dIdt_Array(i-1) = abs(dI/dt);       %Direction of loop current is arbitary
    end
    [dIdt,index] = max(dIdt_Array);                     %Maximum Sol current ramp [A/t]
    dt = time(index+1)-time(index);                     %Timescale of max current ramp [s]
    
    %Compute current surface area and voltage loop path length 
    IAreaLoop = pi*RSol^2;          %Current surface area set by outer sol radius
    VLengthLoop = 2*pi*RGeo;        %Voltage loop path length set by RGeo
    SolLength = ZMaxSol-ZMinSol;    %Solenoid length
    
    %Compute induced voltage during solenoid ramp-down and associated E-field at RGeo
    Vloop = IAreaLoop*((mu0*dIdt*coilturns(iSol))/SolLength);  %Loop voltage from solenoid [V] 
    Eloop = Vloop/VLengthLoop;                                 %E-field at plasma centre [V/m]

    %Compute average solenoid magnetic flux swing during ramp-down
    DeltaPhi = Vloop*dt;                                       %Solenoid flux swing [Vs]
    
    %NOTE :: MAXIMUM POSSIBLE SOLENOID MAGNETIC FLUX
    %(mu0*pi*ncoil*RSol^2)/(HeightSol*2*MaxISol)               %Assume one linear ramp [Vs]
    %(((mu0*pi*230*RSolCentreWinding^2)/(ZMaxSol*2))*2*12500)  %Sanity check ~ 0.252   [Vs]
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute breakdown (avalanche) timescale at given pressure
%Assumes pre-ionisation electron density and breakdown fraction
function [TauAvalanche,Pressure,Alpha,Vde]=AvalancheTimescale(Pressure,Eloop,Lc,ne,ne0,minimise)
   
    %%% Notes ::
    %   Lloyd1991, section 9.1, page 2046 has good citations for Paschen/Townsend Coefficients
    %   Lloyd1991, section 9.4, has equations for burnthrough time calculation
    
     %Compute background electron density assuming pre-ionisation
%    Tau_ECR = 0.005                                    %Pre-ionisation timescale          [s]
%    Iota_ECR = 8000.0                                  %Ionisation Rate (1/Tau_Ion)       [s-1]
%    Loss_ECR = 1.00                                    %Electron loss Rate (1/Tau_Loss)   [s-1]
%    ne0_ECR = ne0*exp((Iota_ECR-Loss_ECR)*Tau_ECR);    %ECR Electron Density [m-3]

    %Define fixed or estimated values (these may need explicit calculation)
    Eta = 43;                   %?????                              [?]   Estimation [Lloyd1991]
%   ne0 = 1.0;                  %Background electron density        [m-3] Estimation [Lloyd1991] ~1.0 with no pre-ionisation
    BDFraction = 0.15;          %Breakdown fraction @ Coulumb phase [%]   [Lloyd1991]
    nbd = BDFraction*ne;        %Electron density post-breakdown    [m-3] [Lloyd1991]
    neFraction = nbd/ne0;       %Relative breakdown density         [%]   [Lloyd1991]    

    %Compute expected pressure and avalanche timescale for breakdown
    Vde = (Eta*Eloop)/Pressure;                             %Electron Drift Speed    [m/s]
    Alpha = 510*Pressure*exp((-1.25e4*Pressure)/Eloop);     %First Townsend coeff    [m-1]
    TauAvalanche = log(neFraction)/(Vde*(Alpha-(1/Lc)) );   %Breakdown time          [s]

    %Compute minimum pressure and avalanche timescale for breakdown
    if minimise == true
        PressureLimits = [1e-6, 1e-3]; Resolution = 25000;
        Pressure_Array = linspace(PressureLimits(1),PressureLimits(2),Resolution);

        for i=1:length(Pressure_Array)
            Alpha_Array(i) = 510*Pressure_Array(i)*exp((-1.25e4*Pressure_Array(i))/Eloop);
            Vde_Array(i) = (Eta*Eloop)/Pressure_Array(i);
            TauAvalanche_Array(i) = log(neFraction)/(Vde_Array(i)*(Alpha_Array(i)-(1/Lc)) );
        end
        %Take minimum values of Vde, Pres, Tau :: corresponding to the maximum townsend coefficient
        Alpha = max(Alpha_Array);                                 %Townsend Coeff          [m-1]
        Vde = Vde_Array( Alpha_Array==Alpha );                    %Electron Drift Speed    [m/s]
        Pressure = Pressure_Array( Alpha_Array==Alpha );          %Townsend coeff          [m-1]
        TauAvalanche = TauAvalanche_Array( Alpha_Array==Alpha );  %Breakdown time          [s]
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Write 2D, 1D and 0D equilibria out to files
function [fileID]=WriteEquilibrium(Equilibrium,config,EquilDir,EquilName,VacuumField)

    %Get any required global variables
    global Grid;
    %Extract Grid arrays for reshaping if required
    rGrid = get(Grid,'r'); zGrid = get(Grid,'z');
    
    %If Equil contains plasma then use in-built geqdsk function
    if VacuumField == false
        %Write 2D qeqdsk equilibrium file
        Filename = strcat(EquilDir,'Equil.txt');
        geqdsk_write_BUXTON(config, Equilibrium, Filename);

        %Write 1D equilibrium qprofile parameters file
        qProfile = qprofile(Equilibrium);
        qProfileVariables = fieldnames(qProfile);
        qProfileParams = struct2cell(qProfile(:,1));
        Filename = strcat(EquilDir,'EquilProfiles',EquilName,'.txt');
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
        EquilibriumParams = parameters(Equilibrium);
        ParamVariables = fieldnames(EquilibriumParams);
        ParamValues = struct2cell(EquilibriumParams(:,1));
        Filename = strcat(EquilDir,'EquilParam',EquilName,'.txt');
        fileID=fopen(Filename,'w');
        for i = 1:length(ParamValues)
            try fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); ParamValues(i)]);
            catch fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); 'NaN']);
            end
        end
    
    %If vacuum field equilibrium is supplied, cannot use geqdsk function
    elseif VacuumField == true 
        %Write 2D vacuum equilibrium file
        Filename = strcat(EquilDir,'Null_Equil.txt');
        fileID = fopen(Filename,'w');
        
        %Vacuum equilibria can't use geqdsk format - save as 2D array
        Psi_RZ = struct2cell(get(Equilibrium,'Psi_vac')); 
        Psi_RZ = Psi_RZ(3); Psi_RZ = Psi_RZ{1,1};                   %len(50200) = len(R)*len(Z)
        Psi_RZ = reshape(Psi_RZ,[length(zGrid),length(rGrid)]);     %[rGrid,zGrid] = [len(R),len(Z)]
        fprintf(fileID,'%s %s %s \n', '     Psi_Null', string(length(rGrid)), string(length(zGrid)));
        for i = 1:size(Psi_RZ,1)
            fprintf(fileID,'%g\t',Psi_RZ(i,:));
            fprintf(fileID,'\n');
        end
        
        %Extract the null poloidal and toroidal B-field vector arrays
        [BrData_Null,BzData_Null,BPhiData_Null,BpolData_Null,BtorData_Null] = ...
            ExtractBField(Equilibrium);
        
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
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fileID]=WriteGeometry(Geometry,EquilDir,Filename)

    %Extract geometry for easier reading
    RGeo = Geometry(1);
    ZGeo = Geometry(2);
    rGeo = Geometry(3);
    kappa = Geometry(4);
    delta = Geometry(5);

    %Write initial target geometry, efit geometry and perturbed geometry
    Filename = strcat(EquilDir,Filename);
    fileID = fopen(Filename,'w');
    fprintf(fileID,'%s %s %s %s %s\r\n', 'RGeo','ZGeo','a','kappa','delta');
    fprintf(fileID,'%1.12f %1.12f %1.12f %1.12f %1.12f\r\n',[RGeo'; ZGeo'; rGeo'; kappa'; delta']);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
