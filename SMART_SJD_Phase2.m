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
FIESTASource = "Source/FIESTA V8.8/";
addpath(genpath(FIESTASource));

%Set maximum CPU thread usage
NumThreads = 2;
NumThreads = maxNumCompThreads(NumThreads);

%%%%%%%%%%%%%%%%%%%%%%  DATA OUTPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define figure colourmap and extension
global colourmap; colourmap = Plasma();     %'Plasma()','Gamma_II()'
FigExt = '.png';                            %'.png','.eps','.pdf'

%Define simulation shot name
ShotName = 'S2-000018';		%Define shot name: typically Sx-xxxxxx"

%Create global output folders for saved data and figures
SimDir = strcat(ShotName,'/'); mkdir(SimDir);
ASCIIDir = strcat(SimDir,'RawData/'); mkdir(ASCIIDir);

%Copy any relevant files into simulation folder before startup
copyfile('SMART_SJD_Phase2.m',SimDir)

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

%Define Vessel Internal Geometry - Do not include wall thickness except for inner wall
%Radius R is defined relative to centre of solenoid, Height Z is defined relative to midplane.
VesselRMinInner=+0.15+VWall_Inboard;	% Vessel Minimum Internal Radius R [m] 
VesselRMaxInner=+0.80;					% Vessel Maximum Internal Radius R [m]
VesselZMinInner=-0.80;					% Vessel Minimum Internal Height Z [m]
VesselZMaxInner=+0.80;					% Vessel Maximum Internal Height Z [m]

%Define center points of vessel walls (Inner Geometry + half wall thickness)
ZMinCentre=VesselZMinInner-(VWall_Lower/2);		% Lower Wall 'grows' outwards (-Z direction)
ZMaxCentre=VesselZMaxInner+(VWall_Upper/2);		% Upper Wall 'grows' outwards (+Z direction)
RMinCentre=VesselRMinInner-(VWall_Inboard/2);	% Inboard wall 'grows' inwards (-R direction)
RMaxCentre=VesselRMaxInner+(VWall_Outboard/2);	% Outboard wall 'grows outwards (+R direction)

%%%%%%%%%%%%%%%%%%%%%%%  DEFINE COIL GEOMETRY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Solenoid Geometry and Parameters
nSol = 230;                                 % Number of Axial Solenoid Windings [-]
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

%Define global constants
global epsilon0; epsilon0 = 8.854E-12;      % Vacuum Permissability       [F m^-1]
global mu0; mu0 = 1.2566e-06;               % Vacuum Permeability         [N/A^2] or [H/m]
global e; e = 1.906e-19;                    % Fundamental Charge          [C]
global mass_e; mass_e = 9.109E-31;          % Electron Mass               [kg]
global kB; kB = 1.3806E-23;                 % Boltzmann's Constant        [m^2 kg s^-2 K-1]
global eV_K; eV_K = 1/8.621738E-5;          % eV to Kelvin Conversion     [K/eV]
global BEarth; BEarth = 1.0E-4;             % Earth's B-Field (Max limit) [T]


%Define initial operating conditions (primarily used for Topeol2)
Te = 250;			% Electron Temperature [eV]
Ti = Te*0.30;		% Ion Temperature      [eV]
BT = +0.3;			% Toroidal B-Field     [T] (Defined at Rgeo)
Ip = +100e3;		% Plasma current       [A]
RGeo = 0.420;		% Geometrical Radius   [m] (~0.420 --> 0.480)
ZGeo = 0.000;		% Geometrical Axis     [m] (=0.000)
RSep = 0.700;		% Separatrix Radius    [m] (=0.700)
rGeo = RSep-RGeo;	% Minor Radius         [m] (=0.250)
Aspect = RGeo/rGeo;	% Aspect ratio         [-] (~1.50 --> 1.85)
Kappa = 1.80;		% Elongation           [-] (~1.70 --> 2.00)
delta = 0.20;		% Triangularity        [-] (~0.20)
Z_eff = 2.0;        % Effective Charge     [C] (~1.44 - 2.00)
li2 = 1;			% Inductance	       [-] (??????)

%Compute further operating conditions (primarily used for Topeol2)
Gr_Frac = 0.40;                            % Greenwald Fraction       [-]
Gr_Limit = 1e20*(Ip*1e-6/(pi*rGeo^2));     % Greenwald Limit          [m-3]
ne = abs(Gr_Limit*Gr_Frac);                % Electron Density         [m-3]
Irod = (BT*2*pi*RGeo)/mu0;                 % Central Rod Current      [A]
S = sqrt( (1.0+Kappa^2)/2.0 );             % Shaping factor           [-]
%deltaUp = (ZGeo-Zup)/a;                   % Upper-Triangularity      [-]
%deltaLo = (ZGeo-Zlo)/a;                   % Lower-Triangularity      [-]
%delta = (deltaUp+deltaLo)/2.0;            % Triangularity            [-]
%betaN = (betaT*BT*a)/(Ip*1e-6*mu0)        % Normalised Beta          [%] 
%betaT = (betaN/a*(Ip*1e-6))/BT;           % Beta toroidal            [%]
betaP = 3/2*ne*(Te+Ti)/(mu0*Ip/(2*pi*rGeo))^2*2*mu0*1.6e-19*Kappa;  % Beta Poloidal  [%]
BZ = -mu0*Ip/(4*pi*RGeo)*(log(8*Aspect)+betaP+0.5*li2-3/2);         % Vertical field [T]

%Define efit Equilibrium Operating Conditions (primarily used for efit)
RGeo_efit = 0.420;					% Geometric Radius      [m] (0.420 --> 0.480) ::
ZGeo_efit = 0.000;					% Geometric Height      [m] (Default 0.000)   ::
Aspect_efit = 1.85;                 % Aspect Ratio          [-] (1.850 --> 2.000) :: RGeo/rGeo
rGeo_efit = RGeo_efit/Aspect_efit;  % Minor Radius	        [m] (Default 0.238)   :: RGeo/Aspect
Kappa_efit = 2.00;					% Elongation			[-] (+1.70 -> +2.00)  :: (Zmax-Zmin)/2rGeo
delta_efit = 0.20;					% Triangularity			[-] (-1.00 -> +1.00)  :: (Zmax-Zgeo)/rGeo (max/min)
efitGeometry_Init = [RGeo_efit, ZGeo_efit, rGeo_efit, Kappa_efit, delta_efit];

%Define plasma stability initial perturbations (primarily used for Feedback control)
deltaRGeo = 0.00;	% Small radial perturbation         [m]
deltaZGeo = 0.01;	% Small axial perturbation          [m]
deltaAspect = 0.00;	% Small aspect ratio perturbation   [-]
deltaKappa = 0.00;	% Small elongation perturbation     [-]
deltadelta = 0.00;	% Small triangiularity perturbation [-]
PertGeometry_Init = [deltaRGeo,deltaZGeo,deltaAspect,deltaKappa,deltadelta];

%Compute parallel and perpendicular plasma resistivities employing Spitzer model with impurities
%Typical values: H=1, He=2, Ar=11.85 (Te < 280eV) https://www.webelements.com/argon/atoms.html
%H discharge: Z_eff = 2 allowing for Carbon impurities in the plasma (Wall tiles)
[EtaPerp,EtaPara] = SpitzerResistivity(ne,Te,Z_eff);    %[Ohm m^-1]
PlasmaResistPerp = EtaPerp*(2*pi*RGeo);                 %[Ohm]
PlasmaResistPara = EtaPara*(2*pi*RGeo);                 %[Ohm]

%Define Coil density, temperature, and resistivity
CoilDensity = 1;                       % Relative Coil Density      [Arb]
CoilTemp = 293.0;                      % Initial Coil Temperature   [K]
CoilResistivity = InterpMaterialResistivity(CoilTemp);

%Define Vessel density and resistivity  (Stainless Steel (AISI 316 L) Vessel)
VesselResistivity = 6.9e-7;             % Absolute Vessel Resistivity   [Ohm m^-1]
VesselDensity = 7.8e3;                  % Absolute Vessel Density       [Kg/m3]

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
                                    %TauR1=15ms  %TauR1=15ms    %TauR1=15ms
                                    %RGeo=0.42   %RGeo=0.42     %RGeo=0.475
%Solenoid coil currents [A]         %Phase2      %Phase2-d      %Phase2+d
I_Sol_Null=+4000;					%+4000;      %+4500;        %+5000;
I_Sol_MidRamp=+0000;				%+0000;      %+0000;        %+0000;
I_Sol_Equil=-0600;                  %-0600;      %-0900;        %-1500;
I_Sol_EndEquil=-2200;           	%-2200;      %-2400;        %-3000;

%PF & Div Equilibrium coil currents [A]         (Default equilibrium: time(4,5,6))
I_PF1_Equil=-1100;					%-1100;      %-1100;        %-1100;
I_PF2_Equil=-1100;					%-1100;      %-1100;        %-1100;
I_Div1_Equil=+1000;					%+1000;      %-3500;        %+2500;     (Pos for +d, Neg for -d)
I_Div2_Equil=+0000;                 %+0000;      %+0000;        %+0000;

%Define TimeIndices (vertices) in the Sol, PF & Div coil current waveforms
TauN  = 0.020;			% Null-Field Timescale      [s] Determines null-field decay timescale
TauR1 = 0.015;			% Breakdown Ramp Timescale  [s] Determines max loop voltage
TauR2 = 0.020;			% PF & Div Ramp Timescale   [s] Determines max PF/Div current ramp
TauR  = TauR1+TauR2;    % Total Ramp Timescale      [s] 
TauP  = 0.100;			% Pulse Timescale      		[s] Determines flat-top timescale

%Create time array, containing Sol, PF & Div coil current waveform time vertices
%Time   [Init      PrePulse   InitRampDown  MidRampDown  EndRampDown  MidEquil     Terminate         ];
time =  [-2*TauN   -TauN      0.0           TauR1        TauR         TauR+TauP    TauR+TauP+(2*TauN)];
nTime = length(time);	% Total Coil Waveform Timesteps	[-]

%Fits linear midpoint to any current defined as 'linear' between times: {pre-ramp, mid-ramp, end-ramp}
I_Sol_MidRamp = FitSolenoidRamp({I_Sol_Null,I_Sol_MidRamp,I_Sol_Equil},time);

%Construct Sol, PF/Div coil current waveforms by specifying current at temporal vertices
%Entries containing NaN will be calculated via RZIp to minimise BpolNull during null-field and breakdown
%					              %!Null-Field! %!Breakdown!   %!Efit Icoil!
%Time   	     [1,  2,          3,            4,             5,             6,               7];
ISol_Waveform =  [0,  I_Sol_Null, I_Sol_Null,   I_Sol_MidRamp, I_Sol_Equil,   I_Sol_EndEquil,  0];
IPF1_Waveform =  [0,  NaN,        NaN,          NaN,           I_PF1_Equil,   I_PF1_Equil,     0];
IPF2_Waveform =  [0,  NaN,        NaN,          NaN,           I_PF2_Equil,   I_PF2_Equil,     0];
IDiv1_Waveform = [0,  NaN,        NaN,          NaN,           I_Div1_Equil,  I_Div1_Equil,    0];
IDiv2_Waveform = [0,  NaN,        NaN,          NaN,           I_Div2_Equil,  I_Div2_Equil,    0];
%%%%%
%CoilWaveforms has structure: [CoilNumber][TimeIndex] - both being integers
CoilWaveforms = [ISol_Waveform; IPF1_Waveform; IPF2_Waveform; IDiv1_Waveform; IDiv2_Waveform];

%Define dynamic coils (i.e. which coil currents are fit by efit)
global efitCoils; efitCoils = {'PF1','PF2'};                        % Default PF1, PF2
global feedbackCoils; feedbackCoils = {'Div2'};                     % Default Div2

%Terminal Outputs for sanity checking
disp([ ' ' ]);
disp([ '%===== Discharge Parameters =====%' ]);
disp([ 'IRod = ' num2str(Irod) ' [kA]' ]);
disp([ 'BT = ' num2str(BT) ' [T]' ]);
disp([ 'BZ = ' num2str(BZ) ' [T]' ]);
disp([ 'ne = ' num2str(ne) ' [m-3]' ]);
disp([ 'Te = ' num2str(Te) ' [eV]' ]);
disp([ 'Ti = ' num2str(Ti) ' [eV]' ]);
disp([ 'Ip = ' num2str(Ip/1000) ' [kA]' ]);
disp([ 'TauP = ' num2str(TauP*1000) ' [ms]' ]);
disp([ ' ' ]);
disp([ '%===== Plasma Shaping =====%' ]);
disp([ 'RGeo = ' num2str(RGeo_efit) ' [m]' ]);
disp([ 'ZGeo = ' num2str(ZGeo_efit) ' [m]' ]);
disp([ 'Minor Radius = ' num2str(rGeo_efit) ' [m]' ]);
disp([ 'AspectRatio = ' num2str(Aspect_efit) ' [-]' ]);
disp([ 'Elongation = ' num2str(Kappa_efit) ' [-]' ]);
disp([ 'Triangularity = ' num2str(delta_efit) ' [-]' ]);
disp([ ' ' ]);
disp([ '%===== Coil Currents =====%' ]);
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

%Define vessel corners, thickness and vessel filament cross-sectional area 
VesselDimensions = [RMinCentre, RMaxCentre, ZMinCentre, ZMaxCentre];       %[m]
WallThickness = [VWall_Upper, VWall_Outboard, VWall_Lower, VWall_Inboard]; %[m]
%Lower filament areas give higher passive current resolution (Note :: FilamentArea Greatly affects convergence)
FilamentArea = 1.50e-4; %(2.5e-4 > A > 1.5e-4 or RZIp M,R matrices fail)   %[m^2]

%Construct SMART vessel wall filaments (false = variable fil area, true = normalised (fixed) fil area)
[vessel_filament,R_Fil_Array,Z_Fil_Array] = ... 
    CreateRectilinearVessel(VesselDimensions,WallThickness,FilamentArea,true);

%Construct passive vessel components and arrange into a vessel object
global passive; passive = fiesta_passive('STVesselPas',vessel_filament,'g',VesselResistivity,VesselDensity);
global vessel; vessel = fiesta_vessel( 'STVessel',passive);

%Compute characteristic magnetic field penetration timescale (Amoskov2005)
TauVessel = (mu0*max(WallThickness)^2)/VesselResistivity;       %[s]

%%%%%%%%%%%%%%%%%%%  INITIATE SOL, PF & DIV COILS  %%%%%%%%%%%%%%%%%%%%%%%%

%Define Sol, PF and Div coil sets - Assigning integers to coil structure names.
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

%Create coil set from parameters defined above. (Function made by Carlos Soria)
%Function createVESTPFCircuit creates two PF coils. One in (R, Z) and another in (R, -Z)
Sol = CreateSMARTSolenoidCircuit('Sol',RSolOuter,RSolInner,ZMaxSol,ZMinSol,coilturns(iSol),nSolR,CoilTemp,CoilResistivity,CoilDensity);
PF1  = CreateSMARTCoilCircuit('PF1',R_PF1,Z_PF1,width_PF1,height_PF1,coilturns(iPF1),nZPF1,nRPF1,CoilTemp,CoilResistivity,CoilDensity,true);
PF2  = CreateSMARTCoilCircuit('PF2',R_PF2,Z_PF2,width_PF2,height_PF2,coilturns(iPF2),nZPF2,nRPF2,CoilTemp,CoilResistivity,CoilDensity,true);
Div1 = CreateSMARTCoilCircuit('Div1',R_Div1,Z_Div1,width_Div1,height_Div1,coilturns(iDiv1),nZDiv1,nRDiv1,CoilTemp,CoilResistivity,CoilDensity,true); 
Div2 = CreateSMARTCoilCircuit('Div2',R_Div2,Z_Div2,width_Div2,height_Div2,coilturns(iDiv2),nZDiv2,nRDiv2,CoilTemp,CoilResistivity,CoilDensity,true);

%Collate global coilset containing Solenoid, PF and Div coil circuits (efit expects a row aligned filament array)
R_Fil_Array = transpose(R_Fil_Array); Z_Fil_Array = transpose(Z_Fil_Array);     
global coilset; coilset = fiesta_coilset('SMARTcoilset',[Sol,PF1,PF2,Div1,Div2],false,R_Fil_Array',Z_Fil_Array');
coilset_init = coilset;

%Divide solenoid current by number of radial windings if nested radial solenoid windings are employed
CoilWaveforms(1,:) = CoilWaveforms(1,:)/nSolR;
% ISSUE :: This doesn't appear to be required for the PF,Div coils, why here?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  END SIMULATION INITIAL SET-UP  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TO DO %%%

%%%   FINISH COMMENTS ON ALL FUNCTIONS
%%%   WRITE MANUAL

%%%   TOGGLEABLE UPPER AND LOWER SINGLE NULL (USN, LSN) CONFIGURATION IN CreateSMARTCoilCircuit
%%%   FEEDBACK SYSTEM WORKING FOR VERTICAL AND HORIZONTAL STABILITY
%%%   FIX THE CORNER OF THE DIFF VESSEL WALLS (LAST FILAMENT IS LARGER)

%%%   GET RZIP ABLE TO TAKE BOTH COIL AND VESSEL FILAMENTS

%%%   FIX THE EVIL TWIN LCFS PROBLEM, FIRST FOR COILS AND THEN FOR THE SOLENOID
%%%   findboundary.m function contains rules for LCFS boundary

%%%   IMPLIMENT THESE OUTPUTS
        %fieldnames(Equil)
        %Current =  reshape(get(get(Equil,'I'),'data'),GridCells_Z,GridCells_R);
        %CurrentDensity =  reshape(get(get(Equil,'J'),'data'),GridCells_Z,GridCells_R);
        %PlotEquilibrium({Current},'Title','CbarLabel','Current.png')
        %PlotEquilibrium({CurrentDensity},'Title','CbarLabel','CurrentDensity.png')

        
        
        
        
        
        
        


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           COMPUTE INITIAL TARGET AND NULL-FIELD EQUILIBRIUA             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TimeArray index for target equilibrium 
TimeIndex_Discharge = 5;                          %default time(5)  End of Sol Ramp
TimeIndex_NullField = 3;                          %default time(3)  Prior to Sol Ramp

%Create initial icoil object at requested TimeIndex - icoil represents coil currents fed to efit.
CoilWaveforms_Init = CoilWaveforms;
global icoil_init; icoil_init = fiesta_icoil(coilset);
%Assign equilibrium coil currents to icoil object [kA]
icoil_init.Sol=CoilWaveforms_Init(iSol,TimeIndex_Discharge);   %Solenoid Equilibrium Current
icoil_init.PF1=CoilWaveforms_Init(iPF1,TimeIndex_Discharge);   %PF1 Equilibrium Current
icoil_init.PF2=CoilWaveforms_Init(iPF2,TimeIndex_Discharge);   %PF2 Equilibrium Current
icoil_init.Div1=CoilWaveforms_Init(iDiv1,TimeIndex_Discharge); %Div1 Equilibrium Current
icoil_init.Div2=CoilWaveforms_Init(iDiv2,TimeIndex_Discharge); %Div2 Equilibrium Current


%%%%%%%%%%%%%%%%%%%%  COMPUTE DISCHARGE EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%

%Compute Jprofile from betaP and Ip employing Topeol type2 model - linear model.
jprofile = fiesta_jprofile_topeol2( 'Topeol2', betaP, 1, li2, Ip );

%Compute equilibrium (Psi(R,Z)) from the supplied jprofile, icoil and geometry
%Returns target equilibrium and CoilWaveforms for PF1 and PF2 at requested time_Index
global config; [Equil,EquilParams,CoilWaveforms,efitGeometry,config] = ...
    efit(jprofile,Irod,'config',efitGeometry_Init,CoilWaveforms_Init,[],TimeIndex_Discharge);

%Save discharge coil currents for all coils at TimeIndex_Discharge
CoilCurrentsEfit = transpose(CoilWaveforms(:,TimeIndex_Discharge));
icoil_efit = fiesta_icoil(coilset, CoilCurrentsEfit);

%%%%%%%%%%%%%%%%%%%%%  COMPUTE OPTIMISED NULL-FIELD  %%%%%%%%%%%%%%%%%%%%%%

%Initiate virtual B-field sensors centered on Rgeo
sensor_btheta = InitiateBSensors(EquilParams.r0_geom,EquilParams.z0_geom,R_Null);

%RZIP computes coefficients [A B C D] using the null field sensors
%Output C is used to compute the null-field PF coil currents 
%Outputs curlyM and curlyR are used to compute the plasma and eddy currents
rzip_config = fiesta_rzip_configuration( 'RZIP', config, vessel, {sensor_btheta} );
[RZIp_A, RZIp_B, RZIp_C, RZIp_D, curlyM, curlyR, Gamma, plasma_parameters, index, label_index, state] = ...
    response(rzip_config, Equil, 'rp', PlasmaResistPerp);   

%Update CoilWaveforms array with null-field values (Using NaN Mask)
%Replaces any "NaN" in CoilWaveforms with the null field value
CoilWaveforms = NullFieldWaveforms(CoilWaveforms,RZIp_C,sensor_btheta,TimeIndex_NullField);

%Save null-field coil currents for all coils at TimeIndex_NullField
CoilCurrentsNull = transpose(CoilWaveforms(:,TimeIndex_NullField));
icoil_Null = fiesta_icoil(coilset, CoilCurrentsNull);

%Compute null-field equilibrium using null-field coil currents
Equil_Null = fiesta_equilibrium('SMART-Null', config, Irod, icoil_Null );
EquilParams_Null = parameters(Equil_Null);


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

%Extract Vessel Eddy Currents during discharge and null-field (time='false' for absolute max)
VesselEddyCurrents = ExtractPassiveCurrents(I_Passive,time_adaptive,time(TimeIndex_Discharge));
VesselEddyCurrentsNull = ExtractPassiveCurrents(I_Passive,time_adaptive,time(TimeIndex_NullField));
% ISSUE :: EDDY CURRENTS CHANGE OVER THE AVALANCHE TIMESCALE, TimeIndex_NullField IS SIMPLIFICATION
%       :: NEED A METHOD TO EXAMINE NULL-FIELD EVOLUTION OVER AVALANCHE TIMESCALE 

%Compute vessel eddy current stress magnitude and direction for each vessel filament
[StressR,StressZ,StressR_max,StressZ_max] = ...
    VesselStresses(Equil,VesselEddyCurrents,FilamentArea);


%%%%%%%%%%%%%%%%%%%%%  COMPUTE BREAKDOWN CRITERIA  %%%%%%%%%%%%%%%%%%%%%%%

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
[Vloop,Eloop,DeltaPhiSol,VloopArray,EloopArray,DeltaPhiSolArray]=...
    LoopVoltage(I_PF_output,time_adaptive,RSolCentreWinding,ZMaxSol,ZMinSol,EquilParams.rin);
Eloop_eff = abs(Eloop)*(BtorAvg_Null/BpolAvg_Null); %[V/m]   %Rough estimate threshold Eloop condition for breakdown
%Generally Eloop_eff > 100 [V/m] for startup with ECRH      (An2015)
%Generally Eloop_eff > 1000 [V/m] for solenoid only startup (Lloyd1991)

%Compute breakdown (avalanche) timescale, TauBD, at given pressure
Pressure = EquilParams.P0*(7.5e-7);                          %[Torr]  (typically ~2e-4 Torr)
[TauBD,Pressure,Alpha,Vde]=AvalancheTimescale(Pressure,Eloop,Lc,ne,1.0,true);

%Compute solenoid OH flux swing employing the Ejima-Wesley coefficient, Cew
Cew = 0.40*EquilParams.aspectratio;                           %[-]      (Gryaznevich2006)
DeltaPhiSolew = Cew*mu0*EquilParams.r0_geom*max(Ip_output);   %[Vs]     (Menard2016 pg36)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           PLOT ITER(0) EQUIL, COIL, EDDY AND STRESS FIGURES             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% PLOT VESSEL AND COIL FILAMENT OVERVIEW %%%%%%%%%%%%%%%%%%%% 

%Plot vessel and coil filaments overview
Filename = '_VesselFilaments';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotVesselOverview(SaveString);

%%%%%%%%%%%%%%%%%%%%%%%% PLOT TARGET EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%%%

%Plot target equilibrium following convergence
Title = {'SMART Target Equilibrium iter(0)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_Equilibrium_00';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotEquilibrium({Equil},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%  PLOT VIRTUAL SENSORS ONTO EQUILIBRIUM  %%%%%%%%%%%%%%%%%

Title = {'SMART Virtual Sensors',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_VirtualBSensors';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotEquilibrium({Equil,sensor_btheta},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%  PLOT NULL-FIELD PHI SURFACES  %%%%%%%%%%%%%%%%%%%%%

%Plot the optimised null-field phi
Title = {'SMART Null-field Equilibrium iter(0)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_NullField_00';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotEquilibrium({Equil_Null},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%% PLOT NULL-FIELD BPOL  %%%%%%%%%%%%%%%%%%%%%

%Log poloidal and toroidal magnetic fields to show details (Sol Obscures)
logBpolData_Null = log10(BpolData_Null);
logBtorData_Null = log10(BtorData_Null);
%Plot the optimised null-field phi
Title = {'SMART Null-field iter(0)',' '};
CbarLabel = 'Null-field B_{\theta} log_{10}([T])';
Filename = '_NullBpol_00';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotEquilibrium({logBpolData_Null},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%% PLOT COIL CURRENT WAVEFORMS %%%%%%%%%%%%%%%%%%%%%%%% 

%Plot figure showing dynamic coil currents
figure('units','inch','position',[10 10 12 12]);
subplot(2,1,1); hold on; grid on; box on;
plot(time_adaptive*1000, I_PF_output/1000, 'LineWidth',2);
title(gca,'SMART Coil Current Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Coil Current I [kA]');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
%%%%%
subplot(2,1,2); hold on; grid on; box on;
plot(time_adaptive(1:end-1)*1000,Delta_IPFoutput/1000, 'LineWidth',2)
title(gca,'SMART Delta Coil Current Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Delta Coil Current \Delta I (kA ms^{-1})');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
%%%%%
Filename = '_CurrentWaveforms';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%% PLOT COIL VOLTAGE WAVEFORMS %%%%%%%%%%%%%%%%%%%%%%%% 

%Plot figure showing dynamic coil currents
figure('units','inch','position',[10 10 12 12]);
subplot(2,1,1); hold on; grid on; box on;
plot(time_adaptive*1000, V_PF_output/1000, 'LineWidth',2);
title(gca,'SMART Coil Voltage Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Coil Voltage V [kV]');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
%%%%%
subplot(2,1,2); hold on; grid on; box on;
plot(time_adaptive(1:end-1)*1000,Delta_VPFoutput/1000, 'LineWidth',2)
title(gca,'SMART Delta Coil Voltage Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Delta Coil Voltage \Delta V (kV ms^{-1})');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
%%%%%
Filename = '_VoltageWaveforms';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%% PLOT PLASMA CURRENT %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot plasma current over full timescale
figure('units','inch','position',[10 10 12 12]);
subplot(2,1,1); hold on; grid on; box on;
plot(time_adaptive*1000, Ip_output/1000, 'LineWidth',2)
title(gca,'SMART Plasma Current iter(0)');
legend(gca,'I_{p}', 'FontSize',16); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Plasma Current I_{p} (kA)');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
%%%%%
subplot(2,1,2); hold on; grid on; box on;
plot(time_adaptive(1:end-1)*1000, Delta_Ip_output/1000, 'LineWidth',2)
title(gca,'SMART Plasma Current iter(0)');
legend(gca,'dI_{p}/dt', 'FontSize',16); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Delta Plasma Current (kA/ms)');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
%set(gca,'XLim',[time(3)*1e3, time(5)*1e3]);        %Breakdown Ramp Closeup
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
%%%%%
Filename = '_PlasmaCurrent';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%% PLOT TOTAL EDDY CURRENT %%%%%%%%%%%%%%%%%%%%%%%%% 

%Sum all filaments (row-wise) to get total net passive current
Net_IPassive = sum(I_Passive,2);
%Plot net passive current density over full timescale
close all
figure('units','inch','position',[12 12 8 8]); hold on; grid on; box on;
plot(time_adaptive*1000, Net_IPassive/1000, 'LineWidth',2)
title(gca,'Net SMART Eddy Current iter(0)');
legend(gca,'Net Eddy Current'); legend boxoff;
xlabel(gca,'Time (ms)');
ylabel(gca,'Net Vessel Current I_{Eddy} [kA]');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca,'YLim',[min(Net_IPassive/1000)*1.15 max(Net_IPassive/1000)*1.20]);
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
Filename = '_1DEddyCurrent';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));

%%%%%%%%%%%%%%%%%% PLOT 2D RESOLVED EDDY CURRENTS %%%%%%%%%%%%%%%%%%%%%%   

%Obtain all passive vessel filaments
Passives = get(vessel,'passives');
Filaments = get(Passives,'filaments');
RFil = get(Filaments(:),'r');             %Radial Filaments
ZFil = get(Filaments(:),'z');             %Axial Filaments
%Plot equilibrium eddy currents within a cross-section of the vessel
close all
figure; hold on; grid on; box on; axis equal;
plot(coilset);
scatter3(RFil,ZFil,VesselEddyCurrents/1000,100,VesselEddyCurrents/1000,'filled');
title('SMART Vessel Eddy Currents iter(0)');
view(2) %2D view
colormap(colourmap);
cbar = colorbar;
cbar.Label.String = 'Eddy-Current I_{eddy} [kA]';
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_EddyCurrent';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));   

%%%%%     %%%%%     %%%%%     %%%%%     %%%%%     %%%%%

%Obtain all passive vessel filaments
Passives = get(vessel,'passives');
Filaments = get(Passives,'filaments');
RFil = get(Filaments(:),'r');             %Radial Filaments
ZFil = get(Filaments(:),'z');             %Axial Filaments
%Plot null-field eddy currents within a cross-section of the vessel
close all
figure; hold on; grid on; box on; axis equal;
plot(coilset);
scatter3(RFil,ZFil,VesselEddyCurrentsNull/1000,100,VesselEddyCurrentsNull/1000,'filled');
title('SMART Vessel Null Eddy Currents iter(0)');
view(2) %2D view
colormap(colourmap);
cbar = colorbar;
cbar.Label.String = 'Eddy-Current I_{eddy} [kA]';
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_EddyCurrentNull';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%% PLOT VESSEL EDDY STRESSES %%%%%%%%%%%%%%%%%%%%%%%%

%Plot figure showing vessel eddy stresses
close all
figure; hold on; grid on; box on; axis equal;
plot(coilset);
%plot(vessel);
quiver(R_Fil_Array,Z_Fil_Array,StressR,StressZ,'color',[1 0 0],'AutoScale','on');
title('SMART Vessel Eddy-Stresses iter(1)');
view(2) %2D view
legend(gca,'hide');
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_EddyStresses';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));
close('all')

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
figure('units','inch','position',[12 12 8 8]); hold on; grid on; box on;
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
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
Filename = '_PaschenCurves_00';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));

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
feedbackCoils = {'Div2'};
%[feedback_config, signals, weights, index] = efit_shape_controller(config, feedbackCoils, PertGeometry);
feedback = shape_controller(config, feedbackCoils, RGeo_Pert, ZGeo_Pert, rGeo_Pert, Kappa_Pert, delta_Pert);
Equil_Pert = set(EquilPert, config, 'feedback', feedback);                                 
EquilParams_Pert = parameters(Equil_Pert);
% ISSUE :: Equil_Pert is not shaped the same as Equil  
% ISSUE :: Coil currents are unreliable until Equil_Pert ~ Equil

%Extract coil currents required to offset perturbation
icoil_pert = get(Equil_Pert,'icoil'); 
CoilCurrentsPert = get(icoil_pert,'currents');
%}

%Obtain the Real Vertical Growth Rate from RZIp, there should only be one positive rate
%i.e. Gamma = eig(-curlyM\curlyR), sort for positive value(s) and save.
if length(Gamma(Gamma>0)) > 0; Gamma_Real = Gamma(Gamma>0); else Gamma_Real = 0; end     %[s-1]

%If feedback fails, overwrite with default equilibrium. 
%NOTE: Vertical control and feedback not yet implimented.
Equil_Pert = Equil; EquilParams_Pert = EquilParams;
PertGeometry = efitGeometry;
icoil_pert = icoil_efit; CoilCurrentsPert = CoilCurrentsEfit;
%end

%%%%%%%%%%%%%%%%%%%%%% PLOT PERTURBED EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%

%Plot perturbed equilibrium following convergence
Title = {'SMART Perturbed Equilibrium \Psi(R,Z)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_Equilibrium_Pert';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotEquilibrium({Equil_Pert},Title,CbarLabel,SaveString);
close('all')

%%%%%%%%%%%%%%%%%%%%%% PLOT VERTICAL GROWTH RATES  %%%%%%%%%%%%%%%%%%%%%%

%Plot the vertical growth rates as calculated by RZIp and manually
figure('units','inch','position',[12 12 8 8]); hold on; grid on; box on;
plot(Gamma,'LineWidth',2);
plot(ones(length(Gamma)),'k--','LineWidth',1.5);
title('SMART Vertical Growth Rates');
LegendString = {strcat('RZIp Gamma: ',string(Gamma_Real),' s^{-1}')};
legend(LegendString);
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
xlabel(gca,'Eigenvalue (Sorted) [-]');
ylabel(gca,'Growth Rate \gamma [s^{-1}]');
Filename = '_GammaAxial';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));
close('all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    ADDING EDDIES TO DISCHARGE EQUIL                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%  RE-COMPUTE DISCHARGE EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%

%Book-keeping for the start of each eddy current loop
close all

%Apply alterations to CoilWaveforms before re-computing to avoid numerical instabilities if required.
%Most crahses arise from LCFS in solenoid - findboundary.m function contains rules for LCFS selection
%Initial discharge currents (CoilWaveforms_Init) are more numerically stable than efit currents (CoilWaveforms_EFIT)
CoilWaveforms(:,TimeIndex_Discharge) = CoilWaveforms_Init(:,TimeIndex_Discharge);   
%CoilWaveforms(iDiv1,TimeIndex_Discharge) = I_Div1_Equil+100;               %Common Trick 1: try increasing IDiv1 and retry
%CoilWaveforms(iSol,TimeIndex_Discharge) = I_Sol_Equil-100;                 %Common Trick 2: try decreasing ISol and retry

%Compute equilibrium (Psi(R,Z)) from the supplied jprofile, icoil and geometry
%Returns target equilibrium and CoilWaveforms for PF1 and PF2 at requested time_Index
if isa(coilset,'fiesta_loadassembly') == 0; coilset = fiesta_loadassembly(coilset, vessel); end
[Equil_Passive,EquilParams_Passive,CoilWaveforms_Passive,efitGeometry_Passive,config_passive] = ...
    efit(jprofile,Irod,'nullconfig',efitGeometry_Init,CoilWaveforms,VesselEddyCurrents,TimeIndex_Discharge);
%NOTE :: NEED TO REPLACE NAN NULL-FIELD VALUES FOR APPROPRIATE COILS IN COILWAVEFORMS IF EDDY LOOP IS EMPLOYED

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%  RE-COMPUTE OPTIMISED NULL-FIELD  %%%%%%%%%%%%%%%%%%%%

%Update CoilWaveforms array with null-field values (Using NaN Mask)
%   CoilWaveforms_Passive = ApplyNaNMask(CoilWaveforms_Passive, CoilWaveforms_Init)
%ISSUE  :: SUBROUTINE ApplyNaNMask() DOES NOT YET EXIST, NEEDS WRITTEN.
%
%   RZIP_C = C_Passive;
%   CoilWaveforms = NullFieldWaveforms(CoilWaveforms_Passive, RZIP_C, sensor_btheta_passive, TimeIndex_NullField);
%ISSUE  :: RZIP_C is computed from RZIP - Need recomputed RZIP with eddy currents
%NOTE   :: REMEMBER THAT NULLFIELDWAVEFORMS EXPECTS NaN VALUES FOR COILS TO BE UPDATED

%Extract previously calculated efit coil currents without eddys
CoilCurrentsNull = transpose(CoilWaveforms_Passive(:,TimeIndex_NullField)); %Null-field coil currents from efit (with eddys)
CoilAndVesselCurrents = [CoilCurrentsNull, VesselEddyCurrentsNull];         %n=5+n_fil; coil + vessel filaments
icoil_Null_Passive = fiesta_icoil(coilset, CoilAndVesselCurrents);
%NOTE :: CoilWaveforms_Passive here does NOT have the re-optimised coil currents!

%Compute null-field equilibrium using null-field coil and vessel eddy currents
Equil_Null_Passive = fiesta_equilibrium('SMART-Null', config_passive, Irod, icoil_Null_Passive);
EquilParams_Null_Passive = parameters(Equil_Null_Passive);

%Extract the new coil currents from the null-field equilibrium:
icoil_Null_Passive = get(Equil_Null_Passive,'icoil'); 
CoilCurrentsNull_Passive = get(icoil_Null_Passive,'currents');

%%%%%%%%%%%%%%%%%  RE-COMPUTE BREAKDOWN FOR NULL-FIELD  %%%%%%%%%%%%%%%%%%%

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
[Vloop_Passive,Eloop_Passive,DeltaPhiSol_Passive,VloopArray_Passive,EloopArray_Passive,DeltaPhiSolArray_Passive]=...
    LoopVoltage(I_PF_output,time_adaptive,RSolCentreWinding,ZMaxSol,ZMinSol,EquilParams_Passive.rin);
Eloop_eff_Passive = abs(Eloop_Passive)*(BtorAvg_Null_Passive/BpolAvg_Null_Passive); %[V/m]   %Rough estimate threshold Eloop condition for breakdown
%Generally Eloop_eff > 100 [V/m] for startup with ECRH      (An2015)
%Generally Eloop_eff > 1000 [V/m] for solenoid only startup (Lloyd1991)

%Compute breakdown (avalanche) timescale at given pressure (!!! NEEDS TESTING !!!)
Pressure_Passive = EquilParams_Passive.P0*(7.5e-7);    %[Torr]  (~2e-4 Torr)
[TauBD_Passive,Pressure_Passive,Alpha_Passive,Vde_Passive] = ...
    AvalancheTimescale(Pressure_Passive,Eloop_Passive,Lc_Passive,ne,1.0,true);


%%%%%%%%%%%%%%  RE-COMPUTE DYNAMIC PLASMA & EDDY CURRENTS  %%%%%%%%%%%%%%%%

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
Filename = '_Equilibrium_01';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotEquilibrium({Equil_Passive},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%  PLOT NULL-FIELD PHI SURFACES  %%%%%%%%%%%%%%%%%%%%%

%Plot the optimised null-field phi
Title = {'SMART Null-field Equilibrium iter(1)',' '};
CbarLabel = 'Flux Surface Function \Psi(R,Z)';
Filename = '_NullField_01';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotEquilibrium({Equil_Null_Passive},Title,CbarLabel,SaveString);

%%%%%%%%%%%%%%%%%%%%%% PLOT NULL-FIELD BPOL  %%%%%%%%%%%%%%%%%%%%%

%Log poloidal and toroidal magnetic fields to show details (Sol Obscures)
logBpolData_Null_Passive = log10(BpolData_Null_Passive);
logBtorData_Null_Passive = log10(BtorData_Null_Passive);
%Plot the optimised null-field phi
Title = {'SMART Null-field iter(1)',' '};
CbarLabel = 'Null-field B_{\theta} log_{10}([T])';
Filename = '_NullBpol_01';
SaveString = strcat(SimDir,ShotName,Filename,FigExt);
PlotEquilibrium({logBpolData_Null_Passive},Title,CbarLabel,SaveString);

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
figure('units','inch','position',[12 12 8 8]); hold on; grid on; box on;
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
set(gca, 'FontSize', 18, 'LineWidth', 0.75);
Filename = '_PaschenCurves_01';
saveas(gcf, strcat(SimDir,ShotName,Filename,FigExt));
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

%Write final 2D, 1D and 0D equilibrium values to text files in a 'geqdsk-like' format
[fileID] = WriteEquilibrium(Equil, config, EquilDir, '_iter0', false);
[fileID] = WriteEquilibrium(Equil_Passive, config_passive, EquilDir, '', false);
[fileID] = WriteEquilibrium(Equil_Null, config, EquilDir, '_iter0', true);
[fileID] = WriteEquilibrium(Equil_Null_Passive, config_passive, EquilDir, '', true);
[fileID] = WriteEquilibrium(Equil_Pert, config, EquilDir, '_Pert', false);

%Write initial target geometry, efit geometry and perturbed geometry
[fileID] = WriteGeometry(efitGeometry_Init, EquilDir, 'efit_Geometry_Init.txt');
[fileID] = WriteGeometry(efitGeometry, EquilDir, 'efit_Geometry_Equil.txt');
[fileID] = WriteGeometry(PertGeometry, EquilDir, 'efit_Geometry_Pert.txt');

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Create subdirectory for coil current related data
icoilDir = strcat(ASCIIDir,'Coil_Data/'); mkdir(icoilDir);

Filename = strcat(icoilDir,'icoil_position.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n',             'Coil','R [m]  ',   'Z [m]  ',  'dR [m] ',  'dZ [m]'      );
fprintf(fileID,'%s %0.5f %0.5f %0.5f %0.5f\r\n', 'Sol ', RSolCentre, ZMaxSol,  Width_Sol/2,  Height_Sol    );
fprintf(fileID,'%s %0.5f %0.5f %0.5f %0.5f\r\n', 'PF1 ', R_PF1,      Z_PF1,    width_PF1/2,  height_PF1/2  );
fprintf(fileID,'%s %0.5f %0.5f %0.5f %0.5f\r\n', 'PF2 ', R_PF2,      Z_PF2,    width_PF2/2,  height_PF2/2  );
fprintf(fileID,'%s %0.5f %0.5f %0.5f %0.5f\r\n', 'Div1', R_Div1,     Z_Div1,   width_Div1/2, height_Div1/2 );
fprintf(fileID,'%s %0.5f %0.5f %0.5f %0.5f\r\n', 'Div2', R_Div2,     Z_Div2,   width_Div2/2, height_Div2/2 );

Filename = strcat(icoilDir,'IRod.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s\r\n', 'IRod [A]');
fprintf(fileID,'%1.12f\r\n', EquilParams_Passive.irod');

Filename = strcat(icoilDir,'icoil_init.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol [A]','PF1 [A]','PF2 [A]','Div1 [A]','Div2 [A]');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_init.Sol'; icoil_init.PF1'; icoil_init.PF2'; icoil_init.Div1'; icoil_init.Div2']);

Filename = strcat(icoilDir,'icoil_efit.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol [A]','PF1 [A]','PF2 [A]','Div1 [A]','Div2 [A]');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_efit_passive.Sol'; icoil_efit_passive.PF1'; icoil_efit_passive.PF2'; icoil_efit_passive.Div1'; icoil_efit_passive.Div2']);

Filename = strcat(icoilDir,'icoil_pert.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol [A]','PF1 [A]','PF2 [A]','Div1 [A]','Div2 [A]');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[icoil_pert.Sol'; icoil_pert.PF1'; icoil_pert.PF2'; icoil_pert.Div1'; icoil_pert.Div2']);

Filename = strcat(icoilDir,'icoil_null.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s %s\r\n', 'ISol [A]','PF1 [A]','PF2 [A]','Div1 [A]','Div2 [A]');
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
DynamicDir = strcat(ASCIIDir,'Current_Data/'); mkdir(DynamicDir);

Filename = strcat(DynamicDir,'Time.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s\r\n', 'time_adaptive [ms]');         %, 'time_linear [ms]');     %Add time_linear here, but array sizes are different
fprintf(fileID,'%1.12f\r\n', time_adaptive'*1000);      %, time_linear'*1000);      %Writes column-aligned data of (len(adapt)+len(linear))/2

Filename = strcat(DynamicDir,'IPlasma.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'time_adaptive [ms]','Ip_output [A]');
fprintf(fileID,'%1.12f %1.12f\r\n', [time_adaptive'*1000; Ip_output']);

Filename = strcat(DynamicDir,'IPassive1D.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'time_adaptive [ms]','Net_I_Passive [A]');
fprintf(fileID,'%1.12f %1.12f\r\n', [time_adaptive'*1000; Net_IPassive']);

Filename = strcat(DynamicDir,'IPassive2D.txt');
[fileID] = WriteMatrixCVS(I_Passive,Filename);            %ISSUE: Header and first line of data share same row, needs fixing!

Filename = strcat(DynamicDir,'VLoop.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s %s\r\n', 'time_adaptive [ms]','Vloop [V]','Eloop [V/m]','DeltaPhi [Vs]');
fprintf(fileID,'%1.12f %1.12f %1.12f %1.12f\r\n', [time_adaptive'*1000, VloopArray_Passive', EloopArray_Passive', DeltaPhiSolArray_Passive']);

%%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%          %%%%%%%%%%

%Misc outputs, save unordered in main RawData directory

Filename = strcat(ASCIIDir,'Eta.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'EtaPerp [Ohm m-1]', 'EtaPara [Ohm m-1]');
fprintf(fileID,'%1.12f %1.12f\r\n', EtaPerp', EtaPara');

Filename = strcat(ASCIIDir,'Bpol.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'Null_Bpol [T]', 'Null_Btor [T]');
fprintf(fileID,'%1.12f %1.12f\r\n', BpolAvg_Null_Passive', BtorAvg_Null_Passive');

Filename = strcat(ASCIIDir,'TauVessel.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s\r\n', 'TauVessel [ms]');
fprintf(fileID,'%1.12f\r\n', TauVessel*1000');

Filename = strcat(ASCIIDir,'betaP.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s\r\n', 'betaP [%]', 'betaP_Pert [%]');
fprintf(fileID,'%1.12f %1.12f\r\n', EquilParams_Passive.betap', EquilParams_Pert.betap');

Filename = strcat(ASCIIDir,'BDMetrics.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%s %s %s\r\n', 'Lc [m]', 'TauBD [ms]', 'Eloop_eff [V/m]');
fprintf(fileID,'%1.12f %1.12f %1.12f\r\n', Lc_Passive, TauBD_Passive*1000', Eloop_eff_Passive');

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

%{
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING ZONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Transpose initial coil currents
CoilCurrents = transpose(CoilWaveforms_Init(:,5));
icoil = fiesta_icoil(coilset, CoilCurrents); 

%Efit outputs coil currents resulting from the supplied jprofile, icoil and geometry
%Returns new currents for the requested coils: {'Coil1, {...}, 'Coiln'}
config = fiesta_configuration('SMART_config', Grid, coilset);
control = fiesta_control('diagnose',true, 'quiet',false, 'convergence',1e-4, 'boundary_method',2);
[efit_config, signals, weights, index] = efit_shape_controller(config, efitCoils, efitGeometry_Init);

%Inverse equilibrium, outputs coil currents resulting in the supplied jprofile and icoil config
Equilibrium = fiesta_equilibrium('SMART', config, -Irod, jprofile, control, efit_config, icoil, signals, weights);
EquilParams = parameters(Equilibrium)

plot(Equilibrium)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Equilibrium,EquilParams,efitCoilWaveforms,efitGeometry,config]= ...
    efit(Jprofile,Irod,InputConfig,InputGeometry,InputCoilWaveforms,InputVesselWaveforms,TimeIndex)
%   Computes equilibrium Psi(R,Z) topology as well as required PF, Div coil currents to achieve target shape
%   Boundary conditions are set by the input Jprofile(R), IRod, CoilCurrents, and target Geometry.
%   A user-defined sub-set of coil currents are allowed to vary, while other currents remain fixed
%   If a stable equilibrium is found, the equilibrium and coil currents required to achieve it are returned.
%   Definitions:
%INPUTS:  
%       Jprofile:               1D array of current density across the midplane (Z=0.0) from TOPEOL         [A/m^2]
%       Irod:                   0D central toroidal field rod current required to achieve Btoroidal         [A]
%       InputConfig:            FIESTA efit configuration object - MUST use same config throught run        [-]
%       InputGeometry:          1D array of target plasma shaping parameters (Rgeo,Zgeo,a,kappa,delta)      [-]
%       InputCoilWaveforms:     2D array (ncoil,timeindex) containing inpit Sol, PF, and Div coil currents  [A] 
%       InputVesselWaveforms:   1D array of vessel eddy currents at time "TimeIndex" - OPTIONAL             [A]
%                               (Origin filament at top left of vessel, proceeding clockwise)
%       TimeIndex:              0D integer of desired current waveform time index to use                    [-]
%
%OUTPUTS:
%       Equilibrium:            2D FIESTA object containing the equilibrium Psi(R,Z) and other data         [-]
%       EquilParams:            1D array of key equilibrium parameters                                      [-]
%       efitCoilWaveforms:      2D array (ncoil,timeindex) containing efit Sol, PF, and Div coil currents   [A]
%                               This represents the coil waveforms required to achieve the Equilibrium
%       efitGeometry:           1D array of achieved plasma shaping parameters (Rgeo,Zgeo,a,kappa,delta)    [-]
%       Config:                 FIESTA efit configuration object - For use if equilibrium is recomputed     [-]
%


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

function [Time_Linear,Time_Adaptive,I_PF_output,V_PF_output,Ip_output,Vp_output,I_Passive]= ...
    DynamicCurrents(CoilWaveforms,TimeIndices,CurlyM,CurlyR)
%   Compute time_resolved plasma and vessel eddy currents employing state_space_including_passive_elements_v4()
%   Computation interpolates instantanious coil waveform currents and employs vessel Inductance/Resistance
%   Definitions:
%INPUTS:
%       CoilWaveforms:      2D array (ncoil,timeindex) containing Sol, PF, and Div coil currents     [A] 
%       TimeIndices:        1D array containing coil waveform time indices                           [-]
%                           (must align with CoilWaveform indices)
%       CurlyM:             2D Inductance Matrix output from RZIp 
%       CurlyR:             2D Resistance Matrix output from RZIp
%
%OUTPUTS:
%       Time_Linear:        1D array of time points with linear spacing                 [s]
%       Time_Adaptive:      1D array of time points with non-linear spacing             [s]
%       I_PF_output:        1D array of temporally resolved PF, Div coil Currents       [A]
%       V_PF_output:        1D array of temporally resolved PF, Div coil Voltages       [V]
%       Ip_output:          1D array of temporally resolved Plasma Current              [A]
%       Vp_output:          1D array of temporally resolved Plasma Voltage              [V]
%       I_Passive:          1D array of temporally resolved vessel filament currents    [A]
%

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
	%CurlyM and CurlyR are the vessel and coil inductance and resistance matrices.
    [V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, Time_Adaptive ] = ...
        state_space_including_passive_elements_v4( CurlyM, CurlyR, Time_Linear, IPFinput_Continous, VPFinput_Continous, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names',coil_names, 'show_plot',false, 'turns',coilturns, 'currentScale',1e3, 'PF_colors',PF_colors );
    
    %Set breakdown time and prepare Ip_long and Vp_long for voltage driven Ip
    %NOTE: Time_Breakdown should really be set to avalanche timescale, although zero is probably close enough...
    Time_Breakdown = 0;                                          %Set time for plasma breakdown (default 0)
    Time_Plasma = Time_Adaptive > Time_Breakdown;                %Set times for which plasma exists
    Vp_output(Time_Plasma) = 0;                                  %Set voltage to zero following breakdown
    Vp_long = interp1(Time_Adaptive, Vp_output, Time_Linear);    %Sets Vp_long = 0 when Time_Linear > 0.
    Ip_long = NaN*Vp_long;                                       %Sets Ip_long to 'NaN' array (Voltage Driven)
    
    %Compute dynamic coil currents employing Voltage Driven Ip
    [ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, Time_Adaptive ] = ...
        state_space_including_passive_elements_v4( CurlyM, CurlyR, Time_Linear, IPFinput_Continous, VPFinput_Continous, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names',coil_names, 'show_plot',false, 'turns',coilturns, 'currentScale',1e3, 'PF_colors',PF_colors );

    %Clean up before returning to main code
    close all
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CoilWaveformsOutput]=...
    NullFieldWaveforms(CoilWaveformsInput,RZIP_C,sensor_btheta,TimeIndex)
   
    %Obtain required global variables
    global iPF1; global iPF2;
    global iDiv1; global iDiv2;
    global iSol; 
    
    %Determine number of coil waveforms and time points within each
    SizeCoilArrays = size(CoilWaveformsInput);
    nCoils = SizeCoilArrays(1);             %Number of PF/Div coils
    nTime = SizeCoilArrays(2);              %Number of TimeVertics

    %Extract scaling factors for null-field coil currents
    %Cn is the part of the matrix C related to the sensors (see response.m)
    C_temp = RZIP_C(end-get(sensor_btheta,'n')+1:end,1:nCoils);
    C1 = C_temp(:,1);          %Elements of C_temp(Cn) for Sol coil
    
    D1_PF1 = C_temp(:,iPF1);   %Elements of C_temp(Cn) for PF1 coil
    D1_PF2 = C_temp(:,iPF2);   %Elements of C_temp(Cn) for PF2 coil
    D1_Div1 = C_temp(:,iDiv1); %Elements of C_temp(Cn) for Div1 coil
    D1_Div2 = C_temp(:,iDiv2); %Elements of C_temp(Cn) for Div2 coil
    
    %Determine if Div1 is in series with solenoid or not and optimise null-field accordingly
    if isnan(CoilWaveformsInput(iDiv1,TimeIndex)) == true                   %If Div1 IS NOT in series with Sol
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


function SensorBTheta=InitiateBSensors(RGeo,ZGeo,R_null)
%   Initiate virtual sensors within null-field region

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

function VesselEddyCurrents=...
    ExtractPassiveCurrents(I_Passive,time_adaptive,DischargeTime)
%   Extract Passive Vessel currents at desired time [s] during the pulse
%   If discharge time is supplied as 'false' then the absolute maximum currents are extracted
%   Definitions:
%INPUTS:  
%       I_Passive:          1D array of temporally resolved vessel filament currents    [A]
%       Time_Adaptive:      1D array of time points with non-linear spacing             [s]
%       DischargeTime:      0D float of desired time to extract Eddy Currents from      [s]
%
%OUTPUTS:
%       VesselEddyCurrents: 1D array of vessel eddy currents at time DischargeTime      [A]
%                           (Origin filament at top left of vessel, proceeding clockwise)

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

function [StressR,StressZ,StressR_max,StressZ_max] = ...
    VesselStresses(Equilibrium,VesselFilaments,FilamentArea)

    %Obtain required global variables
    global Grid;
    global vessel;
    
    %Extract grid sizes and number of filaments in R and Z
    rGrid=get(Grid,'r'); zGrid=get(Grid,'z');
    ptmp = get(vessel,'passives'); ftmp = get(ptmp,'filaments');
    RFil = get(ftmp(:),'r'); %dim 1*number of filaments
    ZFil = get(ftmp(:),'z'); %dim 1*number of filaments
    
    %Obtain the equilibrium B-field in R,Z and Phi
    [BrData,BzData,BPhiData,BpolData,BtorData] = ExtractBField(Equilibrium);

    %Interpolate the B-fields onto a grid that aligns with the vessel grid 
    %These are the values of the B-field at the vessel grid points
    %Br_interp aligns with the vessel filament cells (RR, ZZ) from before
    Br_interp = @(r,z) interpn(zGrid,rGrid,BrData,ZFil,RFil);
    Bz_interp = @(r,z) interpn(zGrid,rGrid,BzData,ZFil,RFil);
    Bphi_interp = @(r,z) interpn(zGrid,rGrid,BPhiData,ZFil,RFil);

    %Extract B-field at vessel walls - meshes are aligned so indexes are the same
    Br_vessel = Br_interp(RFil,ZFil);
    Bz_vessel = Bz_interp(RFil,ZFil);
    Bphi_vessel = Bphi_interp(RFil,ZFil);

    %Combine Bfield values in [R,Phi,Z] into a single array for the vessel
    %size(number of filaments*3), each row is the vector field on one filament
    B_vessel=[Br_vessel' Bphi_vessel' Bz_vessel']; 
    %The maximum current on each vessel filament is I_Passive_fil (size 1*number of filaments)
    %Current vector is in the phi direction so only take magnitude [0, 1*I_Passive(phi), 0]
    VesselEddyCurrentVec=VesselFilaments'*[0 1 0]; 		%size [number of filaments*3]

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
    PressureR = abs(Force_max(1))/(FilamentArea);	%[Pa]
    PressureZ = abs(Force_max(3))/(FilamentArea);	%[Pa]

    %Stress acting on vessel wall is the combined force divided by the unit filiment area
    %These are directional, some are negative and some are positive
    StressR = (Force_fil(:,1))/(FilamentArea);      %[Pa]
    StressZ = (Force_fil(:,3))/(FilamentArea);      %[Pa]
    %Obtain maximum radial and axial stresses - either positive or negative
    StressR_max = max(abs(StressR));				%[Pa]
    StressZ_max = max(abs(StressZ));				%[Pa]
    
    %Clean up before returning to main code
    close all
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILITY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vessel_filaments,R_Fil_Array,Z_Fil_Array]=...
    CreateRectilinearVessel(VesselDimensions,WallThickness,FilamentArea,WallNorm)
%   Creates a rectilinear vessel comprised of passive toroidal current carrying filaments
%   Vessel filaments may employ a fixed cross-sectional area or vary by wall thickness
%   Preferred operation: WallNorm=true scales filaments to maintain constant cross-sectional area
%   
%   Definitions:
%   INPUT  :: VesselDimensions:     1D array containing vessel inner SI dimensions      [m]
%               [RMinCentre, RMaxCentre, ZMinCentre, ZMaxCentre]
%   INPUT  :: WallThickness:        1D array containing vessel SI wall thicknesses      [m]
%               [VWall_Upper, VWall_Outboard, VWall_Lower, VWall_Inboard]
%   INPUT  :: FilamentArea:         0D scalar for vessel filament cross-sectional area  [m]
%               Set FilamentArea to 0.0 to autocompute filament area
%   INPUT  :: WallNorm:             Boolian defining if filament area is scaled or fixed:
%              	'true' for normalised (fixed) filament area as supplied by FilamentArea
%               'false' for variable filament area dependant upon wall thickness/length
%
%   OUTPUT :: vessel_filaments:    (1, nfil) array of fiesta_filament structures        [m]
%               [(R1,Z1), (R2,Z2)... (Rn-1,Zn-1),(Rn,Zn)] 
%   OUTPUT :: R_Fil_Array:         1D array of filament radial coordinates              [m] 
%               (Origin filament at top left of vessel, proceeding clockwise)
%               [R1, R2, R3, R4... Rn-1 Rn]
%   OUTPUT :: Z_Fil_Array:         1D array of filament axial coordinates               [m] 
%               (Origin filament at top left of vessel, proceeding clockwise)
%               [Z1, Z2, Z3, Z4... Zn-1 Zn]
%   NOTES:
%          :: To avoid duplicate filaments, the corner vertices are 'owned' by the axial filaments.
%

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

    %If WallNorm is false, filament area is not normalised and will vary between walls (not recommended)
    if WallNorm == false
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

            %Top and bottom walls 'own' their vertices - i.e. axial walls posess all (NumFil) filaments
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

            %Radial walls don't own the vertices - i.e radial walls possess only (NumFil-2) filaments
            %Remove first and last filaments from each radial array to avoid duplicate filaments at the corners.
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

function [Coil_Circuit]=...
CreateSMARTCoilCircuit(Label,Rc,Zc,DR,DZ,nt,nZ,nR,CoilTemp,CoilResistivity,CoilDensity,symmetry)
%   Creates FIESTA coil object from provided coordinates, size and turns
%   Applies symmetry if required and returns a FIESTA coil circuit object
%   Definitions:
	% Label             Coil object label                                  [-]
	% Rc, Zc            Radial and Axial coordinates of coil center        [m]
	% DR, DZ            Coil total width and height                        [m]
	% nt                Total number of coil windings (turns)              [-]
	% nZ, nR            Coils are distributed in an array of nZ x nR coils [-]
    % CoilTemp          Coil temperature (Time Independent)                [K]
    % CoilResistivity   Coil Material Resistivity                          [Ohm m-1]
    % CoilDensity       Coil Material Density                              [kg m-3]

    %Confirm valid coil configuration has been supplied - return if false
	if rem(nt, nZ*nR) > 0
		warning('number of turns is not a multiple of nZ*nR');
		Coil_Circuit = false;
		return
    end

    %Create coil object of size DR x DZ at position Rc,Zc and label
	Coil1 = create_coil( Rc,Zc,DR,DZ,nt, nZ, nR, CoilTemp, CoilResistivity, CoilDensity);
	Coil1 = set(Coil1,'label','unique');
    
    %Create axially symmetric coil if required and combine into coil circuit
	if symmetry == true
		Coil2 = create_coil( Rc,-Zc,DR,DZ,nt, nZ, nR, CoilTemp, CoilResistivity, CoilDensity);
		Coil1 = set(Coil1,'label','up');
		Coil2 = set(Coil2,'label','down');
		Coil_Circuit = fiesta_circuit(Label,[1 1],[Coil1 Coil2]);
    %Create coil circuit with single coil if symmetry not requested
	else
		Coil_Circuit = fiesta_circuit(Label,[1],[Coil1]);    
	end

end

function c=create_coil(Rc,Zc,DR,DZ,nt,nZ,nR,coil_temperature,resistivity,density)

    % ang1 and ang2 use the DIIID definitions, ang1 is the vertical wall, ang2 horizontal
    % they are tangent angles, and when ==0  the vertical wall is vertical, etc.
    ang1 = 0;
    ang2 = 0;

    %Determine size of each coil filament from the number of turns (nt) and total size (DR,DZ)
	TurnsPerCoil = nt/(nZ*nR); 
	[Zcoils,dz] = divideIntoIntervals( Zc-0.5*DZ, Zc+0.5*DZ, nZ);
	[Rcoils,dr] = divideIntoIntervals( Rc-0.5*DR, Rc+0.5*DR, nR);
	nFilament = nZ*nR;
    
    %For each row of the coil
	for i=nZ:-1:1
        %For each column of the coil
		for j=nR:-1:1
            %Create a fiesta filament at (R,Z) with size (dr,dz) of windings turnsPerCoil
            %No coil angle is applied, all coils are aligned horizontal to the mid-plane
		    filament( nFilament) = fiesta_filament(Rcoils(j),Zcoils(i), dr, dz, TurnsPerCoil, ang1, ang2);
		    nFilament = nFilament - 1 ;
		end
    end
	c = fiesta_coil('',filament, 'Blue', resistivity, density );
end

function [x,dx] = divideIntoIntervals(xi,xf,n)
	% divide the interval [xi,xf] in n parts
	% output: centers and width of the intervals
	x = linspace(xi,xf,2*n+1);
	dx = x(3)-x(1);
	x = x(2:2:end);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EtaPerp,EtaPara] = SpitzerResistivity(ne,Te,Z_eff)
%   Calculates perpendicular and parallel Spitzer Resistivities including effective nuclear charge
%   Citation from: "Self-consistent equilibrium calculation", see definition below eq. 22
%   http://dx.doi.org/10.1088/0741-3335/42/12/304
%
%   Definitions:
%       ne      :: Electron density             [m^-3]
%       Te      :: Electron temperature         [eV]
%       Z_eff   :: Effective Nuclear Charge     [e-]
%       EtaPerp :: Perpendicular Resistivity    [Ohm m^-1]
%       EtaPara :: Parallel Resistivity         [Ohm m^-1]

    %Inport any required global variables
    global e; global mass_e;
    global epsilon0;

    %Compute Electron-ion and Electron-Electron Colomb Logarithms
    Lambda_ee = 23.5-log(sqrt(ne)*Te^(-5/4)-sqrt(1e-5+log(Te)-2)^2/16);
    Lambda_ei = 24.0-log(sqrt(ne)*Te^(-1));
    
    %Compute Relative Resistivity Coefficient from Nuclear Charge
    FZeff = (1.000+1.198*Z_eff+0.222*Z_eff^2)/(1.000+2.966*Z_eff*0.753*Z_eff^2);

    %Compute electron average thermal velocity
    Vth = sqrt(2*e*Te/mass_e);                          %[m s^-1]

    %Compute average electron-electron and electron ion scattering timescales
    %Source, National Plasma Formulary page 32-33
    tau_ee = 12*(pi^(3/2))*(epsilon0^2)*(mass_e^2)*(Vth^3)/(4*ne*(e^4)*Lambda_ee);
    tau_ei = 12*(pi^(3/2))*(epsilon0^2)*(mass_e^2)*(Vth^3)/(4*ne*(e^4)*Lambda_ei);

    %Compute Perpendicular and parallel Spitzer resistivities
    EtaPerp = (Z_eff*mass_e)/(ne*e^2*tau_ee);            %[Ohm m^-1]
    EtaPara = EtaPerp*FZeff;                             %[Ohm m^-1]
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function could be expanded to include other materials if required
function [Resistivity]=InterpMaterialResistivity(InputTemperature) 
%   Resistivity [Ohm m]
%	InputTemperature [K]

    CopperTemperatures  = [20.0                    25.0                   30.0                   35.0                  40.0                  45.0                  50.0                  55.0                  60.0                  70.0                 80.0                 90.0                 100.0                125.0                150.0                175.0                200.0               225.0               250.0               273.150             293.0               300.0               350.0               400.0               500.0               600.0               700.0               800.0               900.0               1000.0              1100.0              1200.0              1300.0              1357.60         ];
    CopperResistivities = [0.000798000000000000    0.00249000000000000    0.00628000000000000    0.0127000000000000    0.0219000000000000    0.0338000000000000    0.0498000000000000    0.0707000000000000    0.0951000000000000    0.152000000000000    0.213000000000000    0.279000000000000    0.346000000000000    0.520000000000000    0.697000000000000    0.872000000000000    1.04400000000000    1.21500000000000    1.38500000000000    1.54100000000000    1.67600000000000    1.72300000000000    2.06100000000000    2.40000000000000    3.08800000000000    3.79000000000000    4.51200000000000    5.26000000000000    6.03900000000000    6.85600000000000    7.71500000000000    8.62400000000000    9.59000000000000    10.1690000000000] * 1.0e-8;

    Resistivity = interp1( CopperTemperatures, CopperResistivities, InputTemperature );
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
function [Vloop,Eloop,DeltaPhi,VloopArray,EloopArray,DeltaPhiArray]=...
    LoopVoltage(CoilCurrents,AdaptiveTime,RSol,ZMaxSol,ZMinSol,RLoop)

    %Obtain required global variables
    global coilturns; global iSol;
    global mu0;
    
    %Compute current surface area and voltage loop path length 
    IAreaLoop = pi*RSol^2;          %Solenoid cross-sectional area      [m^2]
    VLoopLength = 2*pi*RLoop;       %Voltage loop toroidal path length  [m]
    SolLength = ZMaxSol-ZMinSol;    %Solenoid axial length              [m]
    
    %Initiate arrays
    dIdtArray = zeros(1,length(AdaptiveTime));
    VloopArray = zeros(1,length(AdaptiveTime));
    EloopArray = zeros(1,length(AdaptiveTime));
    DeltaPhiArray = zeros(1,length(AdaptiveTime));
    
    %Compute the maximum average change in Sol current over the full discharge
    for i=2:length(AdaptiveTime)
        dI = CoilCurrents(i,iSol)-CoilCurrents(i-1,iSol);       %Change in solenoid current [A]
        dt = AdaptiveTime(i)-AdaptiveTime(i-1);                 %Change in adaptive time [s]
        dIdtArray(i-1) = dI/dt;                                 %Delta solenoid current [A/s]
        
        %Compute induced voltage during solenoid ramp-down and associated E-field at RLoop
        VloopArray(i-1) = IAreaLoop*((mu0*dIdtArray(i-1)*coilturns(iSol))/SolLength);   %Solenoid induced loop voltage [V] 
        EloopArray(i-1) = VloopArray(i-1)/VLoopLength;                                  %Solenoid induced toroidal E-field [V/m]
        DeltaPhiArray(i-1) = VloopArray(i-1)*dt;                                        %Solenoid magnetic flux [Vs]
    end
    VloopArray = transpose(VloopArray);
    EloopArray = transpose(EloopArray);
    DeltaPhiArray = transpose(DeltaPhiArray);
    
    %Extract maximum loop voltage, toroidal E-field and solenoid magnetic flux during solenoid ramp-down
    [dIdt,Index] = min(dIdtArray);          %Maximum negative Sol current ramp [A/t]
    Vloop = VloopArray(Index);              %Loop voltage at maximum negative solenoid ramp [V]
    Eloop = EloopArray(Index);              %Loop E-field at maximum negative solenoid ramp [V/m]
    DeltaPhi = DeltaPhiArray(Index);        %Solenoid flux swing [Vs]
    
    %NOTE :: MAXIMUM POSSIBLE SOLENOID MAGNETIC FLUX
    %(mu0*pi*ncoil*RSol^2)/(HeightSol*2*MaxISol)               %Assume one linear ramp [Vs]
    %(((mu0*pi*nSol*RSolCentreWinding^2)/(ZMaxSol*2))*2*12500)  %Sanity check ~ 0.263   [Vs]
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute breakdown (avalanche) timescale at given pressure
%Assumes fixed pre-ionisation electron density and breakdown fraction
function [TauBD,Pressure,Alpha,Vde]=...
    AvalancheTimescale(Pressure,Eloop,Lc,ne,ne0,minimise)
    %%% Notes ::
    %If tauBD approaches zero from negative value then Alpha^{-1} [m] < Lc [m] and breakdown is impossible
    %   In this case, TauBD will cross the X-axis at ne0 = BreakdownFrac*ne
    %   Notably, the growth rate is negative, so TauBD increases with ne0
    %If TauBD approaches zero from positive value then Alpha^{-1} [m] > Lc [m] and breakdown is possible
    %   In this case, TauBd will approach zero with increasing ne0
    %   Notably, the growth rate is positive, so TauBD decreases with ne0
    %   This data can be used to set a lower limit for the breakdown timescale (TauR1) for a given pre-ionisation density (ne0)
    %%%
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
    TauBD = log(neFraction)/(Vde*(Alpha-(1/Lc)) );   %Breakdown time          [s]

    %Compute minimum pressure and avalanche timescale for breakdown
    if minimise == true
        PressureLimits = [1e-6, 1e-3]; Resolution = 25000;
        Pressure_Array = linspace(PressureLimits(1),PressureLimits(2),Resolution);

        for i=1:length(Pressure_Array)
            Alpha_Array(i) = 510*Pressure_Array(i)*exp((-1.25e4*Pressure_Array(i))/Eloop);
            Vde_Array(i) = (Eta*Eloop)/Pressure_Array(i);
            TauBD_Array(i) = log(neFraction)/(Vde_Array(i)*(Alpha_Array(i)-(1/Lc)) );
        end
        %Take minimum values of Vde, Pres, Tau :: corresponding to the maximum townsend coefficient
        Alpha = max(Alpha_Array);                                 %Townsend Coeff          [m-1]
        Vde = Vde_Array( Alpha_Array==Alpha );                    %Electron Drift Speed    [m/s]
        Pressure = Pressure_Array( Alpha_Array==Alpha );          %Townsend coeff          [m-1]
        TauBD = TauBD_Array( Alpha_Array==Alpha );                %Breakdown time          [s]
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Write 2D, 1D and 0D equilibria out to files
function [fileID]=WriteEquilibrium(Equilibrium,config,EquilDir,EquilName,VacuumField)

    %Get any required global variables
    global Grid;
    %Extract Grid arrays for reshaping if required
    rGrid = get(Grid,'r'); zGrid = get(Grid,'z');
    
    %If Equil contains plasma (i.e. not a vacuum field) then use in-built geqdsk function
    if VacuumField == false
        %Write 2D qeqdsk equilibrium file
        Filename = strcat(EquilDir,'Equil',EquilName,'.txt');
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
        Filename = strcat(EquilDir,'Null_Equil',EquilName,'.txt');
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
        Filename = strcat(EquilDir,'Null_Bpol',EquilName,'.txt');
        fileID=fopen(Filename,'w');
        fprintf(fileID,'%s %s %s \n', '     Bpol_Null', string(length(rGrid)), string(length(zGrid)));
        for i = 1:size(BpolData_Null,1)
            fprintf(fileID,'%g\t',BpolData_Null(i,:));
            fprintf(fileID,'\n');
        end
        
        %Write Null Btor as 2D array
        Filename = strcat(EquilDir,'Null_Btor',EquilName,'.txt');
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

function [fileID]=WriteMatrixCVS(Matrix,Filename)

    %Open file with requested filename
    fileID = fopen(Filename,'w');
    Variable = 'Fil ';              %(hardcoded for now)
    Unit = ' [A]';                  %(hardcoded for now)

    %Write column headers in first row
    for i = 1:length(Matrix(1,:))
          fprintf(fileID,'%s ', strcat(Variable,string(i),Unit) );
    end
    
    %Write data row-wise, containing one data point from each column
    for i = 1:length(Matrix(1,:))        %Length of columns
        for j = 1:length(Matrix)         %Length of data within each column
            fprintf(fileID,'%1.12f ', Matrix(j,i) );
        end
        fprintf(fileID,'\n');            %New line after each row
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig=PlotEquilibrium(Arrays,Title,CbarLabel,SaveString)

    %USEFUL FIESTA CLASSES
    %class(sensor_btheta) = 'fiesta_sensor_btheta'
    %class(equil) = 'fiesta_equilibrium'
    
    %Obtain required global variables
    global colourmap; cmap = colourmap;
    global Grid; global vessel;
    global coilset;
    
    %Extract grid sizes in R and Z
    rGrid=get(Grid,'r'); zGrid=get(Grid,'z');
    GridCells_R = length(rGrid); Rlim = max(rGrid);
    GridCells_Z = length(zGrid); Zlim = max(zGrid);

    %Initiate figure, axes and aspect ratio (fixed for now)
    close all
    figure; hold on; box on; axis equal;
    
    %for each supplied sub-array (Arrays{i})
    for i=1:length(Arrays);
        %If data is a FIESTA class then use in-built function
        if isa(Arrays{i},'fiesta_equilibrium') == true;
            plot(Arrays{i});
        elseif isa(Arrays{i},'fiesta_sensor_btheta') == true;
            plot(Arrays{i});
        %Else plot as a regular contour plot using the supplied grid
        else
            contourf(rGrid,zGrid,Arrays{i});
        end
    end
    
    %Plot Vessel and Coilset
    plot(vessel);
    plot(coilset);
    
    %Colourmap
    colormap(cmap);
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

function fig=PlotVesselOverview(SaveString)

    %Obtain required global variables
    global colourmap; cmap = colourmap;
    global Grid; global vessel;
    global coilset;
    
    %Extract grid sizes in R and Z
    rGrid=get(Grid,'r'); zGrid=get(Grid,'z');
    GridCells_R = length(rGrid); Rlim = max(rGrid);
    GridCells_Z = length(zGrid); Zlim = max(zGrid);

    %Initiate figure, axes and aspect ratio (fixed for now)
    figure('Renderer', 'painters', 'Position', [1,1, 700 1100], 'visible', 'off');
    axes; hold on; box on; AspectRatio = [1,2,1];
    
    %Plot Vessel and Coilset
    ax1 = gca;
    coil = plot(coilset); set(coil, 'EdgeColor', 'k')
    fil = plot(vessel); set(fil, 'EdgeColor', 'k')
    %%%%%
    pbaspect(ax1,AspectRatio)
    set(ax1,'XLim',[0.00 1.0]);     %0 to Rlim
    set(ax1,'YLim',[-1.0 1.0]);     %-Zlim to Zlim
    set(ax1, 'FontSize', 20, 'LineWidth', 0.75);
    xlabel(ax1,'R (m)');
    ylabel(ax1,'Z (m)');
    grid(ax1,'minor')
    %%%%%
    ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
    pbaspect(ax2,AspectRatio)
    set(ax2,'XLim',[0.00 GridCells_R]);
    set(ax2,'YLim',[0.00 GridCells_Z]); 
    set(ax2, 'FontSize', 20, 'LineWidth', 0.75);
    xlabel(ax2,'R_{Grid} (Cells)');
    ylabel(ax2,'Z_{Grid} (Cells)');
    %%%%%
    saveas(gcf, SaveString);

    %Plot a zoomed image of the upper inboard side. 
    %Zoom location fixed for now, ideally would be togglable.
    %{
    figure; axes; hold on; box on;
    set(gca, 'DataAspectRatio', [1,1,1], 'NextPlot', 'add')
    
    %Plot Vessel and Coilset
    coil = plot(coilset); set(coil, 'EdgeColor', 'k');
    fil = plot(vessel); set(fil, 'EdgeColor', 'k');
    %%%%%
    tau=get(vessel, 'tau');
    r=get(vessel, 'r'); z=get(vessel, 'z');
    set(gca,'XLim',[0.05 0.45]);
    set(gca,'YLim',[0.45 1.00]);
    set(gca, 'FontSize', 13, 'LineWidth', 0.75);
    xlabel(gca,'R (m)');
    ylabel(gca,'Z (m)');
    Filename = '_VesselFilaments_Closeup';
    %%%%%
    saveas(gcf, strcat(SaveString,'_Zoomed'));   
%   ISSUE :: need to remove extension first
    %}
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Custom colour map: IDL Gamma-II colormap (PSFT group map)
% INPUTS:
%    - m: Number of points to define the color scale (default = 256)
function cm_data=Gamma_II(m)

    %Base colours in 8 bit format
    T = [0,   0,   0                %// black
         0,   0,   255              %// blue
         255, 0,   0                %// red
         255, 255, 0                %// yellow
         255, 255, 255]./255;       %// white
     %Setting color intervals length
     x = [0                         %// black
         70                         %// blue
         130                        %// red
         200                        %// yellow
         255];                      %// white

     %Interpolation between colors
     if nargin < 1
        map = interp1(x/255,T, linspace(0,1,255));
        cm_data = colormap(map);
     else
        map = interp1(x/255,T, linspace(0,1,m));
        cm_data = colormap(map);
        
        %Test figure: Color bar limits form 0 to 1.0 (black to white)
        I = linspace(0, 1.0, 255);
        imagesc( I(ones(1,10),:)' );
     end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Custom colour map: Python plasma colourmap (SJDoyle's colourmap)
function cm_data=Plasma(m)
cm =  [[  5.03832136e-02,   2.98028976e-02,   5.27974883e-01],
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

    %If no user input, use pre-defined colourmap
    if nargin < 1
        cm_data = cm;
    %User input (integer) determines colour gradient 'spacing' 
    %higher input number --> higher gradients between colours
    else
        hsv = rgb2hsv(cm);
        hsv(153:end,1) = hsv(153:end,1)+1; % hardcoded
        cm_data = interp1( linspace(0,1,size(cm,1)),hsv,linspace(0,1,m) );
        cm_data(cm_data(:,1)>1,1) = cm_data(cm_data(:,1)>1,1)-1;
        cm_data = hsv2rgb(cm_data);
        
         %Test figure: Color bar limits form 0 to 1.0 (black to white)
         I = linspace(0, 1.0, length(cm_data)); colormap(cm_data);
         imagesc( I(ones(1,10),:)' );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% END OF FILE %






































% FIELD LINE INTEGRATOR - IGNORE THIS - %

%% 
%%%%%%%%%%%   L calc by field line integration    %%%%%%%%%%%
%{
         %Lazarus paper 1998, they compute connective length by avergaing on
            %9 lines, the line with Bpol min, and the 8 surroundings. 
         %However can compute the lines in all the VV
      
   IntMethod='Lp';   %Lp faster, Phi for debug!!!      % ''Phi'or 'Lp' to switch the integration mode 

        %Grid
        %I redefine the grid to compute the connection length, for less computer
        %demands (time). Will label it with an L at the end!

        n_pnts_insideL=20;          %100 is the ideal to have good plots of the fields, but the L int failures. 
                %30 gives relative resol of 3.2% (anterior value)
                %20 gives 5.2%, good enough
                
        r_inside_VVL=linspace(VesselRMinInner,VesselRMaxInner,n_pnts_insideL); 
        z_inside_VVL=linspace(VesselZMinInner,VesselZMaxInner,n_pnts_insideL);
        
        %Resolution
        Resol_R=(r_inside_VVL(2)-r_inside_VVL(1))   %[m] R resolution  
        Resol_Z=(z_inside_VVL(2)-z_inside_VVL(1))   %[m] Z resolution        
        
        %Make new mesh
        %Resolution relative to VV size
        Resol_rel_R=Resol_R/(VesselRMaxInner-VesselRMinInner)*100  %[%] relative R resolution
        Resol_rel_Z=Resol_Z/(VesselZMaxInner-VesselZMinInner)*100  %[%] relative Z resolution
        
         %Lets do a meshgrid, will be needed
        [r_ins_VVL,z_ins_VVL]=meshgrid(r_inside_VVL,z_inside_VVL);
        
        %Plot
            figure;
            plot(r_ins_VVL,z_ins_VVL,'r.')
            hold on
            plot(vessel)
            axis equal
            xlabel('R (m)')
            ylabel('Z (m)')
            title(sprintf('meshgrid for the integration with %d^2 points',n_pnts_insideL))

        %To remove the points at the VV:
         %r_inside_VVL=r_inside_VVL(2:end-1);
         %z_inside_VVL=z_inside_VVL(2:end-1); 
        
       %Now, if the coils are inside, the grid points there have to be removed because ode 'se raya', and spent too long time
       global Rin_coils Rout_coils Zdown_coils Zup_coils
       
       %Identifying the coil sets for removal (how removed?)
       Rin_coils=[R_PF2 R_Div1 R_Div2]-[width_PF2 width_Div1 width_Div2]/2; %[m] inner R of coilset (no SOl and PF1 (outside))
       Rout_coils=[R_PF2 R_Div1 R_Div2]+[width_PF2 width_Div1 width_Div2]/2;   %[m]outer R of coilset (no SOl and PF1 (outside))
       Zdown_coils=[Z_PF2 Z_Div1 Z_Div2]-[height_PF2 height_Div1 height_Div2]/2;   %[m]lowerZ of coilset (no SOl and PF1 (outside))
       Zup_coils=[Z_PF2 Z_Div1 Z_Div2]+[height_PF2 height_Div1 height_Div2]/2;    %[m]higher Z of coilset (no SOl and PF1 (outside))    
            %They are good
        
       %Lets rehape the meshgrid to do the loop to remove points
       r_ins_VVL=reshape(r_ins_VVL,length(r_ins_VVL)^2,1);
       z_ins_VVL=reshape(z_ins_VVL,length(z_ins_VVL)^2,1);
       
    
        for co=1:length(Rin_coils) %at each iter, removes the grid points inside the coils
              StoreRZ=[0 0]; %initialization of stored grid points
            for i=1:length(r_ins_VVL) %Have to check each point
                Point=[r_ins_VVL(i) z_ins_VVL(i)]; %grid point to test            
                
                switch sign(Point(2)) %First lets check if Z><0
                    
                    case 1 %z>0
                
                    if Point(1)<Rin_coils(co) | Point(1)>Rout_coils(co) %R out of the coil
                                                                                                   %All Z are good
                               StoreRZ=[StoreRZ; Point]; %store of good points                                                                                           
                    
                    elseif  Point(1)>Rin_coils(co) | Point(1)<Rout_coils(co) %R inside of the coil
                            if Point(2)<Zdown_coils(co) | Point(2)>Zup_coils(co)  %Z out coil
                                StoreRZ=[StoreRZ; Point]; %store of good points
                            end
                    end
                        
                    case -1 %z<0
                
                    if Point(1)<Rin_coils(co) | Point(1)>Rout_coils(co) %R out of the coil
                                                                                                   %All Z are good
                               StoreRZ=[StoreRZ; Point]; %store of good points                                                                                        
                    
                    elseif  Point(1)>Rin_coils(co) | Point(1)<Rout_coils(co) %R inside of the coil
                            if Point(2)>-Zdown_coils(co) | Point(2)<-Zup_coils(co)  %Z out coil
                                StoreRZ=[StoreRZ; Point]; %store of good points
                            end
                    end  
                        
                end             
            end  
            StoreRZ=StoreRZ(2:end,:); %remove first row, the initialization one
            r_ins_VVL=StoreRZ(:,1); %to use the grid created on the following coil loop
            z_ins_VVL=StoreRZ(:,2); %to use the grid created on the following coil loop          
        end             
     
     %Lets reshape it again to do the contour plots later (both are
     %vectors)
     r_ins_VVL_contour=reshape(r_ins_VVL,floor(length(r_ins_VVL)/2),[]); %arbitrary reshape
     z_ins_VVL_contour=reshape(z_ins_VVL,floor(length(z_ins_VVL)/2),[]); %arbitrary reshape
                
        %Plot
                figure;
                subplot(1,2,1)
                plot(StoreRZ(:,1),StoreRZ(:,2),'r.')
                hold on
                plot(vessel)
                axis equal
                title(sprintf('iter %d',co))   
                subplot(1,2,2)
                plot(StoreRZ(:,1),StoreRZ(:,2),'r.')
                hold on
                plot(vessel)
                plot(coilset)
                axis equal
                title(sprintf('iter %d',co))    
 
    %%%End grid
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%         Begin integration         %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    L_max=10000;                             %[m] max L value for the integration; when L achieves
                                             %this value, the integration stops. 
                                             %Iter 55,85,86,81 achieves about 10,000, spending 4h
                                     
    event_colission=@(t,y) Colission(t,y,VesselRMaxInner,...
            VesselRMinInner,VesselZMaxInner,VesselZMinInner,L_max);   
   
     options = odeset('OutputFcn',@ode_progress_bar,'Events',event_colission,'AbsTol',1e-10,'RelTol',1e-6); 
     %I include a Fiesta funciton to show the progress of the ode
                               
   %Both integrators are programmed, so to swich between them will use a switch (xD)
   tic          %to measure time the intergration takes
    switch IntMethod %switch for the starting points (change the number of inputs)
            
         case  'Lp' %Phi as integration method
              %1) Starting points y0
                 y0=[0 0 0 0 0];                    %(r0,z0,L0,phi0,U0) starting points
                    %note it has to be r z L for using the same event function

                    for i=1:length(r_ins_VVL)        
                         points=[r_ins_VVL(i) z_ins_VVL(i) 0 0 0];      %r z L phi U        
                            %U(0)=0 (arbitrary)  
                         y0=[ y0; points];                          
                    end 
                %I have the additional point 0 0 0 form the begining, that i can remove
                %easily with
                y0=y0(2:end,:);
                
               %2)Independt variable values, t0
                 t_values=1000; %1000            %Max value of the independent variable
                 n_iter_t=300000000; %3000000 on s1-14                 %Integer, number of values for tSpan
                 tSpan=linspace(0,t_values,n_iter_t);            %the range of values of independant variable
                    %TOO LITTLE FOR s1-19, MOST LINES DO NOT COLLIDE NOR
                    %ACHIEVE LMAX
                    
                %3)Odefun
                [FieldsBreak, FieldsBreakNoEarth]=fields(Equil_Passive); %Extraction of the fields (Earth's field included)!
                odefun= @(Lp, rzLphiU) Field_LineIntegrator_Lp(Lp,rzLphiU,FieldsBreak.interpn.Br,...
                    FieldsBreak.interpn.Bz, FieldsBreak.interpn.Bphi);                
                
                %4) Integration
                        %%%%%%SINGLE FIELD LINE TRACER and plotter

                        %need to find i for the chosen R,Z value in r0_z0_L0_U0.
                        %I= 85 for a line inside, 49 for a max L outside, 152 for the
                        %outward arm (Z>0). 135 for the outward Z<0 line. 64 for the upper
                        %arm
        
                        i=73%652 %looked in the y0
                        [t_fieldline, y_fieldline]=ode45(odefun,tSpan,y0(i,:),options);    
    
                        %To save the last values of R,Z,L
                        L_single=y_fieldline(end,3);         %L      
        
                        %Plot of one of the line
                            figure;
                            plot3(y0(i,1)*cos(y0(i,4)),y0(i,1)*sin(y0(i,4)),y0(i,2),'k*','LineWidth',3)
                            hold on;
                            plot3(vessel);
                            plot3(coilset);
                            plot3(y_fieldline(:,1).*cos(y_fieldline(:,4)),y_fieldline(:,1).*sin(y_fieldline(:,4)),...
                                y_fieldline(:,2),'r.','LineWidth',3)
                                xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');  
                            plot3(y_fieldline(end,1).*cos(y_fieldline(end,4)),y_fieldline(end,1).*sin(y_fieldline(end,4)),...
                                y_fieldline(end,2),'g*','LineWidth',3)
                            title(sprintf('Field line starting at (R=%2.2f,Z=%2.2f)m, L=%3.2fm ',y0(i,1),y0(i,2),L_single))
                            set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
                            %legend('Starting point (Point with less Bpol)','Field line',...
                        %%%END ONE LINE TRACER
                          
                        for i=1:length(y0)
                            fprintf('Iter %d de %d',i,length(y0))
                            [t_fieldline, y_fieldline]=ode45(odefun,tSpan,y0(i,:),options);        %ode15s Carlos
    
                            %To save the last values of the dependants variables
                            %y_end = 5 columns containing R, Z, Lc, Phi and U
                            y_end(i,1)=y_fieldline(end,1);           %R
                            y_end(i,2)=y_fieldline(end,2);           %Z
                            y_end(i,3)=y_fieldline(end,3);           %Lc           
                            y_end(i,4)=y_fieldline(end,4);           %Phi
                            y_end(i,5)=y_fieldline(end,5);           %U
                        end                 
                        U_int=reshape(y_end(:,5),size(r_ins_VVL_contour,1),size(r_ins_VVL_contour,2));            
        
        case 'Phi' %Lp as integration method
              %1) Starting points y0
                y0=[0 0 0 0]; %r0,z0,L0,U0
                %note it has to be r z fro using the same event function

                for i=1:length(r_ins_VVL)        
                        points=[r_ins_VVL(i) z_ins_VVL(i) 0 0];  %r z L phi U        
                            %U(0)=0 (arbitrary)  
                        y0=[ y0; points];                          
                end
                %I have the additional point 0 0 0 form the begining, that i can remove
                %easily with
                y0=y0(2:end,:);
            
               %2)Independt variable values, t0
                 t_values=10000; %10 for debug                           %Cycles(toroidal turns)
                 n_iter_t=3000000;         %3000                         %Integer, number of values for tSpan
                 tSpan=linspace(0,2*pi*t_values,n_iter_t);    %the range of values of independant variable
               
               %3)Odefun 
                 odefun= @(phi, rzLU) Field_LineIntegrator(phi,rzLU,FieldsBreak.interpn.Br,...
                    FieldsBreak.interpn.Bz,FieldsBreak.interpn.Bphi);               
                %4) Integration
                        %%%%%%%%SINGLE FIELD LINE TRACER and plotter

                        %need to find i for the chosen R,Z value in r0_z0_L0_U0.
                        %I= 85 for a line inside, 49 for a max L outside, 152 for the
                        %outward arm (Z>0). 135 for the outward Z<0 line. 64 for the upper
                        %arm
        
                        i=472 %looked in the y0
                                    %652 line collides with lower PF2
                                    %672 for collision upper PF2
                                    %116 for colission with upper Div1
                                    %472 for coliision with lower wall (VV)
                        [t_fieldline, y_fieldline t_event y_event]=ode45(odefun,tSpan,y0(i,:),options);    
    
                        %To save the last values of R,Z,L
                        L_single=y_fieldline(end,3);         %L      
        
                        %Plot of one of the line
                        figure;
                        plot3(y0(i,1)*cos(t_fieldline(1)),y0(i,1)*sin(t_fieldline(1)),y0(i,2),'k*','LineWidth',3)
                        hold on;
                        plot3(vessel);
                        plot3(coilset);
                        plot3(y_fieldline(:,1).*cos(t_fieldline(:)),y_fieldline(:,1).*sin(t_fieldline(:)),...
                            y_fieldline(:,2),'r.','LineWidth',3)
                        xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');  
                        plot3(y_fieldline(end,1).*cos(t_fieldline(end)),y_fieldline(end,1).*sin(t_fieldline(end)),...
                            y_fieldline(end,2),'g*','LineWidth',3)
                        title(sprintf('Single field line integration phi, L=%3.2f m ',L_single))
                        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
                        %%%END ONE LINE TRACER
                
                        for i=1:length(y0)
                            fprintf('Iter %d de %d',i,length(y0))
                            [t_fieldline, y_fieldline]=ode45(odefun,tSpan,y0(i,:),options);        %ode15s Carlos
    
                            %To save the last values of the dependants variables
                            y_end(i,1)=y_fieldline(end,1);           %R
                            y_end(i,2)=y_fieldline(end,2);           %Z
                            y_end(i,3)=y_fieldline(end,3);           %L           
                            y_end(i,4)=y_fieldline(end,4);           %U
                        end 
                        U_int=reshape(y_end(:,4),size(r_ins_VVL_contour,1),size(r_ins_VVL_contour,2));
        
    end         
   Time_Run_Integrator=toc     %Time spent by the integrator
   
   L_int=reshape(y_end(:,3),size(r_ins_VVL_contour,1),size(r_ins_VVL_contour,2));
    
   %Calc of average L on the null region
        index= y0(:,1)<=max(get(sensor_btheta,'r')) & y0(:,1)>=min(get(sensor_btheta,'r')) &...
            y0(:,2)<=max(get(sensor_btheta,'z')) & y0(:,2)>=min(get(sensor_btheta,'z')); 
                    %Index of R,Z inside null
        L_int_row=y_end(:,3); %[m] L in row form
        L_int_null=mean(L_int_row(index))
        
   %%%Storing of the non colliding starting points
        %1) Non colliding with the VV
            %To store start points that do not collide: first I get the index of both R
            %and Z, but together, since they do not collide if oth R and Z are greater
            %than the min value, and lower than the greatest value
    
            RZ_store_index=y_end(:,1)<VesselRMaxInner & ...
                y_end(:,1)>VesselRMinInner & y_end(:,2)<VesselZMaxInner...
                & y_end(:,2)>VesselZMinInner; %100*5==> error, has to ve vector,, not matrix!!!
                        
            RZ_no_collide=[y0(RZ_store_index,1) y0(RZ_store_index,2)];    
            
       %2) Colliding with the VV
            RZ_no_collide_end=[y_end(RZ_store_index,1) y_end(RZ_store_index,2)];
                                %y_end of points that do not collide with
                                %the VV. Same size as RZ_no_collide
               %Those points have to be tested.
                R_no_collide_end=RZ_no_collide_end(:,1); %R of y_end
                Z_no_collide_end=RZ_no_collide_end(:,2); %Z of y_end
                
                %The corresponding starting points are
                R_no_collide=RZ_no_collide(:,1); % R of y0
                Z_no_collide=RZ_no_collide(:,2); % Z of y0
                
        for co=1:length(Rin_coils)              %at each iter, removes the colliding points
              StoreRZ_y0=[0 0];                        %initialization of stored start points
              StoreRZ_end=[0 0]; %initialization of stored ending points DEBUG!!!
              
            for i=1:length(R_no_collide_end)         %Have to check each ending point
                        %this points will be renewed as the coils are
                        %checked, so that it removes progressively the
                        %points
                Point=[R_no_collide_end(i) Z_no_collide_end(i)];        %ending point to test            
                Point_y0=[R_no_collide(i) Z_no_collide(i)];                 %corresponding starting point
                
                switch sign(Point(2)) %First lets check if Z><0
                    
                    case 1 %z>0
                
                    if Point(1)<Rin_coils(co) | Point(1)>Rout_coils(co) %R out of the coil
                                                                                                   %All Z are good
                               StoreRZ_y0=[StoreRZ_y0; Point_y0]; %store of good start points                                                                                           
                               StoreRZ_end=[StoreRZ_end; Point]; %store of good end points DEBUG!!  
                               
                    elseif  Point(1)>Rin_coils(co) | Point(1)<Rout_coils(co) %R inside of the coil
                            if Point(2)<Zdown_coils(co) | Point(2)>Zup_coils(co)  %Z out coil
                                
                                StoreRZ_y0=[StoreRZ_y0; Point_y0]; %store of good start points
                                StoreRZ_end=[StoreRZ_end; Point]; %store of good end points DEBUG!! 
                            end
                    end
                        
                    case -1 %z<0
                
                    if Point(1)<Rin_coils(co) | Point(1)>Rout_coils(co) %R out of the coil
                                                                                                   %All Z are good
                               StoreRZ_y0=[StoreRZ_y0; Point_y0]; %store of good start points                                                                                        
                               StoreRZ_end=[StoreRZ_end; Point]; %store of good end points DEBUG!! 
                               
                    elseif  Point(1)>Rin_coils(co) | Point(1)<Rout_coils(co) %R inside of the coil
                            if Point(2)>-Zdown_coils(co) | Point(2)<-Zup_coils(co)  %Z out coil
                                                              
                                StoreRZ_y0=[StoreRZ_y0; Point_y0]; %store of good start points
                                StoreRZ_end=[StoreRZ_end; Point]; %store of good end points DEBUG!! 
                            end
                    end  
                        
                end             
            end
            %y0_points
            StoreRZ_y0=StoreRZ_y0(2:end,:); %remove first row, the initialization one
            R_no_collide=StoreRZ_y0(:,1); %to store R of y0 whose yend do not collide with coil 
            Z_no_collide=StoreRZ_y0(:,2); %to store  Z of y0 whose yend do not collide with coil 
            %y_end points
            StoreRZ_end=StoreRZ_end(2:end,:); %remove first row, the initialization one
            R_no_collide_end=StoreRZ_end(:,1); %to store R of yend that do not collide with coil 
            Z_no_collide_end=StoreRZ_end(:,2); %to store  Z of yend that do not collide with coil            
            %Note that after the last loop, this values will be the final
            %values!!!
        end                     
           RZ_no_collide=[R_no_collide Z_no_collide];
   
           
            %However, this is not perfect, when including in the grid the points in the wall,
            %something extrange happens, some of them are store in the non colliding points
            %thought they do not collide since the starting point is also the ending points
            %(you get like some stars just in the VV, but not all, only some of them)
            %To remove it:
            Index=RZ_no_collide(:,1)<VesselRMaxInner & ...
                RZ_no_collide(:,1)>VesselRMinInner & RZ_no_collide(:,2)<VesselZMaxInner...
                & RZ_no_collide(:,2)>VesselZMinInner;
            RZ_no_collide=[RZ_no_collide(Index,1) RZ_no_collide(Index,2)];                

     %%
     %%%Contour plots
      %1) L
        figure;
        contourf(r_ins_VVL_contour,z_ins_VVL_contour,L_int,'ShowText','On')
        %surf(r_ins_VVL_contour,z_ins_VVL_contour,L_int,'EdgeColor','none'), shading('interp')
        view(2)
        hold on
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'g.--','LineWidth', 2)
        plot(RZ_no_collide(:,1),RZ_no_collide(:,2),'m*')
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        colormap(Gamma_II)
        c=colorbar; %colorbar
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        ylabel(c, 'L(m)');
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('L','Field null region','Non colliding points')
        %title(sprintf('L  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))   
        %title(sprintf('L at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        title(sprintf('L at t=%d ms ph1',time_loop(loop)*1e3))        
        Filename = 'L_cont';
        %Filename= sprintf('%s_simu_%d',Filename,sen);     
        saveas(gcf, strcat(FigDir,ProjectName,Filename,FigExt));        

        %1D to 2D plot
        figure;
        tri = delaunay(y0(:,1),y0(:,2));
        plot(y0(:,1),y0(:,2),'.')
        [r,c] = size(tri); %number of triangles there
        trisurf(tri, y0(:,1), y0(:,2),y_end(:,3),'FaceAlpha',1,'EdgeColor','none'), shading('interp');
        view(2)
        colorbar
        hold on                     %this is to make the transition between values continuous,                                                        %instedad of discontinuously between pixels
        colormap(Gamma_II)
        %colormap(plasma)
        plot3([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],...
            ones(1,5)*max(y_end(:,3)),'g.--','LineWidth', 2)
        plot3(RZ_no_collide(:,1),RZ_no_collide(:,2),max(y_end(:,3))*ones(length(RZ_no_collide(:,2)),1),'m*')
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        c=colorbar; %colorbar
        %axis equal
        set(gca,'XLim',[0.1 1]);    
        set(gca,'YLim',[-0.9 0.9]);    
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        ylabel(c, 'L(m)');
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('L','Field null region','Non colliding points')
        %title(sprintf('L  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))   
        %title(sprintf('L at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        title(sprintf('L at t=%d ms ph1',time_loop(loop)*1e3))   
        Filename = 'L';        
        saveas(gcf, strcat(FigDir,ProjectName,Filename,FigExt));   
        print(gcf,strcat(FigDir,ProjectName,Filename,'.eps'),'-depsc2') 
        
        figure;
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'g.--','LineWidth', 2)
        hold on
        grid on
        plot(RZ_no_collide(:,1),RZ_no_collide(:,2),'m*')  %non colliding points
        plot(StoreRZ(:,1),StoreRZ(:,2),'r.') %grid points
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);        
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');      
        %axis equal
        set(gca,'XLim',[0.1 1]);    %To plot upper side of the VV
        set(gca,'YLim',[-0.9 0.9]);    %To plot upper side of the VV
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('Field null region','Non colliding points')
        %title(sprintf('L  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))   
        %title(sprintf('L at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        %title(sprintf('Non colliding at t=%d ms ph1',time_loop(loop)*1e3))   
        title('Cross-section')
        Filename = 'Grid_tracing';                
        saveas(gcf, strcat(FigDir,ProjectName,Filename,FigExt));   
        print(gcf,strcat(FigDir,ProjectName,Filename,'.eps'),'-depsc2')
        
      %2) Pseudo potential U/Vloop
        figure;
        contourf(r_ins_VVL_contour,z_ins_VVL_contour,U_int,10)
        %surf(r_insVV_noLimit,z_insVV_noLimit,U_int), shading('interp')
        hold on
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'U/Vloop');
        view(2) %2D view
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'g.--','LineWidth', 2)
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        xlabel('R (m)')
        ylabel('Z (m)')
        %title(sprintf('Pseudo potential  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))          
        %title(sprintf('Pseudo potential at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        title(sprintf('Pseudo potential at t=%d ms ph1',time_loop(loop)*1e3)) 
        Filename = 'Pseudo_contour';
        %Filename= sprintf('%s_simu_%d',Filename,sen);     
        saveas(gcf, strcat(FigDir,ProjectName,Filename,FigExt)); 
        print(gcf,strcat(FigDir,ProjectName,Filename,'.eps'),'-depsc2')
        %1D to 2D plot
                %1D to 2D plot
        figure;
        tri = delaunay(y0(:,1),y0(:,2));
        plot(y0(:,1),y0(:,2),'.')
        [r,c] = size(tri); %number of triangles there
        
        switch IntMethod
            
            case 'Phi'
                trisurf(tri, y0(:,1), y0(:,2),y_end(:,4),'FaceAlpha',1,'EdgeColor','none'), shading('interp');
            
            case 'Lp'
                trisurf(tri, y0(:,1), y0(:,2),y_end(:,5),'FaceAlpha',1,'EdgeColor','none'), shading('interp');
        end
        view(2)
        colorbar
        hold on                     %this is to make the transition between values continuous,                                                        %instedad of discontinuously between pixels
        colormap(Gamma_II)
        plot3([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],...
            ones(1,5)*max(y_end(:,3)),'g.--','LineWidth', 2)
        plot3(RZ_no_collide(:,1),RZ_no_collide(:,2),max(y_end(:,3))*ones(length(RZ_no_collide(:,2)),1),'m*')
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        colormap(Gamma_II)
        c=colorbar; %colorbar
        %axis equal
        set(gca,'XLim',[0.1 1]);    %To plot upper side of the VV
        set(gca,'YLim',[-0.9 0.9]);    %To plot upper side of the VV
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        ylabel(c, 'U/V');
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('U/Vloop','Field null region','Non colliding points')
        %title(sprintf('Pseudo potential  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))          
        %title(sprintf('Pseudo potential at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        title(sprintf('U/V_{loop} at t=%d ms ph1',time_loop(loop)*1e3))
        Filename = 'Pseudo';        
        saveas(gcf, strcat(FigDir,ProjectName,Filename,FigExt));   
        print(gcf,strcat(FigDir,ProjectName,Filename,'.eps'),'-depsc2')         
        
      %%%3)[Experimental] E_rel plot, to predict where the gas breaks down
        
        p_test=1*10^-4; %Tor
        E_RZmin=C_2(1)*p_test./(log(C_1(1)*p_test*L_int)); %E min, Paschen, but 2D        
        E_RZmin(E_RZmin<0)=NaN; %when Emin<0, there is no breakdwon, so NaN not
                    %to plot it
        
        figure;
        contourf(r_ins_VVL_contour,z_ins_VVL_contour,U_int./L_int*V_loop./E_RZmin)
        %surf(r_insVV_noLimit,z_insVV_noLimit,U_int), shading('interp')
        hold on
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        colormap(Gamma_II)
        c=colorbar; %colorbar
        ylabel(c, 'E_{rel}');
        view(2) %2D view
        plot([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],'g.--','LineWidth', 2)
        axis equal
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        xlabel('R (m)')
        ylabel('Z (m)')
        %title(sprintf('Pseudo potential  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))          
        title(sprintf('E_{rel} at t=%dms for p= %1.2d Tor',time_loop(loop)*1e3,p_test))
        Filename = 'U_L_contour';
        %Filename= sprintf('%s_simu_%d',Filename,sen);   
        saveas(gcf, strcat(FigDir,ProjectName,Filename,FigExt)); 
        print(gcf,strcat(FigDir,ProjectName,Filename,'.eps'),'-depsc2')
        %1D to 2D plot
        
                p_test=1*10^-4; %Tor
        E_RZmin=C_2(1)*p_test./(log(C_1(1)*p_test*y_end(:,3))); %E min, Paschen, but 2D        
        E_RZmin(E_RZmin<0)=NaN; %when Emin<0, there is no breakdwon, so NaN not
                    %to plot it
          figure;
        tri = delaunay(y0(:,1),y0(:,2));
        plot(y0(:,1),y0(:,2),'.')
        [r,c] = size(tri); %number of triangles there        
                switch IntMethod
            
            case 'Phi'
                trisurf(tri, y0(:,1), y0(:,2),y_end(:,4)./y_end(:,3)*V_loop./E_RZmin,'FaceAlpha',1,'EdgeColor','none'), shading('interp');
            
            case 'Lp'
                trisurf(tri, y0(:,1), y0(:,2),y_end(:,5)./y_end(:,3)*V_loop./E_RZmin,'FaceAlpha',1,'EdgeColor','none'), shading('interp');
                end                         
        view(2)
        colorbar
        hold on                     %this is to make the transition between values continuous,                                                        %instedad of discontinuously between pixels
        colormap(Gamma_II)
        %colormap(plasma)
        plot3([min(r_sensors) min(r_sensors) max(r_sensors) max(r_sensors) min(r_sensors)],...
            [min(z_sensors) max(z_sensors) max(z_sensors) min(z_sensors) min(z_sensors)],...
            ones(1,5)*max(y_end(:,3)),'g.--','LineWidth', 2)
        hh=plot(vessel);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        hh=plot(coilset);
        set(hh, 'EdgeColor', 'k')
        set(hh,'HandleVisibility','off');
        c=colorbar; %colorbar
        %axis equal
        set(gca,'XLim',[0.1 1]);    %To plot upper side of the VV
        set(gca,'YLim',[-0.9 0.9]);    %To plot upper side of the VV
        set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
        ylabel(c, 'E_{rel}');
        xlabel('R (m)')
        ylabel('Z (m)')
        %legend('L','Field null region','Non colliding points')
        %title(sprintf('L  at t=%d ms (iter %d/%d)',time_loop(loop)*1e3,loop,length(time_loop)))   
        %title(sprintf('L at t=%dms for simu %d',time_loop(loop)*1e3,sen))
        title(sprintf('E_{rel} at t=%d ms ph1',time_loop(loop)*1e3))   
        Filename = 'U_L';        
%}
%%


% Daniel's Functions Pertaining To The Field Line Integrator
%%
%{
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Extraction of the fields@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%1)Field of things
function [Field_Earth Field_NoEarth]=fields(equil)
%function [Field_grid, Field_VV, Field_interpol, Field_sensor]=fields(equil)
    %beware of the sensors, where they are not defined (yet), its field must not
    %be called!
    %Will add the Earths field! Will be useful when the fields are low, for
    %example at breakdown.
    
    global R_in Z_in R_sensor Z_sensor %to obtain global variables

    %Fiesta fields:
    Br=get(equil,'Br'); 
    Bz=get(equil,'Bz'); 
    Bphi=get(equil,'Bphi_vac'); %not alwais
    
    %Fields data (data in all the grid)
    Br_data = get(get(equil,'Br'),'data');     %Note this is 200*251, GridR*GridZ dimension
    Bz_data = get(get(equil,'Bz'),'data');
    Bphi_data = get(get(equil,'Bphi_vac'),'data');
            %Was Bphi, not always is vac, although the plasma field is
            %negligible in comparison to TF coils field.

    %%2D reshape of the data, to interpolate things
    
    rGrid = get(get(get(equil,'Br'),'grid'),'r');
    zGrid = get(get(get(equil,'Br'),'grid'),'z');
    
    Br_data = reshape( Br_data, length(zGrid), length(rGrid)); %251*200
    Bz_data = reshape( Bz_data, length(zGrid), length(rGrid)); %251*200
    Bphi_data = reshape( Bphi_data, length(zGrid), length(rGrid)); %251*200
    Bpol_data=sqrt(Br_data.^2+Bz_data.^2);
   
    
    %Interpolation vectors
    Br_interp = @(r,z) interpn(zGrid,rGrid,Br_data,z,r,'mikama');
    Bz_interp = @(r,z) interpn(zGrid,rGrid,Bz_data,z,r,'mikama');
    Bphi_interp = @(r,z) interpn(zGrid,rGrid,Bphi_data,z,r,'mikama');
    
    %Finally, the fields inside the vessel are
    Br_VV=Br_interp(R_in,Z_in);
    Bz_VV=Bz_interp(R_in,Z_in);
    Bphi_VV=Bphi_interp(R_in,Z_in);
    Bpol_VV=sqrt(Br_VV.^2+Bz_VV.^2);
    
    %And inside the sensor region:
    Br_sens=Br_interp(R_sensor,Z_sensor);
    Bz_sens=Bz_interp(R_sensor,Z_sensor);
    Bphi_sens=Bphi_interp(R_sensor,Z_sensor);
    Bpol_sens=sqrt(Br_sens.^2+Bz_sens.^2);
    
    %To store them, will create several structures, one for the data, other
    %for the VV data and other for the sensor data
    
    Field_grid.Br=Br_data;
    Field_grid.Bz=Bz_data;
    Field_grid.Bphi=Bphi_data;
    
    Field_VV.Br=Br_VV;
    Field_VV.Bz=Bz_VV;
    Field_VV.Bphi=Bphi_VV;
    Field_VV.Bpol=Bpol_VV;    

    Field_sensor.Br=Br_sens;
    Field_sensor.Bz=Bz_sens;
    Field_sensor.Bphi=Bphi_sens;
    Field_sensor.Bpol=Bpol_sens;
    
    Field_interpol.Br=Br_interp;
    Field_interpol.Bz=Bz_interp;
    Field_interpol.Bphi=Bphi_interp;
    
    %Final structure englobating all
    Field_NoEarth.grid=Field_grid;
    Field_NoEarth.VV= Field_VV;
    Field_NoEarth.sensor= Field_sensor;    
    Field_NoEarth.interpn=Field_interpol;
    
        %%%%%Earths field%%%%%%%%%%%%%
    
    %In seville, with coordinates 37º23'N 5°59′00″W (Wikipedia, spanish)=
    %37+23/60ºN 5+59/60º W=37.3833ºN 5.9833ºW, the components are
    
    X_Earth=27316.6e-9;                 %[T] N-S component, >0 for N, <0 for S
    Y_Earth=-423.4e-9;                      %[T] E-W component, <0 for W, >0 for E
    Z_Earth=33833.4e-9;                 %[T] vertical component, <0 for U, >0 for D
    
    %To create the vectors, I have problems for the r and phi directions, since its
    %magnitude vary when moving the toroidal angle, because Fiesta considers
    %axysymmetry, but Earths field is not axysymmetric. What I will do as an
    %approximation is take the average value of the field over all the phi angles. 
    %This value, the mean, will be used for the r and phi components. 
    %THe z components is norproblematic
    
        BrEarth=0; %initialization
        BphiEarth=0; %initialization
        ang_Earth=linspace(0,2*pi,100);
        
        for i=1:length(ang_Earth)
            %Addition at each step (debug)
            add_Br(i)=X_Earth*cos(ang_Earth(i))+(-Y_Earth)*sin(ang_Earth(i));       
            add_Bphi(i)=-X_Earth*sin(ang_Earth(i))+(-Y_Earth)*cos(ang_Earth(i));
                                                %-Y because Seville is in the West (W)
            BrEarth=BrEarth+add_Br(i);
            BphiEarth=BphiEarth+add_Bphi(i);
            
        end
         
        %The r and phi components are the mean:       
        BrEarth=BrEarth/length(ang_Earth); %mean
        BphiEarth=BphiEarth/length(ang_Earth); %mean      
        
        %The vertical component is
        BzEarth=-Z_Earth;
        BpolEarth=sqrt(BrEarth^2+ BzEarth^2);
    
        %Now lets create the a grid with this constant field values;
       
        BzEarth=BzEarth*ones(length(zGrid), length(rGrid));        
        BrEarth= BrEarth*ones(length(zGrid), length(rGrid));
        BphiEarth=BphiEarth*ones(length(zGrid), length(rGrid));    
        BpolEarth=sqrt(BrEarth.^2+ BzEarth.^2);
    
        
    %%%%%%%%END EARTHS FIELD%%%%%%%%
    %Now will create the same as above but with Earths field
    
    %Interpolation vectors 
    Br_interp = @(r,z) interpn(zGrid,rGrid,Br_data+BrEarth,z,r,'mikama');
    Bz_interp = @(r,z) interpn(zGrid,rGrid,Bz_data+BzEarth,z,r,'mikama');
    Bphi_interp = @(r,z) interpn(zGrid,rGrid,Bphi_data+BphiEarth,z,r,'mikama');
    
    %Finally, the fields inside the vessel are
    Br_VV=Br_interp(R_in,Z_in);
    Bz_VV=Bz_interp(R_in,Z_in);
    Bphi_VV=Bphi_interp(R_in,Z_in);
    Bpol_VV=sqrt(Br_VV.^2+Bz_VV.^2);
    
    %And inside the sensor region:
    Br_sens=Br_interp(R_sensor,Z_sensor);
    Bz_sens=Bz_interp(R_sensor,Z_sensor);
    Bphi_sens=Bphi_interp(R_sensor,Z_sensor);
    Bpol_sens=sqrt(Br_sens.^2+Bz_sens.^2);
    
    %To store them, will create several structures, one for the data, other
    %for the VV data and other for the sensor data
    
    Field_grid.Br=Br_data;
    Field_grid.Bz=Bz_data;
    Field_grid.Bphi=Bphi_data;
    
    Field_VV.Br=Br_VV;
    Field_VV.Bz=Bz_VV;
    Field_VV.Bphi=Bphi_VV;
    Field_VV.Bpol=Bpol_VV;    

    Field_sensor.Br=Br_sens;
    Field_sensor.Bz=Bz_sens;
    Field_sensor.Bphi=Bphi_sens;
    Field_sensor.Bpol=Bpol_sens;
    
    Field_interpol.Br=Br_interp;
    Field_interpol.Bz=Bz_interp;
    Field_interpol.Bphi=Bphi_interp;
    
        %Final structure englobating all
    Field_Earth.grid=Field_grid;
    Field_Earth.VV= Field_VV;
    Field_Earth.sensor= Field_sensor;    
    Field_Earth.interpn=Field_interpol;
    
    %Yes, everything is all right, you are not rewritting things
end


%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Odefun for field line tracing LP@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%Field line integrator function with Lp
     %this solves the field line eq, using Lp (poloidal length) 
         %as an independent variable
     %way more faster than with phi (5 min when Lmax=3000m, 15 inside points)
        
    function [results]=Field_LineIntegrator_Lp(Lp,rzLphiU,Br_interpn,Bz_interpn,Bphi_interpn)
    %rzphiLU=[r z phi L U]
    %Lp= poloidal length (have to write capital L so it not apperas as
    %internsity I). just tchange the variables in the phi equations
    
    %First, the field needs to be evaluated at the point (r,phi,z):
    
    Br_eval=Br_interpn(rzLphiU(1),rzLphiU(2));
    Bphi_eval=Bphi_interpn(rzLphiU(1),rzLphiU(2));
    Bz_eval=Bz_interpn(rzLphiU(1),rzLphiU(2)); 
    Bpol_eval=sqrt(Br_eval^2+Bz_eval^2);
    
    %With the field, the eq to solve is:
    
    dr_dLp=Br_eval/Bpol_eval;
    dphi_dLp=rzLphiU(1)*Bphi_eval/Bpol_eval;
    dz_dLp=Bz_eval/Bpol_eval;
    dlength_dLp=sqrt(1+(Bphi_eval/Bpol_eval)^2);
    dU_Vloop_dLp=1/(2*pi*rzLphiU(1))*dlength_dLp; %pseudo potential U/V_loop
    
    results=zeros(5,1); %column vector to group the results
    results(1)=dr_dLp;
    results(4)=dphi_dLp;
    results(2)=dz_dLp;
    results(3)=dlength_dLp;
    results(5)=dU_Vloop_dLp;
    
    end
    
%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Odefun for field line tracing PHI@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
     function [results]=Field_LineIntegrator(phi,rzLU,Br_interpn,Bz_interpn,Bphi_interpn)
    %rzL=[r z L U]
   
    %First, the field needs to be evaluated at the point (r,phi,z):
    
    Br_eval=Br_interpn(rzLU(1),rzLU(2));
    Bphi_eval=Bphi_interpn(rzLU(1),rzLU(2));
    Bz_eval=Bz_interpn(rzLU(1),rzLU(2));    
    Bpol_eval=sqrt(Br_eval^2+Bz_eval^2);
    
    %With the field, the eq to solve is:
    
    dr_dphi=rzLU(1)*Br_eval/Bphi_eval;
    dz_dphi=rzLU(1)*Bz_eval/Bphi_eval;
    length=rzLU(1)*sqrt(Bphi_eval^2+Bpol_eval^2)/Bphi_eval;
    U_Vloop=1/(2*pi*rzLU(1)); %pseudo potential U/V_loop
    
    results=zeros(4,1); %column vector to group the results
    results(1)=dr_dphi;
    results(2)=dz_dphi;
    results(3)=length;
    results(4)=U_Vloop;
    
  end	
    
%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@Event function for the ode@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   %Only one function can be introduced into ode, so this event has to have
   %all the conditions
   
   function [rz_value isterminal direction]=Colission(phi,rzL,...
       VesselRMaxPoint,VesselRMinPoint,VesselZMaxPoint,VesselZMinPoint,L_max) 
   
    global Rin_coils Rout_coils Zdown_coils Zup_coils %limits of coils inside VV
   
    %when rz_value is zero, stop the integration. There are two possibles
    %coliisions, with the VV and with the coils
   
    %1) Wall colission
   
        %To implement the 4
        %possible colission, we could do like the product of each, since when one of
        %them is zero, the product will be zero, and also to define row vectors
        %for isterminal, direction, and rz_value. Will do this second option
      
        up_colission=rzL(2)-VesselZMaxPoint;            %colission with upper wall
        down_colission=rzL(2)-VesselZMinPoint;      %colission with lower vall
        out_colission=rzL(1)-VesselRMaxPoint;           %colission with outer
        in_colission=rzL(1)-VesselRMinPoint;            %colission with inner
   
        %Max condition of L, to stop the integration  
        L_lim=rzL(3)-L_max;                                                  %[m] Maximum L value
   
        rz_value_VV=[up_colission down_colission out_colission in_colission L_lim];
   
        %Have checked that if I dont use the minR condition, it will impige 
        %in the upper wall, which was was happens at the beggining, when
        %I dont have the inner wall condition
    
    %2) Colission with the coils
               Point=[rzL(1) rzL(2)]; %point of the line
          for coo=1:length(Rin_coils)   %Iter for the coils
            switch sign(Point(2)) %z><0
                
                case 1 %z>0
                    rz_value_Coil(coo)= Point(1)>=Rin_coils(coo) & Point(1)<=Rout_coils(coo) & Point(2)>=Zdown_coils(coo) & Point(2)<=Zup_coils(coo); 
                                    %this value is 0 if the point is inside
                                    %the coil
                                    
                case -1        %z< 0
                    rz_value_Coil(coo)= Point(1)>=Rin_coils(coo) & Point(1)<=Rout_coils(coo) & Point(2)<=-Zdown_coils(coo) & Point(2)>=-Zup_coils(coo); 
                                    %this value is 0 if the point is inside
                                    %the coil
            end
          end
          
          rz_value=[rz_value_VV rz_value_Coil]; %contain both conditions, VV and coils
          isterminal=ones(1,length(rz_value));                   %to stop the integration
          direction=zeros(1,length(rz_value));                    %to not worry about the sloop   
        %works for coil colliding with lower PF2 (652)
        %works for coil colliding with upper PF2 (672)
        %Works for upper Div1(116)
        %Works for colission with lower VV wall (472) ==>work correctly
   end
%}
%%


