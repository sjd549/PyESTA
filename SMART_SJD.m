%% SMall Aspect Ratio Tokamak (SMART), V3p1 init.nam

clear 
clc
close all

%%%%%%%%%%

%Add FIESTA trunk path, include path to any extra functions.
FIESTATrunk = "~/Postdoc Seville/FIESTA/Source Code/FIESTA_V8.8";
FunctionsTrunk = "../../Functions";
addpath(genpath(FIESTATrunk),genpath(FunctionsTrunk));

%Set maximum number of concurrent CPU threads in use
NumThreads = 2;
NumThreads = maxNumCompThreads(NumThreads);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        DEFINE REACTOR GEOMETRY                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define global scaling factors (Vessel Walls Are Not Scaled!)
RScaleVessel=1.00;    %Scale all radial vessel dimensions (Relative to 2.0m)
ZScaleVessel=1.00;    %Scale all axial vessel dimensions (Relative to 2.0m)
RScaleCoil=1.00;      %Scale all radial coil positions (Relative to 2.0m)
ZScaleCoil=1.00;      %Scale all axial coil positions (Relative to 2.0m)

%Define Vessel Wall Thickness
ww_R=0.015;		%Radial Wall Thickness [m]
ww_Z=0.015;		%Axial Wall Thickness [m]

%Define Vessel Internal Geometry (Does not include wall thickness)
VesselRMinInner=0.15*RScaleVessel; % R min position [m]
VesselRMaxInner=0.8*RScaleVessel;  % R max position [m]
VesselZMinInner=-0.8*ZScaleVessel; % Z min position [m]
VesselZMaxInner=0.8*ZScaleVessel;  % Z max position [m]

%Define center points of vessel walls (Inner Geometry + half wall thickness)
ZMinCenter=VesselZMinInner-(ww_Z/2);
ZMaxCenter=VesselZMaxInner+(ww_Z/2);
RMinCenter=VesselRMinInner-(ww_R/2);
RMaxCenter=VesselRMaxInner+(ww_R/2);


%%%%%%%%%%%%%%%%%%%%%%%  DEFINE COIL GEOMETRY  %%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Solenoid Geometry and Parameters
nSol=210;					 	% Number of Solenoid Windings
RSol=0.13-(ww_R/2);			 	% R position of the solenoid [m]	%Inner:0.09, Outer:0.13
ZMinSol=ZMinCenter-(ww_Z/2);	% Solenoid Min Z position			%-0.8*ZScaleVessel;
ZMaxSol=ZMaxCenter+(ww_Z/2);	% Solenoid Max Z position			%+0.8*ZScaleVessel;

%Number of Radial (R) and axial (Z) PF coil windings
nZDiv1=6;
nRDiv1=4;
nZDiv2=6;
nRDiv2=4;
nZPF1=6;
nRPF1=4;
nZPF2=6;
nRPF2=4;

%Define coil turn dimensions to enable cross-section calculation
width_PF=0.042;  % Width of a turn (m)
height_PF=0.035; % Height of a turn (m)

%Define central location of coil sets
R_PF1=(0.9*RScaleCoil)+ww_R/2;   %R position of PF1 (m)			%0.90m
Z_PF1=(0.3*ZScaleCoil)+ww_Z/2;   %Z Position of PF1 (m)			%0.30m
R_PF2=(0.9*RScaleCoil)+ww_R/2;   %R Position of PF2 (m)			%0.90m
Z_PF2=(0.6*ZScaleCoil)+ww_Z/2;   %Z Position of PF2 (m)			%0.60m
R_Div1=(0.25*RScaleCoil)+ww_R/2; %R Position of Div1 (m)		%0.25m
Z_Div1=(0.85*ZScaleCoil)+ww_Z/2; %Z Position of Div1 (m)		%0.85m
R_Div2=(0.55*RScaleCoil)+ww_R/2; %R Position of Div2 (m)		%0.55m
Z_Div2=(0.85*ZScaleCoil)+ww_Z/2; %Z Position of Div2 (m)		%0.85m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     DEFINE OPERATIONAL PARAMETERS                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%  DEFINE INITIAL PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%

%Define any required constants
mu0 = 1.2566e-06; % Magnetic Moment      [I/m^2]

%Define Operating Conditions
Te = 250;         % Electron Temperature [eV]
Ti = Te*0.1;      % Ion Temperature      [eV]
BT = 0.1;         % Toroidal B-Field     [T] (Defined at Rgeo)
Ip = 30e3;        % Plasma current       [A]
TauR = 0.020;     % Ramp Timescale       [s] (Also determines tstep for Ip plot)
TauP = 0.030;     % Pulse Timescale      [s]
RGeo = 0.450;     % Geometrical Radius   [m]
ZGeo = 0.000;     % Geometrical Axis     [m]
RSep = 0.700;     % Separatrix Radius    [m]
a = RSep-RGeo;    % Minor Radius         [m] (~0.25)
Epsilon = RGeo/a; % Aspect ratio         [-] (~1.8)
Kappa = 1.8;      % Elongation           [-]
delta = 0.20;     % Triangularity        [-] (~0.2)
li2 = 1;          % Standard Value?      [-]
%q_cyl = 2.821;   % Safety Factor?       [-]
%betaN = 3.529;   % Normalised Beta      [%] (Obtained via VEST Excel - (2X TOO HIGH)

%Compute Further Operating Conditions
Gr_Limit = 1e20*(Ip*1e-6/(pi*a^2*Kappa));  % Greenwald Limit          [m-3]
Gr_Frac = 0.15;                            % Greenwald Fraction       [-]
ne = Gr_Limit*Gr_Frac;                     % Electron Density         [m-3]  ~3E19
Irod = BT*2*pi*RGeo/mu0;                   % Central Rod Current      [A]
S = sqrt( (1.0+Kappa^2)/2.0 );             % Shaping factor           [-]
%deltaUp = (RGe-Rup)/a;                    % Upper-Triangularity      [-]
%deltaLo = (RGe-Rlo)/a;                    % Lower-Triangularity      [-]
%delta = (deltaUp+deltaLo)/2.0;            % Triangularity            [-]
%betaN = (betaT*BT*a)/(Ip*1e-6*mu0)        % Normalised Beta          [%] 
%betaT = (betaN/a*(Ip*1e-6))/BT;           % Beta toroidal            [%]
betaP = 3/2*ne*(Te+Ti)/(mu0*Ip/(2*pi*a))^2*2*mu0*1.6e-19*Kappa; % Beta Poloidal  [%]
BZ = -mu0*Ip/(4*pi*RGeo)*(log(8*Epsilon)+betaP+0.5*li2-3/2);    % Vertical field [T]

%Coil density, temperature and resistivity
coil_density = 1;                       % Relative Coil Density      [Arb]
coil_temp = 293.0;                      % Initial Coil Temperature   [K]
resistivity = copper_resistivity_at_temperature(coil_temp);

%Gas species analouge - H=1, He=2, Ar=11.85 (for Te < 280eV)
%https://www.webelements.com/argon/atoms.html
Z_eff=1.0;                              % Effective Nuclear Charge   [e-]


%%%%%%%%%%%%%%%%%%  DEFINE SOL RAMP & COIL CURRENTS  %%%%%%%%%%%%%%%%%%%%

%Define number of time-steps (vertices) in the current waveforms
nTime = 6;      %[Steps]
tstep = TauR;	%[s]
time =  [-0.10 -0.05 0 tstep tstep+TauP tstep+0.05];    %Phase1_Daniel
%time = [-0.11 -0.05 0 tstep tstep+0.10 tstep+0.11];    %Phase2_JuanJo

%!!!WOULD BE NICE TO IMPLIMENT CURRENT WAVEFORM ARRAYS!!!
%ISol_Waveform = [+1300, 0000,-1300];
%IPF1_Waveform = [-370];
%IPF2_Waveform = [-400];
%IDiv1_Waveform = [+000];
%IDiv2_Waveform = [+900];

%Solenoid coil currents [kA]	%H+			%He2+
I_Sol_Start=+900;				%+0900;		%+1000
I_Sol_Equil=-500;				%-0500;		%-0500
I_Sol_End=-900;					%-0900;		%-1000
%Symmetric ISol is better for power supply

%PF coil currents (For Equilibrium)
I_PF1=-390;						%-390;		%-0000
I_PF2=-385;						%-385;		%-0000
I_Div1=+000;					%+000;		%+0000
I_Div2=+900;					%+900;		%+0000


%%%%%%%%%%%%%%%%%  DISPLAY VARIABLE OUTPUT TO USER  %%%%%%%%%%%%%%%%%%%%%

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
disp([ 'Minor Radius = ' num2str(a) ' [m]' ]);
disp([ 'AspectRatio = ' num2str(Epsilon) ' [-]' ]);
disp([ 'Elongation = ' num2str(Kappa) ' [-]' ]);
disp([ 'Shaping Factor = ' num2str(S) ' [-]' ]);
%disp([ 'Triangularity = ' num2str(delta) ' [-]' ]);

disp([ ' ' ]);
disp([ '%===== Initial Coil Currents =====%' ]);
disp([ 'I_Sol_Start = ' num2str(I_Sol_Start/1000) ' [kA]' ]);
disp([ 'I_Sol_Equil = ' num2str(I_Sol_Equil/1000) ' [kA]' ]);
disp([ 'I_Sol_End = ' num2str(I_Sol_End/1000) ' [kA]' ]);
disp([ 'I_PF1 = ' num2str(I_PF1/1000) ' [kA]' ]);
disp([ 'I_PF2 = ' num2str(I_PF2/1000) ' [kA]' ]);
disp([ 'I_Div1 = ' num2str(I_Div1/1000) ' [kA]' ]);
disp([ 'I_Div2 = ' num2str(I_Div2/1000) ' [kA]' ]);
disp([ ' ' ]);


%%%%%%%%%%%%%%%%%  DEFINE DATA OUTPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%

%Define figure extension
FigExt = '.png'; 		%'.png','.eps','.pdf'

%Define project and series names
ProjectName = 'SMARTxs-P1';			%Define global project name
SeriesName = 'VaryTauP';		%Define parameter scan series name

%Create global output folders for saved data and figures
ASCIIDir = 'RawData/'; mkdir(ASCIIDir);
%FigDir = 'Figures/'; mkdir(FigDir);			%NOT CURRENTLY USED

%Create simulation name based upon relevant run parameters
SimName = 'DefaultSimName';
disp([ 'SimName: ' SimName ]);
disp([ ' ' ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   INITIATE VESSEL AND COIL OBJECTS                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define FIESTA simulation grid limits and resolution
R_simulation_limits = [0.03*RScaleVessel 1.0*RScaleVessel];	%Must be greater than vessel
Z_simulation_limits = [-1.3*ZScaleVessel 1.3*ZScaleVessel];	%Must be greater than vessel
grid_size_R =200; 	%Arbitrary
grid_size_Z =251;	%Arbitrary

%Generate fiesta_grid object over which equilibrum simulation will be performed
Grid = fiesta_grid( R_simulation_limits(1), R_simulation_limits(2), grid_size_R, Z_simulation_limits(1), Z_simulation_limits(2), grid_size_Z );

%Extract vectors of R and Z grid points for use in further diagnostics
rGrid=get(Grid,'r'); %1*200
zGrid=get(Grid,'z'); %1*251
RGrid=get(Grid,'R'); %1*50200, 50200=251*250
ZGrid=get(Grid,'Z'); %1*50200, 50200=251*250


%%%%%%%%%%%%%%%%%%  INITIATE VACUUM VESSEL FILAMENTS  %%%%%%%%%%%%%%%%%%%

%pnts = [ VesselRMinInner VesselZMaxInner VesselRMaxInner VesselZMaxInner];
%pnts2 = [pnts; [pnts(end:-1:1,1),-pnts(end:-1:1,2)]; pnts(1,:)];    

%Approximate vessel geometry in 'pixelated' form through Bresenham algorithm
%ww =1.5e-2;         %Wall Thickness [m]
%xaccum = [];
%yaccum = [];
%for i=1:length(pnts2)-1
%    [xx,yy]=bresenham( pnts2(i,1)/ww,pnts2(i,2)/ww,pnts2(i+1,1)/ww,pnts2(i+1,2)/ww);
%    xaccum = [xaccum;xx*ww];
%   yaccum = [yaccum;yy*ww];
%end

%Define four vertices defined as the centre of each vessel corner
%Not currently used for calculation - would be good to use these to make linspace...
Vertice1=[RMaxCenter ZMaxCenter];	%Top Right
Vertice2=[RMaxCenter ZMinCenter];	%Bottom Right
Vertice3=[RMinCenter ZMinCenter];	%Bottom Left
Vertice4=[RMinCenter ZMaxCenter];	%Top Left

%Calculate number of filaments within vessel walls (Must be an even number)
n_fil_Z=round((ZMaxCenter-ZMinCenter+2*ww_Z)/ww_Z);	%Filaments in axial wall (Typically 8)
n_fil_R=round((RMaxCenter-RMinCenter+2*ww_R)/ww_R);	%Filaments in radial wall (Typically 8)
%Construct filament arrays for each section of vessel wall:
%Axial Outboard Wall, R=Rmax :: Top Right to Bottom Right (Vertice1 to Vertice2)
R_lin1_2=RMaxCenter*ones(1,n_fil_Z);
Z_lin1_2=linspace(ZMaxCenter,ZMinCenter,n_fil_Z);
%Radial Bottom Wall, Z=Zmin :: Bottom Right to Bottom Left (Vertice2 to Vertice3)
Z_lin2_3=ZMinCenter*ones(1,n_fil_R);
R_lin2_3=linspace(RMaxCenter,RMinCenter,n_fil_R);
%Axial Inboard Wall, R=Rmin :: Bottom Left to Top Left (Vertice3 to Vertice4)
R_lin3_4=RMinCenter*ones(1,n_fil_Z);
Z_lin3_4=linspace(ZMinCenter,ZMaxCenter,n_fil_Z);
%Radial Top Wall, Z=Zmax :: Top Left to Top Right (Vertice3 to Vertice4)
R_lin4_1=linspace(RMinCenter,RMaxCenter,n_fil_R);
Z_lin4_1=ZMaxCenter*ones(1,n_fil_R);

%Assemble the vessel wall filament position arrays
xaccum=[R_lin1_2 R_lin2_3 R_lin3_4 R_lin4_1]';
yaccum=[Z_lin1_2 Z_lin2_3 Z_lin3_4 Z_lin4_1]';
%Remove duplicate cells at wall corners (if any)
dup = (abs(diff(xaccum))+abs(diff(yaccum))) > 0;
xaccum = xaccum(dup);
yaccum = yaccum(dup);

%Construct vessel wall FIESTA filaments using position arrays
%Inputs(R,Z,2*r,2*z,1,0,0) where {R=MajorRadius, Z=Height, r=MinorRadius, z=MinorHeight}
for i=length(xaccum):-1:1
    vessel_filament(i) = fiesta_filament(xaccum(i),yaccum(i),ww_R,ww_Z,1,0,0);	%??? ww/3 ???
end
%Enable induced currents in vessel wall filaments - used only to calculate eddy currents
passive = fiesta_passive('STVesselPas',vessel_filament,'g');
vessel = fiesta_vessel( 'STVessel',passive);


%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIATE PF COILS  %%%%%%%%%%%%%%%%%%%%%%%%%%

%Define and initiate PF coils - Arbitrary numbering of coils
iSol = 1;       %Central Inducting Solenoid
iPF1 = 2;       %Upper Plasma Forming Coil
iPF2 = 3;       %Lower Plasma Forming Coil
iDiv1 = 4;      %Inboard Divertor Coil
iDiv2 = 5;      %Outboard Divertor Coil

%Calculate total number of windings in each coil - Windings defined in coil inputs
nDiv1=nZDiv1*nRDiv1;
nDiv2=nZDiv2*nRDiv2;
nPF1=nZPF1*nRPF1;
nPF2=nZPF2*nRPF2;

%Create array containing number of coil windings - Used to generate coil objects
turns=[];
turns(iSol) = nSol; 
turns(iDiv1) = nDiv1;
turns(iDiv2) = nDiv2;
turns(iPF1) = nPF1;
turns(iPF2) = nPF2;
nPF = 5; 				%Total number of coils including solenoid

%Create coil set from parameters defined above. (Function made by Carlos Soria)
%Function createVESTPFCircuit creates two PF coils. One in (R, Z) and another in (R, -Z)
PF1  = createVestPFCircuit('PF1',R_PF1,Z_PF1,width_PF,height_PF,turns(iPF1),nZPF1,nRPF1,true, coil_temp, resistivity, coil_density);
PF2  = createVestPFCircuit('PF2',R_PF2,Z_PF2,width_PF,height_PF,turns(iPF2),nZPF2,nRPF2,true, coil_temp, resistivity, coil_density);
Div1 = createVestPFCircuit('Div1',R_Div1,Z_Div1,width_PF,height_PF,turns(iDiv1),nZDiv1,nRDiv1,true, coil_temp, resistivity, coil_density); 
Div2 = createVestPFCircuit('Div2',R_Div2,Z_Div2,width_PF,height_PF,turns(iDiv2),nZDiv2,nRDiv2,true, coil_temp, resistivity, coil_density);


%%%%%%%%%%%%%%%%%%%%%%  INITIATE CENTRAL SOLENOID  %%%%%%%%%%%%%%%%%%%%%%

%Number of filaments of the inductor (coil = number of turns)
nfil_ind_coil = turns(iSol); 

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
coilset = fiesta_coilset('STcoilset',[Sol_circuit,PF1,PF2,Div1,Div2],false,xaccum',yaccum');
icoil=fiesta_icoil(coilset);

%Assign coil currents to icoil object [kA]
icoil.Sol=I_Sol_Equil;	        %Ensure post-ramp solenoid current is used for equilibrium
icoil.PF1=I_PF1;                %Equilibrium uses post-ramp I_PF1 current
icoil.PF2=I_PF2;                %Equilibrium uses post-ramp I_PF1 current
icoil.Div1=I_Div1;              %Equilibrium uses post-ramp I_PF1 current
icoil.Div2=I_Div2;              %Equilibrium uses post-ramp I_PF1 current

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      COMPUTE TARGET EQUILIBRIUM                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate current profiles for given simulation grid (Grid) and coilset
%Default method involves Topeol type 2 solver for Grad-Sheranov equations
config = fiesta_configuration( 'STV2C2', Grid, coilset);
control = fiesta_control( 'diagnose',true, 'quiet',false, 'convergence', 1e-5, 'boundary_method',2 );
jprofile = fiesta_jprofile_topeol2( 'Topeol2', betaP, 1, li2, Ip );

%%%%%%%%%%

%Define numerical technique applied to fit equilibrium
IEquilMethod = 'efit';         %'standard','efit','feedback'

%Standard equilibrium model (steady state coil currents)
if strcmp(IEquilMethod, 'standard');

	%Forward equilibrium which computes the jprofile for the input Irod and icoil configuration
	%icoil includes solenoid equilibrium current as default - allows for non-zero isol
	equil = fiesta_equilibrium('STV2C2', config, Irod, jprofile, control, [], icoil);   %'SST'???
	EquilParams = parameters(equil);

%%%%%%%%%%           %%%%%%%%%%           %%%%%%%%%%           %%%%%%%%%%

%Standard efit equilibrium (fit coil currents to jprofile)
elseif strcmp(IEquilMethod, 'efit');
    
    %Define efit plasma geometry (N.B. Kappa and delta are optional)
    %efit_Geometry = [RGeo, ZGeo, a, Kappa, delta];
    efit_Geometry = [0.44, 0.0, 0.44/1.85, 1.8, 0.2];

    %Efit outputs coil currents resulting from the supplied jprofile, icoil and geometry
	%Returns new currents for the requested coils: {'Coil1, {...}, 'Coiln'}
	[efit_config, signals, weights, index] = efit_shape_controller(config, {'PF1','PF2'}, efit_Geometry);

	%Inverse equilibrium, outputs coil currents resulting in the supplied jprofile and icoil config
	equil = fiesta_equilibrium('ST', config, Irod, jprofile, control, efit_config, icoil, signals, weights);
	EquilParams = parameters(equil);

	%Extract the new coil currents from the efit-equilibrium:
	icoil = get(equil,'icoil');
	efitCurrents = get(icoil,'currents');
	I_PF1 = efitCurrents(iPF1);
	I_PF2 = efitCurrents(iPF2);
	I_Div1 = efitCurrents(iDiv1);
	I_Div2 = efitCurrents(iDiv2);

%%%%%%%%%%           %%%%%%%%%%           %%%%%%%%%%           %%%%%%%%%%

%Standard efit equilibrium (fit coil currents to jprofile)
elseif strcmp(IEquilMethod, 'feedback');
    
    %Define feedback plasma geometry (N.B. Kappa and delta are optional)
    %efit_Geometry = [RGeo, ZGeo, a, Kappa, delta];
    efit_Geometry = [0.44, 0.0, 0.44/1.85, 1.8, 0.2];

    %Efit outputs coil currents resulting in the supplied jprofile, icoil and geometry
	%Returns new currents for the requested coils: {'Coil1, {...}, 'Coiln'}
	feedback = shape_controller(config, {'PF2','PF3','Div1','Div2'}, RGeo, ZGeo, a, Kappa, delta);
	[efit_config, signals, weights, index] = efit_shape_controller(config, {'PF1','PF2','Div1','Div2'}, efit_Geometry);

	%Calculate equilibrium fitting coil currents to provided jprofile
	equil = set(equil, config, 'feedback',feedback);
	EquilParams=parameters(equil);

	%Extract the new coil currents from the feedback-equilibrium:
	icoil = get(equil,'icoil');
	efitCurrents = get(icoil,'currents');
	I_PF1 = efitCurrents(iPF1);
	I_PF2 = efitCurrents(iPF2);
	I_Div1 = efitCurrents(iDiv1);
	I_Div2 = efitCurrents(iDiv2);
end

%%%%%%%%%%%%%%%%%%%%%%%% PLOT TARGET EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%%%

%Plot target equilibrium following convergence
close all
figure; hold on; axis equal;
plot(vessel);
plot(coilset);
contour( get(equil,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
contour( get(equil,'Psi'),get(equil,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
title(gca,'SMART Target-Equilibrium');
legend(gca,'hide');
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_TargetEquilibrium';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%

%Plot radial cross-section of key parameters through Z=0.0m
%IMAGE NEEDS SCALING TO FIT THE FIGURE SIZE (OR VISA-VERSA)
XSec = section(equil); hold on;
Filename = '_Equilibrium_XSection';
saveas(XSec, strcat(pwd,'/',ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SET UP VIRTUAL Bp SENSORS                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up virtual sensors to detect the null field region prior to breakdown

%Define null field region as defined by the equilibrium RGeo and ZGeo
RNull = 0.05;
BP_virt_R = linspace(EquilParams.RGeo-RNull,EquilParams.RGeo+RNull,10);
BP_virt_Z = linspace(ZGeo-RNull,ZGeo+RNull,10);

%Create null field region grid
[BP_virt_R,BP_virt_Z] = meshgrid(BP_virt_R,BP_virt_Z);
BP_virt_R = BP_virt_R(:)';
BP_virt_Z = BP_virt_Z(:)';

%Create sensors over the null field region
BP_virt_theta = zeros(1,length(BP_virt_R));
nSensors = length(BP_virt_theta);

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
sensor_btheta = fiesta_sensor_btheta( 'sensor', BP_virt_R, BP_virt_Z, BP_virt_theta, BP_virt_names );

%%%%%%%%%%%%%%  INITIATE RZIP - INCLUDING VIRTUAL SENSORS  %%%%%%%%%%%%%%

%Calculate perpendicular and parallel plasma resistivity using Spitzer model
Lambda=(12*pi*((8.854E-12*1.6E-19*Te)^3/(ne*(1.6E-19)^6))^(1/2));
PlasmaResistPerp=(0.74*1.65E-9*Z_eff*log(Lambda))/((Te*1E-3)^(3/2));
PlasmaResistPara=PlasmaResistPerp/1.96;
%PlasmaResistPerpOLD = 5.94e-6		%!!!OLD VALUE!!!

%RZIP computes coefficients [A B C D] using the null field sensors
%Output C is later used to compute the null-field pre-pulse coil currents in I_PF_NULL
rzip_config = fiesta_rzip_configuration( 'RZIP', config, vessel, {sensor_btheta} );
[A, B, C, D, curlyM, curlyR, gamma, plasma_parameters, index, label_index, state] = response(rzip_config, equil, 'rp', PlasmaResistPerp);

%%%%%%%%%%%%%%%%  PLOT VIRTUAL SENSORS ONTO EQUILIBRIUM  %%%%%%%%%%%%%%%%%

%Plot region of virtual sensors ontop of equilibrium
figure; hold on; axis equal;
plot(vessel);
contour( get(equil,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
contour( get(equil,'Psi'),get(equil,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
plot(sensor_btheta);
title(gca,'SMART VirtualSensors');
legend(gca,'hide');
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_VirtualBSensors';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               CALCULATE OPTIMISED NULL-FIELD EQUILBRIUM               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Copied from ST25D Simulation
C_temp = C(end-get(sensor_btheta,'n')+1:end,1:nPF);
C1 = C_temp(:,1);
D1 = C_temp(:,2:end);

%Produces a pre-pulse current for PF/Div coils which creates the null poloidal field
coil_currents_null = zeros(1,nPF);
I_PF_null = -pinv(D1) * (C1*I_Sol_Start);
coil_currents_null(iSol) = I_Sol_Start;
coil_currents_null(2:end) = I_PF_null';

%Generate new icoil with pre-pulse (null field) configuration
icoil_null = fiesta_icoil( coilset, coil_currents_null );
equil_optimised_null = fiesta_equilibrium( 'ST25D optimised null', config, Irod, icoil_null );
EquilNullParams = parameters(equil_optimised_null);


%%%%%%%%%%%%%  CALCULATE POLOIDAL FIELD & CONNECTION LENGTH  %%%%%%%%%%%%%%

%Extract the null poloidal field
Brdata = get(get(equil_optimised_null,'Br'),'data');         %Optimized null
Bzdata = get(get(equil_optimised_null,'Bz'),'data');         %Optimized null
Bphidata = get(get(equil_optimised_null,'Bphi_vac'),'data');     %Optimized null
rgrid = get(get(get(equil_optimised_null,'Br'),'grid'),'r');
zgrid = get(get(get(equil_optimised_null,'Br'),'grid'),'z');

%Reshape into 3D for plotting
Brdata = reshape( Brdata, length(zgrid), length(rgrid) );
Bzdata = reshape( Bzdata, length(zgrid), length(rgrid) );
Bphidata = reshape( Bphidata, length(zgrid), length(rgrid) );

%Find minimum null poloidal field and associated array indices
Bpoldata = sqrt(Bzdata.^2+Brdata.^2);       %Compute BPoloidal vector
Btordata = sqrt(Bphidata.^2);               %Compute BToroidal vector
Bpolmin = min(min(Bpoldata));				%Minimum BPoloidal value
[BpolminIndex_Row, BpolminIndex_Column] = find(Bpoldata==Bpolmin);

%Average null poloidal and toroidal fields over a small region to improve statistics
%Default to the null field sensor region (Range = 0.05[m] / 0.0055[cell/m] = 10)
BpolMinAvg = 0.0;       %Initiate accumulator value to zero
BtorMinAvg = 0.0;       %Initiate accumulator value to zero
Range = 10;             %Radius over which to average null poloidal field
for i =1:Range
    for j=1:Range
        RowMod = i-ceil(Range/2);
        ColMod = j-ceil(Range/2);
        BpolMinAvg = BpolMinAvg + Bpoldata(BpolminIndex_Row+RowMod, BpolminIndex_Column+ColMod);
        BtorMinAvg = BtorMinAvg + Btordata(BpolminIndex_Row+RowMod, BpolminIndex_Column+ColMod);
    end
end
BpolMinAvg = BpolMinAvg/(Range^2);      %Connection length BPoloidal value
BtorMinAvg = BtorMinAvg/(Range^2);      %Connection length BToroidal value

%Hacky methods of increasing the abnormally low Bpoloidal null-field
%Artificially scale the Bpol by an arbitary scale factor
Bpol_Scale_Factor = 1;						%Default 1E3
BpolMinAvg = BpolMinAvg*Bpol_Scale_Factor;
%Enforce lower limit for Bpolmin as ~Earth's B-field (5.0E-5 [T])
if BpolMinAvg < 5.0E-5;
	BpolMinAvg = 5.0E-5;					
end

%Compute the effective connection length between null-field region and wall
%Effective minor radius (a_eff) assumed equal to null-field region radius
a_eff = RNull								%Effective null-field minor radius
LCon = 0.25*a_eff*(BtorMinAvg/BpolMinAvg)	%Connection length from null field to wall


%%%%%%%%%%%%%  CALCULATE FIESTA NULL-FIELD CONNECTION LENGTH  %%%%%%%%%%%%%%

%Extract null-field RGeo and relevent psi surface for connection length
%RGeo_null_line = fiesta_line('RGeo_Line', EquilNullParams.RGeo, ZGeo)
%psi_null_surf = equil_optimised_null.psi

%Calculate null-field connection length from the seperatrix to the wall for breakdown
%Omitting line 2 defaults to a line from half the grid radius to maximimum grid radius at z=0
%[length_3d, length_2d, connection, phi, path_3d, path_2d, seg1, seg2] = connection_length3(equil, psi_null_surf, RGeo_null_line) %,line2,['REVERSE'])




%%%%%%%%%%%%%  PLOT THE OPTIMISED POLOIDAL NULL-FIELD REGION  %%%%%%%%%%%%%%

%Plot the optimised null field region - required for the breakdown diagnostic
figure;
plot(vessel);
plot(coilset);
plot(equil_optimised_null);
colorbar %colorbar
title(gca,'SMART Optimised Null Field');
legend(gca,'hide');
set(gca,'XLim',[0 1]);
set(gca,'YLim',[-1.5 1.5]);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_OptimisedNull';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     DEFINE LOOP VOLTAGE FOR STARTUP                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Copied from ST25D Simulation
C_temp = C(end-get(sensor_btheta,'n')+1:end,1:nPF);
C1 = C_temp(:,1);
D1 = C_temp(:,2:end);

%Produces a pre-pulse current for PF/Div coils which creates the null poloidal field
I_PF_null = -pinv(D1) * (C1*I_Sol_Start);    %Copied From ST25D Simulation

%Definition of time intervals:
%1--> All coils and Sol start at zero current
%2--> Slow initial ramp up to pre-breakdown current
%3--> Flat pre-pulse plateau prior to plasma breakdown
%4--> Begin pulse, sol ramps down to initiate plasma breakdown
%5--> Sustain plasma, flat PF and Div, Sol ramps down to drive Ip
%6--> All coils and Sol end at zero current

%Initiate PF_input arrays to zero for all times
V_PF_input = NaN(nTime,nPF);
I_PF_input = zeros(nTime,nPF);

%All Sol/PF/Div coil currents default to zero unless subsiquently set
I_PF_input(2,:) = 0;
I_PF_input(3,:) = 0;

%Define Solenoid current waveform vertices :: Startup --> Pre-pulse
I_PF_input(2,iSol) = I_Sol_Start;            %Sets Sol to pre-pulse current
I_PF_input(2,2:end) = I_PF_null;             %Sets all PF/Div coils to I_PF_null

%Define Solenoid current waveform vertices :: Pre-pulse --> Pulse
I_PF_input(3,iSol) = I_Sol_Start;            %Sets Sol to pre-pulse current
I_PF_input(3,2:end) = I_PF_null;             %Sets all PF/Div coils to I_PF_null

%Define coilset current waveforms vertices :: Pulse --> Equilibrium
I_PF_input(4,iSol) = I_Sol_Equil;     %Solenoid Equilibrium current
I_PF_input(4,iPF1) = I_PF1;           %
I_PF_input(4,iPF2) = I_PF2;           %
I_PF_input(4,iDiv1) = I_Div1;         %
I_PF_input(4,iDiv2) = I_Div2;         %

%Define coilset current waveforms vertices :: Equilibrium --> Finish
I_PF_input(5,iSol) = -I_Sol_Start;   %I_Sol_End
I_PF_input(5,iPF1) = I_PF1;          %
I_PF_input(5,iPF2) = I_PF2;          %
I_PF_input(5,iDiv1) = I_Div1;        %
I_PF_input(5,iDiv2) = I_Div2;        %

%All coil currents end at zero
I_PF_input(6,iSol)=0;
I_PF_input(6,iPF1)=0;
I_PF_input(6,iPF2)=0;
I_PF_input(6,iDiv1)=0;
I_PF_input(6,iDiv2)=0;

%%%%%%%%%%

%Convert from 'time' to 'long-time' for increased temporal resolution
nTime_long = 1000;
time_long = linspace(min(time),max(time),nTime_long);
I_PF_input_long = NaN(nTime_long,nPF);
for iPF=1:nPF
    I_PF_input_long(:,iPF) = interp1(time,I_PF_input(:,iPF),time_long);
end
V_PF_input_long = NaN*I_PF_input_long;

%Initiate Plasma Currrent and Plasma Potential arrays in 'long-time'
Ip_long = zeros(size(time_long));
Vp_long = NaN(size(time_long));

%Name and colour coils for plotting
coil_names{iSol} = 'Sol';
coil_names{iPF1} = 'PF2';
coil_names{iPF2} = 'PF3';
coil_names{iDiv1} = 'Div1';
coil_names{iDiv2} = 'Div2';
PF_colors{iSol} = 'Red';
PF_colors{iPF1} = 'Magenta';
PF_colors{iPF2} = 'Black';
PF_colors{iDiv1} = 'Cyan';
PF_colors{iDiv2} = 'Green';

%Plot figure showing dynamic coil currents - With Eddy Currents?
figure;
plot(time_long*1000, I_PF_input_long/1000);
title(gca,'SMART Initial Coil Current Waveforms');
LegendString = {'Sol','PF1','PF2','Div1','Div2'};
legend(gca,LegendString);
xlabel(gca,'Time (ms)');
ylabel(gca,'Coil Current (kA)');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
Filename = '_CurrentWaveforms';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%

%Compute dynamic coil currents employing current driven Ip
[ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, time_adaptive ] = ...
    state_space_including_passive_elements_v4( curlyM, curlyR, time_long, I_PF_input_long, V_PF_input_long, Ip_long, Vp_long, 'adaptive_timesteping',true );

%Set plasma voltage to zero from time zero (assume breakdown) and default plasma current to 'NaN'
iTime_plasma = time_adaptive > 0;						%iTime_plasma: 1 for true, 0 for false
Vp_output(iTime_plasma) = 0;							%Sets Vp = 0 when t > 0.
Vp_long = interp1(time_adaptive, Vp_output, time_long); %Sets Vp_long = 0 when t_long > 0.
Ip_long = NaN*Vp_long;									%Sets Ip_long to 'NaN' array

%Compute dynamic coil currents employing voltage driven Ip
[ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, time_adaptive ] = ...
    state_space_including_passive_elements_v4( curlyM, curlyR, time_long, I_PF_input_long, V_PF_input_long, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names', coil_names, 'show_plot',true, 'turns',turns, 'currentScale',1e3, 'PF_colors',PF_colors );


%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT PLASMA CURRENT %%%%%%%%%%%%%%%%%%%%%%%%%%              

%Plot plasma current over full timescale
close all
plot(time_adaptive*1000, Ip_output/1000)
title(gca,'SMART Plasma Current');
legend(gca,'Plasma Current');
xlabel(gca,'Time (ms)');
ylabel(gca,'Plasma Current (kA)');
set(gca,'XLim',[min(time*1e3) max(time*1e3)]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
Filename = '_PlasmaCurrent';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              DETERMINE EDDY CURRENTS WITHIN THE VESSEL                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%I_Passive contains eddy current of the nf filaments at each instant of time.
%len(I_Passive) = 3811*nf == len(time_adaptive) = 3811*1

%Obtain the filament variables r and z
ptmp = get(vessel,'passives');
ftmp = get(ptmp,'filaments');
RR = get(ftmp(:),'r'); %dim 1*number of filaments; number of filaments=nf
ZZ = get(ftmp(:),'z'); %dim 1*number of filaments

%Okay. We do have the eddy current on each filament an on every time on
%I_Passive. I_passive on the figures is sum over each filament of the eddy
%current, to get the total eddy current induced upon each time. We can not
%sum for every time the eddy current, since it varys its sign, it also
%ceases during certain amounts of tim, so it can not be done.  But there is
%no necesssity, since I_passive contains the eddy current at any time of
%the filament, and this will also provide the force upon each instant; i
%only need to pick up the greatest

%Time intervals intersected with number of filaments (time intervals*number of filaments)
sizeIpas=size(I_Passive);
%For each vessel filament extract the greatest absolute current
for i=1:sizeIpas(2)
	%Obtain largest positive and negative for each filament
    positive=max(I_Passive(:,i));
    negative=min(I_Passive(:,i));
	%Keep the largest absolute value
    if abs(positive)> abs(negative)
        I_Passive_fil(i)=positive;
    else
        I_Passive_fil(i)=negative;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% PLOT EDDY CURRENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%   

%Plot eddy currents within a cross-section of the vessel
close all
figure; hold on; axis equal;
plot(coilset);
%plot(vessel);
scatter3(RR,ZZ,I_Passive_fil/1000,100,I_Passive_fil/1000,'filled');
title('SMART Vessel Eddy-Currents');
view(2) %2D view
colorbar %colorbar
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
%zlabel(gca, 'I (A)')
Filename = '_EddyCurrent';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%

%Plot eddy currents at particular location with respect to time
% for i=1:5%length(time_adaptive)
% figure;
% acumm=I_Passive(i,:)+acumm
%  scatter3(RR,ZZ,acumm,100,acumm,'filled')
%  view(2) %para verlo en 2D
% xlabel(' R (m)')
% ylabel('Z (m)')
% zlabel('I (A)')
% axis([0,1,-1.1,1.1]) %for the tfg
% title('sprintf(iter %d,i)')
% set(gca, 'FontSize', 13, 'LineWidth', 0.75); %<- Set properties TFG
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  DETERMINE FORCES UPON THE VESSEL                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Obtain the equilibrium B-field in R,Z and Phi
Br=get(equil,'Br'); 		%fiesta_field structure
Bz=get(equil,'Bz'); 		%fiesta_field structure
Bphi=get(equil,'Bphi'); 	%fiesta_field structure

%Convert from a fiesta_field structure to a set of 2D grids (200*251 GridR*GridZ)
Brdata = get(get(equil,'Br'),'data');
Bzdata = get(get(equil,'Bz'),'data');
Bphidata = get(get(equil,'Bphi'),'data');

%Reshape data such that it can be interpolated:
Brdata = reshape( Brdata, length(zGrid), length(rGrid));
Bzdata = reshape( Bzdata, length(zGrid), length(rGrid));
Bphidata = reshape( Bphidata, length(zGrid), length(rGrid));

%Interpolate the B-fields onto a grid that aligns with the vessel grid 
%These are the values of the B-field at the vessel grid points
%Br_interp aligns with the vessel filament cells (RR, ZZ) from before
Br_interp = @(r,z) interpn(zGrid,rGrid,Brdata,z,r);
Bz_interp = @(r,z) interpn(zGrid,rGrid,Bzdata,z,r);
Bphi_interp = @(r,z) interpn(zGrid,rGrid,Bphidata,z,r);

%Extract B-field at vessel walls - meshes are aligned so indexes are the same
Br_vessel=Br_interp(RR,ZZ);
Bz_vessel=Bz_interp(RR,ZZ);
Bphi_vessel=Bphi_interp(RR,ZZ);

%Combine Bfield values in [R,Phi,Z] into a single array for the vessel
%size(number of filaments*3), each row is the vector field on one filament
B_vessel=[Br_vessel' Bphi_vessel' Bz_vessel']; 
%The maximum current on each vessel filament is I_Passive_fil (size 1*number of filaments)
%Current vector is in the phi direction so only take magnitude [0, 1*I_Passive(phi), 0]
I_passive_fil_vec=I_Passive_fil'*[0 1 0]; 		%size [number of filaments*3]

%The force upon all the filament would be 2piR*Force_fil_cross. R is stores
%in RR, which contains all the R values in a vector form with number of fil components. 
%Force_fil_cross is a vector of 3 components. It would be difficult to
%multiply them, but we do not need to, right now, because to obtain the
%pressure R cancles out, since the areas are 2piR*anchura (or altura). We
%assimilate the 3D filament as a 2D filament, so that it has no width in
%the R axis, s that its surface is 2piR*altura

%Compute J X B force acting on each filament (J X B computed for all directions) 
Force_fil=cross(I_passive_fil_vec,B_vessel);	%[N] %size [number of filaments*3]
%Take magntude of all forces as some will be negative (only care about the maximum force)
[Force_max, index]=max(abs(Force_fil));			%[N] %Also obtain index of each force

%Pressure acting on vessel wall is force over unit filiment area
%These are absolute numbers - don't include any directionality
Pressure_R=abs(Force_max(1))/(height_PF);		%[Pa]
Pressure_Z=abs(Force_max(3))/(height_PF);		%[Pa]

%Stress acting on vessel wall is the combined force divided by the unit filiment area
%These are directional, some are negative and some are positive
stress_R=(Force_fil(:,1))/(height_PF);			%[Pa]
stress_Z=(Force_fil(:,3))/(height_PF);			%[Pa]
%Obtain maximum radial and axial stresses - either positive or negative
stress_R_max=max(abs(stress_R));				%[Pa]
stress_Z_max=max(abs(stress_Z));				%[Pa]


%%%%%%%%%%%%%%%%%%%%%% PLOT VESSEL EDDY STRESSES %%%%%%%%%%%%%%%%%%%%%%%%

%Scale stresses from [Pa] to [Atm]
ScaleFactor=2e5; 

%Plot figure showing vessel eddy stresses
close all
figure; hold on; axis equal;
plot(coilset);
%plot(vessel);
quiver(xaccum,yaccum,stress_R/ScaleFactor,stress_Z/ScaleFactor,'color',[1 0 0],'AutoScale','off');
title('SMART Vessel Eddy-Stresses 3V-Phase1');
view(2) %2D view
colorbar %colorbar
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
%zlabel(gca, 'I (A)')
Filename = '_EddyStresses';
saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%

%Plot figure showing eddy forces - Note: material tolerences are stated as max stress
%close all
%figure; hold on; axis equal;
%plot(coilset);
%plot(vessel);
%quiver(xaccum,yaccum,Force_fil(:,1),Force_fil(:,3),'color',[1 0 0])
%title('SMART Vessel Eddy-Forces 3V-Phase1');
%view(2) %2D view
%colorbar %colorbar
%set(gca,'XLim',[0 1.1]);
%set(gca,'YLim',[-1.1 1.1]);
%set(gca, 'FontSize', 13, 'LineWidth', 0.75);
%xlabel(gca,'R (m)');
%ylabel(gca,'Z (m)');
%zlabel(gca, 'I (A)')
%Filename = '_EddyForces';
%saveas(gcf, strcat(ProjectName,Filename,FigExt));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      DATA OUTPUT TO TEXT FILES                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Write qeqdsk equilibrium file
filename = strcat(ASCIIDir,'Equil.txt');
geqdsk_write_BUXTON(config, equil, filename)

%Write equilibrium parameters file
EquilParam=parameters(equil);
%struct2csv(EquilParam,Filename);           %csv format if required
ParamVariables = fieldnames(EquilParam);
ParamValues = struct2cell(EquilParam(:,1));
Filename = strcat(ASCIIDir,'EquilParam.txt');
fileID=fopen(Filename,'w');
for i = 1:length(ParamValues)
    try fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); ParamValues(i)]);
    catch fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); 'NaN']);
    end
end

%%%%%%%%%%

Filename = strcat(ASCIIDir,'t.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%1.12f\r\n',time_adaptive);

Filename = strcat(ASCIIDir,'Ip.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%1.12f %1.12f\r\n',[time_adaptive'; Ip_output']);

Filename = strcat(ASCIIDir,'IPass.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%1.12f %1.12f\r\n',[time_adaptive'; I_Passive']);

Filename = strcat(ASCIIDir,'LCon.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%1.12f\r\n',LCon);

%Extract coil current time-traces (Multiplied by number of windings)
a=I_PF_output(:,1).*turns(iSol);   %Sol
b=I_PF_output(:,2).*turns(iPF1);   %PF1
c=I_PF_output(:,3).*turns(iPF2);   %PF2
d=I_PF_output(:,4).*turns(iDiv1);  %Div1
e=I_PF_output(:,5).*turns(iDiv2);  %Div2
Filename = strcat(ASCIIDir,'CoilCurrents.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[a'; b'; c'; d'; e']);

%Extract coil voltage time-traces
a=V_PF_output(:,1);     %Sol
b=V_PF_output(:,2);     %PF1
c=V_PF_output(:,3);     %PF2
d=V_PF_output(:,4);     %Div1
e=V_PF_output(:,5);     %Div2
Filename = strcat(ASCIIDir,'CoilVoltages.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[a'; b'; c'; d'; e']);

%%%%%%%%%%

Filename = strcat(ASCIIDir,'MaxStress.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%1.12f %1.12f\r\n',[stress_R_max'; stress_Z_max']);

%%%%%%%%%%

%Print final equilibrium parameters to terminal
disp([ ' ' ]);
disp([ 'Equilibrium Parameters:' ]);
disp([ EquilParam ]);
disp([ ' ' ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






