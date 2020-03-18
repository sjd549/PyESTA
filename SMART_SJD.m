%% SMall Aspect Ratio Tokamak (SMART), V3p1 init.nam

clear 
clc
close all

%%%%%%%%%%

%Add FIESTA trunk path and any other required paths
FIESTATrunk = "~/Postdoc Seville/FIESTA/Source Code/FIESTA_V8.8"
addpath(genpath(FIESTATrunk),"Functions")

%Define any required global variables
RunFolder = 'CurrentRun_SMART-V3p1'
ProjectName = 'SMART-V3p1';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       INITIATE REACTOR GEOMETRY                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define Vessel Outer Geometry
VesselRInnerPoint=0.15; % R min position [m]
VesselROuterPoint=0.8;  % R max position [m]
VesselZMinPoint=-0.8;   % Z min position [m]
VesselZMaxPoint=0.8;    % Z max position [m]

%Define Solenoid Geometry and Parameters
nSol=800;     % number of turns of the solenoid
RSol=0.09;    % R position of the solenoid [m] (Inner Solenoid)
ZMinSol=-0.8; % Min Z position
ZMaxSol=0.8;  % Max Z position

%Define Div and PF coils
nZDiv1=6;
nRDiv1=4;
nZDiv2=6;
nRDiv2=4;
%nZPF1=6;
%nRPF1=4;
nZPF2=6;
nRPF2=4;
nZPF3=6;
nRPF3=4;

%Calculate total number of turns in each coil
nDiv1=nZDiv1*nRDiv1;
nDiv2=nZDiv2*nRDiv2;
%nPF1=nZPF1*nRPF1;
nPF2=nZPF2*nRPF2;
nPF3=nZPF3*nRPF3; 

%Define coil size to enable cross-section calculation
width_PF=0.042;  % Width of a turn (m)
height_PF=0.035; % Height of a turn (m)

%Define central location of coil sets
RScale=1.00;        %Scale all radial coil positions (Relative to 2.0m)
ZScale=1.0;        %Scale all axial coil positions (Relative to 2.0m)
%%%%%
%R_PF1=0.9*RScale;  %R position of PF1 (m)
%Z_PF1=0.3*ZScale;  %Z position of PF1 (m)
R_PF2=0.9*RScale;   %R position of PF2 (m)
Z_PF2=0.5*ZScale;   %Z Position of PF2 (m)
R_PF3=0.9*RScale;   %R Position of PF3 (m)
Z_PF3=0.8*ZScale;   %Z Position of PF3 (m)
R_Div1=0.25*RScale; %R Position of Div1 (m)
Z_Div1=1.05*ZScale; %Z Position of Div1 (m)
R_Div2=0.55*RScale; %R Position of Div2 (m)
Z_Div2=1.05*ZScale; %Z Position of Div2 (m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     DEFINE OPERATIONAL PARAMETERS                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Explicitly defined parameters - Taken from 2.0m SMART Phase 1
% Added here for reference, not for explicit use in code.
% betaT    = 0.1000   [%]
% betaP    = 0.4578   [%]
% li2      = 1        []
% Ip       = 30e3     [kA]
% Irod     = 4.7750e5 [kA]
% B0       = 0.10     [T]  (0.09 [T] in vacuo)
% TauPulse = 0.020    [s]  (Timestep for pulse)

%%%%%%%%%%%%%%%%%%%%  DEFINE INITIAL PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%
disp([ ' ' ]);
disp([ '%===== Initial Operating Parameters =====%' ]);

%Define any required constants
mu0 = 1.2566e-06;

%Define plasma geometry
RGeo = 0.4763;    % Geometrical Radius [m]
epsilon = 1.985;  % Aspect ratio   []
a = RGeo/epsilon; % Minor radius   [m]
kappa = 1.7;      % Elongation     []
Ip = 30e3;        % Plasma current [A]
li2 = 1;          %                []
%q_cyl = 2.821;   % Safety Factor  []
betaN =3.529;     %                [] (Obtained via VEST Excel)

%Define pulse length
TauPulse = 0.020;    %[s]  (Also Timestep for plotting Ip)
disp([ 'TauPulse = ' num2str(TauPulse*1000) ' [ms]' ]);

%Define/Calculate Toroidal B-field at plasma geometric centre;
% BT = Irod*mu0 / (2*pi*RGeo);  % [T]
% disp([ 'BT = ' num2str(BT) ' [T]' ]);
BT=0.1; %[T]
disp([ 'BT = ' num2str(BT) ' [T]' ]);

%Define/Calculate Central Rod Current
% Irod = q_cyl*2*pi / (5e6*mu0) * Ip * A^2 / kappa;
% disp([ 'Irod = ' num2str(Irod/1e3) ' [kA*Turn]' ]);
Irod=BT*2*pi*RGeo/mu0;
disp([ 'IRod = ' num2str(Irod) ' [kA]' ]);

%Define/Calculate Shaping factor
S = sqrt( (1.0+kappa^2)/2.0 );
disp([ 'Shaping Factor = ' num2str(S) ' [-]' ]);

%Define/Calculate Beta toroidal [%]
betaT = betaN/a * (Ip/1e6)/BT;
disp([ 'beta = ' num2str(betaT) ' [%]' ]);

%Define/Calculate Beta Poloidal [%]
betaP = 25 ./(betaT / 100) .* S^2 .* ((betaN / 100).^2);
disp([ 'betaP = ' num2str(betaP) ' [-]' ]);

%Define/Calculate Combined density and temperature
nT = 2.66*1e20*1e3 * betaT * BT^2;
disp([ 'nT = ' num2str(nT) ' [m^-3 eV]' ]);

%Define/Calculate Central plasma density
n = 3e19;  % m^-3
T = nT/n;  % eV
disp([ 'T = ' num2str(T) ' [eV]' ]);

%Define/Calculate Required vertical field
BZ = - mu0 * Ip / (4*pi*RGeo) * ( log(8*epsilon) + betaP + 0.5*li2 - 3/2 );
disp([ 'BZ = ' num2str(BZ) ' [T]' ]);

%Define/Calculate ?????????????
nGW = Ip * 1e14 / (pi*a^2);

%Define/calculate resistivity and density of the coils
coil_density = 1;   %[Arb]
coil_temp = 293.0;  %[K]
resistivity = copper_resistivity_at_temperature(coil_temp);



%%%%%%%%%%%%%%%%%%  DEFINE SOL RAMP & COIL CURRENTS  %%%%%%%%%%%%%%%%%%%%
disp([ ' ' ]);
disp([ '%===== Initial Coil Currents =====%' ]);

%Phase1 coil currents [kA]                %SJDoyle
I_Sol_start=1500;    %+1500 -> +3000;     %+1500;
I_Sol_ramp=-1100;    %-1100 -> -1150;     %-1100;
I_Sol_equil=-1500;   %-1500 -> -1500;     %-1500;
%
I_PF1=0;             %-0    -> -0         %N/A
I_PF2=-1175;%-825;          %-0800 -> -1100;     %-0825;
I_PF3=-900;%-700          %-0700 -> -900;      %-0700;
I_Div1=-000;         %-0000 -> -0000;     %+0000;
I_Div2=+4440;%2400        %+2400 -> +4400;     %+2400;

%Phase2 coil currents [kA]
%I_Sol_start=+6500;    %+4700 -> +5500;     %+4700;     %+6500;
%I_Sol_ramp=+500;      %+500  -> +500       %+500;      %+500;
%I_Sol_equil=-6500;    %-4700 -> -5500;     %-4700;     %-6500;
%
%I_PF1=0;              %-0    -> -0         %N/A
%I_PF2=-3400;          %-3250 -> -3600      %-3000;     %-3400;
%I_PF3=-940;           %-940  -> -940       %-940;      %-940;
%I_Div1=-000;          %-000  -> -000       %-000;      %-000;
%I_Div2=+9090;         %+9090 -> +9090      %+9090;     %+9090;

%Phase3 coil currents [kA]
%I_Sol_start=4700;   %4200;
%I_Sol_ramp=500;     %
%I_Sol_equil=-4700;  %-5200;
%
%I_PF1=0;            %-0    -> -0         %N/A
%I_PF2=-6000;        %
%I_PF3=-3100;        %
%I_Div1=-000;        %
%I_Div2=15880;       %


%%%%%%%%%%%%%%%%%%%%%  Construct Run Folder Name  %%%%%%%%%%%%%%%%%%%%%

%String relevant inputs for data management purposes
Str_RScale = strcat('RScale#',string(RScale),{', '});
Str_ZScale = strcat('ZScale#',string(ZScale),{', '});
Str_R_PF2 = strcat('RPF2#',string(R_PF2),{', '});
Str_Z_PF2 = strcat('ZPF2#',string(Z_PF2),{', '});
Str_R_PF3 = strcat('RPF3#',string(R_PF3),{', '});
Str_Z_PF3 = strcat('ZPF3#',string(Z_PF3),{', '});
Str_R_Div1 = strcat('RDiv1#',string(R_Div1),{', '});
Str_Z_Div1 = strcat('ZDiv1#',string(Z_Div1),{', '});
Str_R_Div2 = strcat('RDiv2#',string(R_Div2),{', '});
Str_Z_Div2 = strcat('ZDiv2#',string(Z_Div2),{', '});
%
Str_BT = strcat('BT#',string(BT),{', '});
Str_Ip =  strcat('Ip#',string(Ip),{', '});
Str_TauPulse = strcat('Tau#',string(TauPulse),{', '});
%
Str_Sol = strcat('Sol#',string(I_Sol_start),{', '});
Str_PF1 = strcat('PF1#',string(I_PF1),{', '});
Str_PF2 = strcat('PF2#',string(I_PF2),{', '});
Str_PF3 = strcat('PF3#',string(I_PF3),{', '});
Str_Div1 = strcat('Div1#',string(I_Div1),{', '});
Str_Div2 = strcat('Div2#',string(I_Div2),{', '});
%

%Create output folder for saved data and figures
RunName = strcat(Str_ZScale,Str_BT,Str_TauPulse,Str_Sol,Str_PF2,Str_Div2);
%RunName = strcat(Str_BT,Str_TauPulse,Str_Sol,Str_PF2,Str_Div2);
%RunName = strcat(Str_Z_PF2,Str_Z_PF3,Str_Z_Div1,Str_Z_Div2);
%
FolderName = convertStringsToChars( strcat(RunFolder,'/',RunName,'/') );
mkdir(FolderName);

%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      INITIATE COIL CONFIGURATION                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pnts = [ VesselRInnerPoint VesselZMaxPoint                   
              VesselROuterPoint VesselZMaxPoint];
pnts2 = [pnts; [pnts(end:-1:1,1),-pnts(end:-1:1,2)]; pnts(1,:)];    

%Approximate vessel geometry in 'pixelated' form through Bresenham algorithm
ww =1.5e-2;         %Wall Thickness [m]
xaccum = [];
yaccum = [];
for i=1:length(pnts2)-1
    [xx,yy]=bresenham( pnts2(i,1)/ww,pnts2(i,2)/ww,pnts2(i+1,1)/ww,pnts2(i+1,2)/ww);
    xaccum = [xaccum;xx*ww];
    yaccum = [yaccum;yy*ww];
end
%Remove duplicate cells (if any)
dup = (abs(diff(xaccum))+abs(diff(yaccum))) > 0;
xaccum = xaccum(dup);
yaccum = yaccum(dup);


%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIATE COILS  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define and initiate PF coils - Arbitrary numbering of coils
iSol = 1;       %Central Inducting Solenoid
iDiv1 = 4;      %Inboard Divertor Coil
iDiv2 = 5;      %Outboard Divertor Coil
%iPF1 = 6;      %Trangularity coil at Z=0.3 [m]
iPF2 = 2;       %Upper Plasma Forming Coil
iPF3 = 3;       %Lower Plasma Forming Coil

%Create array containing coil turns
turns=[];
turns(iSol) = nSol; 
turns(iDiv1) = nDiv1;
turns(iDiv2) = nDiv2;
%turns(iPF1) = nPF1;
turns(iPF2) = nPF2;
turns(iPF3) = nPF3;

nPF = 5; 
%Create coil set from parameters defined above. (Function made by Carlos Soria)
%Function createVESTPFCircuit creates two PF coils. One in (R, Z) and another in (R, -Z)
%PF1  = createVestPFCircuit( 'PF1',R_PF1,Z_PF1,width_PF,height_PF,turns(iPF1),nZPF1,nRPF1,true, coil_temperature, resistivity, density);
PF2  = createVestPFCircuit('PF2',R_PF2,Z_PF2,width_PF,height_PF,turns(iPF2),nZPF2,nRPF2,true, coil_temp, resistivity, coil_density);
PF3  = createVestPFCircuit('PF3',R_PF3,Z_PF3,width_PF,height_PF,turns(iPF3),nZPF3,nRPF3,true, coil_temp, resistivity, coil_density);
Div1 = createVestPFCircuit('Div1', R_Div1, Z_Div1, width_PF,height_PF, turns(iDiv1), nZDiv1,  nRDiv1, true, coil_temp, resistivity, coil_density); 
Div2 = createVestPFCircuit('Div2', R_Div2, Z_Div2, width_PF,height_PF, turns(iDiv2), nZDiv2,  nRDiv2, true, coil_temp, resistivity, coil_density);


%%%%%%%%%%%%%%%%%%%%%%  INITIATE CENTRAL SOLENOID  %%%%%%%%%%%%%%%%%%%%%%

%Number of filaments of the inductor (coil = number of turns)
nfil_ind_coil = turns(iSol); 
clear('coil_filaments');
Z_filament = linspace(ZMinSol,ZMaxSol,nfil_ind_coil);
%Construct central solenoid filaments
for iFilament=1:nfil_ind_coil
    coil_filaments(iFilament) = fiesta_filament( RSol,Z_filament(iFilament), sqrt(70e-6),sqrt(70e-6) ); 
end
coil_1  = fiesta_coil( 'psh_coil', coil_filaments, 'Blue', resistivity, coil_density );
Sol_circuit = fiesta_circuit( 'Sol', [1], [coil_1] );

%Collate completed coilset
coilset = fiesta_coilset('STcoilset',[Sol_circuit,PF2,PF3,Div1,Div2],false,xaccum',yaccum');
% save('Configuration.mat','vessel','coilset')


%%%%%%%%%%%%%%%%%%  INITIATE VACUUM VESSEL FILAMENTS  %%%%%%%%%%%%%%%%%%%

%Construct vessel wall filaments for eddy current calculation
for i=length(xaccum):-1:1
    vessel_filament(i) = fiesta_filament(xaccum(i),yaccum(i),ww,ww/3,1,0,0);
end
%Enable induced currents in vessel wall filaments
passive = fiesta_passive('STVesselPas',vessel_filament,'g');
vessel = fiesta_vessel( 'STVessel',passive);

%Define FIESTA grid limits and resolution
R_simulation_limits = [0.03 1];
Z_simulation_limits = [-1.3 1.3]; 
grid_size_R =251;
grid_size_Z =200; 
Grif = fiesta_grid( R_simulation_limits(1), R_simulation_limits(2), grid_size_R, Z_simulation_limits(1), Z_simulation_limits(2), grid_size_Z );
rgrid=get(Grif,'r'); 
zgrid=get(Grif,'z'); 

%Calculate current profiles for given simulation grid (Grif) and coilset
config = fiesta_configuration( 'STV2C2', Grif, coilset);
control = fiesta_control( 'diagnose',true, 'quiet',false, 'convergence', 1e-5, 'boundary_method',2 );
jprofile = fiesta_jprofile_topeol2( 'Topeol2', betaP, 1, li2, Ip );

%%%%%%%%%%

%Plot coilset for sanity checking purposes
%figure;
%plot(coilset)
%plot(vessel)
%hold on;
%plot(Sol_circuit);
%plot(PF1)
%plot(PF2)
%plot(PF3);
%plot(Div1);
%plot(Div2);
%axis equal;
%title('SMART-V3p1_TargetEquilibrium Cross-section')
%legend(gca,'hide');
%set(gca,'XLim',[0 1]);
%set(gca,'YLim',[-1.5 1.5]);
%xlabel(gca,'R (m)');
%ylabel(gca,'Z (m)');
%Foldername = 'Data/';
%Filename = 'SMART-V3p1_TargetEquilibrium';
%saveas(gcf, strcat(Foldername,Filename,'.png'));
%saveas(gcf, strcat(Foldername,Filename,'.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   DEFINE COIL SET & RUN EQUILIBRIUM                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create icoil object - Required as FIESTA functions use OO-programming
icoil=fiesta_icoil(coilset);

%Assign coil currents to icoil object [kA]
%icoil.PF1=I_PF1;     %
icoil.PF2=I_PF2;      %
icoil.PF3=I_PF3;      %
icoil.Div1=I_Div1;    %
icoil.Div2=I_Div2;    %

%%%%%%%%%%%%%%%%%%%%%%%%  COMPUTE EQUILIBRIUM  %%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate equilibrium for given coil geometry and currents.
equil=fiesta_equilibrium('SST', config, Irod, jprofile, control, [], icoil);

%%%%%%%%%%

%HACKY NOTE!!! This seems to be a slightly more involved equilibrium calc
%HACKY NOTE!!! efit_shape_controller doesn't work - r variable not defined
%[efit_config, signals, weights, index]=efit_shape_controller(config, {'Div2'}, [0.4, 0, 0.25, 1.8, 0.2]) %[0.47, 0, 0.24, 1.7, 0.3]
%equil=fiesta_equilibrium('ST', config, Irod, jprofile, control, efit_config, icoil, signals, weights)

%HACKY NOTE!!! WHAT ARE THESE FOR???
%config = fiesta_configuration( 'STV2C2', Grif, coilset);
%feedback=shape_controller(config, {'PF2','PF3','Div1','Div2'}, 0.43, 0, 0.22, 1.82, 0.1)
%[efit_config, signals, weights, index]=efit_shape_controller(config, {'PF2','PF3','Div1','Div2'}, [0.43, 0, 0.24])

%feedback=shape_controller(config, {'Div2'}, 0.43, 0, 0.22, 1.8);
%equil=set(equil, config, 'feedback',feedback);
%save('equil.mat','equil')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  PLOT TARGET EQUILIBRIUM CONVERGENCE                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure; hold on; axis equal;
plot(vessel);
plot(coilset);
contour( get(equil,'Psi'),60,'Color','Black', 'LineWidth',0.5 );
contour( get(equil,'Psi'),get(equil,'Psi_boundary')*[1 1],'Color','Black', 'LineWidth',1.5 );
title(gca,'SMART Target-Equilibrium, V3-phase1');
legend(gca,'hide');
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[-1.1 1.1]);
set(gca, 'FontSize', 13, 'LineWidth', 0.75);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_TargetEquilibrium';
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.png'));
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.pdf'));

%Write qeqdsk equilibrium file
filename = strcat(FolderName,'ST_Phase1');
geqdsk_write_BUXTON(config, equil, filename)

%%%%%%%%%%

%Plot radial cross-section of key parameters through Z=0.0m
%IMAGE NEEDS SCALING TO FIT THE FIGURE SIZE (OR VISA-VERSA)
XSec = section(equil); hold on;
Filename = '_EquilibriumXSection';
saveas(XSec, strcat(FolderName,ProjectName,Filename,'.png'));
saveas(XSec, strcat(FolderName,ProjectName,Filename,'.pdf'));

%Write equilibrium parameters file
EquilParam=parameters(equil);
%struct2csv(EquilParam,Filename);           %csv format if required
ParamVariables = fieldnames(EquilParam);
ParamValues = struct2cell(EquilParam(:,1));
Filename = strcat(FolderName,'EquilParam_Phase_1.txt');
fileID=fopen(Filename,'w');
for i = 1:length(ParamValues)
    try fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); ParamValues(i)]);
    catch fprintf(fileID,'%s, %0.5f\r\n',[string(ParamVariables(i)); 'NaN']);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Magnetic sensor to check Bp
% % 
% % % Create Virtual Sensors
% % BP_virt_R = linspace(0.4,0.5,10);
% % BP_virt_Z = linspace(-0.05,0.05,10);
% % 
% % [BP_virt_R,BP_virt_Z] = meshgrid(BP_virt_R,BP_virt_Z);
% % BP_virt_R = BP_virt_R(:)';
% % BP_virt_Z = BP_virt_Z(:)';
% % 
% % BP_virt_theta = zeros(1,length(BP_virt_R));
% % nSensors = length(BP_virt_theta);
% % 
% % BP_virt_names = {};
% % for iSensor=1:nSensors
% %     BP_virt_names{iSensor} = ['Radial Bp Virtual Sensor #' num2str(iSensor) ];
% % end
% % 
% % BP_virt_R = [BP_virt_R  BP_virt_R];
% % BP_virt_Z = [BP_virt_Z  BP_virt_Z];
% % 
% % BP_virt_theta = [BP_virt_theta  BP_virt_theta+pi/2];
% % 
% % for iSensor=nSensors+1:2*nSensors
% %     BP_virt_names{iSensor} = ['Vertical Bp Virtual Sensor #' num2str(iSensor) ];
% % end
% % 
% % bth= fiesta_sensor_btheta( 'sensor', BP_virt_R, BP_virt_Z,BP_virt_theta, BP_virt_names );
% % br=fiesta_sensor_br('xpoint', BP_virt_R,BP_virt_Z);
% % bz=fiesta_sensor_bz('xpoint', BP_virt_R,BP_virt_Z);
% 
% %%
% % save('equil_Div1_0.4m_30kA_Div2_0.7m_30kA_PF2_-4.75kA_PF3_-4.75kA.mat','equil')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             INITIATE RZIP                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rzip_config = fiesta_rzip_configuration( 'ST25_RZIP', config, vessel );
plasma_resistance = 5.94e-6;
[epsilon, B, C, D, curlyM, curlyR, gamma, plasma_parameters, index, label_index, state] = response(rzip_config, equil, 'rp',plasma_resistance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%% DEFINE VESSEL TIME CONSTANT %%%%%%%%%%%%%%%%%%%%%%

[~,tau_vessel] = eig(curlyR(1:end-3,1:end-3)\curlyM(1:end-3,1:end-3));
tau_vessel = max(diag(tau_vessel));
disp([ 'tau_vessel=' num2str(tau_vessel*1e3) 'ms' ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 CALCULATE LOOP VOLTAGE FOR STARTUP                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nTime = 6;
I_PF_input = zeros(nTime,nPF);

%All PF coil currents initiate at zero
I_PF_input(2,:) = 0;
I_PF_input(3,:) = 0;

%Define Solenoid current waveform vertices :: Startup - Ramp
%I_Sol_start = 450;  %[kA]  %OVERRIDE IF REQUIRED
I_PF_input([2,3],iSol) = I_Sol_start;
I_PF_input(4,iSol) = 0;
I_PF_input(5,iSol) = 0;

%Define coilset current waveforms vertices :: Ramp - Equilibrium
I_PF_input(4,iPF2) = icoil.PF2;
I_PF_input(4,iPF3) = icoil.PF3;
I_PF_input(4,iDiv1) = icoil.Div1;
I_PF_input(4,iDiv2) = icoil.Div2;

%Define coilset current waveforms vertices :: Equilibrium - Finish
I_PF_input(5,iPF2) = icoil.PF2;
I_PF_input(5,iPF3) = icoil.PF3;
I_PF_input(5,iDiv1) = icoil.Div1;
I_PF_input(5,iDiv2) = icoil.Div2;

%All coil currents end at zero
I_PF_input(6,iSol) = 0;
I_PF_input(6,iPF2) = 0;
I_PF_input(6,iPF3) = 0;
I_PF_input(6,iDiv1) = 0;
I_PF_input(6,iDiv2) = 0;

%%%%%%%%%%

%Convert from 'time' to 'long-time' for increased temporal resolution
time = [ 0 0.04 tau_vessel tau_vessel+0.025 tau_vessel+0.125 tau_vessel+0.130];
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
coil_names{iPF2} = 'PF2';
coil_names{iPF3} = 'PF3';
coil_names{iDiv1} = 'Div1';
coil_names{iDiv2} = 'Div2';
PF_colors{iSol} = 'Red';
PF_colors{iPF2} = 'Magenta';
PF_colors{iPF3} = 'Black';
PF_colors{iDiv1} = 'Cyan';
PF_colors{iDiv2} = 'Green';

%Plot figure showing coil startup waveforms - No eddy currents?
figure;
plot(time_long, I_PF_input_long);
title(gca,'SMART Startup Coil Current Waveforms, V3-phase1');
LegendString = {'Sol','PF2','PF3','Div1','Div2'};
legend(gca,LegendString);
xlabel(gca,'Time (s)');
ylabel(gca,'Coil Current (kA)');
Filename = '_StartupCurrentWaveforms';
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.png'));
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.pdf'));

%Compute dynamic coil currents employing current driven Ip - state_space_including_passive_elements_v4
[ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, time_adaptive ] = ...
    state_space_including_passive_elements_v4( curlyM, curlyR, time_long, I_PF_input_long, V_PF_input_long, Ip_long, Vp_long, 'adaptive_timesteping',true );

iTime_plasma = time_adaptive>0;
Vp_output(iTime_plasma) = 0;
Vp_long = interp1( time_adaptive,Vp_output, time_long);
Ip_long = NaN*Vp_long;

%Compute dynamic coil currents employing voltage driven Ip - state_space_including_passive_elements_v4
[ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, time_adaptive ] = ...
    state_space_including_passive_elements_v4( curlyM, curlyR, time_long, I_PF_input_long, V_PF_input_long, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names', coil_names, 'show_plot',true, 'turns',turns, 'currentScale',1e3, 'PF_colors',PF_colors );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plot non optimised - No Eddy Currents?  %%%%%%%%%%%%%%%%%%%%%%%%%%%
coil_currents = zeros(1,nPF);
coil_currents(iSol) = I_Sol_start;

icoil = fiesta_icoil( coilset, coil_currents );
equil_non_optimised = fiesta_equilibrium( 'ST25D non optimised', config, Irod, icoil );

figure;
plot(equil_non_optimised);
title(gca,'SMART Initial-Equilibrium, V3-phase1');
legend(gca,'hide');
set(gca,'XLim',[0 1]);
set(gca,'YLim',[-1.5 1.5]);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_InitEquilibrium';
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.png'));
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Make virtual sensors (Breakdown or Bp?)  %%%%%%%%%%%%%%%%%%%%%

%Create Virtual Sensors to check Bp
BP_virt_R = linspace(0.45,0.55,10);
BP_virt_Z = linspace(-0.05,0.05,10);

[BP_virt_R,BP_virt_Z] = meshgrid(BP_virt_R,BP_virt_Z);
BP_virt_R = BP_virt_R(:)';
BP_virt_Z = BP_virt_Z(:)';

BP_virt_theta = zeros(1,length(BP_virt_R));
nSensors = length(BP_virt_theta);

BP_virt_names = {};
for iSensor=1:nSensors
    BP_virt_names{iSensor} = ['Radial Bp Virtual Sensor #' num2str(iSensor) ];
end

BP_virt_R = [BP_virt_R  BP_virt_R];
BP_virt_Z = [BP_virt_Z  BP_virt_Z];

BP_virt_theta = [BP_virt_theta  BP_virt_theta+pi/2];

for iSensor=nSensors+1:2*nSensors
    BP_virt_names{iSensor} = ['Vertical Bp Virtual Sensor #' num2str(iSensor) ];
end

sensor_btheta = fiesta_sensor_btheta( 'sensor', BP_virt_R, BP_virt_Z,BP_virt_theta, BP_virt_names );

%Virtual Sensors gives 'Red Square' breakdown image.
figure; hold on; axis equal;
plot(vessel);
plot(sensor_btheta);
title(gca,'SMART VirtualSensors, V3-phase1');
legend(gca,'hide');
set(gca,'XLim',[0 1]);
set(gca,'YLim',[-1.5 1.5]);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_VirtualBSensors';
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.png'));
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  RZIP with virtual sensors  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rzip_config = fiesta_rzip_configuration( 'ST25_RZIP', config, vessel, {sensor_btheta} );
plasma_resistance = 5.94e-6;
[epsilon, B, C, D, curlyM, curlyR, gamma, plasma_parameters, index, label_index, state] = response(rzip_config, equil, 'rp',plasma_resistance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Optimised null - With Eddy Currents?  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_temp = C(end-get(sensor_btheta,'n')+1:end,1:nPF);
C1 = C_temp(:,1);
D1 = C_temp(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Vessel time constant  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,tau_vessel] = eig(curlyR(1:end-3,1:end-3)\curlyM(1:end-3,1:end-3));
tau_vessel = max(diag(tau_vessel));
disp([ 'tau_vessel=' num2str(tau_vessel*1e3) 'ms' ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plot optimised null - With Eddy Currents?  %%%%%%%%%%%%%%%%%%%%%%%%
coil_currents = zeros(1,nPF);
%I_Sol_start = 200  %[kA]  %OVERRIDE IF REQUIRED
I_PF_null = -pinv(D1) * (C1*I_Sol_start);
coil_currents(iSol) = I_Sol_start;
coil_currents(2:end) = I_PF_null';

icoil = fiesta_icoil( coilset, coil_currents );
equil_optimised_null = fiesta_equilibrium( 'ST25D optimised null', config, Irod, icoil );

figure;
plot(equil_optimised_null);
title(gca,'SMART Optimised-Equilibrium, V3-phase1');
legend(gca,'hide');
set(gca,'XLim',[0 1]);
set(gca,'YLim',[-1.5 1.5]);
xlabel(gca,'R (m)');
ylabel(gca,'Z (m)');
Filename = '_OptimisedEquilibrium';
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.png'));
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Plot Final Current Waveforms  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nTime = 6;
V_PF_input = NaN(nTime,nPF);
I_PF_input = zeros(nTime,nPF);

%All PF coil currents initiate at zero
I_PF_input(2,:) = 0;
I_PF_input(3,:) = 0;

%Define Solenoid current waveform vertices :: Startup --> Ramp
%I_Sol_start = 1500;  %[kA]  %OVERRIDE IF REQUIRED
I_PF_null = -pinv(D1) * (C1*I_Sol_start);
I_PF_input([2,3],iSol) = I_Sol_start;
I_PF_input(2,2:end) = I_PF_null;
I_PF_input(3,2:end) = I_PF_null;

%Define coilset current waveforms vertices :: Ramp --> Equilibrium
I_PF_input(4,iSol) = I_Sol_ramp;    %-1100
I_PF_input(4,iPF2) = icoil.PF2;     %-1100;
I_PF_input(4,iPF3) = icoil.PF3;;
I_PF_input(4,iDiv1) = icoil.Div1;
I_PF_input(4,iDiv2) = icoil.Div2;

%Define coilset current waveforms vertices :: Equilibrium --> Finish
I_PF_input(5,iSol) = I_Sol_equil;   %-1500
I_PF_input(5,iPF2) = icoil.PF2;     %-1100;
I_PF_input(5,iPF3) = icoil.PF3;
I_PF_input(5,iDiv1) = icoil.Div1;
I_PF_input(5,iDiv2) = icoil.Div2;

%All coil currents end at zero
I_PF_input(6,iSol)=0;
I_PF_input(6,iPF2)=0;
I_PF_input(6,iPF3)=0;
I_PF_input(6,iDiv1)=0;
I_PF_input(6,iDiv2)=0;

%%%%%%%%%%

%Convert from 'time' to 'long-time' for increased temporal resolution
tstep = TauPulse  %[s]  %OVERRIDE IF REQUIRED
time = [-0.05 -0.03 0 tstep tstep+0.03 tstep+0.05];
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
coil_names{iPF2} = 'PF2';
coil_names{iPF3} = 'PF3';
coil_names{iDiv1} = 'Div1';
coil_names{iDiv2} = 'Div2';
PF_colors{iSol} = 'Red';
PF_colors{iPF2} = 'Magenta';
PF_colors{iPF3} = 'Black';
PF_colors{iDiv1} = 'Cyan';
PF_colors{iDiv2} = 'Green';

%Plot figure showing dynamic coil currents - With Eddy Currents?
figure;
plot(time, I_PF_input );
title(gca,'SMART Final Coil Current Waveforms, V3-phase1');
LegendString = {'Sol','PF2','PF3','Div1','Div2'};
legend(gca,LegendString);
xlabel(gca,'Time (s)');
ylabel(gca,'Coil Current (kA)');
Filename = '_FinalCurrentWaveforms';
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.png'));
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.pdf'));

%Compute dynamic coil currents employing current driven Ip - state_space_including_passive_elements_v4
[ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, time_adaptive ] = ...
    state_space_including_passive_elements_v4( curlyM, curlyR, time_long, I_PF_input_long, V_PF_input_long, Ip_long, Vp_long, 'adaptive_timesteping',true );

iTime_plasma = time_adaptive>0;
Vp_output(iTime_plasma) = 0;
Vp_long = interp1( time_adaptive,Vp_output, time_long);
Ip_long = NaN*Vp_long;

%Compute dynamic coil currents employing voltage driven Ip - state_space_including_passive_elements_v4
[ V_PF_output, I_PF_output, I_Passive, Vp_output, Ip_output, figure_handle, matlab2tikz_extraAxisOptions, uFinal, time_adaptive ] = ...
    state_space_including_passive_elements_v4( curlyM, curlyR, time_long, I_PF_input_long, V_PF_input_long, Ip_long, Vp_long, 'adaptive_timesteping',true, 'coil_names', coil_names, 'show_plot',true, 'turns',turns, 'currentScale',1e3, 'PF_colors',PF_colors );

%% Inner Solenoid
% Sol Start= 4700 A
% Sol 4 =500
% Sol 5= -4700
% t= 0.025
% time = [-0.1 -0.05 0 t 0.1+t 0.11+t];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Outer Solenoid
% Sol Start= 2165.6 A
% Sol 4 =240.625
% Sol 5= -2165.6
% t= 0.03
% time = [-0.1 -0.05 0 t 0.1+t 0.105+t];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Plasma Current Trace
close all
plot(time_adaptive,Ip_output./1000)
title(gca,'SMART Plasma Current, V3-phase1');
legend(gca,'Plasma Current');
xlabel(gca,'Time (s)');
ylabel(gca,'Plasma Current (kA)');
Filename = '_FinalPlasmaCurrent';
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.png'));
saveas(gcf, strcat(FolderName,ProjectName,Filename,'.pdf'));

%% Write Coil Waveforms, Plasma Current and Passive Current to textfiles

%Scale currents by ???FACTORS??? for ???REASONS???
a=I_PF_output(:,1).*800;    %Sol?
b=I_PF_output(:,2).*24;     %PF2?
c=I_PF_output(:,3).*24;     %PF3?
d=I_PF_output(:,4).*24;     %Div1?
e=I_PF_output(:,5).*24;     %Div2?

Filename = strcat(FolderName,'CoilCurrents_Phase_1.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%0.5f %0.5f %0.5f %0.5f %0.5f\r\n',[a'; b'; c'; d'; e']);

Filename = strcat(FolderName,'t_Phase_1.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%1.12f\r\n',time_adaptive);

Filename = strcat(FolderName,'Ip_Phase_1.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%1.12f\r\n',Ip_output);

Filename = strcat(FolderName,'I-Passive_Phase_1.txt');
fileID=fopen(Filename,'w');
fprintf(fileID,'%1.12f\r\n',I_Passive);
