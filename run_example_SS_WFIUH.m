%   Copyright:  Gabriel Perez
%   Repository : Sanitary Sewer - WFIUH
%   Email:   perezmesagj@ornl.gov
%	Last update: 07/16/2023,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%   Perez, G., Gomez-Velez, J. D., & Grant, S. B. (2023). 
%   The sanitary sewer unit hydrograph model: A comprehensive tool for wastewater flow modeling and inflow-infiltration simulations. 
%   Water Research, 120997. https://doi.org/https://doi.org/10.1016/j.watres.2023.120997

%% Example to run the SS-WFIUH
% 
% The SS-WFIUH model consists on three functions to calculate 
% Base Wastewater Flow (BWF), Groundwater Infiltration (GWI), 
% and Rainfall-Derived Infiltration and Inflow (RDII):
% 
% computeBWF_WIUH.m --> It calculates the BFW
% computeGWI_WIUH.m --> It calculates the GWI
% computeRDII_WIUH.m --> It calculates the RDII
% 
% This script shows examples to calculate the three components of sewer flow
% for different "events" in the Cub Run sewer network located in Fairfax County, 
% Virginia, U.S.A. 



%% Input Data: 

% The model only requires a sewer network (polyline shapefile) 
% with the below attributes and a demand factor representing the diurnal cycle 
% of BWF contribution by the water users. 

% --Shapefile-- of the sewer network with the below attributes 
% ID_Link :     Unique indentifier of the sewer link
% Next_ID :     Identifier of the downstream sewer link
% Acum_L_m :    Total (or accumulated) sewer pipe length to the sewershed outlet [m]
% L_m :         Sewer pipe length [m]
% Ah_km2 :      Hillslope area or contributing area to the sewer pipe [km2]
% isHead :      Binary value to identify if the sewer link is a header
% isOutlet :     Binary value to identify if the sewer link is a sewershed outlet
% BWF_up_cms:   Daily mean Base Wastewater Flow calculated from the users upstream of the sewer pipe [m3/s]

% --Demand factor--
% A time series representing the demand diurnal factors (weekdays and weekends) of BWF.
% we provided those on the folder "Input_data/Diurnal_Cycles/"

% -- Rainfall--
% A time series with the rainfall intensity within the sewershed.

% -- Total sewege flow observations-- This is not needed to run the model. 
% We used this for validation and comparison purposes.

%% Define general parameters
Delta_t = 60;                       % Time step for convolution [seconds]

%% Select any period between '01-Jan-2010' and '31-Dec-2019'
Start_Period = '18-Jul-2018';       % Use format 'DD-MMM-YYYY'
End_Period = '31-Jul-2018';         % Use format 'DD-MMM-YYYY'

% More examples
% Start_Period = '07-Sep-2011';       % Use format 'DD-MMM-YYYY'
% End_Period = '17-Sep-2011';         % Use format 'DD-MMM-YYYY'
% % 
% Start_Period = '09-Jul-2013';       % Use format 'DD-MMM-YYYY'
% End_Period = '20-Jul-2013';         % Use format 'DD-MMM-YYYY'
% % 
% Start_Period = '10-May-2014';       % Use format 'DD-MMM-YYYY'
% End_Period = '23-May-2014';         % Use format 'DD-MMM-YYYY'


%% Load sewer network shapefile
Path=pwd;
S = shaperead([Path '\Input_data\GIS_Data\Sewer_Network_UOSA.shp']);
% Define sewer link to compute flows
ID_outlet = 14972; % This is the LinkID where is located the Cub Run pump station 

%% Create table with sewer network attributes
Sewer_T = table(cell2mat({S(:).ID_Link})',cell2mat({S(:).Next_ID})',cell2mat({S(:).Acum_L_m})', cell2mat({S(:).L_m})',...
    cell2mat({S(:).Ah_km2})', cell2mat({S(:).isHead})', cell2mat({S(:).isOutlet})',  cell2mat({S(:).BWF_up_cms})', ...
    'VariableNames',{'ID_Link','Next_ID','Acum_L_m','L_m','Ah_km2','isHead','isOutlet','BWF_up_cms'});

%% Load the diurnal demand factor for the system. 
% This can be obtained from local information or literature
% This series was estimated based on the method presented by Perez et al.,2023
% This time series must have a time step equal than Delta_t
Table_DF = load('Input_data\Diurnal_Cycle\Table_DF.mat'); Table_DF = Table_DF.Table_DF;
Q_DF = Table_DF.Q_DF( Table_DF.Date>Start_Period &  Table_DF.Date<End_Period);

%% Load hourly Rainfall and Sewage discharge observations
% The hourly rainfall represents the rainfall intensity within the sewershed of analysis.
% The hourly sewage discharge represents the total sewage discharge measured at the outlet of the system.
% Note that the sewage flow observations can have measurement errors and/or missing data (eg. Aug 26, 2019)
Data_P_Q = load('Input_data\Observations\Data_P_Q.mat'); 
Data_P_Q = Data_P_Q.Data_P_Q;
QT_obs = Data_P_Q.Q_obs_cms(Data_P_Q.Date_obs>Start_Period & Data_P_Q.Date_obs<End_Period); % [m3/s]
Rainfall_obs = Data_P_Q.P_mm_hr(Data_P_Q.Date_obs>Start_Period & Data_P_Q.Date_obs<End_Period); % [mm/hr]
Date_obs = Data_P_Q.Date_obs(Data_P_Q.Date_obs>Start_Period & Data_P_Q.Date_obs<End_Period); % 

Area_sewershed = 80 * 1000^2; %[m2] Approximated area for Cubrun sewershed
Rainfall = Rainfall_obs; % [mm/hr]
Rainfall = Rainfall./(1000*60*60); % [m/s].
Rainfall = Rainfall.*Area_sewershed; % Rainfall volume [m3/s]
% Interpolate rainfall for the same Delta_t
Rainfall=interp1((1:length(Rainfall)).*60*60,Rainfall,1:Delta_t:(length(Rainfall).*60*60),'spline');
Rainfall(Rainfall<0)=0; % Remove artifact of negative values generated during the interpolation

%% Model parameters 
% The below parameters were estimated using observations during dry and wet periods,
% by minimizing the sum of the square difference between the observed and the predicted 
% flows at the outlet of the system. (See Perez et al., 2023 for more details) 
% Parameters estimated by calibration using observations during dry periods)
uc = 0.9314;            % Celerity [m/s]    
D = 1.1993;             % Dispersion coefficient [m2/s]  
Qd_mean = 54.19;        % Mean water consume per per person [gallon/day] 
Type_GWI='Pipe_Length';
GWI = 0.0376;           % Mean accumulated graoundwater flow at the outlet [m3/s] 

% Parameters estimated by calibration during 100 wet events 
Type_RDII='Pipe_Length';
uc_f = 0.93;            % Celerity [m/s]    
D_f = 1.1993;           % Dispersion coefficient [m2/s]  
p_fraction_f = 0.0034; 
uc_s =  0.0827;            % Celerity [m/s]    
D_s = 110;            % Dispersion coefficient [m2/s]  
p_fraction_s = 0.2144*sum(Rainfall_obs)^(-0.6914);  % mean value 0.0168


%% Calculate BWF flow at the outlet
[QT_BWF_Mod_CR,wf_t,Weight_Nodes,Qin_BWF] = computeBWF_WIUH(uc,D,Delta_t,Sewer_T,ID_outlet,Q_DF);
% Check mass balance
Qin = sum(Weight_Nodes).*sum((Qin_BWF.*Delta_t)); % [m3]
Qout = sum(QT_BWF_Mod_CR.*Delta_t); % [m3],
ErrorBWF = 100.*(Qout-Qin)/Qin; % This must be close 0, if not, decrease Delta_t

%% Calculate GWI flow at the outlet
[QT_GWI_Mod_CR, wf_t, Weight_Nodes_GWI,Qin_GWI_Mod] = computeGWI_WIUH(uc,D,Delta_t,Sewer_T,ID_outlet,GWI,Type_GWI,Qin_BWF);
% Check mass balance
Qin = sum(Weight_Nodes_GWI).*sum((Qin_GWI_Mod.*Delta_t)); % [m3]
Qout = sum(QT_GWI_Mod_CR.*Delta_t); % [m3],
ErrorGWI = 100.*(Qout-Qin)/Qin; % This must be close to 0, if not, decrease Delta_t

%% Calculate RDII at the outlet 
[QT_RDII_Mod_CR, wf_t_f,wf_t_s,Weight_Network_qii,Qin_RDII_f,Qin_RDII_s] = computeRDII_WIUH(uc_f,D_f,p_fraction_f,uc_s,D_s,p_fraction_s,Delta_t,Sewer_T,ID_outlet,Rainfall,Type_RDII);
% Check mass balance
Qin_f = sum(Weight_Network_qii).*sum((Qin_RDII_f.*Delta_t)); % [m3]
Qin_s = sum(Weight_Network_qii).*sum((Qin_RDII_s.*Delta_t)); % [m3]
Qout = sum(QT_RDII_Mod_CR.*Delta_t); % [m3],
ErrorRDII = 100.*(Qout-(Qin_f+Qin_s))/(Qin_f+Qin_s); % This must be close to 0, if not, decrease Delta_t
% Keep the same vector length of the BWF and GWI
QT_RDII_Mod_CR = QT_RDII_Mod_CR(1:length(QT_BWF_Mod_CR));

%% Create a final matrix with model results
flow_components_WFIUH = [QT_GWI_Mod_CR, QT_BWF_Mod_CR, QT_RDII_Mod_CR];
Date_QT=datetime((datenum(Start_Period)+(0:Delta_t:(length(QT_BWF_Mod_CR)*Delta_t-Delta_t))./(1*24*60*60)),'ConvertFrom','datenum');

%% Do a nice figure showing flow components and observations
fig=figure; set(gcf,'color','white');
h1=subplot(1,1,1);
yyaxis right
bar(Date_obs,Rainfall_obs,'b'); hold on
ax = gca; ax.YAxis(2).Direction = 'reverse';
ylim([0 50]); xlim([Date_obs(1) Date_obs(end)]);
ylabel('Rain [mm/hr]');
set(h1, 'xlim', [datetime(Start_Period)+1 datetime(End_Period)-1]); 
set(gca,'FontSize',14);

yyaxis left
ax = gca; ax.YAxis(1).Direction = 'normal';
A=area(Date_QT,flow_components_WFIUH); hold on;
A(1).FaceColor = [254,178,76]./255;
A(2).FaceColor = [217,95,14]./255;
A(3).FaceColor = [49,130,189]./255;
plot(Date_obs,QT_obs ,'k','LineWidth',2); 
xlabel('Time'); ylabel('Sewage Discharge at Cub Run [m^3/s]');
legend({'GWI','BWF','RDII','Observations'})
xlim([datetime(Start_Period)+1 datetime(End_Period)-1])
set(gca,'FontSize',14); alpha 0.5
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';

