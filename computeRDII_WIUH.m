%   Copyright:  Gabriel Perez
%   Repository : Sanitary Sewer - WFIUH
%   Email:   perezmesagj@ornl.gov
%	Last update: 07/16/2023,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%   Perez, G., Gomez-Velez, J. D., & Grant, S. B. (2023). 
%   The sanitary sewer unit hydrograph model: A comprehensive tool for wastewater flow modeling and inflow-infiltration simulations. 
%   Water Research, 120997. https://doi.org/https://doi.org/10.1016/j.watres.2023.120997


function [QT_RDII,wf_t_f,wf_t_s,Weight_Network_qii,Qin_RDII_f,Qin_RDII_s] = computeRDII_WIUH(uc_f,D_f,p_fraction_f,uc_s,D_s,p_fraction_s,Delta_t,Sewer_T,ID_outlet,Rainfall,Type_RDII)

    % The model accounts for a fast (f) and slow (s) flow response
    %% Inputs
    % This function computes the RDII with the WFIUH using the following inputs
    % -- Sewer_T -- : Table the below fields
    % ID_Link :     Unique indentifier of the sewer link
    % Next_ID :     Identifier of the downstream sewer link
    % Acum_L_m :    Total (or accumulated) sewer pipe length to the sewershed outlet [m]
    % L_m :         Sewer pipe length [m]
    % Ah_km2 :      Hillslope area or contributing area to the sewer pipe [km2]
    % isHead :      Binary value to identify if the sewer link is a header
    % isOutlet :    Binary value to identify if the sewer link is a sewershed outlet
    % BWF_up_cms:   Daily mean Base Wastewater Flow calculated from the users upstream of the sewer pipe [m3/s]
    % ID_outlet:    ID of the outlet which the accumualted flow is estimated
    % Vector_Qnodes:Reference flow for a single node and a single user [m3/s per Qmin]
    % Weight_Nodes: Weight representing the number of times of Qmin at each node [Count of Qmin]

    % uc_s: celerity related to slow flow [m/s]
    % D_s: Dispersion coefficient related to slow flow [m2/s]
    % uc_f: celerity related to fast flow [m/s]
    % D_f: Dispersion coefficient related to fast flow [m2/s]
    % Delta_t: Time step [s]
    % ID_outlet: ID of the outlet which the accumulated RDII flow is estimated
    % Rainfall:  Rainfall volume time series in the sewershed  [m3/s] 
    % p_fraction_s: Fraction of precipitation that becomes RDII for slow flow [-]
    % p_fraction_f: Fraction of precipitation that becomes RDII for fast flow [-]
    % Type_RDII: String to set how to assign the RDII at each pipe. 
    % "Pipe_Length": The RDII is define by pipe lenght [m3/s per m]
    % "Pipe_Area": The RDII is defined by pipe area (pipe perimeter * pipe length) [m3/s per m2]


    %% Outputs
    % QT_RDII: RDII flow [m3/s]
    % wf_t_f:  Weighted Width Function for fast flow (Transfer function) []. 
    % wf_t_s:  Weighted Width Function for slow flow (Transfer function) []. 
    % Weight_Network_qii: Weight representing the number of times of qii_min at each node [Count of qii_min]
    % Qin_RDII_f: Input Q_RDII fast flow  [m3/s]
    % Qin_RDII_s: Input Q_RDII fast flow  [m3/s]


    %% Function
    % Define Input and weights to create transfer function
    RDII_f =  Rainfall*p_fraction_f; % [m3/s]
    RDII_s =  Rainfall*p_fraction_s; % [m3/s]

    switch Type_RDII
        case 'Pipe_Length'
            q0=1; % Reference RDII contribution by unit pipe lenght [m3/s/m]
            qii=Sewer_T.L_m.*q0;  %[m3/s]
            qii_min=min(qii(qii~=0)); %[m3/s]
            Weight_Network=qii./qii_min; % Weight representing the number of times of qii_min at each node [Count of qii_min]
        case 'Pipe_Area'
            q0=1; % Reference RDII contribution by unit area [m3/s/m2]
            qii=(Sewer_T.Ah_km2.*(1000^2)).*q0;  %[m3/s]
            qii_min=min(qii(qii~=0)); %[m3/s]
            Weight_Network=qii./qii_min; % Weight representing the number of times of qii_min at each node [Count of qii_min]
        otherwise
            error('Error: Wrong ID Type_RDII')
    end

    % Get sub-sewershed
    LinkID_Subshed=getSubshed(ID_outlet,Sewer_T.ID_Link,Sewer_T.Next_ID); 
    Vector_L_Network=Sewer_T.Acum_L_m;
    % This is the distance from outlet of interest to the outlet of the entire network sewer network
    L_acum_outlet=min(Vector_L_Network(LinkID_Subshed==1)); 
    % This is the pipe lenght related to the pipe connected to the outlet of interest
    L_pipe=Sewer_T.L_m(Sewer_T.ID_Link==ID_outlet);  
    L_translation=L_acum_outlet-L_pipe;
    Vector_L=Vector_L_Network(Sewer_T.Ah_km2~=0 & LinkID_Subshed==1);
    Vector_L=Vector_L-L_translation; % Translate distances to the outlet


    Weight_Network_qii=Weight_Network(Sewer_T.Ah_km2~=0 & LinkID_Subshed==1);

    Qin_RDII_f=RDII_f./sum(Weight_Network_qii); % Input Q_RDII fast flow [m3/s per qii_min]

    Qin_RDII_s=RDII_s./sum(Weight_Network_qii); % Input Q_RDII slow flow [m3/s per qii_min]

    % Compute the RDII
    if isempty(Weight_Network_qii)
        QT_RDII_f = []; wf_t_f =[]; t_max_f = [];
        QT_RDII_s = []; wf_t_s =[]; t_max_s = [];
    else
        [QT_RDII_f,wf_t_f,t_max_f] = computeWIUH(uc_f,D_f,Delta_t,Qin_RDII_f,Vector_L,Weight_Network_qii);
        if uc_s==0
            QT_RDII_s = 0; wf_t_s =0; t_max_s = 0;
        else
        [QT_RDII_s,wf_t_s,t_max_s] = computeWIUH(uc_s,D_s,Delta_t,Qin_RDII_s,Vector_L,Weight_Network_qii);
        end
        
    end
    % Convert in vector column
    QT_RDII_s = QT_RDII_s(:); % [m3/s]
    QT_RDII_f = QT_RDII_f(:); % [m3/s]

    % Sum of responses
    if length(QT_RDII_s)>length(QT_RDII_f)
        QT_RDII = QT_RDII_s; % [m3/s]
        QT_RDII(1:length(QT_RDII_f)) = QT_RDII_f + QT_RDII(1:length(QT_RDII_f));
    else
        QT_RDII = QT_RDII_f; % [m3/s]
        QT_RDII(1:length(QT_RDII_s)) = QT_RDII_s + QT_RDII(1:length(QT_RDII_s));
    end
        
end













