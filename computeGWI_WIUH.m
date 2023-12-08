%   Copyright:  Gabriel Perez
%   Repository : Sanitary Sewer - WFIUH
%   Email:   perezmesagj@ornl.gov
%	Last update: 07/16/2023,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%   Perez, G., Gomez-Velez, J. D., & Grant, S. B. (2023). 
%   The sanitary sewer unit hydrograph model: A comprehensive tool for wastewater flow modeling and inflow-infiltration simulations. 
%   Water Research, 120997. https://doi.org/https://doi.org/10.1016/j.watres.2023.120997


function [QT_GWI,wf_t,Weight_Network_qii,Qin_GWI] = computeGWI_WIUH(uc,D,Delta_t,Sewer_T,ID_outlet,GWI,Type_GWI,Vector_Qnodes)

    % Inputs:
    % This function computes the GWI with the WFIUH using the following  inputs
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

    % ID_outlet:    ID of the link outlet which the accumulated flow is estimated
    % GWI:          Total Groundwater flow at the outlet of the system [m3/s]. 
    % Note that  GWI input can also be a time series represneting time variation of GWI 
    % Type_GWI string to set how to assign the GWI at each pipe. 
    % "Pipe_Length":The GWI is define by pipe length [m3/s per m]
    % "Pipe_Area":  The GWI is defined by pipe area (pipe perimeter * pipe length) [m3/s per m2]
    % Vector_Qnodes: This is just used to get the same size of Vector_Qnodes

    %% Output
    % QT_GWI:       Total Groundwater at the outlet of the sewer network [m3/s]
    % wf_t:         Weighted Width Function (Transfer function) []. 
    % Weight_Network_qii: Weight representing the number of times of qii_min at each node [Count of qii_min]
    % Qin_GWI:      Input Q_GWI [m3/s per qii_min]

    %% Function
    % Define Input and weights to create transfer function
    switch Type_GWI
        case 'Pipe_Length'
            q0=1; % Reference GWI contribution by unit pipe lenght [m3/s/m]
            qii=Sewer_T.L_m.*q0;  % [m3/s]
            qii_min=min(qii(qii~=0)); % [m3/s]
            Weight_Network=qii./qii_min; % Weight representing the number of times of qii_min at each node [Count of qii_min]
        case 'Pipe_Area'
            q0=1; % Reference GWI contribution by unit area [m3/s/m2]
            qii=(Sewer_T.Ah_km2.*(1000^2)).*q0; % [m3/s]
            qii_min=min(qii(qii~=0)); % [m3/s]
            Weight_Network=qii./qii_min; % Weight representing the number of times of qii_min at each node [Count of qii_min]
        otherwise
            error('Error: Wrong ID Type_GWI')
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
    Q_ref_GWI=GWI./sum(Weight_Network_qii); % Input Q_GWI [m3/s per qii_min]
    if length(Q_ref_GWI)>1
        Qin_GWI=Q_ref_GWI; % Input Q_GWI [m3/s per qii_min]
    elseif length(Q_ref_GWI)==1
        Qin_GWI=Q_ref_GWI.*ones(1,length(Vector_Qnodes)); % Input Q_GWI [m3/s per qii_min]
    end

    % Compute the GWI
    if isempty(Weight_Network_qii)
        QT_GWI = []; wf_t =[]; t_max = [];
    else
        [QT_GWI,wf_t,t_max] = computeWIUH(uc,D,Delta_t,Qin_GWI,Vector_L,Weight_Network_qii);
    end

    QT_GWI =QT_GWI';

end













