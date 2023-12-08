%   Copyright:  Gabriel Perez
%   Repository : Sanitary Sewer - WFIUH
%   Email:   perezmesagj@ornl.gov
%	Last update: 07/16/2023,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%   Perez, G., Gomez-Velez, J. D., & Grant, S. B. (2023). 
%   The sanitary sewer unit hydrograph model: A comprehensive tool for wastewater flow modeling and inflow-infiltration simulations. 
%   Water Research, 120997. https://doi.org/https://doi.org/10.1016/j.watres.2023.120997


function [QT_BWF,wf_t,Weight_Nodes,Qin_BWF] = computeBWF_WIUH(uc,D,Delta_t,Sewer_T,ID_outlet,Q_DF)

    %% Inputs
    % This function computes the BWF with the WFIUH using the following inputs:
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

    % uc : Celerity [m/s]    
    % D : Dispersion coefficient [m2/s]  
    % Delta_t : Time step [s]
    % ID_outlet : ID of the outlet which the accumualted flow is estimated
    % Q_DF :Diurnal demand factor [-]


    %% Outputs
    % QT_BWF :  Total BWF at the outlet [m3/s]
    % wf_t :    Weighted Width Function (Transfer function) []. 
    % Weight_Nodes: Weight representing the number of times of Qmin at each node [Count of Qmin]
    % Qin_BWF:  Input Q_BWF [m3/s per Qmin]

    %% Function
    % Define the Qmin in the system and the weights for each node
    Qmin=Sewer_T.BWF_up_cms(Sewer_T.isHead==1);
    Qmin=min(Qmin(Qmin~=0)); %[m3/s] % 
    % Weight representing the number of times of Qmin at each node [Count of Qmin]
    Weight_Network=Sewer_T.BWF_up_cms./Qmin;
    Qin_BWF=Q_DF.*Qmin; % Input Q_BWF [m3/s per Qmin]
    % Note that Qmin can be replaced by the BWF generated from a single users,
    % doing that Weight_Network will have units of "equivalent users" at each node
    % and Qin_BWF will be m3/s per "equivalent user". We used the term of "equivalent user"
    % because the Sewer_T.BWF_up_cms is calculated from different water 
    % users including residential users and employees.


    % Get sub-sewershed and identify the headers
    LinkID_Subshed=getSubshed(ID_outlet,Sewer_T.ID_Link,Sewer_T.Next_ID); 
    LinkID_Headers=LinkID_Subshed==1 & Sewer_T.isHead==1;
    % Define Parameters
    Vector_L_Network=Sewer_T.Acum_L_m;
    Vector_L=Vector_L_Network(LinkID_Headers);
    % This is the distance from outlet of interest to the outlet of the entire network sewer network
    L_acum_outlet=min(Vector_L_Network(LinkID_Subshed==1)); 
    % This is the pipe lenght related to the pipe connected to the outlet of interest
    L_pipe=Sewer_T.L_m(Sewer_T.ID_Link==ID_outlet);  
    L_translation=L_acum_outlet-L_pipe;
    Vector_L=Vector_L-L_translation; % Translate distances to the outlet of analysis
    Weight_Nodes=Weight_Network(LinkID_Headers);

    % Compute the BWF
    [QT_BWF,wf_t,t_max] = computeWIUH(uc,D,Delta_t,Qin_BWF,Vector_L,Weight_Nodes);

end