%   Copyright:  Gabriel Perez
%   Repository : Sanitary Sewer - WFIUH
%   Email:   perezmesagj@ornl.gov
%	Last update: 07/16/2023,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%   Perez, G., Gomez-Velez, J. D., & Grant, S. B. (2023). 
%   The sanitary sewer unit hydrograph model: A comprehensive tool for wastewater flow modeling and inflow-infiltration simulations. 
%   Water Research, 120997. https://doi.org/https://doi.org/10.1016/j.watres.2023.120997


function [Q_conv, wf_conv, t_max] = computeWIUH(uc,D,Delta_t,Vector_Qnodes,Vector_L,Weight_Nodes)

    % This function computes the flow response using the transfer
    % function estimated from the Weighted Width Function 
    % calculated for a sewer network. This function is intended to be used
    % with the functions computeBFW_WIUH.m computeGWI_WIUH.m, and computeRDII_WIUH.m

    %% Input variables
    % uc: Flow celerity [m/s]
    % D: Dispersion coefficient [m2/s]
    % Delta_t: Time step [s]
    % Vector_Qnodes: Reference flow for a single node [m3/s per Qmin]
    % Vector_L: Distance of each node to the outlet of the system [m]
    % Weight_Nodes: Weight representing the number of times of Qmin at each node [Count of Qmin]

    %% Output variables
    % Q_conv: Flow at the outlet of the system [m3/s]. This is a vector with Q printed every Delta_t
    % wf_conv: Weighted Width Function (Transfer function) []. 
    % t_max: Max time after the convolution is done. 

    %% Function
    % Check ranges of input variables
    if uc <= 0 || D < 0
        error('Error: Check uc and D values')
    end
    
    if D==0 % Kinematic Case      
        wf_x = Vector_L;                      % This is the vector to construct the width function w(x)
        wf_t = wf_x./uc;                      % Convert distance to time, Kinematic case
        t_max = max(wf_t);                    % Find the max travel time in the system to define wf_t histogram
        Edges_histogram = 0:Delta_t:t_max;    % Find the intervals to do the analysis
        t_max = max(Edges_histogram)+Delta_t; % Define max time
        % Find the weighted nodes at each Histogram Interval
        Pos = 1;
        N_Weight_Nodes = zeros(1,size(0:Delta_t:(t_max-Delta_t),2)); 
        for t = 0:Delta_t:(t_max-Delta_t)
            N_Weight_Nodes(Pos) = sum(Weight_Nodes(wf_t>=t & wf_t<(t+Delta_t)));
            Pos = Pos+1;
        end
        NW_Total = sum(N_Weight_Nodes);       % Total Number accounting weighted nodes
        wf_conv = N_Weight_Nodes./NW_Total;   % Transfer function. This is normalized and the sum must be equal to 1
        Q_conv = NW_Total.*conv(wf_conv,Vector_Qnodes(1:end-1));  % Flow at the outlet of the system
    
    else % Including attenuation with dispersion coefficient
        Delta_x = uc*Delta_t;
        wf_x = Vector_L;                      % This is the vector to construct the width function w(x)
        x_max = max(wf_x);                    % Find the max travel time in the system to define wf_t histogram
        Edges_histogram = 0:Delta_x:x_max;    % Find the intervals to do the analysis
        x_max = max(Edges_histogram)+Delta_x; % Find max distance
        Pos = 1;
        N_Weight_Nodes = zeros(1,size(0:Delta_x:(x_max-Delta_x),2)); 
        for t = 0:Delta_x:(x_max-Delta_x)
            N_Weight_Nodes(Pos) = sum(Weight_Nodes(wf_x>=t & wf_x<(t+Delta_x)));
            Pos = Pos+1;
        end
        NW_Total = sum(N_Weight_Nodes);       % Total Number accounting weighted nodes
        wf_x = N_Weight_Nodes./NW_Total;      % This is the normalized width function    
        t_max = x_max/uc;                     % Define max time
        Distance_X = [(Edges_histogram(2:end)+Edges_histogram(1:(end-1)))/2, Edges_histogram(end) + Delta_x/2];
        wf_conv = zeros(size(0:Delta_t:round(t_max),2),1);
        Pos = 1;
        for ii = 0:Delta_t:t_max  
            wf_conv(Pos,1)=sum(((Distance_X.*wf_x)./sqrt(4*pi*D*(ii^3))).*(exp(((Distance_X-uc.*ii).^2)./(-4*D*ii))), 'omitnan');      
            Pos=1+Pos;
        end
        wf_conv=wf_conv./(sum(wf_conv.*Delta_t));        % Transform to a probability distribution             
        wf_conv=wf_conv.*Delta_t;                        % Must be equal to 1  
        Q_conv=NW_Total.*conv(wf_conv,Vector_Qnodes(1:end-1)); % Flow at the outlet of the system
    end

end