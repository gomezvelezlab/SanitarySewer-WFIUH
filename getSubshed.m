%   Copyright:  Gabriel Perez
%   Repository : Sanitary Sewer - WFIUH
%   Email:   perezmesagj@ornl.gov
%	Last update: 07/16/2023,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%   Perez, G., Gomez-Velez, J. D., & Grant, S. B. (2023). 
%   The sanitary sewer unit hydrograph model: A comprehensive tool for wastewater flow modeling and inflow-infiltration simulations. 
%   Water Research, 120997. https://doi.org/https://doi.org/10.1016/j.watres.2023.120997


function [LinkID_Subshed] = getSubshed(ID_outlet,ID_Link,Next_ID)

% Function to get the IDs that belong to a specific subshed defined at the outlet identified
% with ID_outlet

Flag=1;
LinkID_Subshed=ismember(ID_Link,ID_outlet)./1; % Set the outlet as part of the subshed
while Flag==1
   up_LinkID=ID_Link(ismember(Next_ID,ID_outlet)); % Get the upstream LinkID
   LinkID_Subshed=ismember(Next_ID,ID_outlet)./1+LinkID_Subshed;
   if isempty(up_LinkID)==1
       Flag=0;
   else
   ID_outlet=up_LinkID; % Update the link IDs for a new iteration
   end
end

end