# SanitarySewer-WFIUH

This repo contains the functions to implement the Sanitary Sewer Width Function Instantaneous Unit Hydrograph (SS-WFIUH) 
described by Perez et al., 2023.

The model consists on three functions to calculate each sanitary sewer flow component:

SS_WFIUH_BWF.m --> It computes the base wastewater flow (BFW)
SS_WFIUH_GWI.m --> It computes the groundwater infiltration (GWI)
SS_WFIUH_RDII.m --> It computes the Rainfall-derived infiltration and inflow (RDII)

The code run_exampe_SS_WFIUH.m has an example for Cub Run sewer network is 
located in Fairfax County, Virginia, U.S.A. The sewershed in comprised of
around 10,000 sewer conduits and its outlet drains to the UOSA wastewater
treatment plant.

Copyright:  Gabriel Perez
Repository : SanitarySewer-WFIUH (SS-WFIUH)
Last update: 06/07/2023,   MATLAB	2019b  version
IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
Perez, G., Gomez-Velez, J. D., & Grant, S. B. (2023). 
The sanitary sewer unit hydrograph model: A comprehensive tool for wastewater flow modeling and inflow-infiltration simulations. 
Water Research, 120997. https://doi.org/https://doi.org/10.1016/j.watres.2023.120997 
