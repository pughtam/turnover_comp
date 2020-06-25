README_forest_masks.txt

Author: Sarah L. Shafer, U.S. Geological Survey (sshafer@usgs.gov)Last update: 2020-02-6
----------SUMMARY----------
This folder contains two Fortran 90 programs used to develop a forest mask and assign forest types
for analyses of terrestrial biosphere model (TBM) output as described in Pugh et al. (in review, 
Biogeosciences Discuss.). See Pugh et al. (in review) for additional details.


1. pft_to_forest_v2.f90
This program creates a forest mask from monthly or annual LAI simulated by TBMs. For each grid 
cell, forest is assigned if 1) annual total tree PFT LAI is >2.5 or 2) annual total tree PFT 
LAI is >0.5 and the PFT with the maximum LAI is a boreal tree PFT or is a tree PFT that could 
be boreal (i.e., needleleaved evergreen, needleleaved deciduous, broadleaved deciduous). The 
program creates a netCDF output file with annual data indicating whether a grid cell was 
classified as forest in each year during the period 1901-2014 (variable name: forest). 
A second variable in the netCDF output file indicates whether a grid cell was classified as 
forest for any 10 years during the period 1985-2014 (variable name: forest_30yr_any_10_years). 

2. pft_phenology_v2.f90
This program classifies TBM simulated plant functional types (PFTs) into forest types using 
simulated LAI (1985-2014 30-year mean) and PFT phenological traits. The program creates a 
netCDF output file with data indicating for each grid cell the phenology class with the 
maximum LAI for 1985-2014 (30-year mean; variable name: phen_max_lai_phen_number). A second 
netCDF output file is created containing additional LAI statistics for PFTs and phenology 
classes. 


----------REQUIREMENTS----------
These Fortran 90 programs have been compiled and executed on workstations running Microsoft Windows 7 Enterprise, using Intel Visual Studio Compiler Integration for Microsoft Visual 
Studio 2008 and NetCDF version 4.1.3 (Unidata 2011).  Note that these programs were run on a workstation with 192 GB of memory and may not run on systems with less memory.The netCDF software used by these programs for reading and writing netCDF files is available from: https://www.unidata.ucar.edu/software/netcdf/.


----------DISCLAIMER----------
USGS Disclaimer: This software is being provided to meet the need for timely best science.
No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that neither the
USGS nor the U.S. Government shall be held liable for any damages resulting from the
authorized or unauthorized use of the software. The USGS or the U.S. Government shall not
be held liable for improper or incorrect use of the data described and/or contained
herein. Any use of trade, firm, or product names is for descriptive purposes only and
does not imply endorsement by the U.S. Government.


----------REFERENCES----------
Pugh, T.A.M., Rademacher, T., Shafer, S.L., Steinkamp, J., Barichivich, J., Beckage, B., 
Haverd, V., Harper, A., Heinke, J., Nishina, K., Rammig, A., Sato, H., Arneth, A., Hantson, 
S., Hickler, T., Kautz, M., Quesada, B., Smith, B., Thonicke, K. Understanding 
the uncertainty in global forest carbon turnover. Biogeosciences, 2020. 


