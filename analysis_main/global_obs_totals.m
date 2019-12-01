%Calculate global totals for NPP, Cveg and tau_NPP.
%NPP is from Zhao & Running (2010, Science 329, 940-943) for the period 2000 to 2012
%Cveg is from Carvalhais et al. (2014, Nature 514, 213-217)
%(Both pre-processed by T. Rademacher)
%
%T. Pugh
%25.11.19

%Location of *.mat files containing preprocessed observational data
data_obs='/Users/pughtam/Documents/GAP_and_other_work/Mortality/Tim_Dec_plots/raw_data/raw_data_v2/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/Users/pughtam/data/turnover/';

%---
load([data_obs,'obs_NPP_0perc.mat']);
load([data_obs,'obs_Cveg_0perc.mat']);

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[fmask,~,ffrac,~]=get_closed_can_mask(data_mask);

garea=global_grid_area();

farea=ffrac.*garea'.*fmask;

%Calculate global totals for NPP and Cveg
Cveg_farea=Cveg.*garea'; %Multiply by total grid-cell area, as the Cveg data is a biomass density per grid cell and we are making assumption that biomass outside closed-canopy forests is negligible
NPP_farea=NPP.*farea; %Multiply by closed-canopy forest area only. We are assuming that NPP is uniform across the grid cell and we only want the component relating to closed-canopy forests

Cveg_farea_tot=nansum(Cveg_farea(:))/1e12; %Pg C
NPP_farea_tot=nansum(NPP_farea(:))/1e12; %Pg C y-1
tau_NPP=Cveg_farea_tot./NPP_farea_tot; %years

fprintf('Cveg %7.1f\n',Cveg_farea_tot)
fprintf('NPP %7.1f\n',NPP_farea_tot)
fprintf('tau_NPP %7.1f\n',tau_NPP)