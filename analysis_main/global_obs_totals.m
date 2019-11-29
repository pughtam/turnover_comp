%Calculate global totals for NPP, Cveg and tau_NPP.
%NPP is from Zhao & Running (2010, Science 329, 940-943) for the period 2000 to 2012
%Cveg is from Carvalhais et al. (2014, Nature 514, 213-217)
%(Both pre-processed by T. Rademacher)
%
%T. Pugh
%25.11.19

load obs_NPP_0perc.mat
load obs_Cveg_0perc.mat

%Get the forest area
garea=global_grid_area();
ffrac=ncread('/Users/pughtam/Documents/GAP_and_other_work/Mortality/hansen_forested_frac_05.nc4','forested_50_percent');
ffrac=(double(ffrac)/100);
%Make a mask to exclude all grid cells with less than 10% closed-canopy
%forest cover to avoid results being biased by cells with almost no
%closed-canopy forest (although effect is very minimal)
fmask=NaN(size(ffrac));
fmask(ffrac>=0.1)=1;

farea=ffrac.*garea'.*fmask;

%Calculate global totals for NPP and Cveg
Cveg_farea=Cveg.*garea'; %Multiply by total grid-cell area, as the Cveg data is a biomass density per grid cell and we are making assumption that biomass outside closed-canopy forests is negligible
NPP_farea=NPP.*farea; %Multiply by closed-canopy forest area only. We are assuming that NPP is uniform across the grid cell and we only want the component relating to closed-canopy forests

Cveg_farea_tot=nansum(Cveg_farea(:))/1e12;
NPP_farea_tot=nansum(NPP_farea(:))/1e12;

%Calculate a Cveg value assuming that the whole grid cell is closed canopy
%forest. Implicity assumes that biomass outside closed-canopy forests is
%negligible
Cveg_pot=Cveg./ffrac;
tau_NPP=Cveg_pot./NPP;

%Calculate a global turnover time mean weighted by closed-canopy forest
%area in each grid cell.
tau_NPP_farea_mean=nansum(tau_NPP(:).*farea(:))/nansum(farea(:));
