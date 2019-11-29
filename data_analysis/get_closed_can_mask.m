function [fmask,fmask_jules,ffrac,ffrac_jules]=get_closed_can_mask(data_mask)
%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
%A version of the mask was aggregated to the JULES model grid using CDO
%
%T. Pugh
%27.11.19

ffrac_in=ncread([data_mask,'hansen_forested_frac_05.nc4'],'forested_50_percent'); %Percentage of gridcell coverage with at least 50% forest.
%Area fraction
ffrac=(double(ffrac_in)/100); %Convert to fraction
ffrac(ffrac==0)=NaN;
%Mask for grid cells with >10% closed-canopy forest cover
fmask=NaN(size(ffrac_in));
fmask(ffrac_in>10)=1;
clear ffrac_in

ffrac_jules_in=fliplr(ncread([data_mask,'hansen_forested_frac_julesgrid.nc4'],'forested_50_percent')); %Percentage of gridcell coverage with at least 50% forest.
%Area fraction
ffrac_jules=(double(ffrac_jules_in)/100); %Convert to fraction
ffrac_jules(ffrac_jules==0)=NaN;
ffrac_jules(:,1)=[];
%Mask for grid cells with >10% closed-canopy forest cover
fmask_jules=NaN(size(ffrac_jules_in));
fmask_jules(ffrac_jules_in>10)=1;
fmask_jules(:,1)=[];
clear ffrac_jules_in