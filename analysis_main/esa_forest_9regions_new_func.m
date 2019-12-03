function [esa_sel,regions,nregion]=esa_forest_9regions_new_func(makenetcdf)
%Function to read ESA CCI landcover data as aggregated to a 0.5 x 0.5 degree grid using esa_hires_region_mask.mat,
%adjust the grouping to 9 forest type classes and then (if required) write these to a netcdf file.
%
%T. Pugh
%10.11.17

%Load the data preprocessed by esa_hires_region_mask.m
load /data/ESA_landcover/esa_hires_region_mask.mat

%Split Broadleaved deciduous forest types by tropical/temperate (23 degrees split)
[lons,lats]=meshgrid(-179.75:0.5:179.75,-89.75:0.5:89.75);
esa_sel=NaN(360,720); 
%Assign other values to 9
esa_sel(mask_broad_ever & abs(lats)<=23)=1; 
esa_sel(mask_broad_dec & abs(lats)<=23)=2;
esa_sel(mask_broad_ever & abs(lats)>23)=4;
esa_sel(mask_broad_dec & abs(lats)>23)=5;
esa_sel(mask_needl_ever)=6;
esa_sel(mask_needl_dec)=7;
esa_sel(mask_tree_mixed)=8;
esa_sel(mask_tree_shrub_mosaic | mask_tree_flooded)=9;
esa_sel(esa_sel==9 & abs(lats)<=23)=3; %Assign other tropical values to 2
regions={'TrBE','TrBD','OTr','TeBE','TeBD','NE','ND','MX','Other'};

if (makenetcdf)
    %Make netcdf
    outfile='ESA_forest_9regions_v2.nc';
    nccreate(outfile,'region_mask','Dimensions',{'longitude',720,'latitude',360})
    nccreate(outfile,'latitude','Dimensions',{'latitude',360})
    nccreate(outfile,'longitude','Dimensions',{'longitude',720})
    ncwrite(outfile,'region_mask',esa_sel');
    ncwrite(outfile,'latitude',-89.75:0.5:89.75);
    ncwrite(outfile,'longitude',-179.75:0.5:179.75);
    ncwriteatt(outfile,'region_mask','Comment','1=Tropical broadleaved evergreen, 2=Tropical broadleaved deciduous, 3=Other Tropical (<23 degrees), 4=Temperate broadleaved evergreen, 5=Temperate broadleaved deciduous, 6=Needleleaved evergreen, 7=Needleleaved Deciduous, 8=Broadleaved/Needleaved mixed forest, 9=Other')
    ncwriteatt(outfile,'longitude','Units','degrees_east')
    ncwriteatt(outfile,'latitude','Units','degrees_north')
    ncwriteatt(outfile,'/','Comment','Created based on ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif using esa_forest_9regions_new_func.m')
    ncwriteatt(outfile,'/','Version','2')
    ncwriteatt(outfile,'/','Contact','T. Pugh, t.a.m.pugh@bham.ac.uk')
end
