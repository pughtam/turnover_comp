function [phen,phen_jules,phen_label,nphen]=get_forest_type(phen_dir)
%Read in forest type masks
%Forest type files calculated using scripts in 'model_masks' directory as created by Sarah Shafer
%
%T. Pugh
%27.11.19

nmod=6; %Number of models

phen=NaN(360,720,nmod);
phen(:,:,1)=flipud(ncread([phen_dir,'/CABLE-POP_cruncep_lai_annual_1901_2015_phenology_mask.nc'],'phen_max_lai_phen_number')');
phen_julesa=ncread([phen_dir,'/JULESC2_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')';
phen(:,:,3)=ncread([phen_dir,'/lpj-guess_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')';
phen(69:347,:,4)=ncread([phen_dir,'/lpjml_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')';
phen(:,:,5)=flipud(ncread([phen_dir,'/ORCHIDEE_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')');
phen(:,:,6)=flipud(ncread([phen_dir,'/seib_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')');
phen_jules=NaN(192,144);
phen_jules(1:96,29:140)=phen_julesa(:,97:192)';
phen_jules(97:192,29:140)=phen_julesa(:,1:96)';
clear phen_julesa

phen_label={'needleleaved evergreen','needleleaved deciduous','boreal broadleaved deciduous',...
    'temperate broadleaved evergreen','temperate broadleaved deciduous','tropical broadleaved evergreen',...
    'tropical broadleaved raingreen'};
nphen=length(phen_label);