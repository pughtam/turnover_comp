%Read C pool sizes and turnover fluxes from all models.
%The script saves the output in standardised *.mat files which are the
%basis for further analysis in other scripts.
%
%Note that this script requires annual inputs. Where monthly model data was
%provided, this first needs to be converted to annual data using the
%scripts in the 'data_conversion' folder.
%
%Forest masks used are developed using scripts in folder 'model_masks'
%
%Dependencies:
%-global_grid_area.m
%
%T. Pugh
%25.06.19

ifcruncep=false; %Read CRU-NCEP data (true) or IPSL data (false)

ifmask=true; %Apply forest masks for individual models

ifformask=true; %Apply closed-canopy forest mask (based on Hansen et al., 2013, Science) to totals

%---
%In this section the input information is set

%Set indices for each model corresponding to the 'models' array
models={'LPJ-GUESS','LPJmL','SEIB-DGVM','JULES','CABLE-POP','ORCHIDEE'};
nmod=length(models);
ilpjg=find(strcmp(models,'LPJ-GUESS'));
ilpjml=find(strcmp(models,'LPJmL'));
iseib=find(strcmp(models,'SEIB-DGVM'));
ijules=find(strcmp(models,'JULES'));
icable=find(strcmp(models,'CABLE-POP'));
iorchidee=find(strcmp(models,'ORCHIDEE'));

%Set locations and filenames for model simulation data
if ifcruncep
    nyear=114;
    nyear_cable=115; %CABLE-POP ran for one extra year
    climdata='cruncep';
    climdata_lpjml='cruncep';
    climdata_jules='CRUNCEP';
    climdata_cable='cruncep';
    climdata_orchidee='cruncep';
    lyear='2014';
    lyear_jules='1901_2014';
    lyear_cable='2015';
    lpjg_dir='/data/turnover/LPJG/CRUNCEP_v3/netcdfs/';
    lpj_dir='/data/turnover/LPJ/CRUNCEP_v2/';
    lpjml_dir='/data/turnover/LPJmL/CRUNCEP/';
    seib_dir='/data/turnover/SEIB/CRUNCEP/netcdf_v6/';
    jules_dir='/data/turnover/JULES/CRUNCEP_v2/';
    cable_dir='/data/turnover/CABLE/CRUNCEP_v2/';
    orchidee_dir='/data/turnover/ORCHIDEE/CRUNCEP_v2/';
else %IPSL-CM5A-LR bias-corrected RCP 8.5 climate
    nyear=199;
    nyear_cable=199;
    climdata='ipsl';
    climdata_lpjml='ipsl-cm5a_lr';
    climdata_jules='IPSL';
    climdata_cable='ipsl-cm5a-lr';
    climdata_orchidee='ipslcm5alr';
    lyear='2099';
    lyear_jules='1901_2099';
    lyear_cable='2099';
    lpjg_dir='/data/turnover/LPJG/IPSL_v3/netcdfs/';
    lpj_dir='/data/turnover/LPJ/IPSL_v2/';
    lpjml_dir='/data/turnover/LPJmL/IPSL/';
    seib_dir='/data/turnover/SEIB/IPSL/netcdf_v6/';
    jules_dir='/data/turnover/JULES/IPSL_v2/';
    cable_dir='/data/turnover/CABLE/IPSL_v2/';
    orchidee_dir='/data/turnover/ORCHIDEE/IPSL_v2/';
end

%Get the area of the grid cells
area_05=global_grid_area()';
%JULES is the only model not run on a 0.5 x 0.5 degree grid
area_jules=ncread('/data/turnover/JULES/JULES-LandMask.nc','arealand');
area_jules=area_jules(:,29:140); %Only need section for -55 to 83.75, i.e. 29:140

%Read in forest masks for each model (developed using scripts in folder 'model_masks')
mask_dir='/data/turnover/masks/v4/';
mask_lpjg=ncread([mask_dir,'/lpj-guess_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years');
mask_lpj=ncread([mask_dir,'/lpj-wsl_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years');
mask_lpjml=ncread([mask_dir,'/lpjml_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years');
mask_seib=ncread([mask_dir,'/seib_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years');
mask_jules=ncread([mask_dir,'/JULESC2_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years');

mask_cable=ncread([mask_dir,'/CABLE-POP_cruncep_lai_annual_1901_2015_forest_mask_v4.nc'],'forest_30yr_any_10_years');
mask_orchidee=ncread([mask_dir,'/ORCHIDEE_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years');

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
fmask_dir='/data/turnover/';
[~,~,fmask,fmask_jules_in]=get_closed_can_mask(fmask_dir);
%The JULES mask needs to be rearranged to match the JULES data which is based on a different meridian
fmask_jules=NaN(size(fmask_jules_in));
fmask_jules(1:96,:)=fmask_jules_in(97:192,:);
fmask_jules(97:192,:)=fmask_jules_in(1:96,:);

%Set some constants
sec_per_day=86400;
kg_per_Pg=1e12;
days_per_year=365;
days_per_year_jules=360; %Only 360 days per year in JULES

%---
%Now we start the processing

%If masking by closed-canopy forest area, then modify the grid-cell area array accordingly
if ifformask
    area_05=area_05.*fmask;
    area_jules=area_jules.*fmask_jules(:,29:140);
end

%Initialise global total arrays for variables of interest. Note that not every model has
%(or should have) output for each of these
cwood=NaN(nyear,nmod); %C mass in wood
cleaf=NaN(nyear,nmod); %C mass in leaves
croot=NaN(nyear,nmod); %C mass in fine roots
cveg=NaN(nyear,nmod); %C mass in vegetation overall
csoil=NaN(nyear,nmod); %C mass in soil and litter
tleaf=NaN(nyear,nmod); %Annual turnover flux through leaves
troot=NaN(nyear,nmod); %Annual turnover flux through fine roots
tfire=NaN(nyear,nmod); %Annual turnover flux through fire mortality
tgroweff=NaN(nyear,nmod); %Annual turnover flux through growth efficiency mortality
tdist=NaN(nyear,nmod); %Annual turnover flux through background disturbance mortality
tother=NaN(nyear,nmod); %Annual turnover flux through 'other' mortality (sum of all mortality mechanisms except for fire, growth efficiency and background disturbance)
tmort=NaN(nyear,nmod); %Annual turnover flux through all mortality causes
theat=NaN(nyear,nmod); %Annual turnover flux through heat stress mortality
tbioclim=NaN(nyear,nmod); %Annual turnover flux through bioclimatic limit mortality
tshade=NaN(nyear,nmod); %Annual turnover flux through shading mortality (i.e. prescribed self-thinning)
tback=NaN(nyear,nmod); %Annual turnover flux through background mortality
tnbio=NaN(nyear,nmod); %Annual turnover flux through negative biomass mortality
tallom=NaN(nyear,nmod); %Annual turnover flux through bad allometry mortality
allfire=NaN(nyear,nmod); %Annual total fire turnover flux (inc. litter)
npp=NaN(nyear,nmod); %Annual NPP
gpp=NaN(nyear,nmod); %Annual GPP

%-----
%LPJ-GUESS

cwood_lpjg=squeeze(ncread([lpjg_dir,'lpj-guess_',climdata,'_cwood_annual_1901_',lyear,'.nc4'],'cwood',[1 1 12 1],[Inf Inf 1 Inf]));
cleaf_lpjg=squeeze(ncread([lpjg_dir,'lpj-guess_',climdata,'_cleaf_annual_1901_',lyear,'.nc4'],'cleaf',[1 1 12 1],[Inf Inf 1 Inf]));
croot_lpjg=squeeze(ncread([lpjg_dir,'lpj-guess_',climdata,'_croot_annual_1901_',lyear,'.nc4'],'croot',[1 1 12 1],[Inf Inf 1 Inf]));
cveg_lpjg=squeeze(ncread([lpjg_dir,'lpj-guess_',climdata,'_cveg_annual_1901_',lyear,'.nc4'],'cveg',[1 1 12 1],[Inf Inf 1 Inf]));
csoil_lpjg=squeeze(ncread([lpjg_dir,'lpj-guess_',climdata,'_csoil_annual_1901_',lyear,'.nc4'],'csoil'))+...
    squeeze(ncread([lpjg_dir,'lpj-guess_',climdata,'_clitter_annual_1901_',lyear,'.nc4'],'clitter'));

tleaf_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_turnover_leaf_annual_1901_',lyear,'.nc4'],'cmort_turnover_leaf'),3));
troot_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_turnover_root_annual_1901_',lyear,'.nc4'],'cmort_turnover_root'),3));
tfire_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_fire_annual_1901_',lyear,'.nc4'],'cmort_fire'),3));
tgroweff_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_groweff_annual_1901_',lyear,'.nc4'],'cmort_groweff'),3));
tdist_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_dist_annual_1901_',lyear,'.nc4'],'cmort_dist'),3));
tother_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_badallom_annual_1901_',lyear,'.nc4'],'cmort_badallom'),3))+...
    squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_bioclim_annual_1901_',lyear,'.nc4'],'cmort_bioclim'),3))+...
    squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_mmin_annual_1901_',lyear,'.nc4'],'cmort_mmin'),3))+...
    squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_nbiom_annual_1901_',lyear,'.nc4'],'cmort_nbiom'),3));
tbioclim_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_bioclim_annual_1901_',lyear,'.nc4'],'cmort_bioclim'),3));
tback_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_mmin_annual_1901_',lyear,'.nc4'],'cmort_mmin'),3));
tnbio_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_nbiom_annual_1901_',lyear,'.nc4'],'cmort_nbiom'),3));
tallom_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_cmort_badallom_annual_1901_',lyear,'.nc4'],'cmort_badallom'),3));

npp_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_npp_annual_1901_',lyear,'.nc4'],'npp'),3));
gpp_lpjg=squeeze(sum(ncread([lpjg_dir,'lpj-guess_',climdata,'_gpp_annual_1901_',lyear,'.nc4'],'gpp'),3));

allfire_lpjg=squeeze(ncread([lpjg_dir,'lpj-guess_',climdata,'_ffire_annual_1901_',lyear,'.nc4'],'ffire'));

%Remove the last year of each array, which is empty
if ifcruncep
    npp_lpjg(:,:,115)=[];
    gpp_lpjg(:,:,115)=[];
else %IPSL
    npp_lpjg(:,:,200)=[];
    gpp_lpjg(:,:,200)=[];
end

%Apply the LPJ-GUESS-specific mask
if ifmask
    cwood_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    cleaf_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    croot_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    cveg_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    csoil_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tleaf_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    troot_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tfire_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tgroweff_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tdist_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tother_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tbioclim_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tback_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tnbio_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    tallom_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    allfire_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    npp_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
    gpp_lpjg(repmat(mask_lpjg,[1 1 nyear])~=1)=NaN;
end

%Make global sums
cwood(:,ilpjg)=nansum(nansum(cwood_lpjg.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cleaf(:,ilpjg)=nansum(nansum(cleaf_lpjg.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
croot(:,ilpjg)=nansum(nansum(croot_lpjg.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cveg(:,ilpjg)=nansum(nansum(cveg_lpjg.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
csoil(:,ilpjg)=nansum(nansum(csoil_lpjg.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C

tleaf(:,ilpjg)=nansum(nansum(tleaf_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
troot(:,ilpjg)=nansum(nansum(troot_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tfire(:,ilpjg)=nansum(nansum(tfire_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tgroweff(:,ilpjg)=nansum(nansum(tgroweff_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tdist(:,ilpjg)=nansum(nansum(tdist_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tother(:,ilpjg)=nansum(nansum(tother_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tbioclim(:,ilpjg)=nansum(nansum(tbioclim_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tback(:,ilpjg)=nansum(nansum(tback_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tnbio(:,ilpjg)=nansum(nansum(tnbio_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tallom(:,ilpjg)=nansum(nansum(tallom_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tmort(:,ilpjg)=tfire(:,ilpjg)+tgroweff(:,ilpjg)+tdist(:,ilpjg)+tother(:,ilpjg);

allfire(:,ilpjg)=nansum(nansum(allfire_lpjg.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C a-1

npp(:,ilpjg)=nansum(nansum(npp_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
gpp(:,ilpjg)=nansum(nansum(gpp_lpjg.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1

if ifcruncep
    save lpjg_map_arrays_cruncep.mat allfire_lpjg cwood_lpjg cleaf_lpjg croot_lpjg cveg_lpjg csoil_lpjg tleaf_lpjg troot_lpjg tfire_lpjg tgroweff_lpjg tdist_lpjg tother_lpjg npp_lpjg gpp_lpjg tbioclim_lpjg tback_lpjg tnbio_lpjg tallom_lpjg
else
    save lpjg_map_arrays_ipsl.mat allfire_lpjg cwood_lpjg cleaf_lpjg croot_lpjg cveg_lpjg csoil_lpjg tleaf_lpjg troot_lpjg tfire_lpjg tgroweff_lpjg tdist_lpjg tother_lpjg npp_lpjg gpp_lpjg tbioclim_lpjg tback_lpjg tnbio_lpjg tallom_lpjg
end

clear allfire_lpjg cwood_lpjg cleaf_lpjg croot_lpjg cveg_lpjg csoil_lpjg tleaf_lpjg troot_lpjg tfire_lpjg tgroweff_lpjg tdist_lpjg tother_lpjg npp_lpjg gpp_lpjg
clear tbioclim_lpjg tback_lpjg tnbio_lpjg tallom_lpjg

fprintf('Read LPJ-GUESS data\n')

%-----
%LPJmL

cwood_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cwood_annual_1901_',lyear,'.nc'],'cwood'),3));
cleaf_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cleaf_annual_1901_',lyear,'.nc'],'cleaf'),3));
croot_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_croot_annual_1901_',lyear,'.nc'],'croot'),3));
cveg_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cveg_annual_1901_',lyear,'.nc'],'cveg'),3));
csoil_lpjmla=squeeze(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_csoil_fast_annual_1901_',lyear,'.nc'],'csoil_fast'))+...
    squeeze(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_csoil_slow_annual_1901_',lyear,'.nc'],'csoil_slow'))+...
    squeeze(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_clitter_aboveground_leaf_annual_1901_',lyear,'.nc'],'clitter_aboveground_leaf'))+...
    squeeze(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_clitter_aboveground_wood_annual_1901_',lyear,'.nc'],'clitter_aboveground_wood'))+...
    squeeze(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_clitter_belowground_annual_1901_',lyear,'.nc'],'clitter_belowground'));

tleaf_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_turnover_leaf_annual_1901_',lyear,'.nc'],'cmort_turnover_leaf'),3));
troot_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_turnover_root_annual_1901_',lyear,'.nc'],'cmort_turnover_root'),3));
tfire_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_fire_annual_1901_',lyear,'.nc'],'cmort_fire'),3));
tgroweff_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_growth_efficiency_annual_1901_',lyear,'.nc'],'cmort_growth_efficiency'),3));
%no tdist for LPJmL
tother_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_neg_allocation_annual_1901_',lyear,'.nc'],'cmort_neg_allocation'),3))+...
    squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_neg_biomass_annual_1901_',lyear,'.nc'],'cmort_neg_biomass'),3))+...
    squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_heatstress_annual_1901_',lyear,'.nc'],'cmort_heatstress'),3))+...
    squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_shade_annual_1901_',lyear,'.nc'],'cmort_shade'),3));
tnbio_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_neg_biomass_annual_1901_',lyear,'.nc'],'cmort_neg_biomass'),3));
tallom_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_neg_allocation_annual_1901_',lyear,'.nc'],'cmort_neg_allocation'),3));
theat_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_heatstress_annual_1901_',lyear,'.nc'],'cmort_heatstress'),3));
tshade_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_cmort_shade_annual_1901_',lyear,'.nc'],'cmort_shade'),3));

allfire_lpjmla=squeeze(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_ffire_annual_1901_',lyear,'.nc'],'ffire'));

npp_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_npp_annual_1901_',lyear,'.nc'],'npp'),3));
gpp_lpjmla=squeeze(sum(ncread([lpjml_dir,'lpjml_',climdata_lpjml,'_gpp_annual_1901_',lyear,'.nc'],'gpp'),3));

%Apply the LPJmL-specific mask
if ifmask
    cwood_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    cleaf_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    croot_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    cveg_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    csoil_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    tleaf_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    troot_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    tfire_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    tgroweff_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    tother_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    tnbio_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    tallom_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    theat_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    tshade_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    allfire_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    npp_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
    gpp_lpjmla(repmat(mask_lpjml,[1 1 nyear])~=1)=NaN;
end

%Pad out arrays to -89.75:89.75 (note that in LPJmL outputs the latitude dimension goes from -55.75 to 83.25, i.e. 69:347)
cwood_lpjml=NaN(720,360,nyear); cwood_lpjml(:,69:347,:)=cwood_lpjmla;
cleaf_lpjml=NaN(720,360,nyear); cleaf_lpjml(:,69:347,:)=cleaf_lpjmla;
croot_lpjml=NaN(720,360,nyear); croot_lpjml(:,69:347,:)=croot_lpjmla;
cveg_lpjml=NaN(720,360,nyear); cveg_lpjml(:,69:347,:)=cveg_lpjmla;
csoil_lpjml=NaN(720,360,nyear); csoil_lpjml(:,69:347,:)=csoil_lpjmla;
tleaf_lpjml=NaN(720,360,nyear); tleaf_lpjml(:,69:347,:)=tleaf_lpjmla;
troot_lpjml=NaN(720,360,nyear); troot_lpjml(:,69:347,:)=troot_lpjmla;
tfire_lpjml=NaN(720,360,nyear); tfire_lpjml(:,69:347,:)=tfire_lpjmla;
tgroweff_lpjml=NaN(720,360,nyear); tgroweff_lpjml(:,69:347,:)=tgroweff_lpjmla;
tother_lpjml=NaN(720,360,nyear); tother_lpjml(:,69:347,:)=tother_lpjmla;
tnbio_lpjml=NaN(720,360,nyear); tnbio_lpjml(:,69:347,:)=tnbio_lpjmla;
tallom_lpjml=NaN(720,360,nyear); tallom_lpjml(:,69:347,:)=tallom_lpjmla;
theat_lpjml=NaN(720,360,nyear); theat_lpjml(:,69:347,:)=theat_lpjmla;
tshade_lpjml=NaN(720,360,nyear); tshade_lpjml(:,69:347,:)=tshade_lpjmla;
allfire_lpjml=NaN(720,360,nyear); allfire_lpjml(:,69:347,:)=allfire_lpjmla;
npp_lpjml=NaN(720,360,nyear); npp_lpjml(:,69:347,:)=npp_lpjmla;
gpp_lpjml=NaN(720,360,nyear); gpp_lpjml(:,69:347,:)=gpp_lpjmla;
clear cwood_lpjmla cleaf_lpjmla croot_lpjmla cveg_lpjmla csoil_lpjmla tleaf_lpjmla troot_lpjmla tfire_lpjmla tgroweff_lpjmla tother_lpjmla npp_lpjmla gpp_lpjmla
clear tnbio_lpjmla tallom_lpjmla theat_lpjmla tshade_lpjmla allfire_lpjmla

%Make global sums
cwood(:,ilpjml)=nansum(nansum(cwood_lpjml.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cleaf(:,ilpjml)=nansum(nansum(cleaf_lpjml.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
croot(:,ilpjml)=nansum(nansum(croot_lpjml.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cveg(:,ilpjml)=nansum(nansum(cveg_lpjml.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
csoil(:,ilpjml)=nansum(nansum(csoil_lpjml.*repmat(area_05,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C

tleaf(:,ilpjml)=nansum(nansum(tleaf_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
troot(:,ilpjml)=nansum(nansum(troot_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tfire(:,ilpjml)=nansum(nansum(tfire_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tgroweff(:,ilpjml)=nansum(nansum(tgroweff_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tother(:,ilpjml)=nansum(nansum(tother_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tnbio(:,ilpjml)=nansum(nansum(tnbio_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tallom(:,ilpjml)=nansum(nansum(tallom_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
theat(:,ilpjml)=nansum(nansum(theat_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tshade(:,ilpjml)=nansum(nansum(tshade_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tmort(:,ilpjml)=tfire(:,ilpjml)+tgroweff(:,ilpjml)+tother(:,ilpjml);

allfire(:,ilpjml)=nansum(nansum(allfire_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1

npp(:,ilpjml)=nansum(nansum(npp_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
gpp(:,ilpjml)=nansum(nansum(gpp_lpjml.*repmat(area_05,[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1

if ifcruncep
    save lpjml_map_arrays_cruncep.mat allfire_lpjml cwood_lpjml cleaf_lpjml croot_lpjml cveg_lpjml csoil_lpjml tleaf_lpjml troot_lpjml tfire_lpjml tgroweff_lpjml tother_lpjml npp_lpjml gpp_lpjml tnbio_lpjml tallom_lpjml theat_lpjml tshade_lpjml
else
    save lpjml_map_arrays_ipsl.mat allfire_lpjml cwood_lpjml cleaf_lpjml croot_lpjml cveg_lpjml csoil_lpjml tleaf_lpjml troot_lpjml tfire_lpjml tgroweff_lpjml tother_lpjml npp_lpjml gpp_lpjml tnbio_lpjml tallom_lpjml theat_lpjml tshade_lpjml
end

clear cwood_lpjml cleaf_lpjml croot_lpjml cveg_lpjml csoil_lpjml tleaf_lpjml troot_lpjml tfire_lpjml tgroweff_lpjml tother_lpjml npp_lpjml gpp_lpjml
clear tnbio_lpjml tallom_lpjml theat_lpjml tshade_lpjml

fprintf('Read LPJmL data\n')

%-----
%SEIB-DGVM

cwood_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cwood_annual_1901_',lyear,'.nc4'],'cwood'),4));
cleaf_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cleaf_annual_1901_',lyear,'.nc4'],'cleaf'),4));
croot_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_croot_annual_1901_',lyear,'.nc4'],'croot'),4));
cveg_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cveg_annual_1901_',lyear,'.nc4'],'cveg'),4));
csoil_seib=squeeze(ncread([seib_dir,'seib_',climdata,'_csoil_annual_1901_',lyear,'.nc4'],'csoil'))+...
    squeeze(ncread([seib_dir,'seib_',climdata,'_clitter_annual_1901_',lyear,'.nc4'],'clitter'));

tleaf_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_turnover_leaf_annual_1901_',lyear,'.nc4'],'turnover_leaf'),4));
troot_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_turnover_root_annual_1901_',lyear,'.nc4'],'turnover_root'),4));
tfire_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_fire_atmosp_annual_1901_',lyear,'.nc4'],'cmort_fire_atmosp'),4))+...
    squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_fire_litter_annual_1901_',lyear,'.nc4'],'cmort_fire_litter'),4));
tgroweff_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_greff_annual_1901_',lyear,'.nc4'],'cmort_greff'),4));
tdist_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_gap_annual_1901_',lyear,'.nc4'],'cmort_gap'),4));
tother_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_heat_annual_1901_',lyear,'.nc4'],'cmort_heat'),4))+...
    squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_lim_annual_1901_',lyear,'.nc4'],'cmort_lim'),4))+...
    squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_etc_annual_1901_',lyear,'.nc4'],'cmort_etc'),4));
theat_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_heat_annual_1901_',lyear,'.nc4'],'cmort_heat'),4));
tbioclim_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_lim_annual_1901_',lyear,'.nc4'],'cmort_lim'),4));
tback_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_cmort_etc_annual_1901_',lyear,'.nc4'],'cmort_etc'),4));

allfire_seib=squeeze(ncread([seib_dir,'seib_',climdata,'_ffire_annual_1901_',lyear,'.nc4'],'ffire'));

npp_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_npp_annual_1901_',lyear,'.nc4'],'npp'),4));
gpp_seib=squeeze(sum(ncread([seib_dir,'seib_',climdata,'_gpp_annual_1901_',lyear,'.nc4'],'gpp'),4));

%Apply the SEIB-DGVM-specific mask
if ifmask
    cwood_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    cleaf_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    croot_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    cveg_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    csoil_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    tleaf_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    troot_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    tfire_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    tgroweff_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    tdist_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    tother_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    theat_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    tbioclim_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    tback_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    allfire_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    npp_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
    gpp_seib(repmat(mask_seib,[1 1 nyear])~=1)=NaN;
end

%Make global sums
cwood(:,iseib)=nansum(nansum(cwood_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cleaf(:,iseib)=nansum(nansum(cleaf_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
croot(:,iseib)=nansum(nansum(croot_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cveg(:,iseib)=nansum(nansum(cveg_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
csoil(:,iseib)=nansum(nansum(csoil_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C

tleaf(:,iseib)=nansum(nansum(tleaf_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
troot(:,iseib)=nansum(nansum(troot_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tfire(:,iseib)=nansum(nansum(tfire_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tgroweff(:,iseib)=nansum(nansum(tgroweff_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tdist(:,iseib)=nansum(nansum(tdist_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tother(:,iseib)=nansum(nansum(tother_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
theat(:,iseib)=nansum(nansum(theat_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tbioclim(:,iseib)=nansum(nansum(tbioclim_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tback(:,iseib)=nansum(nansum(tback_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tmort(:,iseib)=tfire(:,iseib)+tgroweff(:,iseib)+tdist(:,iseib)+tother(:,iseib);

allfire(:,iseib)=nansum(nansum(allfire_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1

npp(:,iseib)=nansum(nansum(npp_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
gpp(:,iseib)=nansum(nansum(gpp_seib.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1

if ifcruncep
    save seib_map_arrays_cruncep.mat allfire_seib cwood_seib cleaf_seib croot_seib cveg_seib csoil_seib tleaf_seib troot_seib tfire_seib tgroweff_seib tdist_seib tother_seib npp_seib gpp_seib theat_seib tbioclim_seib tback_seib
else
    save seib_map_arrays_ipsl.mat allfire_seib cwood_seib cleaf_seib croot_seib cveg_seib csoil_seib tleaf_seib troot_seib tfire_seib tgroweff_seib tdist_seib tother_seib npp_seib gpp_seib theat_seib tbioclim_seib tback_seib
end

clear cwood_seib cleaf_seib croot_seib cveg_seib csoil_seib tleaf_seib troot_seib tfire_seib tgroweff_seib tdist_seib tother_seib npp_seib gpp_seib
clear theat_seib tbioclim_seib tback_seib

fprintf('Read SEIB-DGVM data\n')

%-----
%JULES

cwood_jules=squeeze(ncread([jules_dir,'JULESC2_',climdata_jules,'_cwood_gb_Annual_',lyear_jules,'.nc'],'cwood_gb'));
cleaf_jules=squeeze(ncread([jules_dir,'JULESC2_',climdata_jules,'_cleaf_gb_Annual_',lyear_jules,'.nc'],'cleaf_gb'));
croot_jules=squeeze(ncread([jules_dir,'JULESC2_',climdata_jules,'_croot_gb_Annual_',lyear_jules,'.nc'],'croot_gb'));
cveg_jules=squeeze(ncread([jules_dir,'JULESC2_',climdata_jules,'_cveg_gb_Annual_',lyear_jules,'.nc'],'cveg_gb'));
csoil_jules=squeeze(ncread([jules_dir,'JULESC2_',climdata_jules,'_csoil_Annual_',lyear_jules,'.nc'],'csoil'));

tleaf_jules=squeeze(sum(ncread([jules_dir,'JULESC2_',climdata_jules,'_leaf_litC_Annual_',lyear_jules,'.nc'],'leaf_litC',[1 1 1 1],[Inf Inf Inf nyear]),3));
troot_jules=squeeze(sum(ncread([jules_dir,'JULESC2_',climdata_jules,'_root_litC_Annual_',lyear_jules,'.nc'],'root_litC',[1 1 1 1],[Inf Inf Inf nyear]),3));

tback_jules=squeeze(sum(ncread([jules_dir,'JULESC2_',climdata_jules,'_wood_litC_Annual_',lyear_jules,'.nc'],'wood_litC',[1 1 1 1],[Inf Inf Inf nyear]),3));
tother_jules=squeeze(sum(ncread([jules_dir,'JULESC2_',climdata_jules,'_lit_c_dyn_wood_Annual_',lyear_jules,'.nc'],'lit_c_dyn_wood',[1 1 1 1],[Inf Inf Inf nyear]),3))+...
    squeeze(sum(ncread([jules_dir,'JULESC2_',climdata_jules,'_lit_c_dyn_leaf_Annual_',lyear_jules,'.nc'],'lit_c_dyn_leaf',[1 1 1 1],[Inf Inf Inf nyear]),3))+...
    squeeze(sum(ncread([jules_dir,'JULESC2_',climdata_jules,'_lit_c_dyn_root_Annual_',lyear_jules,'.nc'],'lit_c_dyn_root',[1 1 1 1],[Inf Inf Inf nyear]),3));

%no tfire for JULES
%no tgroweff for JULES

npp_jules=squeeze(ncread([jules_dir,'JULESC2_',climdata_jules,'_npp_gb_Annual_',lyear_jules,'.nc'],'npp_gb'));
gpp_jules=squeeze(ncread([jules_dir,'JULESC2_',climdata_jules,'_gpp_gb_Annual_',lyear_jules,'.nc'],'gpp_gb'));

%Apply the JULES-specific mask
if ifmask
    cwood_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    cleaf_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    croot_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    cveg_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    csoil_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    tleaf_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    troot_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    tback_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    tother_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    npp_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
    gpp_jules(repmat(mask_jules,[1 1 nyear])~=1)=NaN;
end

%Make global sums
cwood(:,ijules)=nansum(nansum(cwood_jules.*repmat(area_jules,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cleaf(:,ijules)=nansum(nansum(cleaf_jules.*repmat(area_jules,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
croot(:,ijules)=nansum(nansum(croot_jules.*repmat(area_jules,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cveg(:,ijules)=nansum(nansum(cveg_jules.*repmat(area_jules,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
csoil(:,ijules)=nansum(nansum(csoil_jules.*repmat(area_jules,[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C

tleaf(:,ijules)=nansum(nansum(tleaf_jules.*repmat(area_jules,[1 1 nyear]),2),1)*sec_per_day*days_per_year_jules/kg_per_Pg; %Give number in Pg C a-1
troot(:,ijules)=nansum(nansum(troot_jules.*repmat(area_jules,[1 1 nyear]),2),1)*sec_per_day*days_per_year_jules/kg_per_Pg; %Give number in Pg C a-1
tback(:,ijules)=nansum(nansum(tback_jules.*repmat(area_jules,[1 1 nyear]),2),1)*sec_per_day*days_per_year_jules/kg_per_Pg; %Give number in Pg C a-1
tother(:,ijules)=nansum(nansum(tother_jules.*repmat(area_jules,[1 1 nyear]),2),1)*sec_per_day*days_per_year_jules/kg_per_Pg; %Give number in Pg C a-1
tmort(:,ijules)=tback(:,ijules)+tother(:,ijules);

npp(:,ijules)=nansum(nansum(npp_jules.*repmat(area_jules,[1 1 nyear]),2),1)*sec_per_day*days_per_year_jules/kg_per_Pg; %Give number in Pg C a-1
gpp(:,ijules)=nansum(nansum(gpp_jules.*repmat(area_jules,[1 1 nyear]),2),1)*sec_per_day*days_per_year_jules/kg_per_Pg; %Give number in Pg C a-1

if ifcruncep
    save jules_map_arrays_cruncep.mat cwood_jules cleaf_jules croot_jules cveg_jules csoil_jules tback_jules tother_jules npp_jules gpp_jules tleaf_jules troot_jules
else
    save jules_map_arrays_ipsl.mat cwood_jules cleaf_jules croot_jules cveg_jules csoil_jules tback_jules tother_jules npp_jules gpp_jules tleaf_jules troot_jules
end

clear cwood_jules cleaf_jules croot_jules cveg_jules csoil_jules tleaf_jules troot_jules tback_jules tother_jules npp_jules gpp_jules

fprintf('Read JULES data\n')

%-----
%CABLE-POP

cwood_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cwood_year_1901_',lyear_cable,'.nc4'],'cwood'),3));
cleaf_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cleaf_year_1901_',lyear_cable,'.nc4'],'cleaf'),3));
croot_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_croot_year_1901_',lyear_cable,'.nc4'],'croot'),3));
cveg_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cveg_year_1901_',lyear_cable,'.nc4'],'cveg'),3));
csoil_cable=squeeze(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_csoil_year_1901_',lyear_cable,'.nc4'],'csoil'))+...
    squeeze(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_clitter_year_1901_',lyear_cable,'.nc4'],'clitter'));

tleaf_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cmortleaf_year_1901_',lyear_cable,'.nc4'],'cmortleaf'),3));
troot_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cmortfineroot_year_1901_',lyear_cable,'.nc4'],'cmortfineroot'),3));
%No tfire for CABLE-POP
tgroweff_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cmortwoodresourcelim_year_1901_',lyear_cable,'.nc4'],'cmortwoodresourcelim'),3));
tdist_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cmortwooddist_year_1901_',lyear_cable,'.nc4'],'cmortwooddist'),3));
tother_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cmortwoodresourcecrowding_year_1901_',lyear_cable,'.nc4'],'cmortwoodresourcecrowding'),3));
%tmortall_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_cmort_annual_1901_',lyear_cable,'.nc4'],'cmort'),4));

npp_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_npp_year_1901_',lyear_cable,'.nc4'],'npp'),3));
gpp_cable=squeeze(sum(ncread([cable_dir,'CABLE-POP_',climdata_cable,'_gpp_year_1901_',lyear_cable,'.nc4'],'gpp'),3));

%Apply the CABLE-POP-specific mask
if ifmask
    cwood_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    cleaf_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    croot_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    cveg_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    csoil_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    tleaf_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    troot_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    tfire_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    tgroweff_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    tdist_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    tother_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    npp_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
    gpp_cable(repmat(mask_cable,[1 1 nyear_cable])~=1)=NaN;
end

%Make global sums
cwood(:,icable)=nansum(nansum(cwood_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cleaf(:,icable)=nansum(nansum(cleaf_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
croot(:,icable)=nansum(nansum(croot_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cveg(:,icable)=nansum(nansum(cveg_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
csoil(:,icable)=nansum(nansum(csoil_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C

tleaf(:,icable)=nansum(nansum(tleaf_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
troot(:,icable)=nansum(nansum(troot_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tgroweff(:,icable)=nansum(nansum(tgroweff_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tdist(:,icable)=nansum(nansum(tdist_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tother(:,icable)=nansum(nansum(tother_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
tmort(:,icable)=tgroweff(:,icable)+tdist(:,icable)+tother(:,icable);

npp(:,icable)=nansum(nansum(npp_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
gpp(:,icable)=nansum(nansum(gpp_cable(:,:,1:nyear).*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1

if ifcruncep
    save cable_map_arrays_cruncep.mat cwood_cable cleaf_cable croot_cable cveg_cable csoil_cable tleaf_cable troot_cable tgroweff_cable tdist_cable tother_cable npp_cable gpp_cable
else
    save cable_map_arrays_ipsl.mat cwood_cable cleaf_cable croot_cable cveg_cable csoil_cable tleaf_cable troot_cable tgroweff_cable tdist_cable tother_cable npp_cable gpp_cable
end

clear cwood_cable cleaf_cable croot_cable cveg_cable csoil_cable tleaf_cable troot_cable tgroweff_cable tdist_cable tother_cable npp_cable gpp_cable

fprintf('Read CABLE-POP data\n')


%-----
%ORCHIDEE

vegfrac=ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_veget_max_13pft_year_1901_',lyear,'.nc'],'veget_max');

cwood_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_cWood_13pft_year_1901_',lyear,'.nc'],'cWood').*vegfrac,3));
cleaf_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_cLeaf_13pft_year_1901_',lyear,'.nc'],'cLeaf').*vegfrac,3));
croot_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_cRootf_13pft_year_1901_',lyear,'.nc'],'ROOT_M').*vegfrac,3));
cveg_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_cVeg_13pft_year_1901_',lyear,'.nc'],'cVeg').*vegfrac,3));
csoil_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_cSoil_13pft_year_1901_',lyear,'.nc'],'cSoil').*vegfrac,3))+...
    squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_cLitter_13pft_year_1901_',lyear,'.nc'],'cLitter').*vegfrac,3));

tleaf_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_leaf_turn_13pft_year_1901_',lyear,'.nc'],'leaf_turn').*vegfrac,3));
troot_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_root_turn_13pft_year_1901_',lyear,'.nc'],'root_turn').*vegfrac,3));
tmort_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_total_turn_13pft_year_1901_',lyear,'.nc'],'total_turn').*vegfrac,3));
%No breakdown of mortality fluxes for ORCHIDEE

npp_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_npp_13pft_year_1901_',lyear,'.nc'],'npp').*vegfrac,3));
gpp_orchidee=squeeze(sum(ncread([orchidee_dir,'ORCHIDEEr3085_',climdata_orchidee,'_gpp_13pft_year_1901_',lyear,'.nc'],'gpp').*vegfrac,3));
clear vegfrac

%Remove no-data values
maxval=1e6;
cwood_orchidee(cwood_orchidee>maxval)=NaN;
cleaf_orchidee(cleaf_orchidee>maxval)=NaN;
croot_orchidee(croot_orchidee>maxval)=NaN;
cveg_orchidee(cveg_orchidee>maxval)=NaN;
csoil_orchidee(csoil_orchidee>maxval)=NaN;
tleaf_orchidee(tleaf_orchidee>maxval)=NaN;
troot_orchidee(troot_orchidee>maxval)=NaN;
tmort_orchidee(tmort_orchidee>maxval)=NaN;
npp_orchidee(npp_orchidee>maxval)=NaN;
gpp_orchidee(gpp_orchidee>maxval)=NaN;

%Apply masks
if ifmask
    cwood_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    cleaf_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    croot_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    cveg_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    csoil_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    tleaf_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    troot_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    tmort_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    npp_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
    gpp_orchidee(repmat(mask_orchidee,[1 1 nyear])~=1)=NaN;
end

%Make global sums
cwood(:,iorchidee)=nansum(nansum(cwood_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cleaf(:,iorchidee)=nansum(nansum(cleaf_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
croot(:,iorchidee)=nansum(nansum(croot_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
cveg(:,iorchidee)=nansum(nansum(cveg_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C
csoil(:,iorchidee)=nansum(nansum(csoil_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)/kg_per_Pg; %Give number in Pg C

tleaf(:,iorchidee)=nansum(nansum(tleaf_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*days_per_year/1e15; %Give number in Pg C a-1
troot(:,iorchidee)=nansum(nansum(troot_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*days_per_year/1e15; %Give number in Pg C a-1
tmort(:,iorchidee)=nansum(nansum(tmort_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*days_per_year/1e15; %Give number in Pg C a-1

npp(:,iorchidee)=nansum(nansum(npp_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1
gpp(:,iorchidee)=nansum(nansum(gpp_orchidee.*repmat(fliplr(area_05),[1 1 nyear]),2),1)*sec_per_day*days_per_year/kg_per_Pg; %Give number in Pg C a-1

if ifcruncep
    save orchidee_map_arrays_cruncep.mat cwood_orchidee cleaf_orchidee croot_orchidee cveg_orchidee csoil_orchidee tleaf_orchidee troot_orchidee tmort_orchidee npp_orchidee gpp_orchidee
else
    save orchidee_map_arrays_ipsl.mat cwood_orchidee cleaf_orchidee croot_orchidee cveg_orchidee csoil_orchidee tleaf_orchidee troot_orchidee tmort_orchidee npp_orchidee gpp_orchidee
end

clear cwood_orchidee cleaf_orchidee croot_orchidee cveg_orchidee csoil_orchidee tleaf_orchidee troot_orchidee tmort_orchidee npp_orchidee gpp_orchidee

fprintf('Read ORCHIDEE data\n')

%---
%Save summary arrays
if ifcruncep
    save all_cruncep_global_arrays.mat cwood cleaf croot cveg csoil tleaf troot tfire tgroweff tdist tother theat tbioclim tback tmort tshade tnbio tallom npp gpp ifcruncep ifmask ifformask
else
    save all_ipsl_global_arrays.mat cwood cleaf croot cveg csoil tleaf troot tfire tgroweff tdist tother theat tbioclim tback tmort tshade tnbio tallom npp gpp ifcruncep ifmask ifformask
end
