%Script to compare forest cover data from the models and the closed-canopy forest mask based on
%Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
%
%Forest cover is read from mask files calculated using scripts in 'model_masks' directory 
%(provided by Sarah Shafer)
%
%T. Pugh
%08.03.18

fthres=0.1; %Fractional forest cover threshold for comparison with observations

%Location of netcdf files containing preprocessed forest mask data for the models
mask_dir='/Users/pughtam/data/turnover/masks/v4/';
%Location of *.mat files containing preprocessed ESA landcover data
data_esa='/Users/pughtam/data/ESA_landcover/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/Users/pughtam/data/turnover/';

%---
models={'(a) CABLE-POP','(b) JULES','(c) LPJ-GUESS','(d) LPJmL','(e) ORCHIDEE','(f) SEIB-DGVM',};
nmod=length(models);

%Read in forest masks for each model
mmask=NaN(720,360,nmod);
mmask(:,:,1)=ncread([mask_dir,'/CABLE-POP_cruncep_lai_annual_1901_2015_forest_mask_v4.nc'],'forest_30yr_any_10_years');
mmask(:,:,3)=fliplr(ncread([mask_dir,'/lpj-guess_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years'));
mmask(:,14:292,4)=fliplr(ncread([mask_dir,'/lpjml_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years'));
mmask(:,:,5)=ncread([mask_dir,'/ORCHIDEE_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years');
mmask(:,:,6)=ncread([mask_dir,'/seib_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years');
mmask=flip(mmask,2);
mmask_julesa=ncread([mask_dir,'/JULESC2_cruncep_lai_annual_1901_2014_forest_mask_v4.nc'],'forest_30yr_any_10_years')';
mmask_jules=NaN(192,144);
%JULES data is on a diffeent meridian and needs to be converted
mmask_jules(1:96,29:140)=mmask_julesa(:,97:192)';
mmask_jules(97:192,29:140)=mmask_julesa(:,1:96)';
%mmask_jules=fliplr(mmask_jules);
clear mmask_julesa

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[~,~,ffrac,ffrac_jules]=get_closed_can_mask(data_mask);

%Create array for plotting
plotarray=NaN(720,360,nmod);
plotarray_jules=NaN(192,144);
for nn=1:nmod
    if nn==2 %JULES is on different grid
        plotarray_jules(mmask_jules==1 & ffrac_jules>fthres)=1; %Forest simulated and in closed-canopy observations
        plotarray_jules(mmask_jules==2 & ffrac_jules>fthres)=2; %Forest not simulated, but in closed-canopy observations
        plotarray_jules(mmask_jules==1 & (isnan(ffrac_jules) | ffrac_jules<fthres))=3; %Forest simulated, but no forest in closed-canopy observations
        continue
    end
    mmask_nn=squeeze(mmask(:,:,nn));
    plotarray_nn=squeeze(plotarray(:,:,nn));
    plotarray_nn(mmask_nn==1 & ffrac>fthres)=1; %Forest simulated and in closed-canopy observations
    plotarray_nn(mmask_nn==2 & ffrac>fthres)=2; %Forest not simulated, but in closed-canopy observations
    plotarray_nn(mmask_nn==1 & (isnan(ffrac) | ffrac<fthres))=3; %Forest simulated, but no forest in closed-canopy observations
    plotarray(:,:,nn)=plotarray_nn;
    clear mmask_nn plotarray_nn
end

%Load water mask data (from ESA landcover)
load([data_esa,'esa_05_landcover']); %Output from esa_lu_read.m
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
load([data_esa,'esa_jules_landcover']); %Output from esa_lu_read_julesgrid.m
oceanm_jules=NaN(192,144);
oceanm_jules(esa_jules'>200 & esa_jules'<220)=-1;

%Lat and lon arrays for plotting
lats=-89.75:0.5:89.75;
lons=-179.75:0.5:179.75;
latincj=180/144;
lats_jules=-90+(latincj/2):latincj:90-(latincj/2);
lons_jules=-179.0625:1.875:179.0625;

%Set colorbar limits
cmin=-1;
cmax=3;

figure
cmask=[0.9 0.9 0.9; 1 1 1; 0.2 0.2 0.7; 0.7 0 0; 0.5 0.5 0.5]; %Define a custom colormap
colormap(cmask)
for nn=1:nmod
    if nn==2; continue; end %JULES is a special case because of resolution
    s(nn)=subplot(3,2,nn);
    hold on
    l(nn)=pcolor(lons,lats,oceanm');
    set(l(nn),'linestyle','none')
    p(nn)=pcolor(lons,lats,squeeze(plotarray(:,:,nn)'));
    set(p(nn),'linestyle','none')
    box on
    set(gca,'XTick',[],'XTickLabel','')
    set(gca,'YTick',[],'YTickLabel','')
    set(gca,'YLim',[-60 80])
    set(gca,'XLim',[-180 180])
    caxis([cmin-(cmax/100) cmax])
    t(nn)=text(-175,-50,models{nn});
    set(t(nn),'FontSize',12,'FontWeight','Bold')
end
%JULES is a special case because of resolution
s(2)=subplot(3,2,2);
hold on
l(2)=pcolor(lons_jules,lats_jules,oceanm_jules');
set(l(2),'linestyle','none')
p(2)=pcolor(lons_jules,lats_jules,plotarray_jules');
set(p(2),'linestyle','none')
set(gca,'XTick',[],'XTickLabel','')
set(gca,'YTick',[],'YTickLabel','')
set(gca,'YLim',[-60 80])
set(gca,'XLim',[-180 180])
box on
caxis([cmin-(cmax/100) cmax])
t(2)=text(-175,-50,models{2});
set(t(2),'FontSize',12,'FontWeight','Bold')

set(s(1),'Position',[0.05 0.66 0.43 0.28])
set(s(2),'Position',[0.51 0.66 0.43 0.28])
set(s(3),'Position',[0.05 0.37 0.43 0.28])
set(s(4),'Position',[0.51 0.37 0.43 0.28])
set(s(5),'Position',[0.05 0.08 0.43 0.28])
set(s(6),'Position',[0.51 0.08 0.43 0.28])
