%Script to make maps of mean total LAI for a specified time period
%LAI is read from forest-type files calculated using scripts in 'model_masks' directory (provided by Sarah
%Shafer)
%
%Dependencies:
%-get_closed_can_mask.m
%
%T. Pugh
%27.11.19

lai_dir='/data/turnover/masks/phen/multiple_variables_files/';

models={'(a) CABLE-POP','(b) JULES','(c) LPJ-GUESS','(d) LPJmL','(e) ORCHIDEE','(f) SEIB-DGVM'};
nmod=length(models);

%Index of last year to average over in the data (LAI data in these files is the average of the previous 30
%years)
y1=114; %2014
startyear=30; %30 year offset in file begin

y1s=y1-startyear+1;

lai=NaN(360,720,nmod);

%Read the LAI data (note some rotation of arrays to get all on same dimensions)
lai(:,:,1)=flipud(ncread([lai_dir,'/CABLE-POP_cruncep_lai_annual_1901_2015_phenology.nc'],'lai_total',[1 1 y1s],[Inf Inf 1])');
lai_julesa=ncread([lai_dir,'/JULESC2_cruncep_lai_annual_1901_2014_phenology.nc'],'lai_total',[1 1 y1s],[Inf Inf 1])';
lai(:,:,3)=ncread([lai_dir,'/lpj-guess_cruncep_lai_annual_1901_2014_phenology.nc'],'lai_total',[1 1 y1s],[Inf Inf 1])';
lai(69:347,:,4)=ncread([lai_dir,'/lpjml_cruncep_lai_annual_1901_2014_phenology.nc'],'lai_total',[1 1 y1s],[Inf Inf 1])';
lai(:,:,5)=flipud(ncread([lai_dir,'/ORCHIDEE_cruncep_lai_annual_1901_2014_phenology.nc'],'lai_total',[1 1 y1s],[Inf Inf 1])');
lai(:,:,6)=flipud(ncread([lai_dir,'/seib_cruncep_lai_annual_1901_2014_phenology.nc'],'lai_total',[1 1 y1s],[Inf Inf 1])');

%JULES data needs to be reordered because of a different meridian
lai_jules=NaN(192,144);
lai_jules(1:96,29:140)=lai_julesa(:,97:192)';
lai_jules(97:192,29:140)=lai_julesa(:,1:96)';

%Load water mask data (from ESA landcover)
load /data/ESA_landcover/esa_05_landcover %Output from esa_lu_read.m
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
load /data/ESA_landcover/esa_jules_landcover %Output from esa_lu_read_julesgrid.m
oceanm_jules=NaN(192,144);
oceanm_jules(esa_jules'>200 & esa_jules'<220)=-1;

lats=-89.75:0.5:89.75;
lons=-179.75:0.5:179.75;
latincj=180/144;
lats_jules=-90+(latincj/2):latincj:90-(latincj/2);
lons_jules=-179.0625:1.875:179.0625;

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
data_mask='/data/turnover/';
[fmask,fmask_jules,~,~]=get_closed_can_mask(data_mask);

%Colourbar range
cmin=0;
cmax=12;

%Make figure
figure
cmap=colormap(parula(200));
cmap=[0.9 0.9 0.9; cmap];
colormap(cmap)
for nn=1:nmod
    if nn==2; continue; end %JULES is a special case because of resolution
    s(nn)=subplot(3,2,nn);
    hold on
    l(nn)=pcolor(lons,lats,oceanm');
    set(l(nn),'linestyle','none')
    p(nn)=pcolor(lons,lats,lai(:,:,nn).*fmask');
    set(p(nn),'linestyle','none')
    box on
    set(gca,'XTick',[],'XTickLabel','')
    set(gca,'YTick',[],'YTickLabel','')
    set(gca,'YLim',[-60 80])
    set(gca,'XLim',[-180 180])
    caxis([[cmin-(cmax/100) cmax]])
    t(nn)=text(-175,-50,models{nn});
    set(t(nn),'FontSize',12,'FontWeight','Bold')
end
c1=colorbar;
set(c1,'FontSize',12,'FontWeight','Bold')
set(c1,'Limits',[cmin cmax])
set(c1,'Location','southoutside')

%JULES is a special case because of resolution
s(2)=subplot(3,2,2);
hold on
l(2)=pcolor(lons_jules,lats_jules,oceanm_jules');
set(l(2),'linestyle','none')
p(2)=pcolor(lons_jules,lats_jules,lai_jules'.*fmask_jules');
set(p(2),'linestyle','none')
set(gca,'XTick',[],'XTickLabel','')
set(gca,'YTick',[],'YTickLabel','')
set(gca,'YLim',[-60 80])
set(gca,'XLim',[-180 180])
box on
caxis([[cmin-(cmax/100) cmax]])
t(2)=text(-175,-50,models{2});
set(t(2),'FontSize',12,'FontWeight','Bold')

set(s(1),'Position',[0.05 0.66 0.43 0.28])
set(s(2),'Position',[0.51 0.66 0.43 0.28])
set(s(3),'Position',[0.05 0.37 0.43 0.28])
set(s(4),'Position',[0.51 0.37 0.43 0.28])
set(s(5),'Position',[0.05 0.08 0.43 0.28])
set(s(6),'Position',[0.51 0.08 0.43 0.28])
set(c1,'Position',[0.05 0.04 0.89 0.01])
