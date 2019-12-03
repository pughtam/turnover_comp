%Script to make maps of forest type for the CRU-NCEP simulations based on masks created using scripts in
%"forest_masks" directory (by Sarah Shafer) 
%
%The indexes to the forest-types are as follows:
%1=needleleaved evergreen, 
%2=needleleaved deciduous, 
%3=boreal broadleaved deciduous, 
%4=temperate broadleaved evergreen, 
%5=temperate broadleaved deciduous, 
%6=tropical broadleaved evergreen,
%7=tropical broadleaved deciduous
%
%T. Pugh
%22.02.18

%Location of forest mask files created using file in folder "forest_masks"
phen_dir='/data/turnover/masks/phen/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/Users/pughtam/data/turnover/';
%Location of *.mat files containing preprocessed ESA landcover data
data_esa='/Users/pughtam/data/ESA_landcover/';

%---
phen_label={'needleleaved evergreen','needleleaved deciduous','boreal broadleaved deciduous',...
    'temperate broadleaved evergreen','temperate broadleaved deciduous','tropical broadleaved evergreen',...
    'tropical broadleaved raingreen'};

models={'(a) CABLE-POP','(b) JULES','(c) LPJ-GUESS','(d) LPJmL','(e) LPJ-wsl','(f) ORCHIDEE','(g) SEIB-DGVM'};
nmod=length(models);

phen=NaN(360,720,nmod);

%Load in the forest mask files
phen(:,:,1)=flipud(ncread([phen_dir,'/CABLE-POP_cruncep_lai_annual_1901_2015_phenology_mask.nc'],'phen_max_lai_phen_number')');
phen_julesa=ncread([phen_dir,'/JULESC2_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')';
phen(:,:,3)=ncread([phen_dir,'/lpj-guess_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')';
phen(69:347,:,4)=ncread([phen_dir,'/lpjml_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')';
phen(:,:,5)=flipud(ncread([phen_dir,'/lpj-wsl_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')');
phen(:,:,6)=flipud(ncread([phen_dir,'/ORCHIDEE_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')');
phen(:,:,7)=flipud(ncread([phen_dir,'/seib_cruncep_lai_annual_1901_2014_phenology_mask.nc'],'phen_max_lai_phen_number')');

phen_jules=NaN(192,144);
phen_jules(1:96,29:140)=phen_julesa(:,97:192)';
phen_jules(97:192,29:140)=phen_julesa(:,1:96)';

%Load water mask data (from ESA landcover)
load([data_esa,'esa_05_landcover']); %Output from esa_hires_region_mask.m	
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
load([data_esa,'esa_jules_landcover']); %Output from esa_hires_region_mask_jules.m	
oceanm_jules=NaN(192,144);
oceanm_jules(esa_jules'>200 & esa_jules'<220)=-1;

lats=-89.75:0.5:89.75;
lons=-179.75:0.5:179.75;
latincj=180/144;
lats_jules=-90+(latincj/2):latincj:90-(latincj/2);
lons_jules=-179.0625:1.875:179.0625;

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[fmask,fmask_jules,~,~]=get_closed_can_mask(data_mask);

%Set limits for colourbar
cmin=1;
cmax=8;

%Make the figure
figure
cmap=colormap(parula(7));
cmap=[0.9 0.9 0.9; cmap];
colormap(cmap)
for nn=1:nmod
    if nn==2; continue; end %JULES is a special case because of resolution
    s(nn)=subplot(4,2,nn);
    hold on
    l(nn)=pcolor(lons,lats,oceanm');
    set(l(nn),'linestyle','none')
    p(nn)=pcolor(lons,lats,phen(:,:,nn).*fmask');
    set(p(nn),'linestyle','none')
    box on
    set(gca,'XTick',[],'XTickLabel','')
    set(gca,'YTick',[],'YTickLabel','')
    set(gca,'YLim',[-60 80])
    set(gca,'XLim',[-180 180])
    caxis([cmin-1 cmax])
    t(nn)=text(-175,-50,models{nn});
    set(t(nn),'FontSize',12,'FontWeight','Bold')
end
c1=colorbar;
set(c1,'FontSize',12,'FontWeight','Bold')
set(c1,'Limits',[cmin cmax])
set(c1,'Location','eastoutside')
set(c1,'Ticks',1.5:1:8,'TickLabels',phen_label)

%JULES is a special case because of resolution
s(2)=subplot(4,2,2);
hold on
l(2)=pcolor(lons_jules,lats_jules,oceanm_jules');
set(l(2),'linestyle','none')
p(2)=pcolor(lons_jules,lats_jules,phen_jules'.*fmask_jules');
set(p(2),'linestyle','none')
set(gca,'XTick',[],'XTickLabel','')
set(gca,'YTick',[],'YTickLabel','')
set(gca,'YLim',[-60 80])
set(gca,'XLim',[-180 180])
box on
caxis([cmin-1 cmax])
t(2)=text(-175,-50,models{2});
set(t(2),'FontSize',12,'FontWeight','Bold')

set(s(1),'Position',[0.05 0.76 0.44 0.22])
set(s(2),'Position',[0.51 0.76 0.44 0.22])
set(s(3),'Position',[0.05 0.53 0.44 0.22])
set(s(4),'Position',[0.51 0.53 0.44 0.22])
set(s(5),'Position',[0.05 0.30 0.44 0.22])
set(s(6),'Position',[0.51 0.30 0.44 0.22])
set(s(7),'Position',[0.05 0.07 0.44 0.22])
