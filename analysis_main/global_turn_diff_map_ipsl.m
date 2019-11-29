%Script to make difference maps of turnover time between a baseline and future period for the IPSL simulation
%Need to run turnover_pool_flux_read.m first to make the *.mat files used by this script.
%
%Dependencies
%-get_stocks_fluxes.m
%
%T. Pugh
%28.11.19

%Location of *.mat files containing preprocessed model data
data_models='/data/turnover/';
%Location of *.mat files containing preprocessed ESA landcover data
data_esa='/data/ESA_landcover/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/data/turnover/';

%---

%Index of first year for the baseline and future period
y1=[85 170]; %1985, 2070
%Index of last year for the baseline and future period
y2=[114 199]; %2014, 2099

nt=length(y1);

models={'(a) CABLE-POP','(b) JULES','(c) LPJ-GUESS','(d) LPJmL','(f) ORCHIDEE','(g) SEIB-DGVM'};
nmod=length(models);

vstock=NaN(360,720,nmod,nt); %Vegetation carbon stock
mflux=NaN(360,720,nmod,nt); %Annual turnover flux through all mortality causes
lrflux=NaN(360,720,nmod,nt); %Annual turnover flux through leaves and fine roots
nppflux=NaN(360,720,nmod,nt); %Annual NPP
reproflux=NaN(360,720,nmod,nt); %Annual turnover flux through reproduction
vstock_jules=NaN(144,192,nt); %Vegetation carbon stock (JULES only)
mflux_jules=NaN(144,192,nt); %Annual turnover flux through all mortality causes (JULES only)
lrflux_jules=NaN(144,192,nt); %Annual turnover flux through leaves and fine roots (JULES only)
nppflux_jules=NaN(144,192,nt); %Annual NPP (JULES only)

%Get the stock and flux data for the baseline and future periods
ifcruncep=false; %We want to load IPSL data
for nn=1:nt
    [vstock(:,:,:,nn),mflux(:,:,:,nn),lrflux(:,:,:,nn),nppflux(:,:,:,nn),reproflux(:,:,:,nn),...
        vstock_jules(:,:,nn),mflux_jules(:,:,nn),lrflux_jules(:,:,nn),nppflux_jules(:,:,nn),...
        models,nmod]=get_stocks_fluxes(data_models,ifcruncep,y1(nn),y2(nn));
end
clear nn

%Add all the flux components together into a total turnover flux
totturnflux=nansum(cat(5,mflux,lrflux,reproflux),5);
totturnflux_jules=nansum(cat(4,mflux_jules,lrflux_jules),4);

%Calculate turnover time due to mortality
tau_mort_base=vstock(:,:,:,1)./mflux(:,:,:,1);
tau_mort_fut=vstock(:,:,:,2)./mflux(:,:,:,2);
tau_mort_jules_base=vstock_jules(:,:,1)./mflux_jules(:,:,1);
tau_mort_jules_fut=vstock_jules(:,:,2)./mflux_jules(:,:,2);

%Difference in turnover time due to mortality between baseline and future periods (%)
diff_tau_mort=((tau_mort_fut./tau_mort_base)-1)*100;
diff_tau_mort_jules=((tau_mort_jules_fut./tau_mort_jules_base)-1)*100;

%---
%Make plots

%Load water mask data (from ESA landcover)
load([data_esa,'esa_05_landcover']); %Output from esa_lu_read.m
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1000;
load([data_esa,'esa_jules_landcover']); %Output from esa_lu_read_julesgrid.m
oceanm_jules=NaN(192,144);
oceanm_jules(esa_jules'>200 & esa_jules'<220)=-1000;

lats=-89.75:0.5:89.75;
lons=-179.75:0.5:179.75;
latincj=180/144;
lats_jules=-90+(latincj/2):latincj:90-(latincj/2);
lons_jules=-179.0625:1.875:179.0625;

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[fmask,fmask_jules,~,~]=get_closed_can_mask(data_mask);

%Apply forest cover mask
plotarray=diff_tau_mort.*repmat(fmask',1,1,nmod);
plotarray_jules=diff_tau_mort_jules.*fmask_jules';

%Settings for colour scale
cmin=-100;
cmax=100;

%Make plot
figure
cmap=flipud(redblue(200));
cmap=[0.9 0.9 0.9; cmap];
colormap(cmap)
for nn=1:nmod
    if nn==2; continue; end %JULES is a special case because of resolution
    s(nn)=subplot(3,2,nn);
    hold on
    l(nn)=pcolor(lons,lats,oceanm');
    set(l(nn),'linestyle','none')
    p(nn)=pcolor(lons,lats,plotarray(:,:,nn));
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
set(c1,'Location','southoutside')
set(c1,'Position',[0.05 0.04 0.89 0.01])
%JULES is a special case because of resolution
s(2)=subplot(3,2,2);
hold on
l(2)=pcolor(lons_jules,lats_jules,oceanm_jules');
set(l(2),'linestyle','none')
p(2)=pcolor(lons_jules,lats_jules,plotarray_jules);
set(p(2),'linestyle','none')
set(gca,'XTick',[],'XTickLabel','')
set(gca,'YTick',[],'YTickLabel','')
set(gca,'YLim',[-60 80])
set(gca,'XLim',[-180 180])
box on
caxis([cmin-1 cmax])
t(2)=text(-175,-50,models{2});
set(t(2),'FontSize',12,'FontWeight','Bold')

set(s(1),'Position',[0.05 0.66 0.43 0.28])
set(s(2),'Position',[0.51 0.66 0.43 0.28])
set(s(3),'Position',[0.05 0.37 0.43 0.28])
set(s(4),'Position',[0.51 0.37 0.43 0.28])
set(s(5),'Position',[0.05 0.08 0.43 0.28])
set(s(6),'Position',[0.51 0.08 0.43 0.28])
