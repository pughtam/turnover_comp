%Make global maps of various productivity and turnover-related variables
%
%Requires *.mat files of preprocessed model data from turnover_pool_flux_read.m
%Requires *mat files of preprocessed observational data (see notes in global_obs_totals.m)
%
%T. Pugh
%27.11.19

%Location of *.mat files containing preprocessed model data
data_models='/Users/pughtam/Documents/GAP_and_other_work/Mortality/mat_files/';
%Location of *.mat files containing preprocessed observational data
data_obs='/Users/pughtam/Documents/GAP_and_other_work/Mortality/Tim_Dec_plots/raw_data/raw_data_v2/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/Users/pughtam/data/turnover/';
%Path to file containing grid cell areas for JULES (m2)
jules_gridcellarea_file='/data/turnover/JULES/JULES-LandMask.nc';
%Location of *.mat files containing preprocessed ESA landcover data
data_esa='/Users/pughtam/data/ESA_landcover/';

%Choose which variables to plot. 
%1=Vegetation carbon stock
%2=Fraction of total turnover due to mortality
%3=Annual NPP
%4=tau_turn
%5=tau_NPP
%6=Ratio of tau_turn to tau_NPP (%)
plotvar=5;

inc_obs=false; %Also include observational data? (only for NPP, tau_turn, Cveg)

%Index of first and last year to average over in the data
y1=85; %1985
y2=114; %2014

%---
if inc_obs && (plotvar~=1 && plotvar~=3 && plotvar~=5)
    error('inc_obs is true but plotvar is not a compatible value')
end

%---
%Load in the pre-processed model data

ifcruncep=true;
[vstock,mflux,lrflux,nppflux,reproflux,...
    vstock_jules,mflux_jules,lrflux_jules,nppflux_jules,models,nmod]...
    =get_stocks_fluxes(data_models,ifcruncep,y1,y2);

%---
%If necessary, read in the observational data

if (plotvar==1 || plotvar==3 || plotvar==5)
    load([data_obs,'obs_NPP_0perc.mat']);
    load([data_obs,'obs_Cveg_0perc.mat']);
    Cveg=Cveg';
    NPP=NPP';
end

%---
%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[fmask,fmask_jules,ffrac,ffrac_jules]=get_closed_can_mask(data_mask);

%---
%Make plots

lats=-89.75:0.5:89.75;
lons=-179.75:0.5:179.75;
latincj=180/144;
lats_jules=-90+(latincj/2):latincj:90-(latincj/2);
lons_jules=-179.0625:1.875:179.0625;

%Load water mask data (from ESA landcover)
load([data_esa,'esa_05_landcover.mat']); %Output from esa_hires_region_mask.m	
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
load([data_esa,'esa_jules_landcover.mat']); %Output from esa_hires_region_mask_jules.m	
oceanm_jules=NaN(192,144);
oceanm_jules(esa_jules'>200 & esa_jules'<220)=-1;

%Choose data to plot
if plotvar==1
    %Vegetation carbon stock
    plotarray=vstock;
    plotarray_jules=vstock_jules;
    if inc_obs
        plotarray_obs=Cveg;
    end
    cmin=0; %Settings for colour scale
    cmax=40;
elseif plotvar==2
    %Fraction of total turnover due to mortality
    plotarray=mflux./nansum(cat(4,mflux,lrflux,reproflux),4);
    plotarray_jules=mflux_jules./(mflux_jules+lrflux_jules);
    cmin=0; %Settings for colour scale
    cmax=1;
elseif plotvar==3
    %Annual NPP
    plotarray=nppflux;
    plotarray_jules=nppflux_jules;
    if inc_obs
        plotarray_obs=NPP;
    end
    cmin=0; %Settings for colour scale
    cmax=2.5;
elseif plotvar==4
    %tau_turn
    plotarray=vstock./nansum(cat(4,mflux,lrflux,reproflux),4);
    plotarray_jules=vstock_jules./(mflux_jules+lrflux_jules);
    cmin=0; %Settings for colour scale
    cmax=40;
elseif plotvar==5
    %tau_NPP
    plotarray=vstock./nppflux;
    plotarray_jules=vstock_jules./nppflux_jules;
    if inc_obs
        plotarray_obs=Cveg./NPP;
    end
    cmin=0; %Settings for colour scale
    cmax=30;
elseif plotvar==6
    %Ratio of tau_turn to tau_NPP (%)
    plotarray=(((vstock./nansum(cat(4,mflux,lrflux,reproflux),4)) ./ (vstock./nppflux))-1)*100;
    plotarray_jules=(((vstock_jules./(mflux_jules+lrflux_jules)) ./ (vstock_jules./nppflux_jules))-1)*100;
    cmin=0; %Settings for colour scale
    cmax=100;
end

plotarray=plotarray.*repmat(fmask',1,1,nmod);
plotarray_jules=plotarray_jules.*fmask_jules';
if inc_obs
    plotarray_obs=plotarray_obs.*fmask';
end

figure
cmap=colormap(parula(200));
cmap=[0.9 0.9 0.9; cmap];
colormap(cmap)
for nn=1:nmod
    if nn==2 || nn>nmod; continue; end %JULES is a special case because of resolution
    if inc_obs
        s(nn)=subplot(4,2,nn);
    else
        s(nn)=subplot(3,2,nn);
    end
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
    caxis([cmin-(cmax/100) cmax])
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
caxis([cmin-(cmax/100) cmax])
t(2)=text(-175,-50,models{2});
set(t(2),'FontSize',12,'FontWeight','Bold')

if inc_obs
    %Observations also added separately
    s(7)=subplot(4,2,7);
    hold on
    l(2)=pcolor(lons,lats,oceanm');
    set(l(2),'linestyle','none')
    p(2)=pcolor(lons,lats,plotarray_obs);
    set(p(2),'linestyle','none')
    set(gca,'XTick',[],'XTickLabel','')
    set(gca,'YTick',[],'YTickLabel','')
    set(gca,'YLim',[-60 80])
    set(gca,'XLim',[-180 180])
    box on
    caxis([cmin-(cmax/100) cmax])
    t(2)=text(-175,-50,'(g) Satellite-based');
    set(t(2),'FontSize',12,'FontWeight','Bold')
end

if inc_obs
    set(s(1),'Position',[0.05 0.76 0.43 0.22])
    set(s(2),'Position',[0.51 0.76 0.43 0.22])
    set(s(3),'Position',[0.05 0.53 0.43 0.22])
    set(s(4),'Position',[0.51 0.53 0.43 0.22])
    set(s(5),'Position',[0.05 0.30 0.43 0.22])
    set(s(6),'Position',[0.51 0.30 0.43 0.22])
    set(s(7),'Position',[0.05 0.07 0.43 0.22])
else
    set(s(1),'Position',[0.05 0.66 0.43 0.28])
    set(s(2),'Position',[0.51 0.66 0.43 0.28])
    set(s(3),'Position',[0.05 0.37 0.43 0.28])
    set(s(4),'Position',[0.51 0.37 0.43 0.28])
    set(s(5),'Position',[0.05 0.08 0.43 0.28])
    set(s(6),'Position',[0.51 0.08 0.43 0.28])
end
