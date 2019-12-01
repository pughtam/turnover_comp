%Make maps of dominant mortality process from CRU-NCEP or IPSL simulations
%The mortality process with the highest flux in a given grid cell is taken to be the dominant process.
%
%Groups causes of mortality from the models into the following conceptual groupings:
%1=Disturbance (including fire)
%2=Vitality
%3=Background
%4=Heat
%5=Other
%
%Requires *.mat files of preprocessed model data from turnover_pool_flux_read.m
%
%T. Pugh
%22.02.18

%Location of *.mat files containing preprocessed model data
data_models='/Users/pughtam/Documents/GAP_and_other_work/Mortality/mat_files/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/Users/pughtam/data/turnover/';
%Grid cell area file for JULES
jules_landarea_file='/Users/pughtam/data/turnover/JULES/JULES-LandMask.nc';
%Location of *.mat files containing preprocessed ESA landcover data
data_esa='/Users/pughtam/data/ESA_landcover/';

ifcruncep=true; %true=cruncep, false=ipsl

%Set year indexes for averaging
if ifcruncep
    y1=85; %1985
    y2=114; %2014
else %IPSL
    y1=170; %2070
    y2=199; %2099
end

%---
models={'(a) CABLE-POP','(b) JULES','(c) LPJ-GUESS','(d) LPJmL','(e) ORCHIDEE','(f) SEIB-DGVM'};
nmod=length(models);

lats=-90:0.5:89.5;
lons=-180:0.5:179.5;
latincj=180/144;
lats_jules=-90:latincj:90-latincj;
lons_jules=-180:1.875:180-1.875;

%Set some constants
sec_per_day=86400;
g_per_kg=1000;
days_per_year=365;
days_per_year_jules=360; %Only 360 days per year in JULES

%---
%Load data

dom_mort=NaN(360,720,nmod);
mort_dist=NaN(360,720,nmod);
mort_vital=NaN(360,720,nmod);
mort_back=NaN(360,720,nmod);
mort_heat=NaN(360,720,nmod);
mort_other=NaN(360,720,nmod);

%CABLE-POP
if ifcruncep
    load([data_models,'cable_map_arrays_cruncep.mat']);
else
    load([data_models,'cable_map_arrays_ipsl.mat']);
end
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
cable_dom_mort=zeros(720,360);
tdist_cable_mean=mean(tdist_cable(:,:,y1:y2),3)*uconv;
tvital_cable_mean=(mean(tgroweff_cable(:,:,y1:y2),3)+mean(tother_cable(:,:,y1:y2),3))*uconv; %Crowding is a form of vitality mortality
cable_dom_mort(tdist_cable_mean>tvital_cable_mean)=1; %Disturbance
cable_dom_mort(tvital_cable_mean>tdist_cable_mean)=2; %Vitality
cable_dom_mort(cable_dom_mort==0)=NaN;
dom_mort(:,:,1)=flipud(cable_dom_mort');
mort_dist(:,:,1)=flipud(tdist_cable_mean');
mort_vital(:,:,1)=flipud(tvital_cable_mean');
clear cwood_cable cleaf_cable croot_cable cveg_cable csoil_cable tleaf_cable troot_cable tgroweff_cable tdist_cable tother_cable npp_cable gpp_cable

%JULES
if ifcruncep
    load([data_models,'jules_map_arrays_cruncep.mat']);
else
    load([data_models,'jules_map_arrays_ipsl.mat']);
end
uconv=sec_per_day*days_per_year_jules; %Conversion factor for fluxes from seconds to years
jules_dom_mort=zeros(144,192);
tback_jules_mean=NaN(144,192);
tother_jules_mean=NaN(144,192);
tback_jules_mean(29:140,1:96)=mean(tback_jules(97:192,:,y1:y2),3)'*uconv;
tback_jules_mean(29:140,97:192)=mean(tback_jules(1:96,:,y1:y2),3)'*uconv;
tother_jules_mean(29:140,1:96)=mean(tother_jules(97:192,:,y1:y2),3)'*uconv;
tother_jules_mean(29:140,97:192)=mean(tother_jules(1:96,:,y1:y2),3)'*uconv;
jules_dom_mort(tback_jules_mean>tother_jules_mean)=3; %Background (Wood turnover)
jules_dom_mort(tother_jules_mean>tback_jules_mean)=2; %Vitality (Competition mortality is essentially vitality based)
jules_dom_mort(jules_dom_mort==0)=NaN;
dom_mort_jules=jules_dom_mort;
mort_back_jules=tback_jules_mean;
mort_vital_jules=tother_jules_mean;
clear cwood_jules cleaf_jules croot_jules cveg_jules csoil_jules tleaf_jules troot_jules tback_jules tother_jules npp_jules gpp_jules

%LPJ-GUESS
if ifcruncep
    load([data_models,'lpjg_map_arrays_cruncep.mat']);
else
    load([data_models,'lpjg_map_arrays_ipsl.mat']);
end
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
lpjg_dom_mort=zeros(720,360);
tdist_lpjg_mean=(mean(tdist_lpjg(:,:,y1:y2),3)+mean(tfire_lpjg(:,:,y1:y2),3))*uconv;
tvital_lpjg_mean=mean(tgroweff_lpjg(:,:,y1:y2),3)*uconv;
tback_lpjg_mean=mean(tback_lpjg(:,:,y1:y2),3)*uconv;
tother_lpjg_mean=(mean(tother_lpjg(:,:,y1:y2),3)-mean(tback_lpjg(:,:,y1:y2),3))*uconv;
lpjg_dom_mort(tdist_lpjg_mean>tvital_lpjg_mean & tdist_lpjg_mean>tother_lpjg_mean & tdist_lpjg_mean>tback_lpjg_mean)=1; %Disturbance
lpjg_dom_mort(tvital_lpjg_mean>tdist_lpjg_mean & tvital_lpjg_mean>tother_lpjg_mean & tvital_lpjg_mean>tback_lpjg_mean)=2; %Vitality
lpjg_dom_mort(tback_lpjg_mean>tdist_lpjg_mean & tback_lpjg_mean>tother_lpjg_mean & tback_lpjg_mean>tvital_lpjg_mean)=3; %Background
lpjg_dom_mort(tother_lpjg_mean>tdist_lpjg_mean & tother_lpjg_mean>tvital_lpjg_mean & tother_lpjg_mean>tback_lpjg_mean)=5; %Other
lpjg_dom_mort(lpjg_dom_mort==0)=NaN;
dom_mort(:,:,3)=lpjg_dom_mort';
mort_dist(:,:,3)=tdist_lpjg_mean';
mort_vital(:,:,3)=tvital_lpjg_mean';
mort_back(:,:,3)=tback_lpjg_mean';
mort_other(:,:,3)=tother_lpjg_mean';
clear allfire_lpjg cwood_lpjg cleaf_lpjg croot_lpjg cveg_lpjg csoil_lpjg tleaf_lpjg troot_lpjg tfire_lpjg tgroweff_lpjg tdist_lpjg tother_lpjg npp_lpjg gpp_lpjg
clear tbioclim_lpjg tback_lpjg tnbio_lpjg tallom_lpjg

%LPJmL
if ifcruncep
    load([data_models,'lpjml_map_arrays_cruncep.mat']);
else
    load([data_models,'lpjml_map_arrays_ipsl.mat']);
end
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
lpjml_dom_mort=zeros(720,360);
tdist_lpjml_mean=mean(tfire_lpjml(:,:,y1:y2),3)*uconv;
tvital_lpjml_mean=(mean(tgroweff_lpjml(:,:,y1:y2),3)+mean(tshade_lpjml(:,:,y1:y2),3))*uconv;
theat_lpjml_mean=mean(theat_lpjml(:,:,y1:y2),3)*uconv;
tother_lpjml_mean=(mean(tother_lpjml(:,:,y1:y2),3)-mean(theat_lpjml(:,:,y1:y2),3)-mean(tshade_lpjml(:,:,y1:y2),3))*uconv;
lpjml_dom_mort(tdist_lpjml_mean>tvital_lpjml_mean & tdist_lpjml_mean>tother_lpjml_mean & tdist_lpjml_mean>theat_lpjml_mean)=1; %Disturbance
lpjml_dom_mort(tvital_lpjml_mean>tdist_lpjml_mean & tvital_lpjml_mean>tother_lpjml_mean & tvital_lpjml_mean>theat_lpjml_mean)=2; %Vitality
lpjml_dom_mort(theat_lpjml_mean>tdist_lpjml_mean & theat_lpjml_mean>tother_lpjml_mean & theat_lpjml_mean>tvital_lpjml_mean)=4; %Heat
lpjml_dom_mort(tother_lpjml_mean>tvital_lpjml_mean & tother_lpjml_mean>tdist_lpjml_mean & tother_lpjml_mean>theat_lpjml_mean)=5; %Other
lpjml_dom_mort(lpjml_dom_mort==0)=NaN;
dom_mort(:,:,4)=lpjml_dom_mort';
mort_dist(:,:,4)=tdist_lpjml_mean';
mort_vital(:,:,4)=tvital_lpjml_mean';
mort_heat(:,:,4)=theat_lpjml_mean';
mort_other(:,:,4)=tother_lpjml_mean';
clear cwood_lpjml cleaf_lpjml croot_lpjml cveg_lpjml csoil_lpjml tleaf_lpjml troot_lpjml tfire_lpjml tgroweff_lpjml tother_lpjml npp_lpjml gpp_lpjml
clear tnbio_lpjml tallom_lpjml theat_lpjml tshade_lpjml

%ORCHIDEE
if ifcruncep
    load([data_models,'orchidee_map_arrays_cruncep.mat']);
else
    load([data_models,'orchidee_map_arrays_ipsl.mat']);
end
uconv=days_per_year/g_per_kg; %Conversion factor for fluxes from seconds to years
tback_orchidee_mean=mean(tmort_orchidee(:,:,y1:y2),3)*uconv;
orchidee_dom_mort=ones(720,360)*6; %ORCHIDEE doesn't provide a breakdown
dom_mort(:,:,5)=flipud(orchidee_dom_mort');
mort_back(:,:,5)=flipud(tback_orchidee_mean');
clear cwood_orchidee cleaf_orchidee croot_orchidee cveg_orchidee csoil_orchidee tleaf_orchidee troot_orchidee tfire_orchidee tgroweff_orchidee tother_orchidee npp_orchidee gpp_orchidee
clear tnbio_orchidee tallom_orchidee theat_orchidee tshade_orchidee

%SEIB-DGVM
if ifcruncep
    load([data_models,'seib_map_arrays_cruncep.mat']);
else
    load([data_models,'seib_map_arrays_ipsl.mat']);
end
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
seib_dom_mort=zeros(720,360);
tdist_seib_mean=(mean(tdist_seib(:,:,y1:y2),3)+mean(tfire_seib(:,:,y1:y2),3))*uconv;
tvital_seib_mean=mean(tgroweff_seib(:,:,y1:y2),3)*uconv;
theat_seib_mean=mean(theat_seib(:,:,y1:y2),3)*uconv;
tback_seib_mean=mean(tback_seib(:,:,y1:y2),3)*uconv;
tother_seib_mean=(mean(tother_seib(:,:,y1:y2),3)-mean(theat_seib(:,:,y1:y2),3)-mean(tback_seib(:,:,y1:y2),3))*uconv;
seib_dom_mort(tdist_seib_mean>tvital_seib_mean & tdist_seib_mean>tother_seib_mean & tdist_seib_mean>theat_seib_mean & tdist_seib_mean>tback_seib_mean)=1; %Disturbance
seib_dom_mort(tvital_seib_mean>tdist_seib_mean & tvital_seib_mean>tother_seib_mean & tvital_seib_mean>theat_seib_mean & tvital_seib_mean>tback_seib_mean)=2; %Vitality
seib_dom_mort(tback_seib_mean>tdist_seib_mean & tback_seib_mean>tother_seib_mean & tback_seib_mean>theat_seib_mean & tback_seib_mean>tvital_seib_mean)=3; %Background
seib_dom_mort(theat_seib_mean>tdist_seib_mean & theat_seib_mean>tother_seib_mean & theat_seib_mean>tback_seib_mean & theat_seib_mean>tvital_seib_mean)=4; %Heat (note, SEIB-DGVM did not have heat mortality enabled in these simulations)
seib_dom_mort(tother_seib_mean>tdist_seib_mean & tother_seib_mean>tvital_seib_mean & tother_seib_mean>tback_seib_mean & tother_seib_mean>theat_seib_mean)=5; %Other
seib_dom_mort(seib_dom_mort==0)=NaN;
dom_mort(:,:,6)=flipud(seib_dom_mort');
mort_dist(:,:,6)=flipud(tdist_seib_mean');
mort_vital(:,:,6)=flipud(tvital_seib_mean');
mort_heat(:,:,6)=flipud(theat_seib_mean');
mort_back(:,:,6)=flipud(tback_seib_mean');
mort_other(:,:,6)=flipud(tother_seib_mean');
clear cwood_seib cleaf_seib croot_seib cveg_seib csoil_seib tleaf_seib troot_seib tfire_seib tgroweff_seib tdist_seib tother_seib npp_seib gpp_seib
clear theat_seib tbioclim_seib tback_seib

%---
%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[fmask,fmask_jules,fmaskscale,fmaskscale_jules]=get_closed_can_mask(data_mask);

dom_mort=dom_mort.*repmat(fmask',1,1,nmod);
jules_dom_mort=jules_dom_mort.*fmask_jules';

%Calculate totals by mortality class
garea=global_grid_area();
garea_jules_in=ncread(jules_landarea_file,'arealand');
garea_jules=NaN(size(garea_jules_in));
garea_jules(1:96,:)=garea_jules_in(97:192,:);
garea_jules(97:192,:)=garea_jules_in(1:96,:);
garea_jules=garea_jules';
garea_jules(1,:)=[];

%Total flux per grid cell for closed-canopy forest area
mort_dist=mort_dist.*repmat(fmaskscale'.*garea,1,1,nmod);
mort_vital=mort_vital.*repmat(fmaskscale'.*garea,1,1,nmod);
mort_back=mort_back.*repmat(fmaskscale'.*garea,1,1,nmod);
mort_heat=mort_heat.*repmat(fmaskscale'.*garea,1,1,nmod);
mort_other=mort_other.*repmat(fmaskscale'.*garea,1,1,nmod);
mort_back_jules=mort_back_jules.*fmaskscale_jules'.*garea_jules;
mort_vital_jules=mort_vital_jules.*fmaskscale_jules'.*garea_jules;

%Sum per grid cell
mort_dist_sum=squeeze(nansum(nansum(mort_dist,2),1));
mort_vital_sum=squeeze(nansum(nansum(mort_vital,2),1));
mort_back_sum=squeeze(nansum(nansum(mort_back,2),1));
mort_heat_sum=squeeze(nansum(nansum(mort_heat,2),1));
mort_other_sum=squeeze(nansum(nansum(mort_other,2),1));
mort_vital_sum(2)=squeeze(nansum(nansum(mort_vital_jules,2),1));
mort_back_sum(2)=squeeze(nansum(nansum(mort_back_jules,2),1));

%Group arrays and calculate the percentage contribution of each process to the total flux
mort_all_sum=cat(2,mort_dist_sum,mort_vital_sum,mort_back_sum,mort_heat_sum,mort_other_sum)/1e12;
mort_all_sum_rel=mort_all_sum./repmat(sum(mort_all_sum,2),1,5);

%Calculate zonal sums of fluxes
mort_dist_zonsum=squeeze(nansum(mort_dist,2));
mort_vital_zonsum=squeeze(nansum(mort_vital,2));
mort_back_zonsum=squeeze(nansum(mort_back,2));
mort_heat_zonsum=squeeze(nansum(mort_heat,2));
mort_other_zonsum=squeeze(nansum(mort_other,2));
mort_vital_jules_zonsum=squeeze(nansum(mort_vital_jules,2));
mort_back_jules_zonsum=squeeze(nansum(mort_back_jules,2));

%Group arrays and calculate the percentage contribution of each process to the total zonal flux
mort_all_zonsum=cat(3,mort_dist_zonsum,mort_vital_zonsum,mort_back_zonsum,mort_heat_zonsum,mort_other_zonsum)/1e9;
mort_all_zonsum_rel=mort_all_zonsum./repmat(nansum(mort_all_zonsum,3),1,1,5);

mort_all_jules_zonsum=cat(2,NaN(144,1),mort_vital_jules_zonsum,mort_back_jules_zonsum,NaN(144,1),NaN(144,1))/1e9;
mort_all_jules_zonsum_rel=mort_all_jules_zonsum./repmat(nansum(mort_all_jules_zonsum,2),1,5);

mort_all_zonsum(isnan(mort_all_zonsum))=0;
mort_all_zonsum_rel(isnan(mort_all_zonsum_rel))=0;
mort_all_jules_zonsum(isnan(mort_all_jules_zonsum))=0;
mort_all_jules_zonsum_rel(isnan(mort_all_jules_zonsum_rel))=0;

%Load water mask data (from ESA landcover)
load([data_esa,'esa_05_landcover']); %Output from esa_lu_read.m
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
load([data_esa,'esa_jules_landcover']); %Output from esa_lu_read_julesgrid.m
oceanm_jules=NaN(192,144);
oceanm_jules(esa_jules'>200 & esa_jules'<220)=-1;

%---
%Combined figure

figure
cmap=colormap(lines(5));
cmap=[0.9 0.9 0.9; cmap; 0.5 0.5 0.5];
subplotsmap=[2,3,6,7,10,11,14,15];
subplotsbar=[1,4,5,8,9,12,13,16];
for nn=1:nmod
    if nn==2; continue; end %JULES is a special case because of resolution
    s(nn)=subplot(3,4,subplotsmap(nn));
    hold on
    l(nn)=pcolor(lons,lats,oceanm');
    set(l(nn),'linestyle','none')
    p(nn)=pcolor(lons,lats,dom_mort(:,:,nn));
    set(p(nn),'linestyle','none')
    box on
    set(gca,'XTick',[],'XTickLabel','')
    set(gca,'YTick',[],'YTickLabel','')
    set(gca,'YLim',[-60 80])
    set(gca,'XLim',[-180 180])
    caxis([0 6])
    t(nn)=text(-175,-50,models{nn});
    set(t(nn),'FontSize',12,'FontWeight','Bold')
end
clear nn
c1=colorbar;
colormap(cmap)
set(c1,'Limits',[1 5])
set(c1,'Ticks',[1.35 2.1 2.9 3.7 4.55])
set(c1,'TickLabels',{'Dist.','Vitality','Background','Heat','Other'})

%JULES is a special case because of resolution
s(2)=subplot(3,4,3);
hold on
l(2)=pcolor(lons_jules,lats_jules,oceanm_jules');
set(l(2),'linestyle','none')
p(2)=pcolor(lons_jules,lats_jules,jules_dom_mort);
set(p(2),'linestyle','none')
set(gca,'XTick',[],'XTickLabel','')
set(gca,'YTick',[],'YTickLabel','')
set(gca,'YLim',[-60 80])
set(gca,'XLim',[-180 180])
box on
caxis([0 5])
t(2)=text(-175,-50,models{2});
set(t(2),'FontSize',12,'FontWeight','Bold')

%Relative bars
aaa=NaN(nmod,6);
for nn=1:nmod
    ss(nn)=subplot(3,4,subplotsbar(nn));
    if nn~=5
        plotarray=cat(2,squeeze(mort_all_zonsum_rel(:,nn,:)),zeros(360,1)); %Because first and last colors in colormap are grey
    else
        plotarray=cat(2,zeros(360,5),ones(360,1)); %Just plot ones in the last column to give grey output as no breakdown for ORCHIDEE
    end
    aaa(nn,:)=area(lats,plotarray);
    set(gca,'YLim',[0 1],'XLim',[-60 80])
    view(90,270)
    set(gca,'XTickLabel','','YTickLabel','')
    set(gca,'XTick','','YTick','')
end
clear nn
ss(2)=subplot(3,4,4);
plotarray=cat(2,mort_all_jules_zonsum_rel,zeros(144,1));
aaa(2,:)=area(lats_jules,plotarray);
set(gca,'YLim',[0 1],'XLim',[-60 80])
view(90,270)
set(gca,'XTickLabel','','YTickLabel','')
set(gca,'XTick','','YTick','')
for nn=1:nmod
    for mm=1:6
        set(aaa(nn,mm),'LineStyle','none')
    end
end
clear nn mm plotarray

set(s(1),'Position',[0.08 0.67 0.36 0.31])
set(s(2),'Position',[0.45 0.67 0.36 0.31])
set(s(3),'Position',[0.08 0.35 0.36 0.31])
set(s(4),'Position',[0.45 0.35 0.36 0.31])
set(s(5),'Position',[0.08 0.03 0.36 0.31])
set(s(6),'Position',[0.45 0.03 0.36 0.31])

set(ss(1),'Position',[0.02 0.67 0.05 0.31])
set(ss(2),'Position',[0.82 0.67 0.05 0.31])
set(ss(3),'Position',[0.02 0.35 0.05 0.31])
set(ss(4),'Position',[0.82 0.35 0.05 0.31])
set(ss(5),'Position',[0.02 0.03 0.05 0.31])
set(ss(6),'Position',[0.82 0.03 0.05 0.31])

set(c1,'Position',[0.88 0.35 0.025 0.29])

for nn=1:nmod
    subplot(s(nn))
    pp=get(s(nn),'Position');
    aa(nn)=axes('Position',[pp(1)+0.16,pp(2)+0.01,0.18,0.02]);
    if nn~=5
        barh([mort_all_sum_rel(nn,:);NaN(1,5)]*100,'stacked')
    else
        barh([zeros(1,5),1;NaN(1,6)]*100,'stacked')
    end
    set(gca,'YLim',[0.6 1.4],'XLim',[0 100])
    set(gca,'YTick',[],'XTick',[])
    cc=get(gca,'Children');
    cmapflip=flipud(cmap);
    for mm=1:length(cc)
        if nn~=5
            set(cc(mm),'FaceColor',cmapflip(mm+1,:),'EdgeColor',cmapflip(mm+1,:))
        else
            set(cc(mm),'FaceColor',cmapflip(mm,:),'EdgeColor',cmapflip(mm,:)) %Ensure that ORCHIDEE bar is grey
        end
    end
    clear mm cc
    box off
end
clear nn