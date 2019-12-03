%Script to make line plots of standardised mortality rates by forest type
%(standardised to C veg) and line plots of NPP sums by forest type
%
%Dependencies
%-global_grid_area.m
%
%T. Pugh
%25.02.18

%Location of *.mat files containing preprocessed model data
data_models='/data/turnover/';
%Path to file containing grid cell areas for JULES (m2)
jules_gridcellarea_file='/data/turnover/JULES/JULES-LandMask.nc';
%Location of *.mat files containing preprocessed ESA landcover data
data_esa='/data/ESA_landcover/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/data/turnover/';

%---
models={'(a) CABLE-POP','(b) JULES','(c) LPJ-GUESS','(d) LPJmL','(e) ORCHIDEE','(f) SEIB-DGVM'};
modelsonly={'CABLE-POP','JULES','LPJ-GUESS','LPJmL','ORCHIDEE','SEIB-DGVM'};
nmod=length(models);

%Index of first and last year of timeseries
y1=85; %1985
y2=199; %2099
nyear=length(y1:y2);

%Set some constants
sec_per_day=86400;
g_per_kg=1000;
days_per_year=365;
days_per_year_jules=360; %Only 360 days per year in JULES

%---
%Initialise arrays to fill with model data
vstock=NaN(360,720,nyear,nmod); %Vegetation carbon stock
mflux=NaN(360,720,nyear,nmod); %Annual turnover flux through all mortality causes
distflux=NaN(360,720,nyear,nmod); %Annual turnover flux from disturbance mortality (including fire)
vitalflux=NaN(360,720,nyear,nmod); %Annual turnover flux from vitality mortality
backflux=NaN(360,720,nyear,nmod); %Annual turnover flux from background mortality
heatflux=NaN(360,720,nyear,nmod); %Annual turnover flux from heat stress mortality
otherflux=NaN(360,720,nyear,nmod); %Annual turnover flux from other causes of mortality
nppflux=NaN(360,720,nyear,nmod); %Annual NPP
reproflux=NaN(360,720,nyear,nmod); %Annual turnover flux through reproduction

vstock_jules=NaN(144,192,nyear); %Vegetation carbon stock (JULES only)
mflux_jules=NaN(144,192,nyear); %Annual turnover flux through all mortality causes (JULES only)
backflux_jules=NaN(144,192,nyear); %Annual turnover flux from background mortality (JULES only)
vitalflux_jules=NaN(144,192,nyear); %Annual turnover flux from vitality mortality (JULES only)
nppflux_jules=NaN(144,192,nyear); %Annual NPP (JULES only)

%Load the model data

%CABLE-POP
load([data_models,'/cable_map_arrays_ipsl.mat'])
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_cable=tdist_cable+tgroweff_cable+tother_cable;
vstock(:,:,:,1)=flip(squeeze(permute(cveg_cable(:,:,y1:y2),[2 1 3])),1);
mflux(:,:,:,1)=flip(permute(turn_cable(:,:,y1:y2),[2 1 3]).*uconv);
distflux(:,:,:,1)=flip(permute(tdist_cable(:,:,y1:y2),[2 1 3]).*uconv);
vitalflux(:,:,:,1)=flip(permute(tgroweff_cable(:,:,y1:y2)+tother_cable(:,:,y1:y2),[2 1 3]).*uconv);
nppflux(:,:,:,1)=flip(permute(npp_cable(:,:,y1:y2),[2 1 3]).*uconv);
%No reproduction flux
clear turn_cable
clear cleaf_cable croot_cable csoil_cable cveg_cable cwood_cable gpp_cable npp_cable tdist_cable tgroweff_cable tleaf_cable tother_cable troot_cable

%JULES
load([data_models,'/jules_map_arrays_ipsl.mat'])
uconv=sec_per_day*days_per_year_jules; %Conversion factor for fluxes from seconds to years
turn_jules=tback_jules+tother_jules;
vstock_jules(29:140,1:96,:)=permute(cveg_jules(97:192,:,y1:y2),[2 1 3]);
mflux_jules(29:140,1:96,:)=permute(turn_jules(97:192,:,y1:y2),[2 1 3]).*uconv;
backflux_jules(29:140,1:96,:)=permute(tback_jules(97:192,:,y1:y2),[2 1 3]).*uconv;
vitalflux_jules(29:140,1:96,:)=permute(tother_jules(97:192,:,y1:y2),[2 1 3]).*uconv;
nppflux_jules(29:140,1:96,:)=permute(npp_jules(97:192,:,y1:y2),[2 1 3]).*uconv;

vstock_jules(29:140,97:192,:)=permute(cveg_jules(1:96,:,y1:y2),[2 1 3]);
mflux_jules(29:140,97:192,:)=permute(turn_jules(1:96,:,y1:y2),[2 1 3]).*uconv;
backflux_jules(29:140,97:192,:)=permute(tback_jules(1:96,:,y1:y2),[2 1 3]).*uconv;
vitalflux_jules(29:140,97:192,:)=permute(tother_jules(1:96,:,y1:y2),[2 1 3]).*uconv;
nppflux_jules(29:140,97:192,:)=permute(npp_jules(1:96,:,y1:y2),[2 1 3]).*uconv;
%No reproduction flux
clear turn_jules
clear cleaf_jules croot_jules csoil_jules cveg_jules cwood_jules gpp_jules npp_jules tdist_jules tgroweff_jules tleaf_jules tother_jules troot_jules

%LPJ-GUESS
load([data_models,'/lpjg_map_arrays_ipsl.mat'])
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_lpjg=tdist_lpjg+tfire_lpjg+tgroweff_lpjg+tother_lpjg;
vstock(:,:,:,3)=permute(cveg_lpjg(:,:,y1:y2),[2 1 3]);
mflux(:,:,:,3)=permute(turn_lpjg(:,:,y1:y2),[2 1 3]).*uconv;
distflux(:,:,:,3)=permute(tdist_lpjg(:,:,y1:y2)+tfire_lpjg(:,:,y1:y2),[2 1 3]).*uconv;
vitalflux(:,:,:,3)=permute(tgroweff_lpjg(:,:,y1:y2),[2 1 3]).*uconv;
backflux(:,:,:,3)=permute(tback_lpjg(:,:,y1:y2),[2 1 3]).*uconv;
otherflux(:,:,:,3)=permute(tother_lpjg(:,:,y1:y2)-tback_lpjg(:,:,y1:y2),[2 1 3]).*uconv;
nppflux(:,:,:,3)=permute(npp_lpjg(:,:,y1:y2),[2 1 3]).*uconv;
clear turn_lpjg
clear cleaf_lpjg croot_lpjg csoil_lpjg cveg_lpjg cwood_lpjg gpp_lpjg npp_lpjg tdist_lpjg tgroweff_lpjg tleaf_lpjg tother_lpjg troot_lpjg tnbio_lpjg tfire_lpjg
clear tallom_lpjg tbioclim_lpjg tback_lpjg

%LPJmL
load([data_models,'/lpjml_map_arrays_ipsl.mat'])
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_lpjml=tfire_lpjml+tgroweff_lpjml+tother_lpjml;
vstock(:,:,:,4)=permute(cveg_lpjml(:,:,y1:y2),[2 1 3]);
mflux(:,:,:,4)=permute(turn_lpjml(:,:,y1:y2),[2 1 3]).*uconv;
distflux(:,:,:,4)=permute(tfire_lpjml(:,:,y1:y2),[2 1 3]).*uconv;
vitalflux(:,:,:,4)=permute(tgroweff_lpjml(:,:,y1:y2)+tshade_lpjml(:,:,y1:y2),[2 1 3]).*uconv;
heatflux(:,:,:,4)=permute(theat_lpjml(:,:,y1:y2),[2 1 3]).*uconv;
otherflux(:,:,:,4)=permute(tother_lpjml(:,:,y1:y2)-theat_lpjml(:,:,y1:y2)-tshade_lpjml(:,:,y1:y2),[2 1 3]).*uconv;
nppflux(:,:,:,4)=permute(npp_lpjml(:,:,y1:y2),[2 1 3]).*uconv;
clear turn_lpjml
clear cleaf_lpjml croot_lpjml csoil_lpjml cveg_lpjml cwood_lpjml gpp_lpjml npp_lpjml tdist_lpjml tgroweff_lpjml tleaf_lpjml tother_lpjml troot_lpjml
clear tallom_lpjml tfire_lpjml theat_lpjml tnbio_lpjml tshade_lpjml

%ORCHIDEE
load([data_models,'/orchidee_map_arrays_ipsl.mat'])
uconv=days_per_year/g_per_kg; %Conversion factor for fluxes from days to years and from g to kg.
uconv_npp=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_orchidee=tmort_orchidee;
vstock(:,:,:,5)=flip(permute(cveg_orchidee(:,:,y1:y2),[2 1 3]),1);
mflux(:,:,:,5)=flip(permute(turn_orchidee(:,:,y1:y2),[2 1 3]).*uconv,1);
nppflux(:,:,:,5)=flip(permute(npp_orchidee(:,:,y1:y2),[2 1 3]).*uconv_npp,1);
clear turn_orchidee
clear cleaf_orchidee croot_orchidee csoil_orchidee cveg_orchidee cwood_orchidee gpp_orchidee npp_orchidee tdist_orchidee tgroweff_orchidee tleaf_orchidee tother_orchidee troot_orchidee
clear tmort_orchidee

%SEIB-DGVM
load([data_models,'/seib_map_arrays_ipsl.mat'])
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_seib=tdist_seib+tfire_seib+tgroweff_seib+tother_seib;
vstock(:,:,:,6)=flipud(permute(cveg_seib(:,:,y1:y2),[2 1 3]));
mflux(:,:,:,6)=flipud(permute(turn_seib(:,:,y1:y2),[2 1 3]).*uconv);
distflux(:,:,:,6)=flipud(permute(tdist_seib(:,:,y1:y2)+tfire_seib(:,:,y1:y2),[2 1 3]).*uconv);
vitalflux(:,:,:,6)=flipud(permute(tgroweff_seib(:,:,y1:y2),[2 1 3]).*uconv);
heatflux(:,:,:,6)=flipud(permute(theat_seib(:,:,y1:y2),[2 1 3]).*uconv);
backflux(:,:,:,6)=flipud(permute(tback_seib(:,:,y1:y2),[2 1 3]).*uconv);
otherflux(:,:,:,6)=flipud(permute(tother_seib(:,:,y1:y2)-theat_seib(:,:,y1:y2)-tback_seib(:,:,y1:y2),[2 1 3]).*uconv);
nppflux(:,:,:,6)=flipud(permute(npp_seib(:,:,y1:y2),[2 1 3]).*uconv);
clear turn_seib
clear cleaf_seib croot_seib csoil_seib cveg_seib cwood_seib gpp_seib npp_seib tdist_seib tgroweff_seib tleaf_seib tother_seib troot_seib
clear tbacl_seib tbioclim_seib tfire_seib theat_seib

clear uconv uconv_npp

%---
%Make biome-level statistics

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[~,~,ffrac,ffrac_jules]=get_closed_can_mask(data_mask);

%Get the grid cell areas for 0.5 x 0.5 degree grid cells
garea=global_grid_area();
%Get the grid cell areas for JULES grid cells from a file
garea_jules_in=ncread(jules_gridcellarea_file,'arealand');
garea_jules=NaN(size(garea_jules_in));
garea_jules(1:96,:)=garea_jules_in(97:192,:);
garea_jules(97:192,:)=garea_jules_in(1:96,:);
garea_jules=garea_jules';
garea_jules(1,:)=[];

%Total forest area per grid cell (m2)
ffrac=ffrac.*garea';
ffrac_jules=ffrac_jules.*garea_jules';

%Multiple data by forest area (10% min. cover threshold per grid cell)
mflux_mask=mflux.*repmat(ffrac',1,1,nyear,nmod); clear mflux
distflux_mask=distflux.*repmat(ffrac',1,1,nyear,nmod); clear distflux
vitalflux_mask=vitalflux.*repmat(ffrac',1,1,nyear,nmod); clear vitalflux
backflux_mask=backflux.*repmat(ffrac',1,1,nyear,nmod); clear backflux
heatflux_mask=heatflux.*repmat(ffrac',1,1,nyear,nmod); clear heatflux
otherflux_mask=otherflux.*repmat(ffrac',1,1,nyear,nmod); clear otherflux
vstock_mask=vstock.*repmat(ffrac',1,1,nyear,nmod); clear vstock
nppflux_mask=nppflux.*repmat(ffrac',1,1,nyear,nmod); clear nppflux

mflux_jules_mask=mflux_jules.*repmat(ffrac_jules',1,1,nyear); clear mflux_jules
backflux_jules_mask=backflux_jules.*repmat(ffrac_jules',1,1,nyear); clear backflux_jules
vitalflux_jules_mask=vitalflux_jules.*repmat(ffrac_jules',1,1,nyear); clear vitalflux_jules
vstock_jules_mask=vstock_jules.*repmat(ffrac_jules',1,1,nyear); clear vstock_jules
nppflux_jules_mask=nppflux_jules.*repmat(ffrac_jules',1,1,nyear); clear nppflux_jules

%Catch any Inf values
mflux_mask(isinf(mflux_mask))=NaN;
distflux_mask(isinf(distflux_mask))=NaN;
vitalflux_mask(isinf(vitalflux_mask))=NaN;
backflux_mask(isinf(backflux_mask))=NaN;
heatflux_mask(isinf(heatflux_mask))=NaN;
otherflux_mask(isinf(otherflux_mask))=NaN;
vstock_mask(isinf(vstock_mask))=NaN;
mflux_jules_mask(isinf(mflux_jules_mask))=NaN;
backflux_jules_mask(isinf(backflux_jules_mask))=NaN;
vitalflux_jules_mask(isinf(vitalflux_jules_mask))=NaN;
vstock_jules_mask(isinf(vstock_jules_mask))=NaN;
nppflux_jules_mask(isinf(nppflux_jules_mask))=NaN;

%---
%Calculate global totals, normalised to total vegetation stock to give turnover (or growth) rates
mflux_mask_mean=NaN(nmod,nyear);
distflux_mask_mean=NaN(nmod,nyear);
vitalflux_mask_mean=NaN(nmod,nyear);
backflux_mask_mean=NaN(nmod,nyear);
heatflux_mask_mean=NaN(nmod,nyear);
otherflux_mask_mean=NaN(nmod,nyear);
nppflux_mask_mean=NaN(nmod,nyear);
for yy=1:nyear
    for nn=1:nmod
        if nn==2
            vstock_temp=vstock_jules_mask(:,:,yy);
            mflux_temp=mflux_jules_mask(:,:,yy); mflux_mask_mean(2,yy)=nansum(mflux_temp(:))/nansum(vstock_temp(:));
            backflux_temp=backflux_jules_mask(:,:,yy); backflux_mask_mean(2,yy)=nansum(backflux_temp(:))/nansum(vstock_temp(:));
            vitalflux_temp=vitalflux_jules_mask(:,:,yy); vitalflux_mask_mean(2,yy)=nansum(vitalflux_temp(:))/nansum(vstock_temp(:));
            nppflux_temp=nppflux_jules_mask(:,:,yy); nppflux_mask_mean(2,yy)=nansum(nppflux_temp(:))/nansum(vstock_temp(:));
            continue
        end
        vstock_temp=vstock_mask(:,:,yy,nn);
        mflux_temp=mflux_mask(:,:,yy,nn); mflux_mask_mean(nn,yy)=nansum(mflux_temp(:))/nansum(vstock_temp(:));
        distflux_temp=distflux_mask(:,:,yy,nn); distflux_mask_mean(nn,yy)=nansum(distflux_temp(:))/nansum(vstock_temp(:));
        vitalflux_temp=vitalflux_mask(:,:,yy,nn); vitalflux_mask_mean(nn,yy)=nansum(vitalflux_temp(:))/nansum(vstock_temp(:));
        backflux_temp=backflux_mask(:,:,yy,nn); backflux_mask_mean(nn,yy)=nansum(backflux_temp(:))/nansum(vstock_temp(:));
        heatflux_temp=heatflux_mask(:,:,yy,nn); heatflux_mask_mean(nn,yy)=nansum(heatflux_temp(:))/nansum(vstock_temp(:));
        otherflux_temp=otherflux_mask(:,:,yy,nn); otherflux_mask_mean(nn,yy)=nansum(otherflux_temp(:))/nansum(vstock_temp(:));
        nppflux_temp=nppflux_mask(:,:,yy,nn); nppflux_mask_mean(nn,yy)=nansum(nppflux_temp(:))/nansum(vstock_temp(:));
        clear mflux_temp fireflux_temp distflux_temp vitalflux_temp backflux_temp heatflux_temp otherflux_temp lflux_temp rflux_temp vstock_temp
    end
end
clear nn yy

%---
%Calculate forest type values (observed forest-types)

%Read forest type data based on ESA landcover (mask produced in Pugh et al., 2019, PNAS 116(19), 4382-4387 
%and available in the SI of that paper as an ACSII file)
%Produced here using esa_forest_9regions_new_func.m
phen_dir='/data/turnover/masks/phen/';
phen=NaN(360,720,nmod);
ESA_phen=ncread([phen_dir,'ESA_forest_9regions_v2.nc'],'region_mask')';
ESA_phen_label={'(a) Trop. broad. ever.','(b) Trop. broad. dec.','(c) Other trop.','(d) Temp. broad. ever.'...
    ,'(e) Temp. broad. dec.','(f) Needle. ever.','(g) Needle. dec.','(h) Mixed forest','(i) Other'};
nESAphen=length(ESA_phen_label);

ESA_phen_julesa=ncread([phen_dir,'ESA_forest_9regions_v2_julesgrid.nc'],'region_mask')';
ESA_phen_jules=NaN(144,192);
ESA_phen_jules(29:140,:)=ESA_phen_julesa;
clear phen_julesa

%Now make the forest-type-level totals, normalised to total vegetation stock to give turnover (or growth) rates
mflux_mask_mean=NaN(nmod,nyear,nESAphen);
distflux_mask_mean=NaN(nmod,nyear,nESAphen);
vitalflux_mask_mean=NaN(nmod,nyear,nESAphen);
backflux_mask_mean=NaN(nmod,nyear,nESAphen);
heatflux_mask_mean=NaN(nmod,nyear,nESAphen);
otherflux_mask_mean=NaN(nmod,nyear,nESAphen);
nppflux_mask_mean=NaN(nmod,nyear,nESAphen);
nppnotnorm_mask_mean=NaN(nmod,nyear,nESAphen);
for yy=1:nyear
for bb=1:nESAphen
    for nn=1:nmod
        if nn==2
            vstock_temp=vstock_jules_mask(:,:,yy);
            mflux_temp=mflux_jules_mask(:,:,yy); mflux_mask_mean(2,yy,bb)=nansum(mflux_temp(ESA_phen_jules==bb))/nansum(vstock_temp(ESA_phen_jules==bb));
            backflux_temp=backflux_jules_mask(:,:,yy); backflux_mask_mean(2,yy,bb)=nansum(backflux_temp(ESA_phen_jules==bb))/nansum(vstock_temp(ESA_phen_jules==bb));
            vitalflux_temp=vitalflux_jules_mask(:,:,yy); vitalflux_mask_mean(2,yy,bb)=nansum(vitalflux_temp(ESA_phen_jules==bb))/nansum(vstock_temp(ESA_phen_jules==bb));
            nppflux_temp=nppflux_jules_mask(:,:,yy); nppflux_mask_mean(2,yy,bb)=nansum(nppflux_temp(ESA_phen_jules==bb))/nansum(vstock_temp(ESA_phen_jules==bb));
            nppnotnorm_temp=nppflux_jules_mask(:,:,yy); nppnotnorm_mask_mean(2,yy,bb)=nansum(nppflux_temp(ESA_phen_jules==bb));
            continue
        end
        vstock_temp=vstock_mask(:,:,yy,nn);
        mflux_temp=mflux_mask(:,:,yy,nn); mflux_mask_mean(nn,yy,bb)=nansum(mflux_temp(ESA_phen==bb))/nansum(vstock_temp(ESA_phen==bb));
        distflux_temp=distflux_mask(:,:,yy,nn); distflux_mask_mean(nn,yy,bb)=nansum(distflux_temp(ESA_phen==bb))/nansum(vstock_temp(ESA_phen==bb));
        vitalflux_temp=vitalflux_mask(:,:,yy,nn); vitalflux_mask_mean(nn,yy,bb)=nansum(vitalflux_temp(ESA_phen==bb))/nansum(vstock_temp(ESA_phen==bb));
        backflux_temp=backflux_mask(:,:,yy,nn); backflux_mask_mean(nn,yy,bb)=nansum(backflux_temp(ESA_phen==bb))/nansum(vstock_temp(ESA_phen==bb));
        heatflux_temp=heatflux_mask(:,:,yy,nn); heatflux_mask_mean(nn,yy,bb)=nansum(heatflux_temp(ESA_phen==bb))/nansum(vstock_temp(ESA_phen==bb));
        otherflux_temp=otherflux_mask(:,:,yy,nn); otherflux_mask_mean(nn,yy,bb)=nansum(otherflux_temp(ESA_phen==bb))/nansum(vstock_temp(ESA_phen==bb));
        nppflux_temp=nppflux_mask(:,:,yy,nn); nppflux_mask_mean(nn,yy,bb)=nansum(nppflux_temp(ESA_phen==bb))/nansum(vstock_temp(ESA_phen==bb));
        nppnotnorm_temp=nppflux_mask(:,:,yy,nn); nppnotnorm_mask_mean(nn,yy,bb)=nansum(nppflux_temp(ESA_phen==bb));
        clear mflux_temp fireflux_temp distflux_temp vitalflux_temp backflux_temp heatflux_temp otherflux_temp lflux_temp rflux_temp shadeflux_temp nppflux_temp nppnotnorm_temp
    end
end
end
clear nn bb yy


%---
%Make plots
    
%Assign the ORCHIDEE mortality fluxes all to background for plotting purposes
backflux_mask_mean(5,:,:)=mflux_mask_mean(5,:,:);

%Make 31-year running means of the rates of turnover by each of the mortality mechanisms
for nn=1:nmod %Make one plot per model
    figure
    colormap(lines(5))
    for bb=1:nESAphen
        s(bb)=subplot(3,3,bb);
        hold on
        firstyear=1900+y1+15; %Only plot from 15 years into the dataset because using running mean
        lastyear=1900+y2-15;
        
        %Temporary running mean variables
        distflux_runmean=runmean(distflux_mask_mean(nn,:,bb),15);
        vitalflux_runmean=runmean(vitalflux_mask_mean(nn,:,bb),15);
        backflux_runmean=runmean(backflux_mask_mean(nn,:,bb),15);
        heatflux_runmean=runmean(heatflux_mask_mean(nn,:,bb),15);
        otherflux_runmean=runmean(otherflux_mask_mean(nn,:,bb),15);
        
        plot(firstyear:lastyear,distflux_runmean(16:100),'linewidth',2)
        plot(firstyear:lastyear,vitalflux_runmean(16:100),'linewidth',2)
        plot(firstyear:lastyear,backflux_runmean(16:100),'linewidth',2)
        plot(firstyear:lastyear,heatflux_runmean(16:100),'linewidth',2)
        plot(firstyear:lastyear,otherflux_runmean(16:100),'linewidth',2)
        clear distflux_runmean vitalflux_runmean backflux_runmean heatflux_runmean otherflux_runmean
        
        title(ESA_phen_label{bb},'Fontsize',9)
        set(gca,'XLim',[1900+y1+15,1900+y2-15])
        if bb==nESAphen
            legend('Disturbance','Vitality','Background','Heat','Other')
        end
        if bb==1 || bb==4 || bb==7
            ylabel('a^{-1}')
        end
    end
    set(s(1),'Position',[0.07 0.67 0.25 0.23])
    set(s(2),'Position',[0.37 0.67 0.25 0.23])
    set(s(3),'Position',[0.67 0.67 0.25 0.23])
    set(s(4),'Position',[0.07 0.37 0.25 0.23])
    set(s(5),'Position',[0.37 0.37 0.25 0.23])
    set(s(6),'Position',[0.67 0.37 0.25 0.23])
    set(s(7),'Position',[0.07 0.07 0.25 0.23])
    set(s(8),'Position',[0.37 0.07 0.25 0.23])
    set(s(9),'Position',[0.67 0.07 0.25 0.23])
    
    %NPP plots (NOT normalised by Cveg)
    figure
    cmap=[16, 100, 112;... %CABLE
        182, 0, 38;... %JULES
        85, 165, 28;... %LPJ-GUESS
        234, 113, 37;... %LPJmL
        163, 193, 173;... %ORCHIDEE
        104, 172, 229]; %SEIB
    cmap=cmap/256;
    for bb=1:nESAphen
        s(bb)=subplot(3,3,bb);
        hold on
        for mm=1:nmod
            plot(1900+y1:1900+y2,runmean(nppnotnorm_mask_mean(mm,:,bb),15)/1e12,'linewidth',2,'color',cmap(mm,:))
        end
        title(ESA_phen_label{bb},'Fontsize',9)
        set(gca,'XLim',[1900+y1,1900+y2])
        if bb==nESAphen
            legend(modelsonly)
        end
        if bb==1 || bb==4 || bb==7
            ylabel('Pg C a^{-1}')
        end
    end
    set(s(1),'Position',[0.07 0.67 0.25 0.23])
    set(s(2),'Position',[0.37 0.67 0.25 0.23])
    set(s(3),'Position',[0.67 0.67 0.25 0.23])
    set(s(4),'Position',[0.07 0.37 0.25 0.23])
    set(s(5),'Position',[0.37 0.37 0.25 0.23])
    set(s(6),'Position',[0.67 0.37 0.25 0.23])
    set(s(7),'Position',[0.07 0.07 0.25 0.23])
    set(s(8),'Position',[0.37 0.07 0.25 0.23])
    set(s(9),'Position',[0.67 0.07 0.25 0.23])
    clear bb mm
end
clear nn bb

