function [vstock,mflux,lrrflux,nppflux,reproflux,vstock_jules,mflux_jules,lrrflux_jules,...
    nppflux_jules,models,nmod]=get_stocks_fluxes_runmeanstd(data_models,ifcruncep,y1,y2)
%Function to read in *.mat files of preprocessed model data from turnover_pool_flux_read.m
%and extract the standard deviation of the 31-year running mean of fluxes and vegetation 
%stocks for a specified time window.
%
%T. Pugh
%27.11.19

models={'(a) CABLE-POP','(b) JULES','(c) LPJ-GUESS','(d) LPJmL','(e) ORCHIDEE','(f) SEIB-DGVM'};
nmod=length(models);

%Set some constants
sec_per_day=86400;
g_per_kg=1000;
days_per_year=365;
days_per_year_jules=360; %Only 360 days per year in JULES


%Initialise arrays to fill with model data
vstock=NaN(360,720,nmod); %Vegetation carbon stock
mflux=NaN(360,720,nmod); %Annual turnover flux through all mortality causes
lrrflux=NaN(360,720,nmod); %Annual turnover flux through leaves, fine roots and reproduction
nppflux=NaN(360,720,nmod); %Annual NPP
reproflux=NaN(360,720,nmod); %Annual turnover flux through reproduction
vstock_jules=NaN(144,192); %Vegetation carbon stock (JULES only)
mflux_jules=NaN(144,192); %Annual turnover flux through all mortality causes (JULES only)
lrrflux_jules=NaN(144,192); %Annual turnover flux through leaves, fine roots and reproduction (JULES only)
nppflux_jules=NaN(144,192); %Annual NPP (JULES only)

%CABLE-POP
if ifcruncep
    load([data_models,'/cable_map_arrays_cruncep.mat'])
else
    load([data_models,'/cable_map_arrays_ipsl.mat'])
end
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_cable=tdist_cable+tgroweff_cable+tother_cable;
vstock(:,:,1)=flipud(squeeze(std(cveg_cable(:,:,y1:y2),[],3))');
mflux(:,:,1)=flipud(std(runmean(turn_cable(:,:,y1:y2),15,3),[],3)'.*uconv);
nppflux(:,:,1)=flipud(std(npp_cable(:,:,y1:y2),[],3)'.*uconv);
lrrflux(:,:,1)=flipud(std(runmean(tleaf_cable(:,:,y1:y2)+troot_cable(:,:,y1:y2),15,3),[],3)').*uconv;
%No reproduction flux
clear turn_cable
clear cleaf_cable croot_cable csoil_cable cveg_cable cwood_cable gpp_cable npp_cable tdist_cable tgroweff_cable tleaf_cable tother_cable troot_cable

%JULES
if ifcruncep
    load([data_models,'/jules_map_arrays_cruncep.mat'])
else
    load([data_models,'/jules_map_arrays_ipsl.mat'])
end
uconv=sec_per_day*days_per_year_jules; %Conversion factor for fluxes from seconds to years
turn_jules=tback_jules+tother_jules;
vstock_jules(29:140,1:96)=squeeze(std(cveg_jules(97:192,:,y1:y2),[],3))';
mflux_jules(29:140,1:96)=std(runmean(turn_jules(97:192,:,y1:y2),15,3),[],3)'.*uconv;
lrrflux_jules(29:140,1:96)=std(runmean(tleaf_jules(97:192,:,y1:y2)+troot_jules(97:192,:,y1:y2),15,3),[],3)'.*uconv;
nppflux_jules(29:140,1:96)=std(npp_jules(97:192,:,y1:y2),[],3)'.*uconv;
vstock_jules(29:140,97:192)=squeeze(std(cveg_jules(1:96,:,y1:y2),[],3))';
mflux_jules(29:140,97:192)=std(runmean(turn_jules(1:96,:,y1:y2),15,3),[],3)'.*uconv;
lrrflux_jules(29:140,97:192)=std(runmean(tleaf_jules(1:96,:,y1:y2)+troot_jules(1:96,:,y1:y2),15,3),[],3)'.*uconv;
nppflux_jules(29:140,97:192)=std(npp_jules(1:96,:,y1:y2),[],3)'.*uconv;
%No reproduction flux
clear turn_jules
clear cleaf_jules croot_jules csoil_jules cveg_jules cwood_jules gpp_jules npp_jules tdist_jules tgroweff_jules tleaf_jules tother_jules troot_jules

%LPJ-GUESS
if ifcruncep
    load([data_models,'/lpjg_map_arrays_cruncep.mat'])
else
    load([data_models,'/lpjg_map_arrays_ipsl.mat'])
end
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_lpjg=tdist_lpjg+tfire_lpjg+tgroweff_lpjg+tother_lpjg;
vstock(:,:,3)=squeeze(std(cveg_lpjg(:,:,y1:y2),[],3))';
mflux(:,:,3)=std(runmean(turn_lpjg(:,:,y1:y2),15,3),[],3)'.*uconv;
nppflux(:,:,3)=std(npp_lpjg(:,:,y1:y2),[],3)'.*uconv;
reproflux(:,:,3)=nppflux(:,:,3)*0.1;
lrrflux(:,:,3)=std(runmean(tleaf_lpjg(:,:,y1:y2)+troot_lpjg(:,:,y1:y2)+(0.1*npp_lpjg(:,:,y1:y2)),15,3),[],3)'.*uconv;
clear turn_lpjg
clear cleaf_lpjg croot_lpjg csoil_lpjg cveg_lpjg cwood_lpjg gpp_lpjg npp_lpjg tdist_lpjg tgroweff_lpjg tleaf_lpjg tother_lpjg troot_lpjg

%LPJmL
if ifcruncep
    load([data_models,'/lpjml_map_arrays_cruncep.mat'])
else
    load([data_models,'/lpjml_map_arrays_ipsl.mat'])
end
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_lpjml=tfire_lpjml+tgroweff_lpjml+tother_lpjml;
vstock(:,:,4)=squeeze(std(cveg_lpjml(:,:,y1:y2),[],3))';
mflux(:,:,4)=std(runmean(turn_lpjml(:,:,y1:y2),15,3),[],3)'.*uconv;
nppflux(:,:,4)=std(npp_lpjml(:,:,y1:y2),[],3)'.*uconv;
reproflux(:,:,4)=nppflux(:,:,4)*0.1;
lrrflux(:,:,4)=std(runmean(tleaf_lpjml(:,:,y1:y2)+troot_lpjml(:,:,y1:y2)+(0.1*npp_lpjml(:,:,y1:y2)),15,3),[],3)'.*uconv;
clear turn_lpjml
clear cleaf_lpjml croot_lpjml csoil_lpjml cveg_lpjml cwood_lpjml gpp_lpjml npp_lpjml tdist_lpjml tgroweff_lpjml tleaf_lpjml tother_lpjml troot_lpjml

%ORCHIDEE
if ifcruncep
    load([data_models,'/orchidee_map_arrays_cruncep.mat'])
else
    load([data_models,'/orchidee_map_arrays_ipsl.mat'])
end
uconv=days_per_year/g_per_kg; %Conversion factor for fluxes from days to years and from g to kg.
uconv_npp=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_orchidee=tmort_orchidee;
vstock(:,:,5)=flipud(squeeze(std(cveg_orchidee(:,:,y1:y2),[],3))');
mflux(:,:,5)=flipud(std(runmean(turn_orchidee(:,:,y1:y2),15,3),[],3)'.*uconv);
nppflux(:,:,5)=flipud(std(npp_orchidee(:,:,y1:y2),[],3)'.*uconv_npp);
reproflux(:,:,5)=nppflux(:,:,5)*0.1;
lrrflux(:,:,5)=flipud(std(runmean(tleaf_orchidee(:,:,y1:y2)+troot_orchidee(:,:,y1:y2)+(0.1*npp_orchidee(:,:,y1:y2)),15,3),[],3)'.*uconv);
clear turn_orchidee
clear cleaf_orchidee croot_orchidee csoil_orchidee cveg_orchidee cwood_orchidee gpp_orchidee npp_orchidee tdist_orchidee tgroweff_orchidee tleaf_orchidee tother_orchidee troot_orchidee

%SEIB-DGVM
if ifcruncep
    load([data_models,'/seib_map_arrays_cruncep.mat'])
else
    load([data_models,'/seib_map_arrays_ipsl.mat'])
end
uconv=sec_per_day*days_per_year; %Conversion factor for fluxes from seconds to years
turn_seib=tdist_seib+tfire_seib+tgroweff_seib+tother_seib;
vstock(:,:,6)=flipud(squeeze(std(cveg_seib(:,:,y1:y2),[],3))');
mflux(:,:,6)=flipud(std(runmean(turn_seib(:,:,y1:y2),15,3),[],3)'.*uconv);
nppflux(:,:,6)=flipud(std(npp_seib(:,:,y1:y2),[],3)'.*uconv);
reproflux(:,:,6)=nppflux(:,:,6)*0.1; %Repro flux for SEIB is an approximation
lrrflux(:,:,6)=flipud(std(runmean(tleaf_seib(:,:,y1:y2)+troot_seib(:,:,y1:y2)+(0.1*npp_seib(:,:,y1:y2)),15,3),[],3)'.*uconv);
clear turn_seib
clear cleaf_seib croot_seib csoil_seib cveg_seib cwood_seib gpp_seib npp_seib tdist_seib tgroweff_seib tleaf_seib tother_seib troot_seib

clear uconv uconv_npp y1 y2