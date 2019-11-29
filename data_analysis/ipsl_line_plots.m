%Script to make line plots of change in turnover time, mort fractions and Cveg relative to a defined
%reference period, based on the IPSL simulations.
%
%T. Pugh
%24.02.18

data_models='/media/pughtam/rds-2017-pughtam-01/turnover/turnover_comp/data_analysis/'; %Location of *.mat files containing preprocessed model data

%Index of first and last year to average over for the baseline period
y1=85; %1985
y2=114; %2014
nyear=199; %Total number of years in simulation

%---
%Load the globally averaged IPSL data from *.mat file created by turnover_pool_flux_read.m
load([data_models,'all_ipsl_global_arrays']);

models={'CABLE-POP','JULES','LPJ-GUESS','LPJmL','ORCHIDEE','SEIB-DGVM'};
nmod=length(models);

ordind=[5,4,1,2,6,3]; %Order in which to display models (because output from all_ipsl_global_arrays.mat is not in alphabetical order)

%Make a custom colourmap
cmap=[16, 100, 112;... %CABLE
 182, 0, 38;... %JULES
 85, 165, 28;... %LPJ-GUESS
 234, 113, 37;... %LPJmL 
 163, 193, 173;... %ORCHIDEE
 104, 172, 229]; %SEIB
cmap=cmap/256;

%Make reproduction flux array (0.1*NPP for LPJG, LPJ, LPJmL, SEIB and ORCHIDEE, zero for JULES AND CABLE)
trepro=zeros(size(npp));
trepro(:,1:3)=npp(:,1:3)*0.1;
trepro(:,6)=npp(:,6)*0.1;

%Calculate variables to display
tau_npp=cveg./npp;
tau_turn=cveg./(tmort+tleaf+troot+trepro); %Total turnover time defined by output flux
tau_mort=cveg./tmort; %Total turnover time due to mortality
frac_mort=tmort./(tmort+tleaf+troot+trepro); %Fraction of total turnover due to mortality
npp_norm=npp./cveg; %NPP rate per unit vegetation C
tmort_norm=tmort./cveg; %Mortality rate per unit vegetation C

%Calculate mean values for baseline period
tau_npp_mean=mean(tau_npp(y1:y2,:),1);
tau_turn_mean=mean(tau_turn(y1:y2,:),1);
tau_mort_mean=mean(tau_mort(y1:y2,:),1);
frac_mort_mean=mean(frac_mort(y1:y2,:),1);
cveg_mean=mean(cveg(y1:y2,:),1);
tmort_mean=mean(tmort(y1:y2,:),1);
npp_mean=mean(npp(y1:y2,:),1);
tmort_norm_mean=mean(tmort_norm(y1:y2,:),1);
npp_norm_mean=mean(npp_norm(y1:y2,:),1);

%Adjust order of models in array to match ordind
cveg=cveg(:,ordind);
npp=npp(:,ordind);
tau_npp=tau_npp(:,ordind);
tau_turn=tau_turn(:,ordind);
tau_mort=tau_mort(:,ordind);
frac_mort=frac_mort(:,ordind);
npp_norm=npp_norm(:,ordind);
tmort_norm=tmort_norm(:,ordind);

tau_npp_mean=tau_npp_mean(:,ordind);
tau_turn_mean=tau_turn_mean(:,ordind);
tau_mort_mean=tau_mort_mean(:,ordind);
frac_mort_mean=frac_mort_mean(:,ordind);
cveg_mean=cveg_mean(:,ordind);
tmort_mean=tmort_mean(:,ordind);
npp_mean=npp_mean(:,ordind);
tmort_norm_mean=tmort_norm_mean(:,ordind);
npp_norm_mean=npp_norm_mean(:,ordind);

%---
%Make plots

%Change in biomass and tunrover times relative to baseline (10-y running mean)
figure
subplot(4,1,1)
hold on
for nn=1:nmod
    plot(1901:2099,runmean(((cveg(:,nn)./repmat(cveg_mean(:,nn),nyear,1))*100)-100,5),'linewidth',2,'color',cmap(nn,:))
end
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('%')
title('(a) \Delta Biomass')
grid on
set(gca,'XTickLabel','')

subplot(4,1,2)
hold on
for nn=1:nmod
    plot(1901:2099,runmean(((tau_turn(:,nn)./repmat(tau_turn_mean(:,nn),nyear,1))*100)-100,5),'linewidth',2,'color',cmap(nn,:))
end
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('%')
title('(b) \Delta \tau_{turn}')
grid on
set(gca,'XTickLabel','')

subplot(4,1,3)
hold on
for nn=1:nmod
    plot(1901:2099,runmean(((tau_mort(:,nn)./repmat(tau_mort_mean(:,nn),nyear,1))*100)-100,5),'linewidth',2,'color',cmap(nn,:))
end
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('%')
title('(c) \Delta \tau_{mort}')
grid on
set(gca,'XTickLabel','')

subplot(4,1,4)
hold on
for nn=1:nmod
    plot(1901:2099,runmean(((frac_mort(:,nn)./repmat(frac_mort_mean(:,nn),nyear,1))*100)-100,5),'linewidth',2,'color',cmap(nn,:))
end
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('%')
title('(d) Delta fraction as mort.')
grid on

legend(models)

%Change in fluxes relative to baseline (10-y running mean)
figure

subplot(3,1,1)
hold on
for nn=1:nmod
    npp_change=runmean(((npp(:,nn)./repmat(npp_mean(:,nn),nyear,1))*100)-100,5);
    plot(1901:2099,npp_change,'linewidth',2,'color',cmap(nn,:))
end
clear nn npp_change tmort_change
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('\Delta %')
grid on
set(gca,'XTickLabel','')
title('(a) NPP')

subplot(3,1,2)
hold on
for nn=1:nmod
    tmort_change=runmean(((tmort(:,nn)./repmat(tmort_mean(:,nn),nyear,1))*100)-100,5);
    plot(1901:2099,tmort_change,'linewidth',2,'color',cmap(nn,:))
end
clear nn npp_change tmort_change
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('%')
grid on
set(gca,'XTickLabel','')
title('(b) F_{mort}')

subplot(3,1,3)
hold on
for nn=1:nmod
    npp_change=runmean(((npp(:,nn)./repmat(npp_mean(:,nn),nyear,1))*100)-100,5);
    tmort_change=runmean(((tmort(:,nn)./repmat(tmort_mean(:,nn),nyear,1))*100)-100,5);
    plot(1901:2099,npp_change-tmort_change,'linewidth',2,'color',cmap(nn,:))
end
clear nn npp_change tmort_change
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('\Delta (%)')
grid on
set(gca,'XTickLabel','')
title('(c) NPP-F_{mort}')

legend(models)

%Change in normalised fluxes relative to baseline (10-y running mean)
figure

subplot(3,1,1)
hold on
for nn=1:nmod
    npp_change=runmean(((npp_norm(:,nn)./repmat(npp_norm_mean(:,nn),nyear,1))*100)-100,5);
    plot(1901:2099,npp_change,'linewidth',2,'color',cmap(nn,:))
end
clear nn npp_change tmort_change
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('\Delta %')
grid on
set(gca,'XTickLabel','')
title('(a) NPP (norm.)')

subplot(3,1,2)
hold on
for nn=1:nmod
    tmort_change=runmean(((tmort_norm(:,nn)./repmat(tmort_norm_mean(:,nn),nyear,1))*100)-100,5);
    plot(1901:2099,tmort_change,'linewidth',2,'color',cmap(nn,:))
end
clear nn npp_change tmort_change
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('%')
grid on
set(gca,'XTickLabel','')
title('(b) F_{mort} (norm.)')

subplot(3,1,3)
hold on
for nn=1:nmod
    npp_change=runmean(((npp_norm(:,nn)./repmat(npp_norm_mean(:,nn),nyear,1))*100)-100,5);
    tmort_change=runmean(((tmort_norm(:,nn)./repmat(tmort_norm_mean(:,nn),nyear,1))*100)-100,5);
    plot(1901:2099,npp_change-tmort_change,'linewidth',2,'color',cmap(nn,:))
end
clear nn npp_change tmort_change
set(gca,'XLim',[2000 2099])
set(gca,'Fontsize',14)
ylabel('\Delta (%)')
grid on
set(gca,'XTickLabel','')
title('(c) NPP-F_{mort} (norm.)')

legend(models)