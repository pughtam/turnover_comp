%Make global totals of  various productivity and turnover-related variables
%
%Requires *.mat files of preprocessed model data from turnover_pool_flux_read.m
%
%T. Pugh
%30.11.19

%Location of *.mat files containing preprocessed model data
data_models='/Users/pughtam/Documents/GAP_and_other_work/Mortality/mat_files/';

%---
%First do the calculations for CRU-NCEP

ordind=[5,4,1,2,6,3]; %Order in which to display models (because output from all_cruncep_global_arrays.mat is not in alphabetical order)

models={'CABLE-POP','JULES','LPJ-GUESS','LPJmL','ORCHIDEE','SEIB-DGVM'};
nmod=length(models);

%Index of first and last year to average over the selected period
y1=85; %1985
y2=114; %2014

%---
%Load the globally averaged CRU-NCEP data from *.mat file created by turnover_pool_flux_read.m
load([data_models,'all_cruncep_global_arrays']);

%Make reproduction flux array (none for JULES or CABLE)
repro=NaN(size(npp));
repro(:,1:3)=npp(:,1:3)*0.1;
repro(:,6)=npp(:,6)*0.1;

%Calculate overall turnover fluxes
tturn=nansum(cat(3,tmort,tleaf,troot,repro),3);

%Calculate turnover times
tau_npp=cveg./npp;
tau_turn=cveg./tturn;
tau_mort=cveg./tmort;
tau_root=croot./troot;

%Calculate the means over the selected period
npp_mean=mean(npp(y1:y2,:),1);
cveg_mean=mean(cveg(y1:y2,:),1);
tau_npp_mean=mean(tau_npp(y1:y2,:),1);
tau_turn_mean=mean(tau_turn(y1:y2,:),1);
tau_mort_mean=mean(tau_mort(y1:y2,:),1);
tau_root_mean=mean(tau_root(y1:y2,:),1);

%Rearrange the order of the models alphabetically
npp_mean=npp_mean(ordind);
cveg_mean=cveg_mean(ordind);
tau_npp_mean=tau_npp_mean(ordind);
tau_turn_mean=tau_turn_mean(ordind);
tau_mort_mean=tau_mort_mean(ordind);
tau_root_mean=tau_root_mean(ordind);

%---
%Write out a table
fprintf('Variable ')
for nn=1:nmod
    fprintf('%s ',models{nn})
end
fprintf('\n')

fprintf('NPP      ')
for nn=1:nmod
    fprintf('%7.1f ',npp_mean(nn))
end
fprintf('\n')

fprintf('Cveg     ')
for nn=1:nmod
    fprintf('%7.1f ',cveg_mean(nn))
end
fprintf('\n')

fprintf('tau_npp  ')
for nn=1:nmod
    fprintf('%7.1f ',tau_npp_mean(nn))
end
fprintf('\n')

fprintf('tau_turn ')
for nn=1:nmod
    fprintf('%7.1f ',tau_turn_mean(nn))
end
fprintf('\n')

fprintf('tau_mort ')
for nn=1:nmod
    fprintf('%7.1f ',tau_mort_mean(nn))
end
fprintf('\n')

fprintf('tau_root ')
for nn=1:nmod
    fprintf('%7.1f ',tau_root_mean(nn))
end
fprintf('\n')