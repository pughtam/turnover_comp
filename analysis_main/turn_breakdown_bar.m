%Make a stacked bar plot to compare the fractions of turnover through each pathway (phenological turnover, 
%mortality) in the CRU-NCEP simulations.
%Make a second bar plot which shows how these fractions change over the 21st century based on the IPSL
%simulation.
%
%T. Pugh
%28.11.19

%Location of *.mat files containing preprocessed model data
data_models='/Users/pughtam/Documents/GAP_and_other_work/Mortality/mat_files/';

%---
%Some generic settings
ordind=[5,4,1,2,6,3]; %Order in which to display models (because output from all_cruncep_global_arrays.mat is not in alphabetical order)

models={'CABLE-POP','JULES','LPJ-GUESS','LPJmL','ORCHIDEE','SEIB-DGVM'};
nmod=length(models);

%---
%First do the calculations for CRU-NCEP

%Index of first and last year to average over the selected period
y1_cn=85; %1985
y2_cn=114; %2014

%Load the globally averaged CRU-NCEP data from *.mat file created by turnover_pool_flux_read.m
load([data_models,'all_cruncep_global_arrays']);

%Make reproduction flux array (none for JULES or CABLE)
repro=NaN(size(npp));
repro(:,1:3)=npp(:,1:3)*0.1;
repro(:,6)=npp(:,6)*0.1;

%Mean values for turnover flux across the selected period
tmort_cn_mean=mean(tmort(y1_cn:y2_cn,:),1);
tleaf_cn_mean=mean(tleaf(y1_cn:y2_cn,:),1);
troot_cn_mean=mean(troot(y1_cn:y2_cn,:),1);
repro_cn_mean=mean(repro(y1_cn:y2_cn,:),1);
turn_cn_mean=nansum(cat(1,tmort_cn_mean,tleaf_cn_mean,troot_cn_mean,repro_cn_mean),1);

%Clean up
clear cleaf croot csoil cveg cwood gpp ifcruncep ifformask ifmask npp repro tallom tback tbioclim tdist tfire
clear tgroweff theat tleaf tmort tnbio tother troot tshade

%Make array of fraction of turnover by each process, ready for plotting
bararray_cn=cat(1,tmort_cn_mean,tleaf_cn_mean,troot_cn_mean,repro_cn_mean)./repmat(turn_cn_mean,[4 1])*100;
bararray_cn=bararray_cn(:,ordind);

%Make figure
figure
s1=subplot(2,1,1);

bar(bararray_cn','stacked')
set(gca,'XTickLabel',models)
set(gca,'YLim',[0 100])
ylabel('Fraction of F_{turn} (%)')
legend('Mortality','Leaf phen.','Root phen.','Reproduction')
set(gca,'YLim',[0 100])
set(gca,'XLim',[0.5 6.5])
t1=text(0.7,90,'(a)');
set(t1,'FontSize',12,'FontWeight','Bold')
set(gca,'XTickLabel','')

%---
%Now do the calculations for IPSL

%Load the globally averaged IPSL data from *.mat file created by turnover_pool_flux_read.m
load([data_models,'all_ipsl_global_arrays']);

%Make reproduction flux array (none for JULES or CABLE)
repro=NaN(size(npp));
repro(:,1:3)=npp(:,1:3)*0.1;
repro(:,6)=npp(:,6)*0.1;

%Index of first and last year to average over for the baseline period
y1a_ip=85; %1985
y2a_ip=114; %2014
%Index of first and last year to average over for the future period
y1b_ip=170; %2070
y2b_ip=199; %2099

%Mean values for turnover flux across the selected periods
tmort_ip_mean_base=mean(tmort(y1a_ip:y2a_ip,:),1);
tleaf_ip_mean_base=mean(tleaf(y1a_ip:y2a_ip,:),1);
troot_ip_mean_base=mean(troot(y1a_ip:y2a_ip,:),1);
repro_ip_mean_base=mean(repro(y1a_ip:y2a_ip,:),1);
turn_ip_mean_base=nansum(cat(1,tmort_ip_mean_base,tleaf_ip_mean_base,troot_ip_mean_base,repro_ip_mean_base),1);

tmort_ip_mean_fut=mean(tmort(y1b_ip:y2b_ip,:),1);
tleaf_ip_mean_fut=mean(tleaf(y1b_ip:y2b_ip,:),1);
troot_ip_mean_fut=mean(troot(y1b_ip:y2b_ip,:),1);
repro_ip_mean_fut=mean(repro(y1b_ip:y2b_ip,:),1);
turn_ip_mean_fut=nansum(cat(1,tmort_ip_mean_fut,tleaf_ip_mean_fut,troot_ip_mean_fut,repro_ip_mean_fut),1);

%Make array of fraction of turnover by each process, ready for plotting
alloc_frac_ip_base=cat(1,tmort_ip_mean_base,tleaf_ip_mean_base,troot_ip_mean_base,repro_ip_mean_base)./repmat(turn_ip_mean_base,[4 1])*100;
alloc_frac_ip_fut=cat(1,tmort_ip_mean_fut,tleaf_ip_mean_fut,troot_ip_mean_fut,repro_ip_mean_fut)./repmat(turn_ip_mean_fut,[4 1])*100;
alloc_frac_ip_diff=alloc_frac_ip_fut-alloc_frac_ip_base;

bararray_ip=alloc_frac_ip_diff(:,ordind);

%Make figure
s2=subplot(2,1,2);
bar(bararray_ip')
set(gca,'XTickLabel',models)
ylabel('\Delta Fraction of F_{turn} (%)')
set(gca,'YLim',[-5 7])
set(gca,'XLim',[0.5 6.5])
t2=text(0.7,5.5,'(b)');
set(t2,'FontSize',12,'FontWeight','Bold')

%Adjust plot positions
set(s1,'Position',[0.13 0.55 0.775 0.41])
set(s2,'Position',[0.13 0.1 0.775 0.41])
