%Script to assess how much spatial variation of tau is affected by phenology and mortality.
%Calculates:
%tau=vstock/(lrrflux+mortflux)
%tau_fixmort=vstock/(lrrflux+(mean(mortflux))
%tau_fixphen=vstock/(mean(lrrflux)+mortflux)
%Then takes the diffence between tau and tau_fixmort or tau_fixphen and calculates the mean absolute deviation of tau caused
%by fixing lrrflux or mortflux
%
%Dependencies
%-get_stocks_fluxes.m
%-closed_can_mask_read.m
%-get_forest_type.m
%
%T. Pugh
%23.05.20

%Location of *.mat files containing preprocessed model data
data_models='/Users/pughtam/Documents/GAP_and_other_work/Mortality/mat_files/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/Users/pughtam/data/turnover/';
%Forest type files calculated using scripts in 'model_masks' directory
data_phen='/data/turnover/masks/phen/';

modelslabel={'CABLE-POP','JULES','LPJ-GUESS','LPJmL','ORCHIDEE','SEIB-DGVM'};

%Index of first and last year to average over in the data
y1=85; %1985
y2=114; %2014

%Load in the pre-processed model data
ifcruncep=true;
[vstock,mflux,lrflux,~,reproflux,vstock_jules,mflux_jules,lrflux_jules,...
    ~,~,nmod]=get_stocks_fluxes(data_models,ifcruncep,y1,y2);

lrrflux=nansum(cat(4,lrflux,reproflux),4);
lrrflux(isnan(lrflux))=NaN;
lrrflux_jules=lrflux_jules;
clear lrflux reproflux lrflux_jules

%---
%Mask for forest area

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[fmask,fmask_jules,~,~]=get_closed_can_mask(data_mask);

%Mask data by forest cover (10% cover threshold)
lrrflux_mask=lrrflux.*repmat(fmask',1,1,nmod);
mflux_mask=mflux.*repmat(fmask',1,1,nmod);
lrrflux_jules_mask=lrrflux_jules.*fmask_jules';
mflux_jules_mask=mflux_jules.*fmask_jules';
vstock_mask=vstock.*repmat(fmask',1,1,nmod);
vstock_jules_mask=vstock_jules.*fmask_jules';

%---
%Calculate the contributions to tau at the global level

%Mean global values of the fluxes
mflux_mask_mean=NaN(nmod,1);
lrrflux_mask_mean=NaN(nmod,1);
for mm=1:nmod
    tempmflux=mflux_mask(:,:,mm);
    templrrflux=lrrflux_mask(:,:,mm);
    mflux_mask_mean(mm)=nanmean(tempmflux(:));
    lrrflux_mask_mean(mm)=nanmean(templrrflux(:));
end
clear mm tempmflux templrrflux
mflux_mask_mean(2)=nanmean(mflux_jules_mask(:));
lrrflux_mask_mean(2)=nanmean(lrrflux_jules_mask(:));

%Calculate tau and tau using flux means
tau=NaN(360,720,nmod);
tau_fixmort=NaN(360,720,nmod);
tau_fixphen=NaN(360,720,nmod);
for mm=1:nmod
    tempvstock=vstock_mask(:,:,mm);
    tempmflux=mflux_mask(:,:,mm);
    templrrflux=lrrflux_mask(:,:,mm);
    
    tau(:,:,mm)=tempvstock./(templrrflux+tempmflux);
    tau_fixmort(:,:,mm)=tempvstock./(templrrflux+mflux_mask_mean(mm));
    tau_fixphen(:,:,mm)=tempvstock./(lrrflux_mask_mean(mm)+tempmflux);
end
clear mm tempvstock tempmflux templrrflux
tau_jules=vstock_jules_mask./(lrrflux_jules_mask+mflux_jules_mask);
tau_fixmort_jules=vstock_jules_mask./(lrrflux_jules_mask+mflux_mask_mean(2));
tau_fixphen_jules=vstock_jules_mask./(lrrflux_mask_mean(2)+mflux_jules_mask);

%Difference in tau values with and without fixed components
mortcont=tau_fixmort-tau;
phencont=tau_fixphen-tau;
mortcont_jules=tau_fixmort_jules-tau_jules;
phencont_jules=tau_fixphen_jules-tau_jules;

%Mean absolute deviation
mortcont_meanad=NaN(nmod,1);
phencont_meanad=NaN(nmod,1);
for mm=1:nmod
    if mm==2 %JULES is a special case because of resolution
        mortcont_meanad(mm)=nanmean(abs(mortcont_jules(:)));
        phencont_meanad(mm)=nanmean(abs(phencont_jules(:)));
    else
        tempmort=mortcont(:,:,mm);
        tempphen=phencont(:,:,mm);
        mortcont_meanad(mm)=nanmean(abs(tempmort(:)));
        phencont_meanad(mm)=nanmean(abs(tempphen(:)));
    end
end
clear tempmort tempphen mm

%---
%Make plot
    
figure
bb=bar([mortcont_meanad,phencont_meanad]);
set(gca,'XTickLabel',modelslabel)
for nn=1:2
    set(bb(nn),'LineStyle','none')
end
clear nn
legend('Mortality','Phenology')
ylabel('Mean deviation in \tau (years)')
set(gca,'FontSize',12,'FontWeight','normal')
