function [lrrflux_mask_std,mflux_mask_std,lrrflux_mask_std_biom,mflux_mask_std_biom,...
    phen_label,nphen]=global_turn_frac_sd(makeplots)
%Script to make bar plots of standard deviation of turnover flux from phenology and
%mortality in space.
%
%Dependencies
%-get_stocks_fluxes.m
%-closed_can_mask_read.m
%-get_forest_type.m
%
%T. Pugh
%23.02.18

modelslabel={'CABLE-POP','JULES','LPJ-GUESS','LPJmL','ORCHIDEE','SEIB-DGVM'};

%Index of first and last year to average over in the data
y1=85; %1985
y2=114; %2014

%Load in the pre-processed model data
ifcruncep=true;
data_models='/media/pughtam/rds-2017-pughtam-01/turnover/turnover_comp/data_analysis/';
[~,mflux,lrflux,~,reproflux,~,mflux_jules,lrflux_jules,...
    ~,models,nmod]=get_stocks_fluxes(data_models,ifcruncep,y1,y2);

lrrflux=nansum(cat(4,lrflux,reproflux),4);
lrrflux(isnan(lrflux))=NaN;
lrrflux_jules=lrflux_jules;
clear lrflux reproflux lrflux_jules

%---
%Make plots

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
data_mask='/data/turnover/';
[fmask,fmask_jules,~,~]=get_closed_can_mask(data_mask);

%Mask data by forest cover (10% cover threshold)
lrrflux_mask=lrrflux.*repmat(fmask',1,1,nmod);
mflux_mask=mflux.*repmat(fmask',1,1,nmod);
lrrflux_jules_mask=lrrflux_jules.*fmask_jules';
mflux_jules_mask=mflux_jules.*fmask_jules';

%Calculate global values
lrrflux_mask_std=NaN(nmod,1);
mflux_mask_std=NaN(nmod,1);
for mm=1:nmod
    if mm==2 %JULES is a special case because of resolution (NOTE: Comparisons of absolute numbers with other models should not be made because of resolution difference)
        lrrflux_mask_std(mm)=nanstd(lrrflux_jules_mask(:));
        mflux_mask_std(mm)=nanstd(mflux_jules_mask(:));
    else
        temp=lrrflux_mask(:,:,mm);
        temp2=mflux_mask(:,:,mm);
        lrrflux_mask_std(mm)=nanstd(temp(:));
        mflux_mask_std(mm)=nanstd(temp2(:));
        clear temp temp2
    end
end
clear mm

%Calculate values per model forest type
%Forest type files calculated using scripts in 'model_masks' directory
phen_dir='/data/turnover/masks/phen/';
[phen,phen_jules,phen_label,nphen]=get_forest_type(phen_dir);

lrrflux_mask_std_biom=NaN(nmod,nphen);
mflux_mask_std_biom=NaN(nmod,nphen);
lrrflux_mask_mean_biom=NaN(nmod,nphen);
mflux_mask_mean_biom=NaN(nmod,nphen);
for mm=1:nmod
    if mm==2 %JULES is a special case because of resolution
        for nn=1:nphen
            lrrflux_mask_std_biom(mm,nn)=nanstd(lrrflux_jules_mask(phen_jules==nn));
            mflux_mask_std_biom(mm,nn)=nanstd(mflux_jules_mask(phen_jules==nn));
            lrrflux_mask_mean_biom(mm,nn)=nanmean(lrrflux_jules_mask(phen_jules==nn));
            mflux_mask_mean_biom(mm,nn)=nanmean(mflux_jules_mask(phen_jules==nn));
        end
    else
        temp=lrrflux_mask(:,:,mm);
        temp2=mflux_mask(:,:,mm);
        for nn=1:nphen
            lrrflux_mask_std_biom(mm,nn)=nanstd(temp(phen(:,:,mm)==nn));
            mflux_mask_std_biom(mm,nn)=nanstd(temp2(phen(:,:,mm)==nn));
            lrrflux_mask_mean_biom(mm,nn)=nanmean(temp(phen(:,:,mm)==nn));
            mflux_mask_mean_biom(mm,nn)=nanmean(temp2(phen(:,:,mm)==nn));
        end
    end
    clear temp temp2
end
clear mm nn

%Calculate variability within forest types versus variability between forest types

lrrflux_mask_std_biom_withinbiom=sqrt(nanmean(lrrflux_mask_std_biom.^2,2));
lrrflux_mask_std_biom_acrossbiom=nanstd(lrrflux_mask_mean_biom,[],2);
mflux_mask_std_biom_withinbiom=sqrt(nanmean(mflux_mask_std_biom.^2,2));
mflux_mask_std_biom_acrossbiom=nanstd(mflux_mask_mean_biom,[],2);

%---
%Make plots
if makeplots
    
    %Global bar chart
    figure
    bb=bar([mflux_mask_std,lrrflux_mask_std]');
    set(gca,'XTickLabel',{'Mortality','Phenology'})
    for nn=1:nmod
        set(bb(nn),'LineStyle','none')
    end
    clear nn
    legend(models)
    ylabel('\sigma_{space} (kg C m^{-2} a^{-1})')
    
    %Forest-type bar chart
    figure
    bb=NaN(nphen,nmod);
    for mm=1:nmod
        subplot(3,2,mm)
        bb(:,mm)=bar([mflux_mask_std_biom(mm,:);lrrflux_mask_std_biom(mm,:)]);
        set(gca,'YLim',[0 0.45])
        set(gca,'XTickLabel',{'Mortality','Phenology'})
        %Remove lines from bars
        for nn=1:nphen
            set(bb(nn,mm),'LineStyle','none')
        end
        clear nn
        title(models{mm})
    end
    clear mm
    legend(phen_label)
    
    %Within versus across forest-type bar chart
    figure
    bararray=cat(2,mflux_mask_std_biom_withinbiom,mflux_mask_std_biom_acrossbiom,lrrflux_mask_std_biom_withinbiom,lrrflux_mask_std_biom_acrossbiom);
    bar(bararray)
    set(gca,'XTickLabel',modelslabel)
    legend('Mortality (within forest type)','Mortality (across forest types)','Phenology (within forest type)','Phenology (across forest types)')
    ylabel('\sigma_{space} (kg C m^{-2} a^{-1})')
    
end