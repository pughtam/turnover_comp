function [lrrflux_mask_std_mean,mflux_mask_std_mean,lrrflux_mask_std_mean_biom,...
    mflux_mask_std_mean_biom,phen_label,nphen]=global_turn_frac_sd_time(makeplots,data_models,data_mask,data_phen)
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

%Index of first and last year to average over in the data
y1=85; %1985
y2=199; %2099

%Load in the pre-processed model data
ifcruncep=false; %IPSL data
[~,mflux,lrrflux,~,reproflux,~,mflux_jules,lrrflux_jules,...
    ~,models,nmod]=get_stocks_fluxes_runmeanstd(data_models,ifcruncep,y1,y2);

%---
%Make plots

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[fmask,fmask_jules,~,~]=get_closed_can_mask(data_mask);

%Mask data by closed-canopy forest cover (10% threshold per grid cell)
lrrflux_mask=lrrflux.*repmat(fmask',1,1,nmod);
mflux_mask=mflux.*repmat(fmask',1,1,nmod);
lrrflux_jules_mask=lrrflux_jules.*fmask_jules';
mflux_jules_mask=mflux_jules.*fmask_jules';

%Calculate global values
lrrflux_mask_std_mean=NaN(nmod,1);
mflux_mask_std_mean=NaN(nmod,1);
for mm=1:nmod
    if mm==2 %JULES is a special case because of resolution (NOTE: Comparisons of absolute numbers with other models should not be made because of resolution difference)
        lrrflux_mask_std_mean(mm)=sqrt(nanmean(lrrflux_jules_mask(:).^2));
        mflux_mask_std_mean(mm)=sqrt(nanmean(mflux_jules_mask(:).^2));
    else
        temp=lrrflux_mask(:,:,mm);
        temp2=mflux_mask(:,:,mm);
        lrrflux_mask_std_mean(mm)=sqrt(nanmean(temp(:).^2));
        mflux_mask_std_mean(mm)=sqrt(nanmean(temp2(:).^2));
        clear temp temp2
    end
end
clear mm

%Calculate values per model forest type
%Forest type files calculated using scripts in 'model_masks' directory
[phen,phen_jules,phen_label,nphen]=get_forest_type(data_phen);

lrrflux_mask_std_mean_biom=NaN(nmod,nphen);
mflux_mask_std_mean_biom=NaN(nmod,nphen);
for mm=1:nmod
    if mm==2 %JULES is a special case because of resolution
        for nn=1:nphen
            lrrflux_mask_std_mean_biom(mm,nn)=sqrt(nanmean(lrrflux_jules_mask(phen_jules==nn).^2));
            mflux_mask_std_mean_biom(mm,nn)=sqrt(nanmean(mflux_jules_mask(phen_jules==nn).^2));
        end
    else
        temp=lrrflux_mask(:,:,mm);
        temp2=mflux_mask(:,:,mm);
        for nn=1:nphen
            lrrflux_mask_std_mean_biom(mm,nn)=sqrt(nanmean(temp(phen(:,:,mm)==nn).^2));
            mflux_mask_std_mean_biom(mm,nn)=sqrt(nanmean(temp2(phen(:,:,mm)==nn).^2));
        end
    end
    clear temp temp2
end
clear mm nn

%---
%Make plots
if makeplots
    
    %Global bar chart
    figure
    bb=bar([mflux_mask_std_mean,lrrflux_mask_std_mean]');
    set(gca,'XTickLabel',{'Mortality','Phenology'})
    for nn=1:nmod
        set(bb(nn),'LineStyle','none')
    end
    clear nn
    legend(models)
    ylabel('\sigma_{time} (kg C m^{-2} a^{-1})')
    
    %Forest-type bar chart
    figure
    bb=NaN(nphen,nmod);
    for mm=1:nmod
        subplot(3,2,mm)
        bb(:,mm)=bar([mflux_mask_std_mean_biom(mm,:);lrrflux_mask_std_mean_biom(mm,:)]);
        set(gca,'YLim',[0 0.45])
        set(gca,'XTickLabel',{'Mortality','Phenology'})
        %Remove lines from bars
        for nn=1:nphen
            set(bb(nn,mm),'LineStyle','none')
        end
        clear nn
        title(models{mm})
        ylabel('\sigma_{time} (kg C m^{-2} a^{-1})')
    end
    clear mm
    legend(phen_label)
    
end