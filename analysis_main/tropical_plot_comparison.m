%Script to calculate an average biomass turnover time from Brienen et al. (2015, Nature 519, 344-348) and Galbraith et al.
%(2013, Plant Ecol. Divers. 6, 37-41) and then compare this to the turnover intercomparison models.
%
%Brienen et al. data can be downloaded from https://www.forestplots.net/data-packages/brienen-et-al-2015
%Galbraith et al. data obtained from David Galbraith
%
%Note, need to be sure that Brienen et al. data file has unix line endings (run fixascii.sh on the file) in order
%for read to proceed through whole file. Also need to ensure that all NA values in that file are
%replaced by NaN (e.g. using :%s/NA/NaN/g in vim). %In the metadata file, all empty cells in the first column 
%must be filled by some value, in order for the read routine to assign columns correctly.
%
%Requires *.mat files of preprocessed model data from turnover_pool_flux_read.m
%
%T. Pugh
%08.02.18

dir_brienen='/data/turnover/data_package_Long_term_decline_of_Amazon_carbon_sink_2015/'; %Location of Brienen et al. data
data_galbraith='/data/turnover/Galbraith_plot_data/'; %Location of Galbraith et al. plot data
data_models='/data/turnover/'; %Location of *.mat files containing preprocessed model data

firstmodyear=1901;
nmod=6;

%---
%Set some constants
sec_per_day=86400;
g_per_kg=1000;
days_per_year=365;
days_per_year_jules=360; %Only 360 days per year in JULES

%---
%Read Brienen et al. (2015) data

%31 columns:
%,PlotViewName,mid_int_date,AGB_ini,AGB_fin,AGBnetchange,AGBMor_tot,AGB_mort_unobs_yr,AGBGain_tot,
%AGB_Gain_unobs_yr,int_c_length,int_ini,int_fin,plot_area,BAnetchange,Nini,N_rec_yr,N_mor_yr,
%N_change,AGBnetchange_N,AGBGain_N,AGBRec_N,AGBMor_tot_N,BAMor_yr,BAGain_tot_yr,tot_monitor_length,
%N_mort_plus_rec,N_plots,mid_int_since_start,Ncens_interv,AGB comments1,AGB comments2

%Read main data file by census
fid_sa=fopen([dir_brienen,'data_by_cens.csv']);
headers_sa=textscan(fid_sa,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n'...
    ,1,'delimiter',',');
data_sa=textscan(fid_sa,'%s %q %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %q %q\n'...
    ,'delimiter',',');
fclose(fid_sa);

%Get indexes for columns of interest
ncol_sa=length(headers_sa);
for n=1:ncol_sa
    if(strcmp('PlotViewName',headers_sa{n}{1}))
        plot_ind=n;
    end
    if(strcmp('mid_int_date',headers_sa{n}{1}))
        intdate_ind=n;
    end
    if(strcmp('AGB_ini',headers_sa{n}{1}))
        agbini_ind=n;
    end
    if(strcmp('AGB_fin',headers_sa{n}{1}))
        agbfin_ind=n;
    end
    if(strcmp('AGBMor_tot',headers_sa{n}{1}))
        agbmor_ind=n;
    end
    if(strcmp('int_ini',headers_sa{n}{1}))
        intini_ind=n;
    end
    if(strcmp('int_fin',headers_sa{n}{1}))
        intfin_ind=n;
    end
end

%Read metadata file for plots
fid2=fopen([dir_brienen,'metadata.csv']);
meta=textscan(fid2,'%q %q %q %q %q %q %q %q %q\n','delimiter',',','Headerlines',1);
fclose(fid2);

%Extract data that we require
plotn_meta=meta{3};
lats_meta=cellfun(@str2num,meta{6});
lons_meta=cellfun(@str2num,meta{7});
nmeta=length(plotn_meta);

plotn_sa=data_sa{plot_ind};
intdate_sa=cellfun(@str2num,data_sa{intdate_ind});
agbini_sa=cellfun(@str2num,data_sa{agbini_ind});
agbfin_sa=cellfun(@str2num,data_sa{agbfin_ind});
agbmor_sa=cellfun(@str2num,data_sa{agbmor_ind});
intini_sa=cellfun(@str2num,data_sa{intini_ind});
intfin_sa=cellfun(@str2num,data_sa{intfin_ind});

agbmean_sa=mean(cat(2,agbini_sa,agbfin_sa),2);

%Calculate turnover time by census (very peaky)
agbturn_sa=agbmean_sa./agbmor_sa;
agbturn_sa(agbturn_sa>10000)=NaN;

%Make plot index to group from same plots (assumes that plot entries are all sorted together)
ncens_sa=length(plotn_sa);
plotid_sa=NaN(ncens_sa,1);
plotid_sa(1)=1;
for nn=2:ncens_sa
    if strcmp(plotn_sa{nn},plotn_sa{nn-1})
        plotid_sa(nn)=plotid_sa(nn-1);
    else
        plotid_sa(nn)=plotid_sa(nn-1)+1;
    end    
end
clear nn

nplots_sa=max(plotid_sa);

%Perform calculations by plot
agb_p_sa=NaN(nplots_sa,1);
agbmor_p_sa=NaN(nplots_sa,1);
agbturn_p_sa=NaN(nplots_sa,1);
intini_p_sa=NaN(nplots_sa,1);
intfin_p_sa=NaN(nplots_sa,1);
lats_p_sa=NaN(nplots_sa,1);
lons_p_sa=NaN(nplots_sa,1);
for pp=1:nplots_sa
    inds=find(plotid_sa==pp);
    agb_p_sa(pp)=mean(agbmean_sa(inds));
    agbmor_p_sa(pp)=mean(agbmor_sa(inds));
    agbturn_p_sa(pp)=agb_p_sa(pp)./agbmor_p_sa(pp);
    intini_p_sa(pp)=min(intini_sa(inds));
    intfin_p_sa(pp)=max(intfin_sa(inds));
    
    %Assign lat and lon for plot
    for mm=1:nmeta
        if strcmp(plotn_sa{inds(1)},plotn_meta{mm})
            meta_ind=mm;
            break
        end
    end
    clear mm
    lats_p_sa(pp)=lats_meta(meta_ind);
    lons_p_sa(pp)=lons_meta(meta_ind);
    clear meta_ind
end
clear pp inds

%Median turnover across all plots
agbturn_p_median_sa=nanmedian(agbturn_p_sa);

%---
%Read in the Galbraith et al. (2013) data for Africa and Asia

%Asia
fid_asia=fopen([dir_galbraith,'Galbraith2013_Tom_Pugh_Asia_only.csv']);
headers_asia=textscan(fid_asia,'%s %s %s %s %s %s %s %s %s %s %s %s %s\n'...
    ,1,'delimiter',',');
textscan(fid_asia,'%s %s %s %s %s %s %s %s %s %s %s %s %s\n'...
    ,1,'delimiter',',');
data_asia=textscan(fid_asia,'%q %q %s %s %s %s %s %s %s %s %s %s %q\n'...
    ,'delimiter',',');
fclose(fid_asia);

ncol_asia=length(headers_asia);

for n=1:ncol_asia
    if(strcmp('Latitude',headers_asia{n}))
        lat_ind=n;
    end
    if(strcmp('Longitude',headers_asia{n}))
        lon_ind=n;
    end
    if(strcmp('Residence Time (yrs)',headers_asia{n}))
        turn_ind=n;
    end
end

lats_asia=cellfun(@str2num,data_asia{lat_ind});
lons_asia=cellfun(@str2num,data_asia{lon_ind});
agbturn_p_asia=cellfun(@str2num,data_asia{turn_ind});
nplots_asia=length(agbturn_p_asia);

agbturn_p_median_asia=nanmedian(agbturn_p_asia);

%Africa
fid_africa=fopen([dir_galbraith,'Galbraith2013_Tom_Pugh_Africa_only.csv']);
headers_africa=textscan(fid_africa,'%s %s %s %s %s %s %s %s %s %s %s %s %s\n'...
    ,1,'delimiter',',');
textscan(fid_africa,'%s %s %s %s %s %s %s %s %s %s %s %s %s\n'...
    ,1,'delimiter',',');
data_africa=textscan(fid_africa,'%q %q %s %s %s %s %s %s %s %s %s %s %q\n'...
    ,'delimiter',',');
fclose(fid_africa);

ncol_africa=length(headers_africa);

for n=1:ncol_africa
    if(strcmp('Latitude',headers_africa{n}))
        lat_ind=n;
    end
    if(strcmp('Longitude',headers_africa{n}))
        lon_ind=n;
    end
    if(strcmp('Residence Time (yrs)',headers_africa{n}))
        turn_ind=n;
    end
end

lats_africa=cellfun(@str2num,data_africa{lat_ind});
lons_africa=cellfun(@str2num,data_africa{lon_ind});
agbturn_p_africa=cellfun(@str2num,data_africa{turn_ind});
nplots_africa=length(agbturn_p_africa);

agbturn_p_median_africa=nanmedian(agbturn_p_africa);

%-----
%Now load in the model data and calculate a median turnover time relative to mortality for the
%gridcells which have plots (if there is more than one plot in a gridcell, that gridcell should be
%represented more than once, in order to generate a comparable statistic).

load([data_models,'/cable_map_arrays_cruncep.mat'])
tmort_cable=tdist_cable+tgroweff_cable+tother_cable;
tmort_cable=flip(permute(tmort_cable,[2 1 3])*sec_per_day*days_per_year,1);
cveg_cable=flip(permute(cveg_cable,[2 1 3]),1);
clear cleaf_cable croot_cable csoil_cable cwood_cable gpp_cable npp_cable tdist_cable tgroweff_cable tleaf_cable tother_cable troot_cable

load([data_models,'/jules_map_arrays_cruncep.mat'])
tmort_julesa=tback_jules+tother_jules;
tmort_julesa=permute(tmort_julesa,[2 1 3])*sec_per_day*days_per_year_jules;
cveg_julesa=permute(cveg_jules,[2 1 3]);
tmort_jules=NaN(144,192,114); tmort_jules(29:140,:,:)=tmort_julesa;
cveg_jules=NaN(144,192,114); cveg_jules(29:140,:,:)=cveg_julesa;
clear tmort_julesa cveg_julesa
clear cleaf_jules croot_jules csoil_jules cwood_jules gpp_jules npp_jules tdist_jules tgroweff_jules tleaf_jules tother_jules troot_jules

load([data_models,'/lpjg_map_arrays_cruncep.mat'])
tmort_lpjg=tgroweff_lpjg+tfire_lpjg+tdist_lpjg+tother_lpjg;
tmort_lpjg=permute(tmort_lpjg,[2 1 3])*sec_per_day*days_per_year;
cveg_lpjg=permute(cveg_lpjg,[2 1 3]);
clear allfire_lpjg cwood_lpjg cleaf_lpjg croot_lpjg csoil_lpjg tleaf_lpjg troot_lpjg tfire_lpjg tgroweff_lpjg tdist_lpjg tother_lpjg npp_lpjg gpp_lpjg
clear tbioclim_lpjg tback_lpjg tnbio_lpjg tallom_lpjg

load([data_models,'/lpjml_map_arrays_cruncep.mat'])
tmort_lpjml=tfire_lpjml+tgroweff_lpjml+tother_lpjml;
tmort_lpjml=permute(tmort_lpjml,[2 1 3])*sec_per_day*days_per_year;
cveg_lpjml=permute(cveg_lpjml,[2 1 3]);
clear cleaf_lpjml croot_lpjml csoil_lpjml cwood_lpjml gpp_lpjml npp_lpjml tdist_lpjml tgroweff_lpjml tleaf_lpjml tother_lpjml troot_lpjml

load([data_models,'/orchidee_map_arrays_cruncep.mat'])
tmort_orchidee=tmort_orchidee;
tmort_orchidee=flip(permute(tmort_orchidee,[2 1 3])*days_per_year/g_per_kg,1);
cveg_orchidee=flip(permute(cveg_orchidee,[2 1 3]),1);
clear cleaf_orchidee croot_orchidee csoil_orchidee cwood_orchidee gpp_orchidee npp_orchidee tdist_orchidee tgroweff_orchidee tleaf_orchidee tother_orchidee troot_orchidee

load([data_models,'/seib_map_arrays_cruncep.mat'])
tmort_seib=tdist_seib+tfire_seib+tgroweff_seib+tother_seib;
tmort_seib=flip(permute(tmort_seib,[2 1 3])*sec_per_day*days_per_year,1);
cveg_seib=flip(permute(cveg_seib,[2 1 3]),1);
clear cleaf_seib croot_seib csoil_seib cwood_seib gpp_seib npp_seib tdist_seib tgroweff_seib tleaf_seib tother_seib troot_seib

%Cycle through plots and regions
nregion=3;
plotns=cat(1,nplots_sa,nplots_africa,nplots_asia);
plotn_max=max(plotns);

cruncep_corr=0.5; %Adjust latitudes for all except CABLE to correct for CRU-NCEP v5 indexing error (take one gridcell further north)

tmort_cable_p=NaN(plotn_max,nregion);
cveg_cable_p=NaN(plotn_max,nregion);
tmort_jules_p=NaN(plotn_max,nregion);
cveg_jules_p=NaN(plotn_max,nregion);
tmort_lpjg_p=NaN(plotn_max,nregion);
cveg_lpjg_p=NaN(plotn_max,nregion);
tmort_lpjml_p=NaN(plotn_max,nregion);
cveg_lpjml_p=NaN(plotn_max,nregion);
tmort_orchidee_p=NaN(plotn_max,nregion);
cveg_orchidee_p=NaN(plotn_max,nregion);
tmort_seib_p=NaN(plotn_max,nregion);
cveg_seib_p=NaN(plotn_max,nregion);
for nn=1:nregion
    for pp=1:plotns(nn)
        if nn==1 %South America
            lon_ind=floor((lons_p_sa(pp)+180)*2);
            lat_ind_cable=floor((lats_p_sa(pp)+90)*2);
            lat_ind=floor((lats_p_sa(pp)+cruncep_corr+90)*2);
            lons_p_jules=lons_p_sa; lons_p_jules(lons_p_sa<0)=lons_p_jules(lons_p_sa<0)+360; %JULES has different resolution and starts at greenwich meridian, therefore different indices
            lon_jules_ind=floor(lons_p_jules(pp)/1.875);
            lat_jules_ind=floor((lats_p_sa(pp)+cruncep_corr+90)/1.25);
            year_min=floor(intini_p_sa(pp))-firstmodyear+1; %Use actual census intervals
            year_max=floor(intfin_p_sa(pp))-firstmodyear+1;
            %year_min=1985-firstmodyear+1; 
            %year_max=2014-firstmodyear+1;
        elseif nn==2 %Africa
            lon_ind=floor((lons_africa(pp)+180)*2);
            lat_ind_cable=floor((lats_africa(pp)+90)*2);
            lat_ind=floor((lats_africa(pp)+cruncep_corr+90)*2);
            lons_p_jules=lons_africa; lons_p_jules(lons_africa<0)=lons_p_jules(lons_africa<0)+360;
            lon_jules_ind=floor(lons_p_jules(pp)/1.875);
            lat_jules_ind=floor((lats_africa(pp)+cruncep_corr+90)/1.25);
            year_min=1985-firstmodyear+1; %Galbraith et al. (2013) don't give the census intervals, so use standard year range for this study
            year_max=2014-firstmodyear+1;
        elseif nn==3  %Asia
            lon_ind=floor((lons_asia(pp)+180)*2);
            lat_ind_cable=floor((lats_asia(pp)+90)*2);
            lat_ind=floor((lats_asia(pp)+cruncep_corr+90)*2);
            lons_p_jules=lons_asia; lons_p_jules(lons_asia<0)=lons_p_jules(lons_asia<0)+360;
            lon_jules_ind=floor(lons_p_jules(pp)/1.875);
            lat_jules_ind=floor((lats_asia(pp)+cruncep_corr+90)/1.25);
            year_min=1985-firstmodyear+1; %Galbraith et al. (2013) don't give the census intervals, so use standard year range for this study
            year_max=2014-firstmodyear+1;
        else
            error('nregion > 3')
        end
        
        tmort_cable_p(pp,nn)=nanmean(tmort_cable(lat_ind_cable,lon_ind,year_min:year_max));
        cveg_cable_p(pp,nn)=nanmean(cveg_cable(lat_ind_cable,lon_ind,year_min:year_max));
        
        tmort_lpjg_p(pp,nn)=nanmean(tmort_lpjg(lat_ind,lon_ind,year_min:year_max));
        cveg_lpjg_p(pp,nn)=nanmean(cveg_lpjg(lat_ind,lon_ind,year_min:year_max));
        tmort_lpjml_p(pp,nn)=nanmean(tmort_lpjml(lat_ind,lon_ind,year_min:year_max));
        cveg_lpjml_p(pp,nn)=nanmean(cveg_lpjml(lat_ind,lon_ind,year_min:year_max));
        tmort_orchidee_p(pp,nn)=nanmean(tmort_orchidee(lat_ind,lon_ind,year_min:year_max));
        cveg_orchidee_p(pp,nn)=nanmean(cveg_orchidee(lat_ind,lon_ind,year_min:year_max));
        tmort_seib_p(pp,nn)=nanmean(tmort_seib(lat_ind,lon_ind,year_min:year_max));
        cveg_seib_p(pp,nn)=nanmean(cveg_seib(lat_ind,lon_ind,year_min:year_max));
        
        tmort_jules_p(pp,nn)=nanmean(tmort_jules(lat_jules_ind,lon_jules_ind,year_min:year_max));
        cveg_jules_p(pp,nn)=nanmean(cveg_jules(lat_jules_ind,lon_jules_ind,year_min:year_max));
    end
end

%Calculate turnover times for each plot location
turn_cable_p=cveg_cable_p./tmort_cable_p;
turn_jules_p=cveg_jules_p./tmort_jules_p;
turn_lpjg_p=cveg_lpjg_p./tmort_lpjg_p;
turn_lpjml_p=cveg_lpjml_p./tmort_lpjml_p;
turn_orchidee_p=cveg_orchidee_p./tmort_orchidee_p;
turn_seib_p=cveg_seib_p./tmort_seib_p;

%Calculate median turnover times across all plot locations
turn_cable_p_median=nanmedian(turn_cable_p);
turn_jules_p_median=nanmedian(turn_jules_p);
turn_lpjg_p_median=nanmedian(turn_lpjg_p);
turn_lpjml_p_median=nanmedian(turn_lpjml_p);
turn_orchidee_p_median=nanmedian(turn_orchidee_p);
turn_seib_p_median=nanmedian(turn_seib_p);

%-----
%Make a boxplot comparing the models and observations

figure
for nn=1:nregion
    subplot(3,1,nn)
    
    if nn==1 %South America
        obs=agbturn_p_sa;
        titlestr='(a) South America';
    elseif nn==2 %Africa
        obs=agbturn_p_africa;
        titlestr='(b) Africa';
    elseif nn==3 %Asia
        obs=agbturn_p_asia;
        titlestr='(c) Asia/Oceania';
    end
    
    notnan_cable=find(isnan(turn_cable_p(:,nn))==0);
    notnan_jules=find(isnan(turn_jules_p(:,nn))==0);
    notnan_lpjg=find(isnan(turn_lpjg_p(:,nn))==0);
    notnan_lpjml=find(isnan(turn_lpjml_p(:,nn))==0);
    notnan_orchidee=find(isnan(turn_orchidee_p(:,nn))==0);
    notnan_seib=find(isnan(turn_seib_p(:,nn))==0);
    N=[length(notnan_cable),length(notnan_jules),length(notnan_lpjg),...
        length(notnan_lpjml),length(notnan_orchidee),length(notnan_seib)];
    
    plotarray=[obs(notnan_cable); turn_cable_p(notnan_cable,nn); NaN(size(notnan_cable));...
        obs(notnan_jules); turn_jules_p(notnan_jules,nn); NaN(size(notnan_jules));...
        obs(notnan_lpjg); turn_lpjg_p(notnan_lpjg,nn); NaN(size(notnan_lpjg));...
        obs(notnan_lpjml); turn_lpjml_p(notnan_lpjml,nn); NaN(size(notnan_lpjml));...
        obs(notnan_orchidee); turn_orchidee_p(notnan_orchidee,nn); NaN(size(notnan_orchidee));...
        obs(notnan_seib); turn_seib_p(notnan_seib,nn)];
    plotgrouping=[ones(size(notnan_cable));2*ones(size(notnan_cable));3*ones(size(notnan_cable));...
        4*ones(size(notnan_jules));5*ones(size(notnan_jules));6*ones(size(notnan_jules));...
        7*ones(size(notnan_lpjg));8*ones(size(notnan_lpjg));9*ones(size(notnan_lpjg));...
        10*ones(size(notnan_lpjml));11*ones(size(notnan_lpjml));12*ones(size(notnan_lpjml));...
        13*ones(size(notnan_orchidee));14*ones(size(notnan_orchidee));15*ones(size(notnan_orchidee));...
        16*ones(size(notnan_seib));17*ones(size(notnan_seib))];
    
    boxplot(plotarray,plotgrouping,...
        'plotstyle','compact','notch','on','colors',[0.5 0.5 0.5],'symbol','.');
    
    set(gca,'XTick',1.5:3:18)
    set(gca,'XTickLabel',{'CABLE-POP','JULES','LPJ-GUESS','LPJmL','ORCHIDEE','SEIB'})
    set(gca,'YLim',[0 200])
    ylabel('\tau_{mort} (years)')
    title(titlestr)
    cc=0;
    for pp=1:3:17
        cc=cc+1;
        text(pp,190,mat2str(N(cc)))
    end
    clear cc pp
end
clear nn

% set(gca,'XTickLabel',' ')
% set(gca,'Position',[0.1 0.7 0.85 0.24])
% set(gca,'Position',[0.1 0.4 0.85 0.24])
% set(gca,'Position',[0.1 0.1 0.85 0.24])