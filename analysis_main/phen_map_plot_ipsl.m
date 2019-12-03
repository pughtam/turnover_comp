%Script to make maps of forest type for the IPSL simulation for two time
%periods, following the protocol laid down by Sarah Shafer (no ISPL
%phenology masks from Sarah available). See scripts in forest_masks directory
%For more details on that protocol.
%
%The indexes to the forest-types are as follows:
%1=needleleaved evergreen, 
%2=needleleaved deciduous, 
%3=boreal broadleaved deciduous, 
%4=temperate broadleaved evergreen, 
%5=temperate broadleaved deciduous, 
%6=tropical broadleaved evergreen,
%7=tropical broadleaved deciduous
%
%T. Pugh
%24.02.18

%Location of forest mask files created using file in folder "forest_masks"
phen_dir='/data/turnover/';

%Sub-paths for the individual LAI files
lai_file_cable='CABLE/IPSL_v2/CABLE-POP_ipsl-cm5a-lr_lai_month_1901_2099.nc4';
lai_file_jules='JULES/IPSL_v2/JULESC2_IPSL_lai_Monthly_1901_2099.nc';
lai_file_lpjg='LPJG/IPSL_v3/netcdfs/lpj-guess_ipsl_lai_monthly_1901_2099.nc4';
lai_file_lpjml='LPJmL/IPSL/lpjml_ipsl-cm5a_lr_lai_annual_1901_2099.nc';
lai_file_orchidee='ORCHIDEE/IPSL_v2/ORCHIDEEr3085_ipslcm5alr_lai_13pft_year_1901_2099.nc';
lai_file_seib='SEIB/IPSL/netcdf_v6/seib_ipsl_lai_month_1901_2099.nc4';

%Sub-path for the vegetation cover frac file which is needed to correct for grid cells with less
%than 100% land area in the ORCHIDEE data
vegfrac_orchidee='ORCHIDEE/IPSL_v2/ORCHIDEEr3085_ipslcm5alr_veget_max_13pft_year_1901_2099.nc';

%Location of *.mat files containing preprocessed ESA landcover data
data_esa='/data/ESA_landcover/';
%Location of netcdf files containing closed canopy forest mask
data_mask='/data/turnover/';

%---
models={'(a) CABLE-POP','(b) JULES','(c) LPJ-GUESS','(d) LPJmL','(e) ORCHIDEE','(f) SEIB-DGVM'};
nmod=length(models);

phen_label={'needleleaved evergreen','needleleaved deciduous','boreal broadleaved deciduous',...
    'temperate broadleaved evergreen','temperate broadleaved deciduous','tropical broadleaved evergreen',...
    'tropical broadleaved deciduous'};
nphen=length(phen_label);

%Index of first year for the baseline and future period
y1=[85 170]; %1985, 2070
%Index of last year for the baseline and future period
y2=[114 199]; %2014, 2099

nyear=[length(y1(1):y2(1)),length(y1(2):y2(2))];

%Convert to months in order to extract from netcdf files which are monthly
m1=(y1*12)-11;
m2=y2*12;
nmon=[length(m1(1):m2(1)),length(m1(2):m2(2))];

if nyear(1)~=nyear(2)
    error('Different number of years for each period')
end

%---
%Read LAI data

%CABLE (monthly file)
lai_cable_ann=NaN(720,360,10,nyear(1),2);
for tt=1:length(m1)
    lai_cable=ncread([phen_dir,lai_file_cable],'lai',[1 1 1 m1(tt)],[Inf Inf Inf nmon(tt)]);
    for yy=1:nyear(tt)
        m_s=(yy*12)-11;
        m_e=yy*12;
        lai_cable_ann(:,:,:,yy,tt)=nanmax(lai_cable(:,:,:,m_s:m_e),[],4); %Maximum monthly LAI
    end
    clear yy m_s m_e
end
lai_cable_annmean=squeeze(mean(lai_cable_ann,4));
lai_cable_annmean=flip(permute(lai_cable_annmean,[2 1 3 4]),1);
clear lai_cable lai_cable_ann

%JULES (monthly file)
lai_jules_ann=NaN(192,144,9,nyear(1),2);
for tt=1:length(m1)
    lai_jules=ncread([phen_dir,lai_file_jules],'lai',[1 1 1 m1(tt)],[Inf Inf Inf nmon(tt)]);
    for yy=1:nyear(tt)
        m_s=(yy*12)-11;
        m_e=yy*12;
        lai_jules_ann(1:96,29:140,:,yy,tt)=nanmax(lai_jules(97:192,:,:,m_s:m_e),[],4); %Maximum monthly LAI
        lai_jules_ann(97:192,29:140,:,yy,tt)=nanmax(lai_jules(1:96,:,:,m_s:m_e),[],4);
    end
    clear yy m_s m_e
end
lai_jules_annmean=squeeze(mean(lai_jules_ann,4));
lai_jules_annmean=permute(lai_jules_annmean,[2 1 3 4]);
clear lai_jules lai_jules_ann

%LPJ-GUESS (monthly file)
lai_lpjg_ann=NaN(720,360,11,nyear(1),2);
for tt=1:length(m1)
    lai_lpjg=ncread([phen_dir,lai_file_lpjg],'lai',[1 1 1 m1(tt)],[Inf Inf Inf nmon(tt)]);
    for yy=1:nyear(tt)
        m_s=(yy*12)-11;
        m_e=yy*12;
        lai_lpjg_ann(:,:,:,yy,tt)=nanmax(lai_lpjg(:,:,:,m_s:m_e),[],4); %Maximum monthly LAI
    end
    clear yy m_s m_e
end
lai_lpjg_annmean=squeeze(mean(lai_lpjg_ann,4));
lai_lpjg_annmean=permute(lai_lpjg_annmean,[2 1 3 4]);
clear lai_lpjg lai_lpjg_ann

%LPJmL (annual file)
lai_lpjml_ann=NaN(720,360,9,nyear(1),2);
for tt=1:length(m1)
    lai_lpjml_ann(:,69:347,:,:,tt)=ncread([phen_dir,lai_file_lpjml],'lai',[1 1 1 y1(tt)],[Inf Inf Inf nyear(tt)]);
end
lai_lpjml_annmean=squeeze(mean(lai_lpjml_ann,4));
lai_lpjml_annmean=permute(lai_lpjml_annmean,[2 1 3 4]);
clear lai_lpjml_ann

%ORCHIDEE (annual file, needs vegetation cover fraction correction)
lai_orchidee_ann=NaN(720,360,13,nyear(1),2);
for tt=1:length(m1)
    lai_orchidee_ann(:,:,:,:,tt)=ncread([phen_dir,lai_file_orchidee],'lai',[1 1 1 y1(tt)],[Inf Inf Inf nyear(tt)]).*...
        ncread([phen_dir,vegfrac_orchidee],'veget_max',[1 1 1 y1(tt)],[Inf Inf Inf nyear(tt)]);
    lai_orchidee_ann(lai_orchidee_ann>100)=NaN;
end
lai_orchidee_annmean=squeeze(mean(lai_orchidee_ann,4));
lai_orchidee_annmean=flip(permute(lai_orchidee_annmean,[2 1 3 4]),1);
clear lai_orchidee lai_orchidee_ann

%SEIB (monthly file)
lai_seib_ann=NaN(720,360,14,nyear(1),2);
for tt=1:length(m1)
    lai_seib=ncread([phen_dir,lai_file_seib],'lai',[1 1 m1(tt) 1],[Inf Inf nmon(tt) Inf]);
    for yy=1:nyear(tt)
        m_s=(yy*12)-11;
        m_e=yy*12;
        lai_seib_ann(:,:,:,yy,tt)=squeeze(nanmax(lai_seib(:,:,m_s:m_e,:),[],3)); %Maximum monthly LAI
    end
    clear yy m_s m_e
end
lai_seib_annmean=squeeze(mean(lai_seib_ann,4));
lai_seib_annmean=flip(permute(lai_seib_annmean,[2 1 3 4]),1);
clear lai_seib lai_seib_ann

%---
%Group to phenology classes

%Note that for models which don't distinguish boreal or tropical PFTs,
%latitude cut-offs of 55 degrees and 23 degrees respectively are used.
[~,lats]=meshgrid(-179.75:0.5:179.75,-89.75:0.5:89.75);
tropmask=false(360,720);
tropmask(abs(lats)<=23)=true;
tempmask=false(360,720);
tempmask(abs(lats)>23 & abs(lats)<=55)=true;
bormask=false(360,720);
bormask(abs(lats)>55)=true;
[~,lats_jules]=meshgrid(-179.0625:1.875:179.0625,-89.375:1.25:89.375);
tropmask_jules=false(144,192);
tropmask_jules(abs(lats_jules)<=23)=true;
tempmask_jules=false(144,192);
tempmask_jules(abs(lats_jules)>23 & abs(lats_jules)<=55)=true;
bormask_jules=false(144,192);
bormask_jules(abs(lats_jules)>55)=true;

lai_annmean_phen=NaN(nmod,360,720,nphen,2);
lai_annmean_phen_jules=NaN(144,192,nphen,2);
%CABLE
%1=Evergreen Needleleaf Forest, 2=Evergreen Broadleaf Forest, 3=Deciduous Needleleaf Forest, 4=Deciduous Broadleaf Forest, 5=shrub, 6=C3 grass, 7=C4 grass, 8=tundra
lai_annmean_phen(1,:,:,1,:)=lai_cable_annmean(:,:,1,:); %needleleaved evergreen
lai_annmean_phen(1,:,:,2,:)=lai_cable_annmean(:,:,3,:); %needleleaved deciduous
lai_annmean_phen(1,:,:,3,:)=squeeze(lai_cable_annmean(:,:,4,:)).*repmat(bormask,1,1,2); %boreal broadleaved deciduous
lai_annmean_phen(1,:,:,4,:)=squeeze(lai_cable_annmean(:,:,2,:)).*repmat(tempmask,1,1,2); %temperate broadleaved evergreen
lai_annmean_phen(1,:,:,5,:)=squeeze(lai_cable_annmean(:,:,4,:)).*repmat(tempmask,1,1,2); %temperate broadleaved deciduous
lai_annmean_phen(1,:,:,6,:)=squeeze(lai_cable_annmean(:,:,2,:)).*repmat(tropmask,1,1,2); %tropical broadleaved evergreen
lai_annmean_phen(1,:,:,7,:)=squeeze(lai_cable_annmean(:,:,4,:)).*repmat(tropmask,1,1,2); %tropical broadleaved deciduous

%JULES
%1=Tropical broadleaved evergreen tree, 2=Temperate broadleaved evergreen tree, 3=Broadleaf deciduous tree, 4=Needleleaf evergreen tree, 5=Needleleaf deciduous tree,
%6=C3 grass, 7=C4 grass, 8=evergreen shrub, 9=deciduous shrub
lai_annmean_phen_jules(:,:,1,:)=lai_jules_annmean(:,:,4,:); %needleleaved evergreen
lai_annmean_phen_jules(:,:,2,:)=lai_jules_annmean(:,:,5,:); %needleleaved deciduous
lai_annmean_phen_jules(:,:,3,:)=squeeze(lai_jules_annmean(:,:,3,:)).*repmat(bormask_jules,1,1,2); %boreal broadleaved deciduous
lai_annmean_phen_jules(:,:,4,:)=lai_jules_annmean(:,:,2,:); %temperate broadleaved evergreen
lai_annmean_phen_jules(:,:,5,:)=squeeze(lai_jules_annmean(:,:,3,:)).*repmat(tempmask_jules,1,1,2); %temperate broadleaved deciduous
lai_annmean_phen_jules(:,:,6,:)=lai_jules_annmean(:,:,1,:); %tropical broadleaved evergreen
lai_annmean_phen_jules(:,:,7,:)=squeeze(lai_jules_annmean(:,:,3,:)).*repmat(tropmask_jules,1,1,2); %tropical broadleaved deciduous

%LPJ-GUESS
%1=Boreal needleleaf evergreen, 2=Boreal shade-intolerant needleleaf evergreen, 3=Boreal needleleaf summergreen, 4=Temperate broadleaf summergreen,
%5=Temperate shade-intolerant broadleaf summergreen, 6=Temperate broadleaf evergreen, 7=Tropical broadleaf evergreen, 
%8=Tropical shade-intolerant broadleaf evergreen, 9=Tropical broadleaf raingreen, 10=C3 grasses, 11=C4 grasses
lai_annmean_phen(3,:,:,1,:)=lai_lpjg_annmean(:,:,1,:)+lai_lpjg_annmean(:,:,2,:); %needleleaved evergreen
lai_annmean_phen(3,:,:,2,:)=lai_lpjg_annmean(:,:,3,:); %needleleaved deciduous
lai_annmean_phen(3,:,:,3,:)=squeeze(lai_lpjg_annmean(:,:,4,:)).*repmat(bormask,1,1,2)+squeeze(lai_lpjg_annmean(:,:,5,:)).*repmat(bormask,1,1,2); %boreal broadleaved deciduous
lai_annmean_phen(3,:,:,4,:)=lai_lpjg_annmean(:,:,6,:); %temperate broadleaved evergreen
lai_annmean_phen(3,:,:,5,:)=squeeze(lai_lpjg_annmean(:,:,4,:)).*repmat(tempmask,1,1,2)+squeeze(lai_lpjg_annmean(:,:,5,:)).*repmat(tempmask,1,1,2); %temperate broadleaved deciduous
lai_annmean_phen(3,:,:,6,:)=lai_lpjg_annmean(:,:,7,:)+lai_lpjg_annmean(:,:,8,:); %tropical broadleaved evergreen
lai_annmean_phen(3,:,:,7,:)=lai_lpjg_annmean(:,:,9,:); %tropical broadleaved deciduous

%LPJmL
%1=Tropical broad-leaved evergreen tree, 2=Tropical broad-leaved raingreen tree, 3=Temperate needle-leaved evergreen tree, 4=Temperate broad-leaved evergreen tree,
%5=Temperate broad-leaved summergreen tree, 6=Boreal needle-leaved evergreen tree, 7=Boreal broad-leaved summergreen tree, 8=Temperate herbaceous, 9=Tropical herbaceous
lai_annmean_phen(4,:,:,1,:)=lai_lpjml_annmean(:,:,6,:)+lai_lpjml_annmean(:,:,3,:); %needleleaved evergreen
%needleleaved deciduous
lai_annmean_phen(4,:,:,3,:)=lai_lpjml_annmean(:,:,7,:); %boreal broadleaved deciduous
lai_annmean_phen(4,:,:,4,:)=lai_lpjml_annmean(:,:,4,:); %temperate broadleaved evergreen
lai_annmean_phen(4,:,:,5,:)=lai_lpjml_annmean(:,:,5,:); %temperate broadleaved deciduous
lai_annmean_phen(4,:,:,6,:)=lai_lpjml_annmean(:,:,1,:); %tropical broadleaved evergreen
lai_annmean_phen(4,:,:,7,:)=lai_lpjml_annmean(:,:,2,:); %tropical broadleaved deciduous

%ORCHIDEE
%1=Bare ground, 2=tropical broad-leaved evergreen, 3=tropical broad-leaved raingreen, 4=temperate needleleaf evergreen, 5=temperate broad-leaved evergreen,
%6=temperate broad-leaved summergreen, 7=boreal needleleaf evergreeen, 8=boreal broad-leaved summergreen, 9=boreal needleleaf summergreen, 10=C3 grass, 11=C4 grass
lai_annmean_phen(5,:,:,1,:)=lai_orchidee_annmean(:,:,7,:)+lai_orchidee_annmean(:,:,4,:); %needleleaved evergreen
lai_annmean_phen(5,:,:,2,:)=lai_orchidee_annmean(:,:,9,:); %needleleaved deciduous
lai_annmean_phen(5,:,:,3,:)=lai_orchidee_annmean(:,:,8,:); %boreal broadleaved deciduous
lai_annmean_phen(5,:,:,4,:)=lai_orchidee_annmean(:,:,5,:); %temperate broadleaved evergreen
lai_annmean_phen(5,:,:,5,:)=lai_orchidee_annmean(:,:,6,:); %temperate broadleaved deciduous
lai_annmean_phen(5,:,:,6,:)=lai_orchidee_annmean(:,:,2,:); %tropical broadleaved evergreen
lai_annmean_phen(5,:,:,7,:)=lai_orchidee_annmean(:,:,3,:); %tropical broadleaved deciduous

%SEIB
%1=Tropical broad-leaved evergreen tree (Canopy, Shade tolerant), 2=Tropical broad-leaved evergreen tree (sub-Canopy, Shade tolerant),
%3=Tropical broad-leaved evergreen tree (Pioneer, Light demanding), 4=Tropical broad-leaved evergreen tree (Short, Intermediate shade tolerant),
%5=Tropical broad-leaved evergreen tree (Only distributes on Africa), 6=Tropical broad-leaved raingreen tree (Only distributes on Africa),
%7=Temperate needle-leaved evergreen tree, 8=Temperate broad-leaved evergreen tree, 9=Temperate broad-leaved summergreen tree,
%10=Boreal needle-leaved evergreen tree, 11=Boreal needle-leaved summergreen tree, 12=Boreal broad-leaved summergreen tree, 13=Temperate herbaceous (C3),
%14=Temperate herbaceous (C3)
lai_annmean_phen(6,:,:,1,:)=lai_seib_annmean(:,:,7,:)+lai_seib_annmean(:,:,10,:); %needleleaved evergreen
lai_annmean_phen(6,:,:,2,:)=lai_seib_annmean(:,:,11,:); %needleleaved deciduous
lai_annmean_phen(6,:,:,3,:)=lai_seib_annmean(:,:,12,:); %boreal broadleaved deciduous
lai_annmean_phen(6,:,:,4,:)=lai_seib_annmean(:,:,8,:); %temperate broadleaved evergreen
lai_annmean_phen(6,:,:,5,:)=lai_seib_annmean(:,:,9,:); %temperate broadleaved deciduous
lai_annmean_phen(6,:,:,6,:)=lai_seib_annmean(:,:,1,:)+lai_seib_annmean(:,:,2,:)+lai_seib_annmean(:,:,3,:)+lai_seib_annmean(:,:,4,:)+...
    lai_seib_annmean(:,:,5,:); %tropical broadleaved evergreen
lai_annmean_phen(6,:,:,7,:)=lai_seib_annmean(:,:,6,:); %tropical broadleaved deciduous

%Find phenology class with maximum LAI

lai_annmean_phenmax=NaN(nmod,360,720,2);
lai_annmean_phenmax_jules=NaN(144,192,2);
for xx=1:720
    for yy=1:360
        for tt=1:2
            for nn=1:nmod
                temp=lai_annmean_phen(nn,yy,xx,:,tt);
                inds=find(temp==nanmax(temp)); %This formulation is a catch in case multiple with the same LAI, in which case just take first one (assuming this is a rare occurrence)
                if ~isempty(inds)
                    lai_annmean_phenmax(nn,yy,xx,tt)=inds(1);
                end
            end
            if xx<=192 && yy<=144
                temp=lai_annmean_phen_jules(yy,xx,:,tt);
                inds=find(temp==nanmax(temp));
                if ~isempty(inds)
                    lai_annmean_phenmax_jules(yy,xx,tt)=inds(1);
                end
            end
        end
    end
end
clear xx yy tt nn temp inds

%---

%Load water mask data (from ESA landcover)
load([data_esa,'esa_05_landcover']); %Output from esa_hires_region_mask.m	
oceanm=NaN(720,360);
oceanm(esa_05'>200 & esa_05'<220)=-1;
load([data_esa,'esa_jules_landcover']); %Output from esa_hires_region_mask_jules.m	
oceanm_jules=NaN(192,144);
oceanm_jules(esa_jules'>200 & esa_jules'<220)=-1;

lats=-89.75:0.5:89.75;
lons=-179.75:0.5:179.75;
latincj=180/144;
lats_jules=-90+(latincj/2):latincj:90-(latincj/2);
lons_jules=-179.0625:1.875:179.0625;

%Read year 2000 closed-canopy forest mask derived from Hansen et al. (2013) data (Pugh et al, 2019, Nature Geoscience 12, 730-735)
[fmask,fmask_jules,~,~]=get_closed_can_mask(data_mask);

%Set colorbar limits
cmin=1;
cmax=8;

%Make plots
figure
cmap=colormap(parula(7));
cmap=[0.9 0.9 0.9; cmap];
colormap(cmap)

%First do time period 1
mm=0;
for nn=[1 3 5 7 9 11];
    mm=mm+1;
    if mm==2; continue; end %JULES is a special case because of resolution
    s(nn)=subplot(6,2,nn);
    hold on
    l(nn)=pcolor(lons,lats,oceanm');
    set(l(nn),'linestyle','none')
    p(nn)=pcolor(lons,lats,squeeze(lai_annmean_phenmax(mm,:,:,1)).*fmask');
    set(p(nn),'linestyle','none')
    box on
    set(gca,'XTick',[],'XTickLabel','')
    set(gca,'YTick',[],'YTickLabel','')
    set(gca,'YLim',[-60 80])
    set(gca,'XLim',[-180 180])
    caxis([cmin-1 cmax])
    t(nn)=text(-175,-50,models{mm});
    set(t(nn),'FontSize',12,'FontWeight','Bold')
end
clear nn mm
c1=colorbar;
set(c1,'FontSize',12,'FontWeight','Bold')
set(c1,'Limits',[cmin cmax])
set(c1,'Location','eastoutside')
set(c1,'Ticks',1.5:1:8,'TickLabels',phen_label)

%JULES is a special case because of resolution
s(3)=subplot(6,2,3);
hold on
l(3)=pcolor(lons_jules,lats_jules,oceanm_jules');
set(l(3),'linestyle','none')
p(3)=pcolor(lons_jules,lats_jules,squeeze(lai_annmean_phenmax_jules(:,:,1)).*fmask_jules');
set(p(3),'linestyle','none')
set(gca,'XTick',[],'XTickLabel','')
set(gca,'YTick',[],'YTickLabel','')
set(gca,'YLim',[-60 80])
set(gca,'XLim',[-180 180])
box on
caxis([cmin-1 cmax])
t(3)=text(-175,-50,models{2});
set(t(3),'FontSize',12,'FontWeight','Bold')

%Then do time period 2
mm=0;
for nn=[2 4 6 8 10 12];
    mm=mm+1;
    if mm==2; continue; end %JULES is a special case because of resolution
    s(nn)=subplot(6,2,nn);
    hold on
    l(nn)=pcolor(lons,lats,oceanm');
    set(l(nn),'linestyle','none')
    p(nn)=pcolor(lons,lats,squeeze(lai_annmean_phenmax(mm,:,:,2)).*fmask');
    set(p(nn),'linestyle','none')
    box on
    set(gca,'XTick',[],'XTickLabel','')
    set(gca,'YTick',[],'YTickLabel','')
    set(gca,'YLim',[-60 80])
    set(gca,'XLim',[-180 180])
    caxis([cmin-1 cmax])
    t(nn)=text(-175,-50,models{mm});
    set(t(nn),'FontSize',12,'FontWeight','Bold')
end
clear nn mm

%JULES is a special case because of resolution
s(4)=subplot(6,2,4);
hold on
l(4)=pcolor(lons_jules,lats_jules,oceanm_jules');
set(l(4),'linestyle','none')
p(4)=pcolor(lons_jules,lats_jules,squeeze(lai_annmean_phenmax_jules(:,:,2)).*fmask_jules');
set(p(4),'linestyle','none')
set(gca,'XTick',[],'XTickLabel','')
set(gca,'YTick',[],'YTickLabel','')
set(gca,'YLim',[-60 80])
set(gca,'XLim',[-180 180])
box on
caxis([cmin-1 cmax])
t(4)=text(-175,-50,models{2});
set(t(4),'FontSize',12,'FontWeight','Bold')

set(s(1),'Position',[0.05 0.85 0.44 0.12])
set(s(2),'Position',[0.51 0.85 0.44 0.12])
set(s(3),'Position',[0.05 0.72 0.44 0.12])
set(s(4),'Position',[0.51 0.72 0.44 0.12])
set(s(5),'Position',[0.05 0.59 0.44 0.12])
set(s(6),'Position',[0.51 0.59 0.44 0.12])
set(s(7),'Position',[0.05 0.46 0.44 0.12])
set(s(8),'Position',[0.51 0.46 0.44 0.12])
set(s(9),'Position',[0.05 0.33 0.44 0.12])
set(s(10),'Position',[0.51 0.33 0.44 0.12])
set(s(11),'Position',[0.05 0.2 0.44 0.12])
set(s(12),'Position',[0.51 0.2 0.44 0.12])

set(c1,'Position',[0.51 0.03 0.04 0.16])
