%Script to make combined bar plot of standard deviation of turnover fluxes for both space and time
%
%Dependencies
%-global_turn_frac_sd.m
%-global_turn_frac_sd_time.m
%And dependencies within dependencies
%-get_stocks_fluxes.m
%-get_closed_can_mask.m
%-get_forest_type.m
%
%T. Pugh
%23.02.18

makeplots=false;

[lrrflux_mask_std,mflux_mask_std,lrrflux_mask_std_biom,mflux_mask_std_biom...
    ,~,~]=global_turn_frac_sd(makeplots);

[lrrflux_mask_std_mean,mflux_mask_std_mean,lrrflux_mask_std_mean_biom,...
    mflux_mask_std_mean_biom,phen_label,nphen]=...
    global_turn_frac_sd_time(makeplots);

models={'CABLE-POP','JULES','LPJ-GUESS','LPJmL','ORCHIDEE','SEIB-DGVM'};

figure
subplot(2,1,1)
bb=bar([mflux_mask_std,lrrflux_mask_std]);
set(gca,'XTickLabel',models)
for nn=1:2
    set(bb(nn),'LineStyle','none')
end
clear nn
legend('Mortality','Phenology')
ylabel('\sigma_{space} (kg C m^{-2} a^{-1})')
set(gca,'YLim',[0 0.4])
t1=text(0.2,0.35,'(a)');
set(t1,'FontSize',12,'FontWeight','Bold')

subplot(2,1,2)
bb=bar([mflux_mask_std_mean,lrrflux_mask_std_mean]);
set(gca,'XTickLabel',models)
for nn=1:2
    set(bb(nn),'LineStyle','none')
end
clear nn
ylabel('\sigma_{time} (kg C m^{-2} a^{-1})')
set(gca,'YLim',[0 0.4])
t2=text(0.2,0.35,'(b)');
set(t2,'FontSize',12,'FontWeight','Bold')
