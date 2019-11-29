%Plot the leaf cost versus leaf longevity for those models that provide this information
%
%T.Pugh
%13.04.19

figure
hold on

%---
%LPJ-GUESS
%Eq. C21 from Smith et al. (2014, Biogeosciences 11, 2027-2054)

%Needleleaf
aleaf=0.5:0.1:3; %Leaf longevity
SLA = 0.2 * 10.^(2.41 - 0.38 * log10(12.0.* aleaf)); %Specific leaf area

plot(aleaf,(1./(SLA))./aleaf,'k','linewidth',2)
hold on

%Broadleaf
aleaf=0.5:0.1:3;
SLA = 0.2 * 10.^(2.29 - 0.4 * log10(12.0.* aleaf));

plot(aleaf,(1./SLA)./aleaf,'k--','linewidth',2)

%---
%LPJmL
%Eq. 6 in Sitch et al. (2003, Global Change Biology 9, 161-185)

aleaf=0.5:0.1:3;
SLA=2e-4.*(exp(6.15)./((12.*aleaf).^0.46));
SLA=SLA*1000; %Convert from m2 gC-1 to m2 kgC-1

plot(aleaf,(1./SLA)./aleaf,'linewidth',2)

%---
%SEIB-DGVM
%SLA taken from Table B4 in Sato et al. (2007, Ecological Modelling 200, 279-307), units m2 / g Dry-Mass (Hisashi Sato, pers comm.)

SLA_raw=[0.01 0.013 0.004 0.007 0.015 0.004 0.015 0.016]; %Order as in Table B4
aleaf_raw=[1.69 0.63 4.55 2.63 0.46 4.55 0.25 NaN];
SLA_raw=SLA_raw*2; %Convert from m2 gDM-1 to m2 gC-1
SLA_raw=SLA_raw*1000; %Convert from m2 gC-1 to m2 kgC-1
[aleaf,I]=sort(aleaf_raw);
SLA=SLA_raw(I);

plot(aleaf,(1./(SLA))./aleaf,'linewidth',2)

%Plot labelling
ylabel('Leaf cost (kg C m^{-2} y^{-1})')
xlabel('Leaf longevity (y)')
legend('LPJ-GUESS (needleleaf)','LPJ-GUESS (broadleaf)','LPJmL','SEIB-DGVM')