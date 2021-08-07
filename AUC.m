%This script determines the overlap in the area under the curve (AUC) for
%Weibull distributions (from 1st to 99th case percentiles) for URT and LRT
%shedding.

%This is estimates how accurate URT and LRT rVL may potentially be as a
%prognistic indicator of COVID-19 severity.

clear all
close all
clc 

x = 0:0.05:15; %rVLs from 0 to 15 log10 copies/ml

%Weibull distributions
A_NS = 7.06457; %Scale parameter for nonsevere
B_NS = 4.86281; %Shape parameter for nonsevere
A_S = 9.47998; %Scale parameter for severe
B_S = 5.81418; %Shape parameter for severe

W_NS = (B_NS/A_NS).*( (x/A_NS).^(B_NS-1) ).* ...
    exp(-((x/A_NS).^B_NS)); %Weibull pdf for nonsevere
W_S = (B_S/A_S).*( (x/A_S).^(B_S-1) ).* ...
    exp(-((x/A_S).^B_S)); %Weibull pdf for severe

y_d = [W_S(W_S<W_NS) W_NS(W_NS<W_S)]; %Overlapped portions
AUC_overlap  = trapz(y_d) %AUC of overlapped portions
area_int_NS  = trapz(W_NS) %AUC of nonsevere
area_int_S  = trapz(W_S) %AUC of severe

P_nonoverlap = (area_int_NS + area_int_S - 2*AUC_overlap)/...
    (area_int_NS + area_int_S)*100 %Percentage of distribution where no overlap
P_overlap = AUC_overlap / area_int_NS * 100 %Percentage of distribution overlap

%plots in case you want to visualize
plot(x,W_NS)
hold on
plot(x,W_S)
plot(x(2:length(y_d)+1),y_d,'k-o')
legend('nonsevere','severe','overlap')
xlabel('rVL (log_1_0 copies/ml)')
ylabel('Density')