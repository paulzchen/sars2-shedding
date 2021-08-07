%This script determines the sensitivity and specificity of prognostic 
%thresholds based URT or LRT shedding. To do so, it uses the Weibull
%cumulative density function (cdf) and then calculates true positive and
%true negative rates.

clear all
close all
clc 

x = 0:0.05:15; %rVLs from 0 to 15 log10 copies/ml

%Weibull distributions
A_NS = 5.36532; %Scale parameter for nonsevere
B_NS = 2.59011; %Shape parameter for nonsevere
A_S = 5.84757; %Scale parameter for severe
B_S = 2.98665; %Shape parameter for severe

W_NS = (B_NS/A_NS).*( (x/A_NS).^(B_NS-1) ).* ...
    exp(-((x/A_NS).^B_NS)); %Weibull pdf for nonsevere
W_S = (B_S/A_S).*( (x/A_S).^(B_S-1) ).* ...
    exp(-((x/A_S).^B_S)); %Weibull pdf for severe

cdf_NS = 1 - exp(-(x./A_NS).^B_NS); %Weibull cdf for nonsevere
cdf_S = 1 - exp(-(x./A_S).^B_S); %Weibull cdf for severe

%Initialize arrays
TP = zeros(length(x),1); %true positive
TN = zeros(length(x),1); %true negative
FP = zeros(length(x),1); %false positive
FN = zeros(length(x),1); %false negative
Sens = zeros(length(x),1); %sensitivity (true positive rate)
Spec = zeros(length(x),1); %specificity (true negative rate)
overall = zeros(length(x),1); %specificity + sensitivity

for i = 2 : length(x)
    TP(i) = 1 - cdf_S(i);
    FP(i) = 1 - cdf_NS(i);
    TN(i) = cdf_NS(i-1);
    FN(i) = cdf_S(i-1);
    Sens(i) = TP(i) / (TP(i) + FP(i));
    Spec(i) = TN(i) / (TN(i) + FN(i));
    overall(i) = Sens(i) + Spec(i);
end

output = [x(3:end-2)' Sens(3:end-2) Spec(3:end-2)];
output2 = [x(3:end-2)' TP(3:end-2) TN(3:end-2) FP(3:end-2) FN(3:end-2)];
[m,index] = max(overall);
x_max = x(index) %rVL threshold for highest overall accuracy
                 %This may not be the best threshold for prognostication,
                 %based on the relative specificity and sensitivity

%plots in case you want to visualize
subplot(2,2,1)
plot(x,cdf_NS)
hold on
plot(x,cdf_S)
xlabel('rVL (log_1_0 copies/ml)')
ylabel('%')
legend('nonsevere','severe')
hold off

subplot(2,2,2)
plot(x,TP)
hold on
plot(x,TN)
plot(x,FP)
plot(x,FN)
xlabel('rVL (log_1_0 copies/ml)')
ylabel('%')
legend('true pos.','true neg.','false pos.','false neg.')
hold off

subplot(2,2,3)
plot(x,Sens)
hold on
plot(x,Spec)
xlabel('rVL (log_1_0 copies/ml)')
ylabel('%')
legend('sensitivity','specificity')
hold off