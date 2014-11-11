clear all;
close all;
clc;

%load the data
load '../data/SignauxMelange.mat';
load '../data/SignauxReference.mat';


%% performance de ICA, SOBI et SOBI avec correction du bruit en fonction du SNR
%function [errICA, errSOBI, errSOBINoise] = 
%function_noisy_ICA_SOBI(Melange, Signal, SNR)

SNRdebut = 0; 
SNRfin = 5;
SNRpas = 0.01;
vectorSNR = SNRdebut:SNRpas:SNRfin;
nbSNR = length(vectorSNR);

results = zeros(3, nbSNR);

indSNR = 1;
for k=SNRdebut:SNRpas:SNRfin
    res = function_noisy_ICA_SOBI(Melange, Signal, k);
    results(1,indSNR) = res(1);
    results(2,indSNR) = res(2);
    results(3,indSNR) = res(3);
    indSNR = indSNR +1;
end;

impactSOBINoise = 100*((results(3,:) - results(2,:))./results(2,:));
pos = length(find(impactSOBINoise>=0));
neg = length(find(impactSOBINoise<0));
percentagePos = 100*(pos/(pos+neg))
mean(impactSOBINoise(find(impactSOBINoise>=0)))
mean(impactSOBINoise(find(impactSOBINoise<0)))


%plot des courbes EQM = f(SNR) pour ICA, SOBI, SOBI+corrNoise :
cc = hsv(3);
str = cell(1,3);
h = zeros(1,3);

figure;
hold on;
for k=1:3
    h(k) = plot(vectorSNR, results(k,:),'color',cc(k,:));
end;
str{1} = sprintf('ICA');
str{2} = sprintf('SOBI');
str{3} = sprintf('SOBI+corrNoise');

legend(h(1:3),str{1:3});

hold off;