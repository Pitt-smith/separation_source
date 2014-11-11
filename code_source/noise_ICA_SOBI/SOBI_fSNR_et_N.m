clear all;
close all;
clc;

%load the data
load '../data/SignauxMelange.mat';
load '../data/SignauxReference.mat';


%% tracage des courbes EQM = f(SNR) avec plusieurs tailles de fenetres N differentes

load '../data/SignauxMelange.mat';

%on calcule les valeurs pour N dans [10,100] (par pas de 10)

%definition de la plage de SNR
SNRdebut = 0;
SNRfin = 50;
SNRpas = 1;
nbMesSNR = length(SNRdebut:SNRpas:SNRfin);

Ndebut = 1; 
Nfin = 10;
Npas = 1;
nbMesN = length(Ndebut:Npas:Nfin);
matResults = zeros(nbMesN,nbMesSNR);


indSNR = 1;
display(sprintf('calcul des EQM pour SNR %d:%d:%d  et N %d:%d:%d', SNRdebut, SNRpas, SNRfin, Ndebut, Npas, Nfin));
for SNR=SNRdebut:SNRpas:SNRfin
    indN = 1;
    for N=Ndebut:Npas:Nfin
        matResults(indN,indSNR) = SOBI_SNR(Melange, Signal, SNR, N);
        %display(sprintf('SNR = %d', SNR));
        indN = indN +1;
    end;
    indSNR = indSNR +1;
end;

%plot des courbes EQM = f(SNR) pour differentes valeur de la fenetre :
cc = hsv(nbMesN);
str = cell(1,nbMesN);
h = zeros(1,nbMesN);
SNR = SNRdebut:SNRpas:SNRfin;

figure;
hold on;
for k=1:nbMesN
    h(k) = plot(SNR, matResults(k,:),'color',cc(k,:));
    str{k} = sprintf('N = %d',Ndebut+(k-1)*Npas);
end;

legend(h(1:nbMesN),str{1:nbMesN});

hold off;





