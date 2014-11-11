%function [errVect] = SOBI_functionv2(SizeChild, Melange, Signal, varNoise)
%errVect = [SizeChild ; errOiseau ; errGong ; errOiseau + errGong];

%function [errOiseau, errGong] = ICA_function(Melange, Signal)



clear all;
close all;
clc;

%load the data
load '../data/SignauxMelange.mat';
load '../data/SignauxReference.mat';


%% ajout de bruit
SNR = 40;
sigmaNoise = sqrt(var(Melange(1,:))*10^(-SNR/10));
noise = sigmaNoise*randn(size(Melange));
Melange = Melange + noise;


%% performance de SOBI et ICA avec ajout du bruit
display('performance de ICA avec bruit additif:')
err = ICA_function(Melange, Signal); 
err

display('performance de SOBI avec bruit additif, sans correction de bruit:')
errVect = SOBI_functionv2(40, Melange, Signal,0);
errVect(2)
errVect(3)

errVect = SOBI_functionv2(40, Melange, Signal,var(noise(1,:)));
display('performance de SOBI avec correction de bruit:')
errVect(2)
errVect(3)

%% tracage des courbes EQM = f(SNR) avec plusieurs tailles de fenetres N differentes
% 
% load '../data/SignauxMelange.mat';
% 
% %on calcule les valeurs pour N dans [10,100] (par pas de 10)
% 
% %definition de la plage de SNR
% SNRdebut = 0;
% SNRfin = 50;
% SNRpas = 10;
% nbMesSNR = length(SNRdebut:SNRpas:SNRfin);
% 
% Ndebut = 10; 
% Nfin = 30;
% Npas = 10;
% nbMesN = length(Ndebut:Npas:Nfin);
% matResults = zeros(nbMesN,nbMesSNR);
% 
% 
% indSNR = 1;
% display(sprintf('calcul des EQM pour SNR %d:%d:%d  et N %d:%d:%d', SNRdebut, SNRpas, SNRfin, Ndebut, Npas, Nfin));
% for SNR=SNRdebut:SNRpas:SNRfin
%     indN = 1;
%     for N=Ndebut:Npas:Nfin
%         matResults(indN,indSNR) = SOBI_SNR(Melange, Signal, SNR, N);
%         %display(sprintf('SNR = %d', SNR));
%         indN = indN +1;
%     end;
%     indSNR = indSNR +1;
% end;
% 
% %plot des courbes EQM = f(SNR) pour differentes valeur de la fenetre :
% cc = hsv(nbMesN);
% str = cell(1,nbMesN);
% h = zeros(1,nbMesN);
% SNR = SNRdebut:SNRpas:SNRfin;
% 
% figure;
% hold on;
% for k=1:nbMesN
%     h(k) = plot(SNR, matResults(k,:),'color',cc(k,:));
%     str{k} = sprintf('N = %d',Ndebut+(k-1)*Npas);
% end;
% 
% legend(h(1:nbMesN),str{1:nbMesN});
% 
% hold off;





