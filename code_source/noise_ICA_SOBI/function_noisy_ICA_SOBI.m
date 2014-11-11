

function [vectErr] = function_noisy_ICA_SOBI(Melange, Signal, SNR)
%% ajout de bruit
sigmaNoise = sqrt(var(Melange(1,:))*10^(-SNR/10));
noise = sigmaNoise*randn(size(Melange));
Mel = Melange + noise;


%% performance de SOBI et ICA avec ajout du bruit
err = ICA_function(Mel, Signal); 
errICA = err(1) + err(2);

errVect = SOBI_functionv2(40, Mel, Signal,0);
errSOBI = errVect(2)+errVect(3);


errVect = SOBI_functionv2(40, Mel, Signal,var(noise(1,:)));
errSOBINoise = errVect(2)+errVect(3);

vectErr = [errICA, errSOBI, errSOBINoise];

