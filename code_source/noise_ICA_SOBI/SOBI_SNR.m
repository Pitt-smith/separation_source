function [errGlobal] = SOBI_SNR(Melange, Signal, SNR, N)

sigmaNoise = sqrt(var(Melange(1,:))*10^(-SNR/10));
noise = sigmaNoise*randn(size(Melange));
Mel = Melange + noise;
resu = SOBI_functionv2(N, Mel, Signal,0);
errGlobal = resu(4);
