
%-------------------------------------------------------------------------
%
%  S�paration de source en utilisant l'algorithme ICA
%
%  Matrice de m�lange : A=[1,1; -1,2]; 
%
%  M�lange des signaux en ulisant une matrice pleine A=[1,1; -1,2];  
%  Les signaux m�lang�s sont sauvegard� dans des fichiers :
%     - Melange.mat pour Matlab 
%     - signalMelange1.txt & signalMelange2.txt pour le langage C
%  Le format des donn�es �crites est 'double'
%
%-------------------------------------------------------------------------
clear all;
close all;
clc;


%% lecture des signaux m�lang�s 
%------------------------------

load SignauxMelange.mat 
nombreEchantillons = length(Melange);

%% Calcul des valeurs moyennes des m�langes et Centrage des m�langes 
%-------------------------------------------------------------------

Melange(1,:) = Melange(1,:) - (1/nombreEchantillons)*sum(Melange(1,:));
Melange(2,:) = Melange(2,:) - (1/nombreEchantillons)*sum(Melange(2,:));

%% Autocorrelation normalis�e
%-----------------------------

Rx = Melange*Melange'/nombreEchantillons;

%% Blanchiment 
%-------------

[V,D]=eig(Rx);                 % Valeurs propres et vecteurs propres

Xw = V*(diag(diag(D).^(-0.5)))*V'*Melange;   
 
%% Recherche d'une seule composante ind�pendante
%----------------------------------------------
% Vecteur de s�paration al�atoirement initialis�  
    %W=[rand(),rand()];
%  Vecteur de s�paration initialis� de fa�on � comparer les r�sultats Matlab avec le code C  
W = [1.5,-2];
W = W/norm(W);
    
% Initialisation des signaux estim�s
y1=W*Xw;

% Nombre d'iterations 
Nbiterations  = 100;

for i=1:Nbiterations
	% Estimation de l'entropie du signal
	Y = tanh(y1);
    Z = 1./((cosh(y1)).^2);  
	W(1) =  (1/nombreEchantillons)*(sum(Y.*Xw(1,:))-(sum(Z)*W(1)) );
	W(2) =  (1/nombreEchantillons)*(sum(Y.*Xw(2,:))-(sum(Z)*W(2)) );
    W=W/norm(W);
    y1 = W*Xw;
end;
% Recherche de la seconde ligne.

Z = [-W(2), W(1) ];

B=[W ; Z] ;

%% Calcul des signaux estim�s 
%----------------------------

signalEstime = B*Xw;


%% lecture des signaux r�f�rences 
%--------------------------------

load SignauxReference.mat 
nombreEchantillons = min(length(signalEstime),length(Signal));

if nombreEchantillons < length(signalEstime)
    signalEstime = signalEstime(:,1:nombreEchantillons);
else
    Signal = Signal(:,1:nombreEchantillons);
end

%% Calcul des erreurs
%-----------------------

 Erreur1 = 10*log10(1-(Signal(1,:)*signalEstime(1,:)'/(norm(Signal(1,:))*norm(signalEstime(1,:))))^2)
 Erreur2 = 10*log10(1-(Signal(2,:)*signalEstime(2,:)'/(norm(Signal(2,:))*norm(signalEstime(2,:))))^2)


%% Affichage des histogrammes des diff�rents signaux m�lang�s
%------------------------------------------------------------

figure ;hist(Melange(1,:),50); drawnow;grid on;
%figure ;hist(Signal(1,:),50); drawnow;grid on;
figure ;hist(Melange(2,:),50); drawnow;grid on;
%figure ;hist(Signal(2,:),50); drawnow;grid on;

%% Affichage des signaux m�lang�s centr�s
%----------------------------------------
figure()
subplot(2, 1, 1)
plot (Melange(1,:),'r')
grid on
title('Signal m�lang� 1 centr�')
ylabel('Amplitude')
maximum = max(abs(Melange(1,:)));
axis([0 nombreEchantillons -1.1*maximum 1.1*maximum])
subplot(2, 1, 2)
plot (Melange(2,:),'r')
grid on
title('Signal m�lang� 2 centr�')
ylabel('Amplitude')
maximum = max(abs(Melange(2,:)));
axis([0 nombreEchantillons -1.1*maximum 1.1*maximum])
 
%% Affichage des signaux estim�s
%--------------------------------

figure()
subplot(2, 1, 1)
plot (signalEstime(1,:),'r')
grid on
title('Signal estim� 1')
ylabel('Amplitude')
maximum = max(abs(signalEstime(1,:)));
axis([0 nombreEchantillons -1.1*maximum 1.1*maximum])
subplot(2, 1, 2)
plot (signalEstime(2,:),'r')
grid on
title('Signal estim� 2')
ylabel('Amplitude')
maximum = max(abs(signalEstime(2,:)));
axis([0 nombreEchantillons -1.1*maximum 1.1*maximum])
 
%% Affichage des signaux de r�f�rence
%------------------------------------

figure()
subplot(2, 1, 1)
plot (Signal(1,:),'r')
grid on
title('Signal r�f�rence 1')
ylabel('Amplitude')
maximum = max(abs(Signal(1,:)));
axis([0 nombreEchantillons -1.1*maximum 1.1*maximum])
subplot(2, 1, 2)
plot (Signal(2,:),'r')
grid on
title('Signal r�f�rence 2')
ylabel('Amplitude')
maximum = max(abs(Signal(2,:)));
axis([0 nombreEchantillons -1.1*maximum 1.1*maximum])
 

