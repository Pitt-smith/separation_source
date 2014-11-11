%Auteur : Selim RABOUDI
%dans cette version, j'utilise un algo de fenetre glissante pour creer plus
%de vecteurs enfants des signaux melangés, et ainsi mieux estimer la
%matrice d'intercovariance

clear all;
close all;
clc;

tic;
%load les donnees
load '../../data/SignauxMelange.mat';
load '../../data/SignauxReference.mat';

N = length(Melange);

x1 = Melange(1,:)';
x2 = Melange(2,:)';

R11 = x1*x1';
R22 = x2*x2';
R12 = x1*x2';

%je divise chaque xi en plusieurs vecteurs, de taille moindre, pour trois
%raisons :
% 1 - afin de pouvoir estimer l'esperance de chaque x*x' 
%     (en faisant une moyenne lineaire de chaque x*x' pour x de taille plus petite)
% 2 - afin d'eviter de calculer l'intercorrelation entre des echantillons
%     tres eloignes dans le temps, qui seront proches de 0 et alourdirons les
%     calculs sans benefices.
% 3 - l'"off" des matrices d'intercovariance est proche de 0 pour
%     N = 13129. (division par 1/(N*(N-1))). Cela fausse le calcul de
%     Achapeau


%taille des vecteurs enfants ( = Taille de la fenetre glissante)
SizeChild = 40

%nombre de vecteurs enfants
nbX = N - SizeChild + 1;

%matrices contenant les vecteurs enfants
X1 = zeros(SizeChild,nbX);
X2 = zeros(SizeChild,nbX);

%remplissage des matrices contenant les enfants de x1 et x2
for k=1:nbX
    X1(:,k) = x1(k:k+SizeChild-1);
    X2(:,k) = x2(k:k+SizeChild-1); 
end

%matrice contenant les matrices d'intercovariance pour chaque vecteur
%enfant de x1
Xinter11 = zeros(SizeChild,SizeChild);
Xinter22 = zeros(SizeChild,SizeChild);
Xinter12 = zeros(SizeChild,SizeChild);

for k=1:nbX
    Xinter11 = Xinter11 + X1(:,k)*X1(:,k)';
    Xinter22 = Xinter22 + X2(:,k)*X2(:,k)';
    Xinter12 = Xinter12 + X1(:,k)*X2(:,k)';
end;

%l'esperance est approximee par la moyenne lineaire
R11 = 1/(nbX)*Xinter11;
R22 = 1/(nbX)*Xinter22;
R12 = 1/(nbX)*Xinter12;

%calcul des valeurs F et T
T1 = trace(R11);
T2 = trace(R22);
T12 = trace(R12);

F1 = 1/(SizeChild*(SizeChild-1))*(sum(sum(R11)) - T1);
F2 = 1/(SizeChild*(SizeChild-1))*(sum(sum(R22)) - T2);
F12 = 1/(SizeChild*(SizeChild-1))*(sum(sum(R12)) - T12);

%calcul des autres constantes
alpha = 2*F12*T12 - (F1*T2 + F2*T1);
beta = 2*((T12)^2 - T1*T2);
gamma2 = (F1*T2 - F2*T1)^2 + 4*(F12*T2-T12*F2)*(F12*T1 - T12*F1);
gamma = sqrt(gamma2);
d1 = alpha - gamma;
d2 = alpha + gamma;

%finalement Achap

Achap = [beta*F1-T1*d1 , beta*F12-T12*d2 ; beta*F12-T12*d1 , beta*F2-T2*d2];

invAChap = inv(Achap);

Sestime = invAChap*Melange;

s1Estime = Sestime(1,:);
s2Estime = Sestime(2,:);

toc;

%calcul d'erreur :
sOiseau = Signal(1,:);
sGong = Signal(2,:);


errOiseau = 10*log10(1 - (s1Estime*sOiseau'/(norm(s1Estime)*norm(sOiseau)))^2);
errGong = 10*log10(1 - (s2Estime*sGong'/(norm(s2Estime)*norm(sGong)))^2);
sprintf('erreur signal oiseau : %.1f',errOiseau);
sprintf('erreur signal gong : %.1f',errGong);
errOiseau
errGong
errOiseau + errGong 
sound(s1Estime, 8000);
pause(3);
sound(s2Estime, 8000);

