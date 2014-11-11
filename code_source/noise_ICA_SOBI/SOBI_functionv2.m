%fonction implementant l'algorithme de SOBI
%Auteur : Selim RABOUDI

function [errVect] = SOBI_functionv2(SizeChild, Melange, Signal, varNoise)

N = length(Melange);

x1 = Melange(1,:)';
x2 = Melange(2,:)';

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
%SizeChild = 100;

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
R11 = zeros(SizeChild,SizeChild);
R22 = zeros(SizeChild,SizeChild);
R12 = zeros(SizeChild,SizeChild);

for k=1:nbX
    R11 = R11 + X1(:,k)*X1(:,k)';
    R22 = R22 + X2(:,k)*X2(:,k)';
    R12 = R12 + X1(:,k)*X2(:,k)';
end;

%l'esperance est approximee par la moyenne lineaire
R11 = 1/(nbX)*R11;
R22 = 1/(nbX)*R22;
R12 = 1/(nbX)*R12;

%calcul des valeurs F et T
T1 = trace(R11);
T2 = trace(R22);
T12 = trace(R12);

F1 = 1/(SizeChild*(SizeChild-1))*(sum(sum(R11)) - T1);
F2 = 1/(SizeChild*(SizeChild-1))*(sum(sum(R22)) - T2);
F12 = 1/(SizeChild*(SizeChild-1))*(sum(sum(R12)) - T12);

%calcul des autres constantes
alpha = 2*F12*T12 - (F1*(T2-varNoise) + F2*(T1-varNoise));
beta = 2*((T12)^2 - (T1-varNoise)*(T2-varNoise));
gamma2 = (F1*(T2-varNoise) - F2*(T1-varNoise))^2 + 4*(F12*(T2-varNoise)-T12*F2)*(F12*(T1-varNoise) - T12*F1);
gamma = sqrt(gamma2);
d1 = alpha - gamma;
d2 = alpha + gamma;

%finalement Achap

Achap = [beta*F1-(T1-varNoise)*d1 , beta*F12-T12*d2 ; beta*F12-T12*d1 , beta*F2-(T2-varNoise)*d2];

invAChap = inv(Achap);

Sestime = invAChap*Melange;

s1Estime = Sestime(1,:);
s2Estime = Sestime(2,:);

%calcul d'erreur :
sOiseau = Signal(1,:);
sGong = Signal(2,:);


errOiseau = 10*log10(1 - (s1Estime*sOiseau'/(norm(s1Estime)*norm(sOiseau)))^2);
errGong = 10*log10(1 - (s2Estime*sGong'/(norm(s2Estime)*norm(sGong)))^2);
%sprintf('erreur signal oiseau : %.1f',errOiseau)
%sprintf('erreur signal gong : %.1f',errGong)
errVect = [SizeChild ; errOiseau ; errGong ; errOiseau + errGong];



