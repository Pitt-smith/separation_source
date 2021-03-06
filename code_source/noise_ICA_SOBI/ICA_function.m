% implementation d'ICA en fonction
% AUTEUR : Selim RABOUDI

function [err] = ICA_function(Melange, Signal)


%load '../data/SignauxMelange.mat';
%load '../data/SignauxReference.mat';


%ETAPE 1 : blanchiment de Melange
x = Melange;
[l,n] = size(x);    
xCov = (1/n)*(x*x'); %matrice de covariance de x (esperance approximee par la moyenne lineaire)
[E,D] = eig(xCov);
xTild = E*inv(sqrtm(D))*E'*x;
xTildCov = (1/n)*(xTild*xTild'); %la matrice de covariance de xTild est bien egal a In

x = xTild; %x est maintenant blanchi
%%ETAPE 2 : Realisation de l'algorithme ICA
w = rand(1,2); % w initial aleatoire
%mean(bsxfun(@times,x,tanh(w*x)),2)
%mean((1-(tanh(w*x)).^2))
normW = mean((w*x).^2);
for k=1:100
    % expliquer role de bsxfun(.)
    w = mean(bsxfun(@times,x,tanh(w*x)),2)' - mean((1-(tanh(w*x)).^2))*w;
    normW = mean((w*x).^2);
    w = w/normW;

end;

sFiltre = w*x;

%orthogonalisation de Gram-Schmidt (en dimension 2, immediat)
w2 = [-w(2) w(1)];
sFiltre2 = w2*x;

%sound(sFiltre,8000)
sOiseau = Signal(1,:);
sGong = Signal(2,:);

errOiseau1 = 10*log10(1 - (sFiltre*sOiseau'/(norm(sFiltre)*norm(sOiseau)))^2);
errGong1 = 10*log10(1 - (sFiltre*sGong'/(norm(sFiltre)*norm(sGong)))^2);
errOiseau2 = 10*log10(1 - (sFiltre2*sOiseau'/(norm(sFiltre2)*norm(sOiseau)))^2);
errGong2 = 10*log10(1 - (sFiltre2*sGong'/(norm(sFiltre2)*norm(sGong)))^2);
err = [min(errOiseau1, errGong1) min(errOiseau2, errGong2)];
%sprintf('erreur signal oiseau : %.1f',errOiseau);
%sprintf('erreur signal gong : %.1f',errGong);

% sound(sFiltre, 8000);
% pause(3)
% sound(sFiltre2,8000);

