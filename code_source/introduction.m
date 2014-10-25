%% Auteur : Selim RABOUDI
%% introduction du TP de séparation aveugle de source.

y1 = 1/(2*sqrt(3))*rand(1,1000);
y2 = 1/(2*sqrt(3))*rand(1,1000);
Y = [y1 ; y2];
figure;
plot(y1,y2,'+');

%% Generation des x

A = [1 1; -1 2]
X = A*Y;
x1 = X(1,:);
x2 = X(2,:);

figure;
plot(x1,x2,'+');