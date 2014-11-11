clear all
close all;
clc;

%load the data
load '../../data/SignauxMelange.mat';
load '../../data/SignauxReference.mat';

fprintf('length melange : %d\n',length(Melange));
fprintf('length reference : %d\n',length(Signal));

sizeStep = 10;
%nbIter = round(length(Melange)/sizeStep);
nbIter = 10;

Merr = zeros(4,nbIter);
fprintf('beginning of the iteration...\n\n');

for j=1:nbIter
    fprintf('calculating for SizeW=%d....\n',j*sizeStep);
    Merr(:,j) = SOBI_functionv2(j*sizeStep,Melange,Signal);
end;

figure;
subplot(1,2,1);
plot(Merr(1,:),Merr(2,:),'--bs');
title('Signal 1 (Oiseau) : log(EQMN) = f(sizeFenetre)');


subplot(1,2,2);
plot(Merr(1,:),Merr(3,:), '--bs');
title('Signal 2 (Gong) : log(EQMN) = f(sizeFenetre)');

figure;
plot(Merr(1,:),Merr(4,:),'--bs');
title('log(EQMN1) + log(EQMN2) = f(sizeFenetre)');    