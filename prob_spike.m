% Routine to calculate the probability of spiking from a set ofd data
clear all
close all
clc

tic

data_matrix = load('190_60_47_13.dat');
%data_matrix = data_matrix.m_;
% Number of row and columns in data matrix
[r, c] = size(data_matrix);

% Finding index of each peak

% Store number of peaks in each experiment
Npeaks = [];

% Store peak index
index = [];

% Store thrshold values
thr_values = [];


% Peaks above -10 mV.
[Npeaks, index] = findNpeaks(data_matrix, -10);


%%
% data_matrix2 is a cell each line of the cell will store the potencial
% values below the found threhold for a particular expriment
data_matrix2 = [];

t = 0:0.05:10;

% figure(1);
% subplot(1,2,1);
% hold on;
count = 0;
for i = 1:Npeaks-1
    % Separate each peak, to find the threshold
    % P is one peak of a given experimente
    % P is choosed from the peak to 3ms ago
    P = data_matrix(index(i)-60:index(i));
    % dP/dt
    P1 = diff(P);
    % d?P/dt?
    P2 = diff(P1);
    % Method VI
    Kp = P2.*(1+(P1(1:end-1).^2)).^(-3/2);
    % Find the max Kp and its index
    [max_h, aux2] = max(Kp);
    thr_values(i) = P(aux2);
    aux = data_matrix(index(i)-60:index(i+1)-59);
    aux( aux > thr_values(i) ) = NaN;
    data_matrix2 = [data_matrix2 ; aux];
    %Plot the separated peaks, with the found threshold
%     if mod(count, 10000) == 0
%         figure(1);
%         plot(t,data_matrix(index(i)-60:index(i)+140),'black');
%         hold on;
%         plot(t(60-60+aux2), P(aux2),'Ored');
%     end
%      count = count + 1;
end

% ylabel('Membrane Potencial [mV]')
% legend('Action Potential','Threshold values')
% % axis([5 10 -80 60])
% % axis([25 45 0 -20])
% print('potencial_thr_values','-dpng','-r600')

%%
bl = 0.5;            % Length of the bins

% Dividing into bins
v_m = -76:bl:1;
[n1, x1] = hist(data_matrix2,v_m);
[n2, x2] = hist(thr_values,v_m);
phi_v = n2 ./ n1;
figure
% plot(v_m, phi_v)
plot(v_m(~isnan(phi_v))-10,phi_v(~isnan(phi_v))); %cesar
% axis([-50 -44.6 0 1.1])  % simul
% axis([-60 -40 0 1]) %cesar
ylabel('Probability')
xlabel('Membrane Potencial [mV]')
% print('graf_prob_disp','-dpng','-r600')

%%
figure
[ax, h1, h2] = plotyy(x1,n1,x2,n2,@stairs);
% set([h1,h2],'plottype', 'stairs')
% hold(ax)
% h1 = stairs(x1,n1,'Color',[0 0 1]);
% ax(1) = gca;
% axpos = get(ax(1),'position');
% axes('position',axpos);
% h2 = stairs(x2,n2,'color',[0 1 0]);
% ax(2) = gca;
% set(ax(2),'Box','off','Visible','off');
% ylabel('Counts')
% xlabel('Membrane Potential [mV]');
legend([h1, h2],'V_{m}', 'V_{th}');

% print('hist_prob_disp','-dpng','-r600')

toc
