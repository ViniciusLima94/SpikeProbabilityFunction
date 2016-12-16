% Routine to calculate the probability of spiking from a set ofd data
clear all
close all
clc

% Choose whether to read, from .dat or .mat. Use .dat
load_from = 'mat_file'; % Or mat_file
curr = [];
if strcmp(load_from, 'mat_file')
    
    cd data
    
    files = dir('*.mat'); % Read all .mat files in folder
    data_matrix = [];     % Stores all data
    
    for i = 1:length(files)
        load(files(i).name) % Load ith file
        
        aa = who('Trace*');
        b = aa{1};
        c = strfind(b,'_');
        
        indsweep1 = b(c(1)+1:c(2)-1);
        indsweep2 = b(c(2)+1:c(3)-1);
        indsweep3 = b(c(3)+1:c(4)-1);
        indsweep4 = b(c(4)+1:end);
        
        eval(['t = Trace_' indsweep1 '_' indsweep2 '_1_' indsweep4 '(:,1);']);
        raw_signal = [];
        smooth_signal = [];
        for j = 1:length(aa)
            eval(['raw_signal = [raw_signal Trace_' indsweep1 '_' indsweep2 '_' num2str(j) '_' indsweep4 '(:,2)];'])
            smooth_signal = [smooth_signal smooth(raw_signal(:,j),40)];
            curr(end+1) = (j-1)*20;
        end
        
        raw_signal = 1000*raw_signal;
        smooth_signal = 1000*smooth_signal;
        if  i == 1
            data_matrix = smooth_signal;
        else
            data_matrix = horzcat(data_matrix, smooth_signal);
        end
        clear aa b c indsweep1 indsweep2 indsweep3 indsweep4 Trace*
    end
    
    cd ..
    
elseif strcmp(load_from, 'dat_file')
    data_matrix = load('dataMatrix.dat'); data_matrix = data_matrix';
end

% If you wish to plot the data matrix uncomment

%     figure()
%     plot(t,data_matrix,'b'); hold on;
%     hold on
%     xlabel('s')
%     ylabel('mV')

%%
% Number of row and columns in data matrix
[r, c] = size(data_matrix);

% Finding peak index from each experiment (each column from data_matrix is
% an experiment)

% Store number of peaks in each experiment
Npeaks = [];

% Store peak index
index = {};

% Store thrshold values
thr_values = [];

for j = 1:c
    % Peaks above -10 mV.
    % For each experiment it will store the index of all peaks found
    % findNpeks return the number of peaks, and a cell containing its
    % indexes
    [Npeaks(j), index{j}] = findNpeaks(data_matrix(:,j), -10);
end

%%
% data_matrix2 is a cell each line of the cell will store the potencial
% values below the found threhold for a particular expriment
data_matrix2 = {};
for j = 1:c
    % If it have any peaks
    if ~isempty(index{j})
        % Just convert to a matrix
        aux = cell2mat(index{j});
        % Runs towards all peaks
        for i = 1:length(aux)
            % Separate each peak, to find the threshold
            % P is one peak of a given experimente
            % P is choosed from the peak to 3ms ago
            P = data_matrix(aux(i)-60:aux(i),j);
            % dP/dt
            P1 = diff(P);
            % d²P/dt²
            P2 = diff(P1);
            % Method VI
            Kp = P2.*(1+(P1(1:end-1).^2)).^(-3/2);
            % Find the max Kp and its index
            [max_h, aux2] = max(Kp);
            thr_values(end+1) = P(aux2);
            data_matrix2{end+1} = P(P <= P(aux2));
            % Plot the separated peaks, with the found threshold
            figure(1);
            plot(P,'black');
            hold on;
            plot(aux2, P(aux2),'Ored');
        end
    end
end
ylabel('Membrane Potencial [mV]')
legend('Action Potential','Threshold values')
% axis([25 45 0 -20])
% print('potencial_thr_values','-dpng','-r600')
%%
bl = 2;             % Length of the bins
v_m = [];           % Potencial value of the center of each bar.
pot_disp = [];      % Spike potencials.
bins_pot = [];      % Distribution of the potencial values.
bins_potdisp = [];  % Distribution of the potencial values.

% Dividing into bins
for v = -70:bl:1
    v_m(end+1) = v + bl/2;
    aux = 0;
    for j = 1:length(data_matrix2)
        aux = aux + sum( (data_matrix2{j} >= v & data_matrix2{j} < v+bl) );
    end
    bins_pot(end+1) = aux;
    bins_potdisp(end+1) = sum( (thr_values >= v & thr_values < v+bl) );
end
% Remove zeros from bins_pot
[zeros, index] = find(bins_pot > 0);
bins_pot = bins_pot(index);
bins_potdisp = bins_potdisp(index);
% Calculate phi(v)
v_m = v_m(index);
phi_v = bins_potdisp ./ bins_pot;
figure
plot(v_m,phi_v);
ylabel('Probability')
xlabel('Membrane Potencial [mV]')
% print('graf_prob_disp','-dpng','-r600')
%% Histograms
figure
subplot(2,2,1)
bar(v_m, bins_pot)
title('Vm')
ylabel('Counts')
xlabel('Membrane Potencial [mV]')
subplot(2,2,2)
normplot(bins_pot)

subplot(2,2,3)
bar(v_m, bins_potdisp)
title('Threshold Vm')

ylabel('Counts')
xlabel('Membrane Potencial [mV]')

subplot(2,2,4)
normplot(bins_potdisp)
% print('hist','-dpng','-r600')
%% Fits phi_v

idx = find(v_m <= -40);
idx2 = find(v_m > -40);

bseline = mean(phi_v(idx));
sbseline = std(phi_v(idx));
f = fit((v_m(idx2)+max(abs(v_m(idx2))))', (phi_v(idx2)-bseline+0.0115)','exp1');
phi_v_fit = [];
for v = -55:.01:-34
    if v <= -40
        phi_v_fit(end+1) = 0;%bseline;
    elseif v > -40 && v < -34
        phi_v_fit(end+1) = f.a*exp(f.b *(v + 40)) -f.a;%+ bseline;
    else
        phi_v_fit(end+1) = 1;
    end
end

v = -55:.01:-34;
