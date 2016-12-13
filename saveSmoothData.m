% Routine to create smooth_signal of N set of data
clear all
clc

data_matrix = load('dataMatrix.dat');
data_matrix = data_matrix';

% If you wish to plot the data matrix uncomment

% figure()
% plot(t,data_matrix,'b'); hold on;
% hold on
% xlabel('s')
% ylabel('mV')

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
    [Npeaks(j), index{j}] = findNpeaks(data_matrix(:,j), -10);
end

for j = 1:c
    if ~isempty(index{j})
        aux = cell2mat(index{j});
        for i = 1:length(aux)
            % Separate each peak, to find the threshold
            P = data_matrix(aux(i)-60:aux(i),j);
            % dP/dt
            P1 = diff(P);
            % d²P/dt²
            P2 = diff(diff(P));
            Kp = P2.*(1+(P1(1:end-1).^2)).^(-3/2);
            [max_h, aux2] = max(Kp);
            thr_values(end+1) = P(aux2);
        end
    end
end

%% Saving all potencial values in a vector
v = [];
for i = 1:r
    for j = 1:c
        v(end+1) = data_matrix(i,j);
    end
end

%%

nbins = 50;
% Create potentials histrogram
[counts, center] = hist(v,nbins);
% Bin length
bin = diff(center);
bin = bin(1);
% Creating threshold potential values
counts_1 = [];
for i = 1:length(center)
    counts_1(end+1) = sum(thr_values >= center(i)-bin/2 & thr_values < center(i)+bin/2);
end
% Calculating probability density
gamma = sum(counts_1./counts)^-1;
ro = gamma .* counts_1./counts;

% Calculating probability
phi_v = cumsum(ro);
plot(center, phi_v);
title('Spike probability')
ylabel('Probability')
xlabel('Membrane Potential [mV]')
print('graf_prob_disp','-dpng','-r600')
