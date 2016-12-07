% Routine to create smooth_signal of N set of data
clear all
clc

load_from = 'dat_file'; % Or mat_file

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
    
    % If you wish to plot the data matrix uncomment
    
%     figure()
%     plot(t,data_matrix,'b'); hold on;
%     hold on
%     xlabel('s')
%     ylabel('mV')
    
elseif strcmp(load_from, 'dat_file')
    data_matrix = load('dataMatrix.dat'); data_matrix = data_matrix';
end

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
    [Npeaks(j), index{j}] = findNpeaks(data_matrix(:,j), -10);
end

%%
close all

data_matrix2 = {};
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
            if P(aux2) < -30
                thr_values(end+1) = P(aux2);
                data_matrix2{end+1} = P(P <= P(aux2));
                figure(1);
                plot(P)
                hold on
                plot(aux2, P(aux2),'O')
            end            
        end
    end
end
%%
bl = 1;           % Length of the bins
v_m = [];           % Potencial value of the center of each bar.
pot_disp = [];      % Spike potencials.
bins_pot = [];      % Distribution of the potencial values.
bins_potdisp = [];      % Distribution of the potencial values.

% Dividing into bins
for v = -70:bl:1
    v_m(end+1) = v + bl/2;
    aux = 0;
    for j = 1:length(data_matrix2)
        %         aux = aux + length( find(data_matrix2{j} >= v & data_matrix2{j} < v+bl) );
        aux = aux + sum( (data_matrix2{j} >= v & data_matrix2{j} < v+bl) );
    end
    bins_pot(end+1) = aux;
    bins_potdisp(end+1) = sum( (thr_values >= v & thr_values < v+bl) );
end

phi_v = bins_potdisp ./ bins_pot;

figure
plot(v_m,phi_v);
ylabel('Probability')
xlabel('Membrane Potencial [mV]')
% print('graf_prob_disp','-dpng','-r600')

