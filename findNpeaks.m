% Count the number of peaks from a membrane potencial signal above a given
% th value.
% signal: signal from neuron

function [Npeaks, index] = findNpeaks(signal, th)
    index = {};
    peaks = findpeaks(signal);
    peaks = peaks(peaks>th);
    Npeaks = length( peaks );
    if Npeaks > 0
        for i = 1:length(peaks)
            index{end+1} = find(signal == peaks(i));
        end
        index = index';
    end
end