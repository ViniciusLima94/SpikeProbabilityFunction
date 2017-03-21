% Count the number of peaks from a membrane potencial signal above a given
% th value.
% signal: signal from neuron

function [Npeaks, index] = findNpeaks(signal, th)
    [peaks, index] = findpeaks(signal);
    index = index( peaks > th );
    Npeaks = length( index );
end
