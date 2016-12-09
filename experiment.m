% Experiment with the experimental phi_v in a GL_neuron

clear all
clc

%

dt = 1; % Integration step [ms]

%% Parameters
tau = 20;  % Membrane time constant [ms]
Vr  = -60; % Rest potential [mV]
R   = 100; % Membrane resistence [Mohm]
Vl  = -50; % Threshold potencial [mV]

%%
Tsim = 1000; % Maximum simulation time [ms].
t = 0:dt:Tsim;
I = 0:1:300; % Currents
% Cell parameters
V = Vr.*ones(1, length(t));

% Probability function parameters

p = 0;
a = 0.0066;
b = 0.835;

single_neuron = 0;

if single_neuron == 1
    
    for i = 1:length(t)-1
        V(i+1) = V(i) + (dt/tau) * (-V(i) + Vr + R*0.3);
        if V(i+1) <= -40
            p = 0;
        elseif V(i+1) > -40 && V(i+1) < -34
            p = a * (exp(b*(V(i+1)+40)) - 1);
        else
            p = 1;
        end
        if rand() < p
            V(i+1) = Vr;
        end
    end
    
elseif single_neuron == 0
    
    counts = zeros(1,length(I));
    for j = 1:length(I)
        for k = 1:20
            V = Vr.*ones(1, length(t));
            for i = 1:length(t)-1
                %             V(i+1) = V(i) + (dt/tau) * (-V(i) + Vr + R*I(j));
                V(i+1) = 0.9*V(i) + 0.1*Vr + 10e-3*I(j);
                if V(i+1) <= -40
                    p = 0;
                elseif V(i+1) > -40 && V(i+1) < -34
                    p = a * (exp(b*(V(i+1)+40)) - 1);
                else
                    p = 1;
                end
                if rand() < p
                    V(i+1) = Vr;
                    counts(j) = counts(j) + 1;
                end
            end
        end
        counts(j) = counts(j)/20;
    end
    
end
