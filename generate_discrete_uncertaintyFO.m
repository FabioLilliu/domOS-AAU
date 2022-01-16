%collect data
T = 8; %time horizon
data = struct();
data.InitialSoC = 7; %initial state of charge
data.L = sqrt(0.9); %loss of the battery (square root of roundtrip efficiency)
data.Qmin = 0; %minimum possible state of charge
data.Qmax = 14;%maximum possible state of charge
data.Pmin = -5; %maximum discharging power
data.Pmax = 5; %maximum charging power
data.decay = 1; %decay at each time unit
data.percentage = 100; %percentage of flexibility exploited by DFO
pt = 0.95; %probability threshold
gr = 100; %granularity

%generate slices
slicestruct = prob_slices(data,T,gr);
slices = slicestruct.slices;
slicespace = slicestruct.space;

%plot graph
for t = 1:T
subplot(ceil(T/4), 4, t);
plot(slices(t,:),slicespace);
xlim([0,1.2]);
title(sprintf('T=%d', t));    
    ylabel(sprintf('Energy(kWh)'));
    xlabel(sprintf('Probability'));    
end