T = 6; %time horizon
MinPrice = 1; %minimum possible price
MaxPrice = 400; %maximum possible price
space = [MinPrice:MaxPrice]; %total space of prices
s = size(space);
prob = repelem(0,T,s(2)); %initialize matrix of probability distributions

%insert predicted mean values (the vector has to have size T)
spotprices = [1.622800000000000e+02,1.967400000000000e+02,1.942900000000000e+02,1.838600000000000e+02,1.841600000000000e+02,1.340600000000000e+02];

for t = 1:T
    for k = space
        prob(t,k) = normpdf(k,spotprices(t),spotprices(t)/5);
    end
end
for t = 1:T
subplot(ceil(T/3), 3, t);
plot(prob(t,:),space);
xlim([0,0.015]);
title(sprintf('T=%d', t));    
    ylabel(sprintf('Price(DKK)'));
    xlabel(sprintf('Probability'));    
end
