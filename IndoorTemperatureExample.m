
TempMin = 20; %Minimum indoor temperature
TempMax = 23; %Maximum indoor temperature
T = 5; %Time units for FO
TempStart = 22; %Initial temperature
MinOT = 0; %Minimum possible value for outdoor temperature
MaxOT = 40; %Maximum possible value for outdoor temperature

% Model is taken from:
%       https://control.me.berkeley.edu/~sanandaji/my_papers/Allerton2013_TCL.pdf
% LTI Representation
%    x(t+Ts) = A*x(t) + B*u(t)+f 
%       y(t) = C*x(t) + D*u(t) 
%

A = 0;
B = 0.41;
f = 13.7;
C = 0;
D = 1;
lti = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'f', f);
lti.u.min = MinOT;
lti.u.max = MaxOT;
lti.x.min = TempMin;
lti.x.max = TempMax;
lti.initialize(TempStart);

 % Outer approximation performance
 lti.instantiate(T);
 F1 = FlexSystem(lti);
 
 ODFO = DFOSystem(F1,-1);
 
ODFO.plot_slices()
