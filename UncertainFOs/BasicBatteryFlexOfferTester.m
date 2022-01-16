SoC0 = 0; %initial state of charge
Qmin = 0; %minimum possible state of charge
Qmax = 14;%maximum possible state of charge
Pmin = -5; %maximum discharging power
Pmax = 5; %maximum charging power
decay = 1; %decay at each time unit
T = 5; %time horizon
L = sqrt(0.9);%loss of the battery (square root of roundtrip efficiency)

% Model is taken from:
%       https://control.me.berkeley.edu/~sanandaji/my_papers/Allerton2013_TCL.pdf
% LTI Representation
%    x(t+Ts) = A*x(t) + B*u(t)+f 
%       y(t) = C*x(t) + D*u(t) 
%
A = [decay];
B = [L,1/L];
f = [0];
C = [0];
D = [1,1];
lti = LTISystem('A', A, 'B', B, 'C', C, 'D', D, 'f', f);

lti.u.min = [0,L*Pmin];
lti.u.max = [Pmax,0];
lti.x.min = Qmin;
lti.x.max = Qmax;
lti.initialize(SoC0);
lti.instantiate(T); 
F1 = FlexSystem(lti);
 
ODFO = DFOSystem(F1,-1);
ODFO.plot_slices()
