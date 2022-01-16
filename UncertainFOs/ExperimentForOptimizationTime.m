format shortG
yalmip('solver','glpk');
sett = sdpsettings('verbose',0);

%Import M and N, the matrices with spot and imbalance prices respectively
load('M.mat');
load('N.mat');

timeposs = 18:24; %which times we want to check
attempts = 5; %how many optimizations we want to try
MRes = repelem(0,3);
MMRes = [];

for time = timeposs

    T = time;
    InitialSoC = 0; %initial state of charge
    L = sqrt(0.9); %loss of the battery (square root of roundtrip efficiency)
    Qmin = 0; %minimum possible state of charge
    Qmax = 14;%maximum possible state of charge
    Pmin = -5; %maximum discharging power
    Pmax = 5; %maximum charging power
    decay = 1; %decay at each time unit
    percentage = 100; %percentage of flexibility exploited by DFO
    pt = 0.95; %probability threshold
    gr = 100; %granularity

    days = 1;%number of days for iteration
    Length = 1; %number of cycles needed

    RMatrix = repelem(0,19,Length); %initializing result matrix
    SoC0 = 0; %initializing "continuity array", 
              %it describes the state of charge at the beginning of next cycle
    SoC0(:,1) = InitialSoC'; %initial state of charge
    daysahead = 0; %how many days ahead we want to start; mainly for debugging
    Unit = eye(T);

    spotprices = M(1:T,1:1)'; %spot prices
  
    A = decay;
    B3 = [L,1/L]; %FOR DFOS ONLY
    f = [0];
    C = [0];
    D1 = [1,1]; %FOR DFOS ONLY
    PminV = [0;L*Pmin]; %FOR DFOS ONLY
    PmaxV = [Pmax;0]; %FOR DFOS ONLY
    lti = LTISystem('A', A, 'B', B3, 'C', C, 'D', D1); %FOR DFOS ONLY
    lti.u.min = PminV;
    lti.u.max = PmaxV;
    lti.x.min = Qmin;
    lti.x.max = Qmax;
    lti.initialize(InitialSoC);
    data1 = struct();
    data1.InitialSoC = InitialSoC;
    data1.L = L;
    data1.Pmin = Pmin;
    data1.Pmax = Pmax;
    data1.Qmin = Qmin;
    data1.Qmax = Qmax;
    data1.percentage = 100;

    P = Polyhedron('lb', [InitialSoC], 'ub', [InitialSoC]);
    lti.x.with('initialSet');
    lti.x.initialSet = P;

    lti.instantiate(T);
    F1 = FlexSystem(lti);
    
    %generate exact constraints
    
    Q = sdpvar(1,T+1);
    psum = F1.getEnergyVars();

    Q(1) = InitialSoC;

    for k = 1:T
       Q(k+1) = decay*Q(k) + L * max(psum(k),0) + 1/L * min(psum(k),0);
    end

    constraints = [repelem(L*Pmin,T) <= psum <= repelem(Pmax,T), repelem(Qmin,T+1) <= Q <= repelem(Qmax,T+1)];
    ccost = spotprices(1:T) * psum';

    %generate SFO constraints
    
    evarsSFO = sdpvar(1,T);
    SFOcost = spotprices(1:T) * transpose(evarsSFO);
    c_AAO = const_AAO(data1,T);
    cAAO = [c_AAO(1,:) <= evarsSFO <= c_AAO(2,:)];
    
    %generate uncertain FO constraints
    
    evarsPFO = sdpvar(1,T);
    PFOcost = spotprices(1:T) * transpose(evarsPFO);
    %measure time for generating uncertain FO constraints
    tic
    probslice1=prob_slices(data1,T,gr);
    cprob = const_prob(probslice1.slices,probslice1.space,pt);
    aa = toc
    PC = [cprob(1,:) <= evarsPFO <= cprob(2,:)];
    
    %optimization
    
    vres = repelem(0,3,attempts);
    for k = 1:attempts
        tic
        %Solution = optimize(constraints, ccost,sett);
        vres(1,k) = toc;
        tic
        SFOSolution = optimize(cAAO,SFOcost,sett);
        vres(2,k) = toc;
        tic
        PFOSolution = optimize(PC,PFOcost,sett);
        vres(3,k) = toc;        
    end
    for k = 1:3
        MRes(k) = sum(vres(k,:)/attempts);
    end
    MMRes = [MMRes; MRes];
end
