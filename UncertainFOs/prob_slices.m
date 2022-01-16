function slicestructure=prob_slices(d,t,gr)
%generate slices for uncertain FOs

Cmax = d.L*d.Pmax;

Q = [];
totalspace = [d.Qmin:1/gr:d.Qmax]; %possible values for SoC
length = size(totalspace);

%Generating probability distribution for SoC

%first slice is Dirac delta
slicezero = repelem(0,length(2));
slicezero(floor(d.InitialSoC*gr)+1) = length(2)/(d.Qmax-d.Qmin);
Q = [Q,slicezero];

s1 = slicezero;
CountMax = floor(Cmax*gr); 
CountMin = floor(d.Pmin*gr);

for j = 1:t-1
    s2 = repelem(0,length(2));
    total = sum(s1);
    for k = 1:length(2)
        ccmax = max(d.Qmin*gr+1,k-CountMax);
        ccmin = min(d.Qmax*gr+1,k-CountMin);
        s2(k) = sum(s1(ccmax:ccmin))/total;
    end
    ss = sum(s2);
    s2 = s2*gr/(ss);
    Q = [Q;s2];
    s1 = s2;
end

si = size(Q);
for k = 1:si(1)
ev = 0;
for j = 1:si(2)
ev = ev + j*Q(k,j)/gr;
end
ev = ev/gr;
end

%Generating slices

slices = [];
slicespace = [d.Pmin:1/gr:Cmax];
slicelength = size(slicespace);
slicezero = find(slicespace == min(abs(slicespace)));
slicespace1 = [d.L*d.Pmin:d.L*abs(d.Pmin)/(slicezero(1)-1):0,d.Pmax/(slicelength(2)-slicezero(1)):d.Pmax/(slicelength(2)-slicezero(1)):d.Pmax];
tt = -d.Pmin*gr + 1;
for j = 1:t
    slice = repelem(0,slicelength(2));
    QS = Q(j,:);
    for k = 1:slicelength(2)
        slice(k) = sum(QS(max(tt-k,1):min(length(2),length(2)+tt-k)))/sum(QS);
    end
    slices=[slices;slice];
end
slicestructure = struct;
slicestructure.slices = slices;
slicestructure.space = slicespace1;
