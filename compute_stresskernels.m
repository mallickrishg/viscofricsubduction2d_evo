function evl =  compute_stresskernels(rcv,shz)

% fault-fault
[Ktau,~] = computetractionkernels(rcv,rcv);

evl = [];
evl.KK = Ktau;

%fault - shearzone
[Ktau,~] = computetractionkernels(rcv,shz);
evl.KL = Ktau;

%shearzone - shearzone
[Ktau,~] = computetractionkernels(shz,shz);
evl.LL = Ktau;

%shearzone - fault
[Ktau,~] = computetractionkernels(shz,rcv);
evl.LK = Ktau;


end