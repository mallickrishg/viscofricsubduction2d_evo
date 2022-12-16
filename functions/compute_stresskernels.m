function evl =  compute_stresskernels(rcv,shz)

% fault-fault
[Ktau,~] = fgeom.computetractionkernels(rcv,rcv);

evl = [];
evl.KK = Ktau;

%fault - shearzone
[Ktau,~] = fgeom.computetractionkernels(rcv,shz);
evl.KL = Ktau;

%shearzone - shearzone
[Ktau,~] = fgeom.computetractionkernels(shz,shz);
evl.LL = Ktau;

%shearzone - fault
[Ktau,~] = fgeom.computetractionkernels(shz,rcv);
evl.LK = Ktau;


end