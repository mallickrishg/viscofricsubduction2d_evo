function slip_fullfault = compute_earthquakeslip(ss,evl,maxslip,taumax)

vw = ss.a<ss.b;
vs = ~vw;

slip_fullfault = zeros(ss.N,1);
slip_fullfault(vw) = maxslip;
slipco = slip_fullfault(vw);

Nco = length(slipco);
Naf = ss.N-Nco;

% set bounds on slip amplitude
lb = zeros(Naf,1);
ub = maxslip.*ones(Naf,1);

% stress kernels (sub-sets)
Kaf = evl.KK(vs,vs);
Kco = evl.KK(vs,vw);

tauco = Kco*slipco;

Ain = [Kaf;-Kaf];
bin = [(taumax.*ones(Naf,1) - tauco);tauco];

slip_af = lsqlin(eye(Naf),zeros(Naf,1),Ain,bin,[],[],lb,ub);
slip_fullfault(vs) = slip_af;

end