function Yp = odeimposedViscoelastic_a_zerowidth(~,Y,ss,shz,evl)
% ODE function for faults and shear zones evolution
% Y = [slip,V - flt
%     [s12,e12] - shz
% We express dY/dt as f(Y) and using RK-4th order to integrate with
% adaptive timesteps
% Rishav Mallick, 2021, EOS

%% FAULTS
damping = ss.earthModel.G./ss.Vs/2;
VS = ss.a>ss.b;

% Slip velocities 
dummy = ss.Vo.*exp(Y(2:ss.dgf:ss.N*ss.dgf));
V = zeros(ss.N,1);
V(VS) = dummy(VS);

% Initiate state derivative
Yp=zeros(size(Y));  

% Slip velocity
Yp(1:ss.dgf:ss.N*ss.dgf) = V;

%% SHEAR ZONES
% for shz -> shz.a is viscosity (MPa-s), and shz.tMax is the power (for now only have 1)

s12 = Y(ss.N*ss.dgf+1:shz.dgf:ss.N*ss.dgf+shz.N*shz.dgf);

if shz.tMax==1  % power  
    e12dot = (s12)./shz.a;
else
    error('Not a valid rheology - change power/combo')
end

% strain rates
Yp(ss.N*ss.dgf+2:shz.dgf:ss.N*ss.dgf+shz.N*shz.dgf) = e12dot;

%% INTERACTIONS
% Acceleration (rate of log(V/Vo))
kv = evl.KK*V + evl.tau0 + evl.LK*e12dot;
kv(~VS) = 0;

Yp(2:ss.dgf:ss.N*ss.dgf) = (kv)./((ss.a-ss.b).*ss.sigma + damping.*V);

% Stressing Rates 
s12dot = evl.sigma0 + evl.LL*e12dot + evl.KL*V;

Yp(ss.N*ss.dgf+1:shz.dgf:ss.N*ss.dgf+shz.N*shz.dgf) = s12dot;

end


