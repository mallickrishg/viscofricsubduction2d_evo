function Yp=interseismicode_fault_shearzones2d(~,Y,o,remoteloading)
% ODE function: call as follows
%   [t,y]=ode45(yp,[0 200],y0,options);
%
% where o is a class that contains all the relevent kernels and geometry
% information (evl class)
%
%        /        s          \            
%        |       tau         |         
%    y = | log(theta Vo / L) |           
%        |   log( V / Vo )   |
%        |       sigma       |
%        |      epsilon_ij   | (3 strains)
%        \      sigma_ij     / (3 stresses)
%
% we use Rate-state friction for the faults
% mu = mu0 + a*log(v/v0) + b*log(v0 theta/L)
% ageing law for theta evolution
%
% and for shear zones we use a maxwell rheology
%
%   d epsilon_ij          sigma_ij
%   ------------      =   --------
%      dt                  etaM
%
% etaM is the viscosity of the maxwell body
% 
%% frictional stresses and interactions
% find velocity weakening region
VW = zeros(o.flt.N,1);
VW(o.flt.b>o.flt.a) = 1;
VW = logical(VW);
% Initiate state derivative
Yp=zeros(size(Y));  

damping = o.flt.earthModel.G./o.flt.Vs/2;

% State variables
th   = Y(3:o.flt.dgf:o.flt.dgf*o.flt.N);
th(VW) = 0;
sigma = Y(5:o.flt.dgf:o.flt.dgf*o.flt.N);

% Slip velocities 
V    = o.flt.Vo.*exp(Y(4:o.flt.dgf:o.flt.dgf*o.flt.N));
% pin the velocities of the VW segments
V(VW) = 0;

% Slip velocity
Yp(1:o.flt.dgf:o.flt.dgf*o.flt.N) = V;

% Rate of state (rate of log(theta/theta0))
dth=(o.flt.Vo.*exp(-th)-V)./o.flt.l;
dth(VW) = 0;
Yp(3:o.flt.dgf:o.flt.dgf*o.flt.N) = dth;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%
%                 S H E A R   Z O N E S
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% isotropic stress

p=(Y(o.flt.dgf*o.flt.N+4:o.shz.dgf:end)+Y(o.flt.dgf*o.flt.N+6:o.shz.dgf:end))/2;

% deviatoric stress components
s22p=Y(o.flt.dgf*o.flt.N+4:o.shz.dgf:end)-p; 
s23p=Y(o.flt.dgf*o.flt.N+5:o.shz.dgf:end);   
s33p=Y(o.flt.dgf*o.flt.N+6:o.shz.dgf:end)-p; 

% Maxwell strain rates d epsilon / dt = sigma' / etaM

e22d=s22p./o.shz.etaM; 
e23d=s23p./o.shz.etaM; 
e33d=s33p./o.shz.etaM; 

% strainrate goes in the time-derivative
Yp(o.flt.dgf*o.flt.N+1:o.shz.dgf:end)=e22d;
Yp(o.flt.dgf*o.flt.N+2:o.shz.dgf:end)=e23d;
Yp(o.flt.dgf*o.flt.N+3:o.shz.dgf:end)=e33d;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%
%          S T R E S S   I N T E R A C T I O N S
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% normal stress change for fault
dsigdt = (remoteloading.dsigmabar_inf + o.KK{2,3}*V + [o.LK{4,3} o.LK{5,3} o.LK{6,3}]*[e22d;e23d;e33d]);
Yp(5:o.flt.dgf:o.flt.dgf*o.flt.N) = 0*dsigdt;


% shear stress rate of change (Stress_kernel*strainrate_inelastic + stresskernal*faultvelocity + remoteloading)
taudot22 = o.LL{4,4}*e22d+o.LL{5,4}*e23d+o.LL{6,4}*e33d + o.KL{2,4}*V + remoteloading.dsigma22_inf;
taudot23 = o.LL{4,5}*e22d+o.LL{5,5}*e23d+o.LL{6,5}*e33d + o.KL{2,5}*V + remoteloading.dsigma23_inf;
taudot33 = o.LL{4,6}*e22d+o.LL{5,6}*e23d+o.LL{6,6}*e33d + o.KL{2,6}*V + remoteloading.dsigma33_inf;

Yp(o.flt.dgf*o.flt.N+4:o.shz.dgf:end)=taudot22;
Yp(o.flt.dgf*o.flt.N+5:o.shz.dgf:end)=taudot23;
Yp(o.flt.dgf*o.flt.N+6:o.shz.dgf:end)=taudot33;

% Acceleration (rate of log(V/Vo))
kv = remoteloading.dtau_inf + o.KK{2,2}*V + [o.LK{4,2} o.LK{5,2} o.LK{6,2}]*[e22d;e23d;e33d];
dzetadt = (kv - o.flt.b.*sigma.*dth - (o.flt.mu0 + o.flt.a.*Y(4:o.flt.dgf:o.flt.dgf*o.flt.N) + o.flt.b.*th).*dsigdt)./(o.flt.a.*sigma + damping.*V);
dzetadt(VW) = 0;
Yp(4:o.flt.dgf:o.flt.dgf*o.flt.N) = dzetadt;

% Shear stress rate on fault due to fault
Yp(2:o.flt.dgf:o.flt.dgf*o.flt.N) = kv - damping.*V.*Yp(4:o.flt.dgf:o.flt.dgf*o.flt.N);

end 