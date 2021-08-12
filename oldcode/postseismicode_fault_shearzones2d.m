function Yp=postseismicode_fault_shearzones2d(~,Y,o,remoteloading)
% ODE function: call as follows
%   [t,y]=ode45(yp,[0 200],y0,options);
%
% where o is a class that contains all the relevent kernels and geometry
% information (evl class)
%        |slip, velocity .   |
%        |       sigma       |
%        |      epsilon_ij   | (3 strains)
%        \      sigma_ij     / (3 stresses)
%
% we use steady-sate Rate-state friction for the faults
% mu = mu0 + a-b*log(v/v0) 
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

% Slip velocities 
V    = o.flt.Vo.*exp(Y(2:o.flt.dgf:o.flt.dgf*o.flt.N));
% pin the velocities of the VW segments
V(VW) = 0;

% Slip velocity
Yp(1:o.flt.dgf:o.flt.dgf*o.flt.N) = V;

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
m = o.shz.m;
e22d=(s22p.*(s22p.^(m-1)))./o.shz.etaM; 
e23d=(s23p.*(s23p.^(m-1)))./o.shz.etaM; 
e33d=(s33p.*(s33p.^(m-1)))./o.shz.etaM; 

% e22d=(s22p)./o.shz.etaM; 
% e23d=(s23p)./o.shz.etaM; 
% e33d=(s33p)./o.shz.etaM; 

zerostrain = 1e-16;
e22d(abs(e22d)<zerostrain | o.shz.sigma0==1) = 0;
e23d(abs(e23d)<zerostrain | o.shz.sigma0==1) = 0;
e33d(abs(e33d)<zerostrain | o.shz.sigma0==1) = 0;

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

% shear stress rate of change (Stress_kernel*strainrate_inelastic + stresskernal*faultvelocity + remoteloading)
% taudot22 = o.LL{4,4}*e22d+o.LL{5,4}*e23d+o.LL{6,4}*e33d + o.KL{2,4}*V + remoteloading.dsigma22_inf;
% taudot23 = o.LL{4,5}*e22d+o.LL{5,5}*e23d+o.LL{6,5}*e33d + o.KL{2,5}*V + remoteloading.dsigma23_inf;
% taudot33 = o.LL{4,6}*e22d+o.LL{5,6}*e23d+o.LL{6,6}*e33d + o.KL{2,6}*V + remoteloading.dsigma33_inf;

taudot22 = o.LL{1,1}*e22d+o.LL{3,1}*e23d+o.LL{6,1}*e33d + 0*o.KL{2,1}*V + remoteloading.dsigma22_inf;
taudot23 = o.LL{1,3}*e22d+o.LL{3,3}*e23d+o.LL{6,3}*e33d + 0*o.KL{2,3}*V + remoteloading.dsigma23_inf;
taudot33 = o.LL{1,6}*e22d+o.LL{3,6}*e23d+o.LL{6,6}*e33d + 0*o.KL{2,6}*V + remoteloading.dsigma33_inf;

Yp(o.flt.dgf*o.flt.N+4:o.shz.dgf:end)=taudot22;
Yp(o.flt.dgf*o.flt.N+5:o.shz.dgf:end)=taudot23;
Yp(o.flt.dgf*o.flt.N+6:o.shz.dgf:end)=taudot33;

% Acceleration (rate of log(V/Vo))
dzetadt = (remoteloading.dtau_inf + o.KK{2,2}*V + 0*[o.LK{1,2} o.LK{3,2} o.LK{6,2}]*[e22d;e23d;e33d])./((o.flt.a - o.flt.b).*o.flt.sigma);
dzetadt(VW) = 0;
Yp(2:o.flt.dgf:o.flt.dgf*o.flt.N) = dzetadt;


end 