% Testing a viscoelastic earthquake cycle model in planestrain (2d)
% With normal stress perturbation and a maxwell mantle
% Rishav Mallick, 2019, Earth Observatory of Singapore

clear

addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*
minmax = @(x) [min(x) max(x)];

s2y=60*60*24*365;
y2s=1./s2y;
G=30e3;
nu=1/4;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%           R E C E I V E R   F A U L T S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

rcv=geometry.receiver('./faults/megathrust_2d.seg',greens.okada92(G,nu));
% define dip-distance in (m)
dipdist = zeros(rcv.N,1);
for i = 1:rcv.N
    if i==1
        dipdist(i) = rcv.W(1)/2;
    else
        dipdist(i) = dipdist(i-1) + rcv.W(i);
    end
end
fault_width = max(dipdist);
% create velocity weakening region
VW = zeros(rcv.N,1);
topvw = floor(20e3/(fault_width/rcv.N));
botvw = ceil(80e3/(fault_width/rcv.N));
VW(topvw:botvw) = 1;
VW = logical(VW);

% effective confining pressure on fault (MPa) 
rcv.sigma = 100.*ones(rcv.N,1);

% frictional parameters 
rcv.a = 1e-2.*ones(rcv.N,1);
% if you want to set the entire thing as VS
rcv.b = rcv.a - 5e-3;
rcv.b(VW) = rcv.a(VW) + 4e-3;
% static friction coefficient
rcv.mu0 = 0.6*ones(rcv.N,1);
% numer of state parameters
rcv.dgf=5;
% characteristic weakening distance (m)
rcv.l = 0.03.*ones(rcv.N,1);
% plate velocity (m/s)
rcv.Vpl = 1e-9*ones(rcv.N,1);

% reference slip rate (m/s)
rcv.Vo = 1e-6*ones(rcv.N,1);
% shear wave speed (m/s)
rcv.Vs = 3e3.*ones(rcv.N,1);

% minimum grid size
fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),G/(1-nu)*rcv.l(topvw+1)/rcv.b(topvw+1)/rcv.sigma(topvw+1)) 
fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(1)/rcv.b(1))
fprintf(1,'a/b in VW regions = %.2f\n',rcv.a(topvw+1)/rcv.b(topvw+1))

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%               S H E A R   Z O N E S                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% load shear zone
shz=unicycle.geometry.shearZoneReceiver('./faults/lowercrust',unicycle.greens.shearZone16(G,nu));

Pas2MPay = 1e-6/3.15e7;
% setup different viscosity on either side of the dauki
shz.etaM= zeros(shz.N,1) + 1e18/1e6; % MPa s
shz.etaM(shz.xc(:,2)<-1.9*shz.xc(:,3) - 120e3) = (1e24/1e6);

%% GEOMETRY REVIEW

if false
    figure(1),clf, set(gcf,'Name','Geometry Review')
    subplot(2,1,1)
    shz.plotShearZoneEdges(log10(shz.etaM*1e6)), hold on
    shz.plotShearZoneById([1:shz.N])
    rcv.plotPatch
    view(0,90)
    axis tight equal
    xlabel('East (m)'),ylabel('North (m)')
    
    subplot(2,1,2)
    shz.plotShearZoneEdges(log10(shz.etaM*1e6)), hold on
    shz.plotShearZoneById([1:shz.N])
    rcv.plotPatch
    view(90,0)
    axis tight equal, grid on
    xlabel('East (m)'),ylabel('North (m)'),zlabel('Z (m)')
    colorbar
    colormap(parula(20))
    caxis([10 30])
end
% return
%% load kernels
tic
if ~exist('evl','var')
    evl=unicycle.ode.rateStrengtheningMaxwell([],rcv,shz,[],'./kernels/');
    evl.flt.dgf = 5;
    evl.shz.dgf = 6;
    dgf = evl.flt.dgf + evl.shz.dgf;
end
disp('Loaded Stress Kernels')
toc
remoteloading = [];
% remote tectonic loading on fault (weird units - MPa/s)
remoteloading.dtau_inf = -evl.KK{2,2}*rcv.Vpl;
remoteloading.dsigmabar_inf = -evl.KK{2,3}*rcv.Vpl;

remoteloading.dsigma22_inf = -evl.KL{2,4}*rcv.Vpl;
remoteloading.dsigma23_inf = -evl.KL{2,5}*rcv.Vpl;
remoteloading.dsigma33_inf = -evl.KL{2,6}*rcv.Vpl;
%% total duration of simulation
tend = 1000*3.15e7;

%% Initialize State Vector
Y0=zeros(evl.flt.N*evl.flt.dgf + evl.shz.N*evl.shz.dgf,1);
% Fault patches
% Y0(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N) = zeros(rcv.N,1);         % slip
Y0(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N) = evl.flt.sigma.*(evl.flt.mu0+(evl.flt.a-evl.flt.b).*log(evl.flt.Vpl./evl.flt.Vo))+G*evl.flt.Vpl./(2*evl.flt.Vs); %stress 
Y0(3:evl.flt.dgf:evl.flt.dgf*evl.flt.N) = log(evl.flt.Vo./evl.flt.Vpl);      % log(theta Vo / L)
Y0(4:evl.flt.dgf:evl.flt.dgf*evl.flt.N) = log(evl.flt.Vpl*0.999./evl.flt.Vo); % log(V/Vo)
Y0(5:evl.flt.dgf:evl.flt.dgf*evl.flt.N) = evl.flt.sigma; % normal stress
% Shear zones
% Y0(evl.flt.dgf*evl.flt.N + 1:evl.shz.dgf:evl.flt.dgf*evl.flt.N+evl.shz.dgf*evl.shz.N) = zeros

% initialize the function handle with
% set constitutive parameters
yp=@(t,y) ode_fault_shearzones2d(t,y,evl,remoteloading);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',1e-7,'InitialStep',1e-5,'MaxStep',1e7); 
[t,Y]=ode45(yp,[0 tend],Y0,options);
disp('Finished running ODE solver')
toc

%% extract important quantities

V = evl.flt.Vo(1).*exp(Y(1:end,4:evl.flt.dgf:evl.flt.N*evl.flt.dgf));
slip = Y(:,1:evl.flt.dgf:evl.flt.N*evl.flt.dgf);
tau = Y(:,2:evl.flt.dgf:evl.flt.N*evl.flt.dgf);
sigmabar = Y(:,5:evl.flt.dgf:evl.flt.N*evl.flt.dgf);
e22 = Y(:,evl.flt.N*evl.flt.dgf+1:evl.shz.dgf:end);
e23 = Y(:,evl.flt.N*evl.flt.dgf+2:evl.shz.dgf:end);
e33 = Y(:,evl.flt.N*evl.flt.dgf+3:evl.shz.dgf:end);
emax = sqrt((0.5*(e22-e33)).^2 + (e23/2).^2);

Yp = zeros(size(Y));
parfor i = 1:length(t)
    Yp(i,:) = yp(t(i),Y(i,:)');
end

e22dot = Yp(:,evl.flt.N*evl.flt.dgf+1:evl.shz.dgf:end);
e23dot = Yp(:,evl.flt.N*evl.flt.dgf+2:evl.shz.dgf:end);
e33dot = Yp(:,evl.flt.N*evl.flt.dgf+3:evl.shz.dgf:end);
edotmax = sqrt((0.5*(e22dot-e33dot)).^2 + (e23dot/2).^2);







