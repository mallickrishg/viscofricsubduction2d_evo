function [rcv,shz,src,results] = RUN_simpleSZ_imposedcycles(faultparams,shzparams,Trecur,ncycles)
% Simulates postseismic and interseismic fault creep in response to imposed
% earthquakes at a recgular interval for a subduction zone in 2-d
% Rishav Mallick, EOS, 2021

% Rigidity (MPa)
G = 30e3;
nu = 0.25;

% starting point
y2i = 0;
y3i = 0;

% Fault parameters
Fwidth = 250e3;
Mf = 100;

% Viscous domain
Vwidth = 200e3;
Mv = 80;
Vwidthbot = 500e3;
Mvbot = 100;

% plate thickness
Tplate = faultparams.Tplate;
%% Create faults and shear zones
earthModel = fgeom.LDhs(G,nu);
[rcv,shz,src] = create_flt_shz(earthModel,[y2i,y3i],faultparams.dip,Fwidth,Mf,Vwidth,Mv,Vwidthbot,Mvbot,Tplate);

% plate velocity
Vpl = faultparams.Vpl; %(m/s)

%% impose earthquake parameters
Teq = Trecur.*3.15e7;% earthquake every Teq years
% ncycles = 5000*3.15e7/Teq; % number of earthquake cycles

%% Stress Kernels and EVL object

evl = compute_stresskernels(rcv,shz);
disp('Loaded Stress Kernels')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% effective confining pressure on fault (MPa)
rcv.sigma = 50*ones(rcv.N,1);

% frictional parameters
rcv.a = 1e-2*ones(rcv.N,1);
rcv.b = rcv.a - 5e-3*ones(rcv.N,1);
% velocity-weakening
vw = abs(rcv.xc(:,2))>=0e3 & abs(rcv.xc(:,2))<=35e3;
rcv.b(vw) = rcv.a(vw) + 5e-3;

% static friction coefficient
rcv.mu0 = 0.6*ones(rcv.N,1);

% plate velocity (m/s)
rcv.Vpl = Vpl.*rcv.Vpl;

% reference slip rate (m/s)
rcv.Vo = 1e-6*ones(rcv.N,1);

% shear wave speed (m/s)
rcv.Vs = 3e3*ones(rcv.N,1);

% Degrees of Freedom [slip, log(v/vo)]
rcv.dgf = 2;

disp('Assigned Frictional Properties')

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   R H E O L O G Y                    %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% long-term strain rate
shz.Vpl = Vpl.*shz.Vpl;

% Viscosity in Pa-s
etaval_arc = shzparams.eta_arc;
etaval_oc = shzparams.eta_oc;
power = 1;

shz.n = power;% store power in 'tmax' parameter
shz.Ainverse = ones(shz.N,1);
shz.Ainverse(shz.Vpl>0) = etaval_arc*1e-6;% store viscosity in 'a' parameter
shz.Ainverse(shz.Vpl<0) = etaval_oc*1e-6;

% degrees of freedom
shz.dgf = 2;

% Loading Stresses - backslip loading
evl.tau0     = - evl.KK*rcv.Vpl - evl.LK*shz.Vpl;
evl.sigma0   = - evl.KL*rcv.Vpl - evl.LL*shz.Vpl;

%% solve ODE
disp('Solving ODE')

% Initialize State Vector
Y0=zeros(rcv.N*rcv.dgf+shz.N*shz.dgf,1);

% Fault patches
taumax = 3;
eqslip = compute_earthquakeslip(rcv,evl,Teq*Vpl,taumax);
taueq = evl.KK*eqslip;
taueq(taueq<0) = 0;

% initialize Y0
Y0(1:rcv.dgf:rcv.N*rcv.dgf) = zeros(rcv.N,1);% slip
Y0(2:rcv.dgf:rcv.N*rcv.dgf) = log(rcv.Vpl*0.99./rcv.Vo); % v

% Shear zones
Y0(rcv.N*rcv.dgf+1:shz.dgf:end) = 1e-9 + zeros(shz.N,1); %stress12
Y0(rcv.N*rcv.dgf+2:shz.dgf:end) = zeros(shz.N,1); %strain12

%% Simulation 

tic
% initialize the function handle with
% set constitutive parameters

yp=@(t,y) odeimposedViscoelastic_a_zerowidth(t,y,rcv,shz,evl);

tic
% Solve the system
options=odeset('Refine',1,'AbsTol',1e-6,'RelTol',1e-6,'InitialStep',1e-6,'MaxStep',3e8); 
for i = 1:ncycles
    if i==1
        [t,Y]=ode45(yp,[0 Teq],Y0,options);
        %t=t';Y=Y';
    else
        % provide new initial conditions (loaded by earthquake)
        Y0 = Y(end,:);
        
        Y0(1:rcv.dgf:rcv.N*rcv.dgf) = eqslip;        
        Y0(2:rcv.dgf:rcv.N*rcv.dgf) = Y0(2:rcv.dgf:rcv.N*rcv.dgf)' + taueq./(rcv.a-rcv.b)./rcv.sigma;
        
        Y0(rcv.N*rcv.dgf+1:shz.dgf:end) = Y0(rcv.N*rcv.dgf+1:shz.dgf:end)' + evl.KL*eqslip; %stress12
        Y0(rcv.N*rcv.dgf+2:shz.dgf:end) = 0;
        
        [tmod,Ymod]=ode45(yp,[0,Teq],Y0,options);
        t = tmod;
        Y = Ymod;        
    end
    disp(['Cycle ' num2str(i) ' in progress'])
end
toc

%% results
V=repmat(rcv.Vo',size(Y,1),1).*exp(Y(:,2:rcv.dgf:rcv.N*rcv.dgf));
V(:,vw) = 0;
% slip = Y(:,1:rcv.dgf:rcv.N*rcv.dgf);

s12 = Y(:,rcv.N*rcv.dgf+1:shz.dgf:rcv.N*rcv.dgf+shz.N*shz.dgf);
% e12 = Y(:,rcv.N*rcv.dgf+2:shz.dgf:rcv.N*rcv.dgf+shz.N*shz.dgf);
% e12d = s12./repmat(shz.a',size(Y,1),1);
e12d = (s12.^repmat(shz.n',size(Y,1),1))./repmat(shz.Ainverse',size(Y,1),1);

results = [];
results.t = t;
results.V = V;
results.e12d = e12d;

toc
disp('ODE Solution')

end