% Simulates postseismic and interseismic fault creep in response to imposed
% earthquakes at a recgular interval for a subduction zone in 2-d
% Rishav Mallick, EOS, 2021

clear
% close all
% addpath ~/Dropbox/scripts/utils/
addpath ODESolving/
addpath odefunction/
addpath greenfunc/
addpath ~/Dropbox/scripts/topotoolbox/colormaps/
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

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
Mv = 50;
Vwidthbot = 500e3;
Mvbot = 100;

% plate thickness
Tplate = 50e3;
%% Create faults and shear zones
earthModel = unicycle.greens.okada92(G,nu);
[rcv,shz,src] = create_flt_shz(earthModel,y2i,y3i,Fwidth,Mf,Vwidth,Mv,Vwidthbot,Mvbot,Tplate);

% plate velocity
Vpl = 1e-9; %(m/s)

%% impose earthquake parameters
Teq = 200.*3.15e7;% earthquake every Teq years
ncycles = 5000*3.15e7/Teq; % number of earthquake cycles

%% Stress Kernels and EVL object

evl = compute_stresskernels(rcv,shz);
disp('Loaded Stress Kernels')


figure(2000),clf
plotpatch2d(rcv,rcv.Vpl), hold on
plotpatch2d(shz,shz.Vpl), shading faceted
plotpatch2d(src,src.Vpl)
xlabel('x_2 (km)'), ylabel('x_3 (km)')
axis tight equal, box on, grid on
xlim([-1 1].*600)
set(gca,'FontSize',20,'LineWidth',2)
colorbar
% return
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
vw = abs(rcv.xc(:,3))>=5e3 & abs(rcv.xc(:,3))<=35e3;
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
% Viscosity in Pa-s
etaval = 2e14;
power = 1;

shz.tMax = power;% store power in 'tmax' parameter
shz.a = etaval*1e-6.*ones(shz.N,1);% store viscosity in 'a' parameter

shz.Vpl = Vpl.*shz.Vpl;

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

if true
    figure(10),clf
    subplot(211)
    plot(rcv.xc(:,1)./1e3,eqslip,'LineWidth',2)
    axis tight
    subplot(212)
    plot(rcv.xc(:,1)./1e3,evl.KK*eqslip,'rx'), hold on
    plot(shz.xc(:,1)./1e3,evl.KL*eqslip,'kx')
    axis tight, grid on
    ylim([-0.1 1].*taumax)
end

% initialize Y0
Y0(1:rcv.dgf:rcv.N*rcv.dgf) = zeros(rcv.N,1);% slip
Y0(2:rcv.dgf:rcv.N*rcv.dgf) = log(rcv.Vpl*0.99./rcv.Vo); % v

% Shear zones
Y0(rcv.N*rcv.dgf+1:shz.dgf:end) = 1e-9 + zeros(shz.N,1); %stress12
Y0(rcv.N*rcv.dgf+2:shz.dgf:end) = zeros(shz.N,1); %strain12
% return
%% Simulation 

tic
% initialize the function handle with
% set constitutive parameters

yp=@(t,y) odeimposedViscoelastic_a_zerowidth(t,y,rcv,shz,evl);

tic
% Solve the system
options=odeset('Refine',1,'AbsTol',1e-6,'RelTol',1e-6,'InitialStep',1e-6,'MaxStep',3e8,'oDir','ode_out/'); 
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
slip = Y(:,1:rcv.dgf:rcv.N*rcv.dgf);

s12 = Y(:,rcv.N*rcv.dgf+1:shz.dgf:rcv.N*rcv.dgf+shz.N*shz.dgf);
e12 = Y(:,rcv.N*rcv.dgf+2:shz.dgf:rcv.N*rcv.dgf+shz.N*shz.dgf);
e12d = s12./repmat(shz.a',size(Y,1),1);

toc
disp('ODE Solution')

%% displacement kernels
ng = 200;
ox = linspace(-100e3,500e3,ng)';

GF_d = compute_displacementkernels([ox 0.*ox 0.*ox],rcv,shz,src);

% surface velocity
gps = [];
gps.rcv.vh = (GF_d.rcv.Gh*V')';
gps.rcv.vz = (GF_d.rcv.Gz*V')';
gps.shz.vh = (GF_d.shz.Gh*e12d')';
gps.shz.vz = (GF_d.shz.Gz*e12d')';
gps.src.vh = Vpl.*(GF_d.src.Gh*src.Vpl)';
gps.src.vz = Vpl.*(GF_d.src.Gz*src.Vpl)';

% plot surface velocities
figure(11),clf
plotindex = [0.05:0.05:0.95].*Teq;
% plotindex = [2,3,5,10,20,30,50,70,80,90,99].*3.15e7;

cspec = cool(length(plotindex));

subplot(211)
plot(ox./1e3,(gps.src.vh + gps.shz.vh(end,:) + gps.rcv.vh(end,:))./Vpl,'r-','LineWidth',3), hold on
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    plot(ox./1e3,(gps.src.vh + gps.shz.vh(tindex,:) + gps.rcv.vh(tindex,:))./Vpl,'-','LineWidth',1,'Color',cspec(i,:))
end
plot(rcv.xc(vw,1)./1e3,0.*rcv.Vpl(vw),'k-','LineWidth',3)
axis tight, grid on
ylim([-1 1])
plot(max(rcv.xc(:,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
plot(max(shz.xc(shz.Vpl>0,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
plot(0.*[1 1],[-1 1],'k-','Linewidth',2)
xlabel('x_2 (km)')
ylabel('v_h/v_{pl}')
set(gca,'Fontsize',20,'LineWidth',2)

subplot(212)
plot(ox./1e3,(gps.src.vz + gps.shz.vz(end,:) + gps.rcv.vz(end,:))./Vpl,'r-','LineWidth',3), hold on
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    plot(ox./1e3,(gps.src.vz + gps.shz.vz(tindex,:) + gps.rcv.vz(tindex,:))./Vpl,'-','LineWidth',1,'Color',cspec(i,:))
end
axis tight, grid on
ylim([-1 1].*0.8)
plot(max(rcv.xc(:,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
plot(max(shz.xc(shz.Vpl>0,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
plot(0.*[1 1],[-1 1],'k-','Linewidth',2)
xlabel('x_2 (km)')
ylabel('v_z/v_{pl}')
set(gca,'Fontsize',20,'Linewidth',2)

%% plot each component's effect on the surface
% horizontal motion

cspec = parula(length(plotindex));

figure(12),clf
set(gcf,'Name','Horizontal velocity contributions')
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    e12dplot = e12d(tindex,:);
    e12dplot(shz.Vpl>0) = 0;
    
    subplot(321)
    toplot = GF_d.shz.Gh*e12dplot';
    plot(ox./1e3,(toplot)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
    subplot(322)
    toplot = GF_d.shz.Gz*e12dplot';
    plot(ox./1e3,(toplot)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
end
subplot(321)
title('Oceanic mantle')
axis tight, grid on
ylim([-1 2])
xlabel('x_2 (km)')
ylabel('v_h/v_{pl}')
set(gca,'Fontsize',20,'LineWidth',2)
subplot(322)
title('Oceanic mantle')
axis tight, grid on
ylim([-1 1])
xlabel('x_2 (km)')
ylabel('v_z/v_{pl}')
set(gca,'Fontsize',20,'LineWidth',2)

for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    e12dplot = e12d(tindex,:);
    e12dplot(shz.Vpl<0) = 0;
    
    subplot(323)
    toplot = GF_d.shz.Gh*e12dplot';
    plot(ox./1e3,(toplot)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
    subplot(324)
    toplot = GF_d.shz.Gz*e12dplot';
    plot(ox./1e3,(toplot)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
end
subplot(323)
title('Arc mantle')
axis tight, grid on
ylim([-1 1])
xlabel('x_2 (km)')
ylabel('v_h/v_{pl}')
set(gca,'Fontsize',20,'LineWidth',2)
subplot(324)
title('Arc mantle')
axis tight, grid on
ylim([-1 1])
xlabel('x_2 (km)')
ylabel('v_z/v_{pl}')
set(gca,'Fontsize',20,'LineWidth',2)

for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);    
    
    subplot(325)
    toplot = GF_d.rcv.Gh*V(tindex,:)';
    plot(ox./1e3,(toplot)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
    subplot(326)
    toplot = GF_d.rcv.Gz*V(tindex,:)';
    plot(ox./1e3,(toplot)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
end
subplot(325)
title('Fault')
axis tight, grid on
ylim([-1 1])
xlabel('x_2 (km)')
ylabel('v_h/v_{pl}')
set(gca,'Fontsize',20,'LineWidth',2)
subplot(326)
title('Fault')
axis tight, grid on
ylim([-1 1])
xlabel('x_2 (km)')
ylabel('v_z/v_{pl}')
set(gca,'Fontsize',20,'LineWidth',2)
% return
%% make plots of fault+shz slip rates
% timeseries snapshots
index = shz.Vpl>0;
figure(3),clf
subplot(211)
plot(rcv.xc(:,1)./1e3,V(end,:)./Vpl,'r-','LineWidth',3), hold on
plot(shz.xc(index,1)./1e3,e12d(end,index)./Vpl,'r-','LineWidth',3)
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    plot(rcv.xc(:,1)./1e3,V(tindex,:)./Vpl,'-','LineWidth',2,'Color',cspec(i,:))
    plot(shz.xc(index,1)./1e3,e12d(tindex,index)./Vpl,'-','LineWidth',2,'Color',cspec(i,:))
end
axis tight, grid on
xlim([min(shz.xc(:,1)) max(shz.xc(:,1))]./1e3)
xlabel('x_2 (km)')
ylabel('V/V_{pl}')
set(gca,'LineWidth',2,'FontSize',20)

subplot(212)
index = shz.Vpl<0;
plot(shz.xc(index,1)./1e3,-e12d(end,index)./Vpl,'s','LineWidth',1,'MarkerFaceColor',rgb('orangered')), hold on
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    [xsort,isort] = sort(shz.xc(index,1)./1e3);
    toplot = -e12d(tindex,index)./Vpl;
    plot(xsort,toplot(isort),'-','LineWidth',2,'Color',cspec(i,:))
end
axis tight, grid on
% ylim([0.6 1.4])
xlim([min(shz.xc(:,1)) max(shz.xc(:,1))]./1e3)
xlabel('x_2 (km)')
ylabel('V/V_{pl}')
set(gca,'LineWidth',2,'FontSize',20)


% late interseismic sliprates
index = shz.Vpl>0;
figure(4),clf
subplot(211)
plot(rcv.xc(:,1)./1e3,V(end,:)./Vpl,'-','LineWidth',3), hold on
plot(shz.xc(index,1)./1e3,e12d(end,index)./Vpl,'-','LineWidth',3)
axis tight, box on, grid on
xlabel('x_2 (km)')
ylabel('V/V_{pl}')
legend('friction','arc viscous')
set(legend,'location','best')
ylim([0 1])
xlim([min(shz.xc(:,1)) max(shz.xc(:,1))]./1e3)
set(gca,'LineWidth',2,'FontSize',15)

index = shz.Vpl<0;
subplot(212)
plot(shz.xc(index,1)./1e3,e12d(end,index)./Vpl,'s','LineWidth',1,'MarkerFaceColor',rgb('steelblue'))
axis tight, box on, grid on
ylim([-1.1 0])
xlim([min(shz.xc(:,1)) max(shz.xc(:,1))]./1e3)
xlabel('x_2 (km)')
ylabel('V/V_{pl}')
legend('oceanic viscous')
set(legend,'location','best')
set(gca,'LineWidth',2,'FontSize',15,'YDir','reverse')


figure(5),clf
plotpatch2d(rcv,V(end,:)./Vpl), hold on
plotpatch2d(shz,abs(e12d(end,:)./Vpl))
plotpatch2d(src,0.5+0.*abs(src.Vpl./Vpl))
plot(src.x(:,1)./1e3,src.x(:,3)./1e3,'rx')
axis tight equal, box on, grid on
xlim([-1 1].*600)
cb= colorbar;
colormap(jet(20))
caxis([0 1])
set(gca,'LineWidth',2,'FontSize',20)