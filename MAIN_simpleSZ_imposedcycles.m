% Simulates postseismic and interseismic fault creep in response to imposed
% earthquakes at a recgular interval for a subduction zone in 2-d
% Rishav Mallick, EOS, 2021

clear

addpath ODESolving/
addpath odefunction/
addpath greenfunc/
addpath ~/Dropbox/scripts/topotoolbox/colormaps/
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

% Rigidity (MPa)
G = 40e3;
nu = 0.25;

% starting point
y2i = 0;
y3i = 0;

% Fault parameters
Fwidth = 200e3;
Mf = 100;% number of patches
dip = 15;% in degrees

% Viscous domain
Vwidth = 250e3;
Mv = 50;% number of patches
Vwidthbot = 500e3;
Mvbot = 50;% number of patches

% plate thickness
Tplate = 50e3;
%% Create faults and shear zones
earthModel = unicycle.greens.okada92(G,nu);
[rcv,shz,src] = create_flt_shz(earthModel,y2i,y3i,dip,Fwidth,Mf,Vwidth,Mv,Vwidthbot,Mvbot,Tplate);

% plate velocity
Vpl = 1.5e-9; %(m/s)
% Vpl = 3.0420564301467792e-09;

% impose earthquake parameters
Teq = 200.*3.15e7;% earthquake every Teq years
ncycles = 10;%20000*3.15e7/Teq; % number of earthquake cycles

%% Stress Kernels and EVL object

evl = compute_stresskernels(rcv,shz);
disp('Loaded Stress Kernels')


figure(2000),clf
plotpatch2d(rcv,rcv.Vpl), hold on
plotpatch2d(shz,shz.Vpl), shading faceted
% plotpatch2d(src,src.Vpl)
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
rcv.sigma = 100*ones(rcv.N,1);
% rcv.sigma(rcv.xc(:,1)<150e3) = 150;
% rcv.sigma = linspace(100,400,rcv.N)';

% frictional parameters
rcv.a = 1e-2*ones(rcv.N,1);
rcv.b = rcv.a - 5e-3*ones(rcv.N,1);
% velocity-weakening
vw = abs(rcv.xc(:,1))>=0e3 & abs(rcv.xc(:,1))<=350e3;
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
% return
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                   R H E O L O G Y                    %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
% long-term strain rate
shz.Vpl = Vpl.*shz.Vpl;

% Viscosity in Pa-s
% etaval_arc = 5e12;
% etaval_arc = logspace(13,12,length(find(shz.Vpl>0)))';
etaval_arc = 1e12.*linspace(10,1,length(find(shz.Vpl>0)))';
etaval_oc = 1e16;
power = 1;

shz.tMax = power;% store power in 'tmax' parameter
shz.a(shz.Vpl>0) = etaval_arc*1e-6;% store viscosity in 'a' parameter
shz.a(shz.Vpl<0) = etaval_oc*1e-6;

% degrees of freedom
shz.dgf = 2;

% Loading Stresses - backslip loading
evl.tau0     = - evl.KK*rcv.Vpl - evl.LK*shz.Vpl;
evl.sigma0   = - evl.KL*rcv.Vpl - evl.LL*shz.Vpl;

%% solve ODE
disp('Solving ODE')

% Initialize State Vector
Y0=zeros(rcv.N*rcv.dgf+shz.N*shz.dgf,1);

% Calculate shear stress change for IVP
taumax = 5;
% eqslip = compute_earthquakeslip(rcv,evl,Teq*Vpl,taumax);
eqslip = Teq*Vpl.*ones(rcv.N,1);
eqslip(~vw) = 0;%eqslip(find(~vw,5,'first')) = Teq*Vpl;

taueq = evl.KK*eqslip;
taueq(taueq>taumax) = taumax;
% np = 20;
% taueq(find(~vw,1):find(~vw,1)+np) = linspace(0,taueq(find(~vw,1)+np),np+1);
% taueq(taueq<0) = 0;

taubot = evl.KL*eqslip;
taubot(taubot>taumax) = taumax;

if true
    figure(10),clf
    subplot(211)
    plot(rcv.xc(:,1)./1e3,eqslip,'LineWidth',2)
    axis tight
    subplot(212)
    plot(rcv.xc(:,1)./1e3,taueq,'r.'), hold on
    plot(shz.xc(:,1)./1e3,taubot,'kx')
    axis tight, grid on
    ylabel('\Delta\tau (MPa)')
    xlabel('x (km)')
    ylim([-0.2 1].*taumax)
end

% initialize Y0
Y0(1:rcv.dgf:rcv.N*rcv.dgf) = zeros(rcv.N,1);% slip
Y0(2:rcv.dgf:rcv.N*rcv.dgf) = log(rcv.Vpl*1./rcv.Vo); % v

% Shear zones
Y0(rcv.N*rcv.dgf+1:shz.dgf:end) = shz.a.*shz.Vpl + zeros(shz.N,1); %stress12
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
        
        Y0(rcv.N*rcv.dgf+1:shz.dgf:end) = Y0(rcv.N*rcv.dgf+1:shz.dgf:end)' + taubot; %stress12
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
ng = 700;
ox = linspace(-100e3,600e3,ng)';

GF_d = compute_displacementkernels([ox 0.*ox 0.*ox],rcv,shz,src);

% surface velocity
gps = [];
gps.rcv.vh = (GF_d.rcv.Gh*V')';
gps.rcv.vz = (GF_d.rcv.Gz*V')';
gps.shz.vh = (GF_d.shz.Gh*e12d')';
gps.shz.vz = (GF_d.shz.Gz*e12d')';
gps.src.vh = Vpl.*(GF_d.src.Gh*src.Vpl)';
gps.src.vz = Vpl.*(GF_d.src.Gz*src.Vpl)';

%% plot surface displacement timeseries
usurf_h = (GF_d.rcv.Gh*(slip-eqslip')')' + (GF_d.shz.Gh*e12')' + 1*(GF_d.src.Gh*(t*(Vpl.*src.Vpl'))')';
usurf_z = (GF_d.rcv.Gz*(slip-eqslip')')' + (GF_d.shz.Gz*e12')' + 1.*(GF_d.src.Gz*(t*(Vpl.*src.Vpl'))')';

% remove late interseismic effects
% usurf_h = usurf_h - t*(gps.src.vh.*1 + gps.shz.vh(end,:) + gps.rcv.vh(end,:));
% usurf_z = usurf_z - t*(gps.src.vz.*1 + gps.shz.vz(end,:) + gps.rcv.vz(end,:));

%% plot cumulative displacements
if false
    figure(1),clf
    subplot(121)
    pcolor(t./3.15e7,ox./1e3,usurf_h'), shading interp
    xlabel('t (yrs)'),ylabel('x_2 (km)')
    cb=colorbar;cb.Location = 'northoutside';
    cb.Label.String = 'u_h (m)';
    caxis([-1 1]*1)%max(abs(usurf_h(:))))
    xlim([8e4 Teq]./3.15e7)
    set(gca,'FontSize',20,'Linewidth',2,'TickDir','out','XScale','log')
    
    subplot(122)
    pcolor(t./3.15e7,ox./1e3,usurf_z'), shading interp
    xlabel('t (yrs)'),ylabel('x_2 (km)')
    cb=colorbar;cb.Location = 'northoutside';
    cb.Label.String = 'u_z (m)';
    caxis([-1 1]*1)%max(abs(usurf_z(:))))
    colormap(ttscm('vik',20))
    xlim([8e4 Teq]./3.15e7)
    set(gca,'FontSize',20,'Linewidth',2,'TickDir','out','XScale','log')
end
% plot lines - time series
figure(100),clf
plot(ox./1e3,usurf_h(end,:),'k-','Linewidth',2), hold on
plot(ox./1e3,usurf_z(end,:),'r-','Linewidth',2)
axis tight, grid on
% ylim([-1 1]*2)
xlabel('x_2 (km)'), ylabel('u (m)')
set(gca,'FontSize',20,'Linewidth',2,'TickDir','out')

figure(102),clf
plotindex = [0.1,1,10].*3.15e7;
for i = 1:length(plotindex)    
    subplot(length(plotindex),1,i)

    tindex = find(t>plotindex(i),1);
    plot(rcv.xc(:,1)./1e3,slip(tindex,:)'-eqslip-V(end,:)'.*t(tindex),'x-','LineWidth',3), hold on
    plot(shz.xc(:,1)./1e3,e12(tindex,:)-e12d(end,:).*t(tindex),'d','Linewidth',2)
    ylabel('Slip (m)'), xlabel('x (km)')
    axis tight, grid on
    xlim([-100 500]), 
    title(['\Delta t = ' num2str(plotindex(i)./3.15e7,'%.1f'), ' yr'])
    set(gca,'XScale','lin','Fontsize',20,'Linewidth',1.5)
end

figure(101),clf
for i = 1:length(plotindex)
    subplot(length(plotindex),1,i)
    tindex = find(t>plotindex(i),1);
    plot(ox./1e3,-usurf_h(tindex,:),'k-','Linewidth',2), hold on
    plot(ox./1e3,usurf_z(tindex,:),'r-','Linewidth',2)
    axis tight, grid on
    ylabel('Displacement (m)'), xlabel('x (km)')
    xlim([-100 500]),ylim([-1 1]*max(abs(get(gca,'YLim'))))
    title(['\Delta t = ' num2str(plotindex(i)./3.15e7,'%.1f'), ' yr'])
    legend('u_h','u_z')
    set(gca,'XScale','lin','Fontsize',20,'Linewidth',1.5)
end

figure(103),clf
plotindex = [200,250,300,350,400].*1e3;
for i = 1:length(plotindex)
    oxindex = find(ox>plotindex(i),1);
    subplot(221)
    semilogx(t,usurf_h(:,oxindex),'-','Linewidth',2), hold on
    axis tight, grid on
    xlim([8e4, max(t)])
    ylabel('u_h (m)')
    xlabel('\Delta t (s)')
    set(gca,'Fontsize',20,'Linewidth',1.5)
    
    subplot(222)
    plot(t./3.15e7,usurf_h(:,oxindex),'-','Linewidth',2), hold on
    axis tight, grid on
    xlim([0, 10])
    ylabel('u_h (m)')
    xlabel('\Delta t (yr)')
    set(gca,'Fontsize',20,'Linewidth',1.5)
    
    subplot(223)
    semilogx(t,usurf_z(:,oxindex),'-','Linewidth',2), hold on
    axis tight, grid on
    xlim([8e4, max(t)])
    ylabel('u_z (m)')
    xlabel('\Delta t (s)')
    set(gca,'Fontsize',20,'Linewidth',1.5)
    
    subplot(224)
    plot(t./3.15e7,usurf_z(:,oxindex),'-','Linewidth',2), hold on
    axis tight, grid on
    xlim([0, 10])
    ylabel('u_z (m)')
    xlabel('\Delta t (yr)')
    set(gca,'Fontsize',20,'Linewidth',1.5)
end
subplot(223)
legend('200 km','250 km','300 km','350 km','400 km','location','best')
return

%% plot surface velocities
figure(11),clf
% plotindex = [0.05:0.05:0.95].*Teq;
% plotindex = [1,3,5,10,20,30,50,70,80,90,100,150,200].*3.15e7;
plotindex = [10,50,100,150,200].*3.15e7;
p = [];
lgd = {};

cspec = cool(length(plotindex));

subplot(211)
plot(ox./1e3,(gps.src.vh + gps.shz.vh(end,:) + gps.rcv.vh(end,:))./Vpl,'r-','LineWidth',3), hold on
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    p(i) = plot(ox./1e3,(gps.src.vh + gps.shz.vh(tindex,:) + gps.rcv.vh(tindex,:))./Vpl,'-','LineWidth',2,'Color',cspec(i,:));
    lgd{i} = [num2str(round(plotindex(i)./3.15e7)) ' yrs'];
end
plot(rcv.xc(vw,1)./1e3,0.*rcv.Vpl(vw),'k-','LineWidth',3)
axis tight, grid on
% ylim([-1 1])
plot(max(rcv.xc(:,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
plot(max(shz.xc(shz.Vpl>0,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
plot(0.*[1 1],[-1 1],'k-','Linewidth',2)
xlabel('x_2 (km)')
ylabel('v_h/v_{pl}')
legend(p,lgd)
set(legend,'Box','off','location','eastoutside')
set(gca,'Fontsize',20,'LineWidth',2)

subplot(212)
plot(ox./1e3,(gps.src.vz + gps.shz.vz(end,:) + gps.rcv.vz(end,:))./Vpl,'r-','LineWidth',3), hold on
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    p(i) = plot(ox./1e3,(gps.src.vz + gps.shz.vz(tindex,:) + gps.rcv.vz(tindex,:))./Vpl,'-','LineWidth',2,'Color',cspec(i,:));
end
axis tight, grid on
plot(max(rcv.xc(:,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
plot(max(shz.xc(shz.Vpl>0,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
plot(0.*[1 1],[-1 1],'k-','Linewidth',2)
xlabel('x_2 (km)')
ylabel('v_z/v_{pl}')
legend(p,lgd)
ylim([-1 1].*0.7)
set(legend,'Box','off','location','eastoutside')
set(gca,'Fontsize',20,'Linewidth',2)



%% plot each component's effect on the surface
% horizontal motion

% cspec = parula(length(plotindex));

figure(12),clf
set(gcf,'Name','Horizontal velocity contributions')
for i = 1:length(plotindex)
    tindex = find(abs(t-plotindex(i))==min(abs(t-plotindex(i))),1);
    e12dplot = e12d(tindex,:);
    e12dplot(shz.Vpl>0) = 0;
    
    subplot(321)
    toplot = GF_d.shz.Gh*e12dplot';
    plot(ox./1e3,(gps.src.vh' + toplot)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
    subplot(322)
    toplot = GF_d.shz.Gz*e12dplot';
    plot(ox./1e3,(gps.src.vz' + toplot)./Vpl,'-','LineWidth',2,'Color',cspec(i,:)), hold on
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
    plot([rcv.xc(:,1);shz.xc(index,1)]./1e3,[V(tindex,:),e12d(tindex,index)]./Vpl,'-','LineWidth',2,'Color',cspec(i,:))
    %plot(shz.xc(index,1)./1e3,e12d(tindex,index)./Vpl,'-','LineWidth',2,'Color',cspec(i,:))
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
plotpatch2d(src,abs(src.Vpl./Vpl))
plot(src.x(:,1)./1e3,src.x(:,3)./1e3,'rx')
axis tight equal, box on, grid on
xlim([-1 1].*600)
cb= colorbar;
colormap(jet(20))
caxis([0 1])
set(gca,'LineWidth',2,'FontSize',15)