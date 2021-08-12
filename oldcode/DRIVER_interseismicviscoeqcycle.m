% Testing a viscoelastic earthquake cycle model in planestrain (2d)
% With normal stress perturbation and a maxwell mantle
% Rishav Mallick, 2019, Earth Observatory of Singapore

clear
addpath ~/Dropbox/scripts/utils/
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
topvw = floor(10e3/(fault_width/rcv.N));
botvw = ceil(35e3/(fault_width/rcv.N));
VW(topvw:botvw) = 1;
VW = logical(VW);

% effective confining pressure on fault (MPa) 
rcv.sigma = 50.*ones(rcv.N,1);

% frictional parameters 
rcv.a = 1e-2.*ones(rcv.N,1);
% if you want to set the entire thing as VS
rcv.b = rcv.a - 8e-3;
rcv.b(VW) = rcv.a(VW) + 4e-3;
% static friction coefficient
rcv.mu0 = 0.6*ones(rcv.N,1);
% numer of state parameters
rcv.dgf=5;
% characteristic weakening distance (m)
rcv.l = 0.03.*ones(rcv.N,1);
% plate velocity (m/s)
rcv.Vpl = 3.3e-9*ones(rcv.N,1);

% reference slip rate (m/s)
rcv.Vo = 1e-6*ones(rcv.N,1);
% shear wave speed (m/s)
rcv.Vs = 3e3.*ones(rcv.N,1);

% minimum grid size
fprintf(1,'grid size = %.2f m, minimum grid size = %.2f m\n',rcv.W(1),G/(1-nu)*rcv.l(topvw+1)/rcv.b(topvw+1)/rcv.sigma(topvw+1)) 
fprintf(1,'a/b in VS regions = %.2f\n',rcv.a(1)/rcv.b(1))
fprintf(1,'a/b in VW regions = %.2f\n',rcv.a(topvw+1)/rcv.b(topvw+1))

%% DEEP DRIVING SOURCE
src = unicycle.geometry.receiver('faults/deep_driver.seg',unicycle.greens.okada92(G,nu));
% dummy=unicycle.ode.rateStrengtheningMaxwell([],src,shz,[],'./deepkernels/');
[~,~,KL22] = grdread('deepkernels/KL_d22.grd');
[~,~,KL23] = grdread('deepkernels/KL_d23.grd');
[~,~,KL33] = grdread('deepkernels/KL_d33.grd');
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%               S H E A R   Z O N E S                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% load shear zone
shz=unicycle.geometry.shearZoneReceiver('./faults/lowercrust',unicycle.greens.shearZone16(G,nu));

% setup different viscosity on either side of the dauki
shz.etaM= zeros(shz.N,1) + 1e18/1e6; % MPa s
shz.etaM = 10.^(17 + (3)./(max(abs(shz.xc(:,3))) - min(abs(shz.xc(:,3)))).*(abs(shz.xc(:,3)) - min(abs(shz.xc(:,3)))))/1e6;
% shz.etaM = 10.^(19 - (3)./(max(abs(shz.xc(:,3))) - min(abs(shz.xc(:,3)))).*(abs(shz.xc(:,3)) - min(abs(shz.xc(:,3)))))/1e6;

% create slab as zone of high viscosity
shz.etaM((shz.xc(:,2)<-1.75*shz.xc(:,3) - 50e3) & (shz.xc(:,2)>-1.75*shz.xc(:,3) - 150e3)) = (1e30/1e6);
% %%
% figure(11),clf
% shz.plotShearZoneEdges(log10(shz.etaM*1e6))
% % shz.plotShearZoneEdges(double(Plotshz))
% hold on
% shz.plotShearZoneEdges
% rcv.plotPatch
% src.plotPatch
% view(90,0)
% axis tight equal
% colorbar
% colormap jet
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
remoteloading.dsigmabar_inf = -0*evl.KK{2,3}*rcv.Vpl;

% loading from fault
remoteloading.dsigma22_inf = KL22*(rcv.Vpl(1)*[1 -1 -1]');%-evl.KL{2,4}*rcv.Vpl;
remoteloading.dsigma23_inf = KL23*(rcv.Vpl(1)*[1 -1 -1]');%-evl.KL{2,5}*rcv.Vpl;
remoteloading.dsigma33_inf = KL33*(rcv.Vpl(1)*[1 -1 -1]');%-evl.KL{2,6}*rcv.Vpl;

% % horizontal strain loading
% e22_0 = ((rcv.Vpl(1)./(100e3))*ones(shz.N,1));
% e22_0(shz.xc(:,2)<0) = 0;
% remoteloading.dsigma22_inf = -1*evl.LL{4,4}*e22_0;
% remoteloading.dsigma23_inf = -0*evl.LL{4,5}*e22_0;
% remoteloading.dsigma33_inf = -0*evl.LL{4,6}*e22_0;
%% total duration of simulation
tend = 100*3.15e7;

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
yp=@(t,y) interseismicode_fault_shearzones2d(t,y,evl,remoteloading);
tic
% Solve the system
options=odeset('Refine',1,'RelTol',1e-7,'InitialStep',1e-5,'MaxStep',1e7); 
[t,Y]=ode45(yp,[0 tend],Y0,options);
disp('Finished running ODE solver')
toc

%% extract important quantities

% V = evl.flt.Vo(1).*exp(Y(1:end,4:evl.flt.dgf:evl.flt.N*evl.flt.dgf));
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
V = Yp(:,1:evl.flt.dgf:evl.flt.N*evl.flt.dgf);
e22dot = Yp(:,evl.flt.N*evl.flt.dgf+1:evl.shz.dgf:end);
e23dot = Yp(:,evl.flt.N*evl.flt.dgf+2:evl.shz.dgf:end);
e33dot = Yp(:,evl.flt.N*evl.flt.dgf+3:evl.shz.dgf:end);
edotmax = sqrt((0.5*(e22dot-e33dot)).^2 + (e23dot/2).^2);
ej2dot = e22dot.*e33dot - e23dot.^2;
%% PLOT megathrust behaviour
down = 10; %temporal downslampling
% Plotshz = (true(1,shz.N));
Plotshz = (shz.xc(:,2)<-1.8*shz.xc(:,3) + 10e3) & (shz.xc(:,2)>-1.8*shz.xc(:,3) - 40e3);
% Plotshz = abs(shz.xc(:,2))<=10e3;

% Plotshz(shz.xc(:,2)<-1.9*shz.xc(:,3) - 100e3 & shz.xc(:,2)>-1.9*shz.xc(:,3) - 120e3) = 1;
shzindex = find(Plotshz==1);

% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
subplot(1,2,1)

pcolor(t(1:down:end)./3.15e7,dipdist/1e3,(V(1:down:end,:)')./evl.flt.Vpl(1)), shading flat, hold on

box on
set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
title(h,'$\frac{V}{V_{pl}}$','interpreter','latex','Fontsize',20),ylabel('\zeta_d (km)');
colormap(jet(40))
caxis([0 1])
set(gca,'FontSize',12)
xlabel('Time (yrs)','Interpreter','latex')

subplot(1,2,2)
plot((V(end,:)')./evl.flt.Vpl(1),dipdist./1e3,'k','LineWidth',1)
xlim([0 1])
grid on
xlabel('$\frac{V}{V_{pl}}$','interpreter','latex','Fontsize',20)
ylabel('$\zeta_d$','Interpreter','latex','Fontsize',20)
set(gca,'YDir','reverse')

figure(4),clf
shz.plotShearZoneEdges((edotmax(length(e33dot(:,1)),:).*s2y*1e9)) %nano strain per year
% shz.plotShearZoneEdges(abs(ej2dot(length(e33dot(:,1)),:).*s2y*1e9)) %nano strain per year
hold on
shz.plotShearZoneEdges
rcv.plotPatch, shading flat
view(90,0)
axis tight equal
grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;
h.FontSize = 15;
h.Location='eastoutside';
set(gca,'ColorScale','log')
h.Label.Interpreter='latex';
h.Label.String='$\frac{\partial \epsilon_{max}}{\partial t} 10^{-9}yr^{-1}$';
caxis(10.^[0 2.5])
colormap(flipud(hot))
%% plot of strain and strain rate
cspec = jet(length(shzindex));

figure(2),clf
for i = 1:length(shzindex)
    subplot(2,1,1)
    plot(t./s2y,emax(:,shzindex(i)),'-','LineWidth',1,'Color',cspec(i,:))
    axis tight, hold on
    set(gca,'YScale','log','XScale','log','YLim',10.^[-10 -4],'FontSize',12)
    grid on
    xlim(10.^[-1 3])
    xlabel('Time (yrs)','Interpreter','latex')
    ylabel('\epsilon_{max}')
    
    subplot(2,1,2)
    plot(t./s2y,edotmax(:,shzindex(i)),'-','LineWidth',1,'Color',cspec(i,:))
    axis tight, grid on, hold on
    set(gca,'YScale','log','XScale','log','YLim',10.^[-18 -14],'FontSize',12)
    xlim(10.^[-1 3])
    xlabel('Time (yrs)','Interpreter','latex')
    ylabel('$\dot{\epsilon}_{max}$','Interpreter','latex')
end
%% GPS deformation

gps = unicycle.manifold.gpsReceiver('./gps/gpsnetwork.dat',evl,3);
ushz = [gps.LO{4} gps.LO{5} gps.LO{6}]*[e22';e23';e33'];
vshz = [gps.LO{4} gps.LO{5} gps.LO{6}]*[e22dot';e23dot';e33dot'];

uflt = [gps.KO{2}]*slip';
vflt = [gps.KO{2}]*(V-1*rcv.Vpl(1))';


figure(3),clf
subplot(2,1,1)
plot(gps.x(:,2)./1e3,vflt(2:3:end,end).*s2y.*1e3,'k','LineWidth',1), hold on
plot(gps.x(:,2)./1e3,vshz(2:3:end,end).*s2y.*1e3,'Color',rgb('peru'),'LineWidth',2)
axis tight
grid on
ylabel('V_{fault} mm/yr')
ylim([-1 1].*max(abs(get(gca,'YLim'))))
yyaxis right
plot(gps.x(:,2)./1e3,vshz(2:3:end,end).*s2y.*1e3,'Color',rgb('peru'),'LineWidth',2)
axis tight
ylabel('V_{shear zone} mm/yr')
title('Horizontal Component','Interpreter','latex')
ylim([-1 1].*max(abs(get(gca,'YLim'))))
set(gca,'FontSize',15,'Color','none')

subplot(2,1,2)
plot(gps.x(:,2)./1e3,vflt(3:3:end,end).*s2y.*1e3,'k','LineWidth',1)
axis tight
hold on
grid on
plot(gps.x(:,2)./1e3,vshz(3:3:end,end).*s2y.*1e3,'Color',rgb('peru'),'LineWidth',2)
ylabel('V_{fault} mm/yr')
ylim([-1 1].*max(abs(get(gca,'YLim'))))
yyaxis right
plot(gps.x(:,2)./1e3,vshz(3:3:end,end).*s2y.*1e3,'Color',rgb('peru'),'LineWidth',2)
axis tight
ylabel('V_{shear zone} mm/yr')
title('Vertical Component','Interpreter','latex')
ylim([-1 1].*max(abs(get(gca,'YLim'))))
xlabel('Distance (km)','Interpreter','latex')
set(gca,'FontSize',15,'Color','none')