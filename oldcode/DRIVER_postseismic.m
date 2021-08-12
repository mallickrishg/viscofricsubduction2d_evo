% Testing a viscoelastic earthquake cycle model in planestrain (2d)
% Post-seismic relaxation with continuous loading superimposed
% Rishav Mallick, 2019, Earth Observatory of Singapore

clear
addpath ~/Dropbox/scripts/eqphysics/ODESolving/
addpath ~/Dropbox/scripts/topotoolbox/colormaps/
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
topvw = ceil(0.1e3/(fault_width/rcv.N));
botvw = ceil(56e3/(fault_width/rcv.N));
VW(topvw:botvw) = 1;
VW = logical(VW);

% effective confining pressure on fault (MPa) 
rcv.sigma = 50.*ones(rcv.N,1);

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
rcv.l = 0.01.*ones(rcv.N,1);
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

% uniform viscosity
shz.etaM= zeros(shz.N,1) + 1e18/1e6; % MPa s
% shz.etaM = 10.^(17 + (7)./(max(abs(shz.xc(:,3))) - min(abs(shz.xc(:,3)))).*(abs(shz.xc(:,3)) - min(abs(shz.xc(:,3)))))/1e6;
% shz.etaM = 10.^(19 - (3)./(max(abs(shz.xc(:,3))) - min(abs(shz.xc(:,3)))).*(abs(shz.xc(:,3)) - min(abs(shz.xc(:,3)))))/1e6;

% create slab as zone of high viscosity
Plotshz = (shz.xc(:,2)<-2*shz.xc(:,3) - 30e3) & (shz.xc(:,2)>-2*shz.xc(:,3) - 71e3);
% shz.etaM((shz.xc(:,2)<-1.9*shz.xc(:,3) - 40e3) & (shz.xc(:,2)>-1.9*shz.xc(:,3) - 120e3)) = (1e30/1e6);
shz.etaM(Plotshz) = (1e38/1e6);

% force one side to be low viscosity
% shz.etaM((shz.xc(:,2)<-1.9*shz.xc(:,3) - 40e3)) = (1e22/1e6);

shz.sigma0 = zeros(shz.N,1);
shz.sigma0(Plotshz) = 1;
shz.m = shz.m.*0 + 1;

figure(11),clf
subplot(211)
shz.plotShearZoneEdges(log10(shz.etaM*1e6))
hold on
shz.plotShearZoneEdges
rcv.plotPatch
view(90,0)
axis tight equal
colorbar
colormap jet
subplot(212)
shz.plotShearZoneEdges(double(Plotshz)), hold on
shz.plotShearZoneEdges
rcv.plotPatch
view(90,0)
axis tight equal
colorbar
colormap jet
drawnow

% return
%% load kernels
tic
if ~exist('evl','var')
    evl=unicycle.ode.rateStrengtheningMaxwell([],rcv,shz,[],'./kernels/');
    evl.flt.dgf = 2; 
    evl.shz.dgf = 6;
    dgf = evl.flt.dgf + evl.shz.dgf;
end
disp('Loaded Stress Kernels')
toc
remoteloading = [];
% remote tectonic loading on fault (weird units - MPa/s)
remoteloading.dtau_inf = -0*evl.KK{2,2}*rcv.Vpl;
% remoteloading.dsigmabar_inf = -evl.KK{2,3}*rcv.Vpl;

remoteloading.dsigma22_inf = -0*evl.KL{2,1}*rcv.Vpl;
remoteloading.dsigma23_inf = -0*evl.KL{2,3}*rcv.Vpl;
remoteloading.dsigma33_inf = -0*evl.KL{2,6}*rcv.Vpl;


%% total duration of simulation

tend = 30*s2y;

%% Initialize State Vector
coseismic_slip = zeros(rcv.N,1);
coseismic_slip(10:12) = 100;
Y0=zeros(evl.flt.N*evl.flt.dgf + evl.shz.N*evl.shz.dgf,1);
% Fault patches
% Y0(1:evl.flt.dgf:evl.flt.dgf*evl.flt.N) = zeros(rcv.N,1);         % slip
tau_cs = (evl.KK{2,2}*coseismic_slip)./(rcv.a.*rcv.sigma);
tau_cs(tau_cs<0) = 0;
tau_cs(tau_cs>10) = 10;
% return
Y0(2:evl.flt.dgf:evl.flt.dgf*evl.flt.N) = log(rcv.Vpl./rcv.Vo) + tau_cs;

% Shear zones
Y0(evl.flt.dgf*evl.flt.N + 4:evl.shz.dgf:evl.flt.dgf*evl.flt.N+evl.shz.dgf*evl.shz.N) = evl.KL{2,1}*coseismic_slip;
Y0(evl.flt.dgf*evl.flt.N + 5:evl.shz.dgf:evl.flt.dgf*evl.flt.N+evl.shz.dgf*evl.shz.N) = evl.KL{2,3}*coseismic_slip;
Y0(evl.flt.dgf*evl.flt.N + 6:evl.shz.dgf:evl.flt.dgf*evl.flt.N+evl.shz.dgf*evl.shz.N) = evl.KL{2,6}*coseismic_slip;

figure(12),clf
shz.plotShearZoneEdges(evl.KL{2,3}*coseismic_slip)
% dummyshz = 0*ones(shz.N,1);
% dummyshz(375) = 1;
% shz.plotShearZoneEdges(evl.LL{3,6}*dummyshz + 0*evl.LL{3,4}*dummyshz + 0*evl.LL{3,6}*dummyshz)
hold on
shz.plotShearZoneEdges
rcv.plotPatch
view(90,0)
axis tight equal
colorbar
colormap bluewhitered
set(gca,'FontSize',15)
drawnow
% return
%% initialize the function handle with
% set constitutive parameters
yp=@(t,y) postseismicode_fault_shearzones2d(t,y,evl,remoteloading);
tic
% Solve the system
% options=odeset('Refine',1,'RelTol',1e-4,'AbsTol',1e-4,'InitialStep',1e-3,'MaxStep',1e7,'ODir','./odeout');
options=odeset('Refine',1,'RelTol',1e-6,'AbsTol',1e-6,'InitialStep',1e-3,'MaxStep',1e7);
tvec = [0,tend];
[t,Y]=ode45(yp,tvec,Y0,options);
% t=t';
% Y=Y';
disp('Finished running ODE solver')
toc


%% extract important quantities

% V = evl.flt.Vo(1).*exp(Y(1:end,4:evl.flt.dgf:evl.flt.N*evl.flt.dgf));
slip = Y(:,1:evl.flt.dgf:evl.flt.N*evl.flt.dgf);

e22 = Y(:,evl.flt.N*evl.flt.dgf+1:evl.shz.dgf:end);
e23 = Y(:,evl.flt.N*evl.flt.dgf+2:evl.shz.dgf:end);
e33 = Y(:,evl.flt.N*evl.flt.dgf+3:evl.shz.dgf:end);
emax = sqrt((0.5*(e22-e33)).^2 + (e23).^2);
% emax = abs(e22+e33);

ej2 = abs(sqrt(e22.*e33 - e23.^2));

Yp = zeros(size(Y));
for i = 1:length(t)
    Yp(i,:) = yp(t(i),Y(i,:)');
end
V = Yp(:,1:evl.flt.dgf:evl.flt.N*evl.flt.dgf);
e22dot = Yp(:,evl.flt.N*evl.flt.dgf+1:evl.shz.dgf:end);
e23dot = Yp(:,evl.flt.N*evl.flt.dgf+2:evl.shz.dgf:end);
e33dot = Yp(:,evl.flt.N*evl.flt.dgf+3:evl.shz.dgf:end);
edotmax = sqrt((0.5*(e22dot-e33dot)).^2 + (e23dot).^2);
ej2dot = sqrt(e22dot.*e33dot - e23dot.^2);
%% PLOT megathrust behaviour
down = 1; %temporal downslampling
% Plotshz = (true(1,shz.N));
% Plotshz = (shz.xc(:,2)<-1.9*shz.xc(:,3) - 30e3) & (shz.xc(:,2)>-1.9*shz.xc(:,3) - 40e3);
% Plotshz = abs(shz.xc(:,2))<=10e3;

% Plotshz(shz.xc(:,2)<-1.9*shz.xc(:,3) - 100e3 & shz.xc(:,2)>-1.9*shz.xc(:,3) - 120e3) = 1;
% shzindex = find(Plotshz==1);
shzindex = find(abs(shz.xc(:,2))<15e3);

% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
subplot(2,1,1)

pcolor(t(1:down:end)./3.15e7,dipdist/1e3,log10(V(1:down:end,:)'./evl.flt.Vpl(1))), shading flat, hold on
box on
set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
title(h,'$\log_{10}\frac{V}{V_{pl}}$','interpreter','latex','Fontsize',20),ylabel('\zeta_d (km)');
colormap(jet(50))
caxis([-1 2])
set(gca,'FontSize',12,'XScale','log')
% set(gca,'ColorScale','log')
xlabel('Time (yrs)','Interpreter','latex')

subplot(2,1,2)
pcolor(t(1:down:end)./3.15e7,shz.xc(shzindex,3)./1e3,log10(edotmax(1:down:end,shzindex)'))
shading flat
% plot((V(end,:)')./evl.flt.Vpl(1),dipdist./1e3,'k','LineWidth',1)
cb=colorbar;
cb.Location = 'northoutside';
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\log_{10}\dot{\epsilon}_{max}$';
cb.Label.FontSize = 15;
% set(gca,'ColorScale','log')
grid on
% caxis(10.^[-17 -13])
% caxis([-17 -13])
% xlabel('$\frac{V}{V_{pl}}$','interpreter','latex','Fontsize',20)
xlabel('Time (yrs)','Interpreter','latex')
ylabel('$\zeta_d$ (km)','Interpreter','latex','Fontsize',12)
set(gca,'FontSize',12,'XScale','log')

figure(10),clf
plot(rcv.x(:,2)./1e3,slip(end,:),'k-','LineWidth',2), hold on
yyaxis right
plot(rcv.x(:,2)./1e3,coseismic_slip,'-','LineWidth',2)
axis tight
xlabel('X (km)');
ylabel('Slip (m)')

%% plot strain at snapshots
figure(4),clf
t1 = 0.1;
t2 = 2;
t3 = 20;

t1ind = find(abs(t-t1*3.15e7) == min(abs(t-t1*3.15e7)));
t2ind = find(abs(t-t2*3.15e7) == min(abs(t-t2*3.15e7)));
t3ind = find(abs(t-t3*3.15e7) == min(abs(t-t3*3.15e7)));

toplot = e33;
toplot(:,Plotshz) = nan;
% strains
subplot(3,2,5)
shz.plotShearZoneEdges((toplot(t3ind,:)))
hold on
%shz.plotShearZoneEdges
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
% set(gca,'ColorScale','log');
h.Label.Interpreter='latex';h.Label.String='$\log_{10} \epsilon_{max} $';
title(['Time = ' num2str(t(t3ind)./s2y) ' yrs'])
caxis([-1 1].*max(abs(get(gca,'CLim'))))

subplot(3,2,1)
shz.plotShearZoneEdges((toplot(t1ind,:)))
hold on
% shz.plotShearZoneEdges
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
% set(gca,'ColorScale','log');
h.Label.Interpreter='latex';h.Label.String='$\log_{10} \epsilon_{max} $';
title(['Time = ' num2str(t(t1ind)./s2y) ' yrs'])
caxis([-1 1].*max(abs(get(gca,'CLim'))))

subplot(3,2,3)
shz.plotShearZoneEdges((toplot(t2ind,:)))
hold on
% shz.plotShearZoneEdges
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
% set(gca,'ColorScale','log');
h.Label.Interpreter='latex';h.Label.String='$\log_{10} \epsilon_{max} $';
title(['Time = ' num2str(t(t2ind)./s2y) ' yrs'])
caxis([-1 1].*max(abs(get(gca,'CLim'))))

toplot = (e33dot);
toplot(:,Plotshz) = nan;
%strain rates
subplot(3,2,2)
shz.plotShearZoneEdges(((toplot(t1ind,:))))
hold on
% shz.plotShearZoneEdges
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';%set(gca,'ColorScale','log');
h.Label.Interpreter='latex';h.Label.String='$\log_{10} \dot{\epsilon}_{max} $';
caxis([-1 1].*max(abs(get(gca,'CLim'))))

subplot(3,2,4)
shz.plotShearZoneEdges((toplot(t2ind,:)))
hold on
% shz.plotShearZoneEdges
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';%set(gca,'ColorScale','log');
h.Label.Interpreter='latex';h.Label.String='$\log_{10} \dot{\epsilon}_{max} $';
caxis([-1 1].*max(abs(get(gca,'CLim'))))

subplot(3,2,6)
shz.plotShearZoneEdges((toplot(t3ind,:)))
hold on
% shz.plotShearZoneEdges
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,
grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';%set(gca,'ColorScale','log');
h.Label.Interpreter='latex';h.Label.String='$\log_{10} \dot{\epsilon}_{max} $';
caxis([-1 1].*max(abs(get(gca,'CLim'))))

colormap(jet)
% colormap(flipud(magmacolor(20)))
set(findall(gcf,'-property','FontSize'),'FontSize',15)

%% plot of strain and strain rate
cspec = jet(length(shzindex));

figure(2),clf
for i = 1:length(shzindex)
    subplot(2,1,1)
    plot(t./s2y,emax(:,shzindex(i)),'-','LineWidth',1,'Color',cspec(i,:))
    axis tight, hold on
    set(gca,'YScale','log','XScale','log','FontSize',12)
    grid on
    xlim(10.^[-3 2])
    xlabel('Time (yrs)','Interpreter','latex')
    ylabel('\epsilon_{max}')
    
    subplot(2,1,2)
    plot(t./s2y,edotmax(:,shzindex(i)),'-','LineWidth',1,'Color',cspec(i,:))
    axis tight, grid on, hold on
    set(gca,'YScale','log','XScale','log','FontSize',12)
    xlim(10.^[-3 2])
    xlabel('Time (yrs)','Interpreter','latex')
    ylabel('$\dot{\epsilon}_{max}$','Interpreter','latex')
end
%% GPS deformation

gps = unicycle.manifold.gpsReceiver('./gps/gpsnetwork.dat',evl,3);

% e22 = e22.*0;
% e33 = e33.*0;
% e23 = e23.*0;
% 
% e22(end,100) = 1;
% e33(end,100) = 1;

ushz = [gps.LO{1} gps.LO{3} gps.LO{6}]*[e22';e23';e33'];
vshz = [gps.LO{1} gps.LO{3} gps.LO{6}]*[e22dot';e23dot';e33dot'];

uflt = [gps.KO{2}]*slip';
% vflt = [gps.KO{2}]*(V-1*rcv.Vpl(1))';
vflt = [gps.KO{2}]*(V)';

tplot = find(abs(t-100*3.15e7) == min(abs(t-100*3.15e7)));

figure(3),clf
subplot(2,1,1)
plot(gps.x(:,2)./1e3,1e3*uflt(2:3:end,tplot),'k','LineWidth',1), hold on
axis tight
grid on
ylabel('u_{x} (mm)')
ylim([-1 1].*max(abs(get(gca,'YLim'))))
xlim([-100 500])

% yyaxis right
plot(gps.x(:,2)./1e3,1e3*ushz(2:3:end,tplot),'Color',rgb('steelblue'),'LineWidth',2)
axis tight
% ylabel('u_{shear zone} (m)')
title('Horizontal Component')
ylim([-1 1].*max(abs(get(gca,'YLim'))))
set(gca,'FontSize',15,'Color','none')
xlim([-100 500])

subplot(2,1,2)
plot(gps.x(:,2)./1e3,1000*uflt(3:3:end,end),'k','LineWidth',1)
axis tight
hold on
grid on
ylabel('u_{z} (mm)')
ylim([-1 1].*max(abs(get(gca,'YLim'))))
xlim([-100 500])

% yyaxis right
plot(gps.x(:,2)./1e3,1000*ushz(3:3:end,tplot),'Color',rgb('steelblue'),'LineWidth',2)
axis tight
% ylabel('u_{shear zone} (m)')
title('Vertical Component')
ylim([-1 1].*max(abs(get(gca,'YLim'))))
xlabel('Distance (km)')
set(gca,'FontSize',15,'Color','none')
xlim([-100 500])

%% do simple matrix inversion
Aeq = [1.*eye(shz.N) eye(shz.N).*0 1.*eye(shz.N)];
beq = zeros(shz.N,1);

einv = lsqlin([evl.LL{1,1} evl.LL{3,1} evl.LL{6,1};evl.LL{1,3} evl.LL{3,3} evl.LL{6,3};evl.LL{1,6} evl.LL{3,6} evl.LL{6,6}],...
    [evl.KL{2,1}*coseismic_slip;evl.KL{2,3}*coseismic_slip;evl.KL{2,6}*coseismic_slip],[],[],...
    Aeq,beq);

% einv = lsqlin([evl.LL{4,4} evl.LL{5,4} evl.LL{6,4};evl.LL{4,5} evl.LL{5,5} evl.LL{6,5};evl.LL{4,6} evl.LL{5,6} evl.LL{6,6}],...
%     [evl.KL{2,4}*coseismic_slip;evl.KL{2,5}*coseismic_slip;evl.KL{2,6}*coseismic_slip],[],[],...
%     Aeq,beq);

figure(15),clf
subplot(3,1,1)
shz.plotShearZoneEdges(einv(1:shz.N))
hold on
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
h.Label.Interpreter='latex';h.Label.String='$\epsilon_{22} $';
colormap jet
caxis([-1 1].*max(abs(get(gca,'CLim'))))

subplot(3,1,2)
shz.plotShearZoneEdges(einv(shz.N+1:2*shz.N))
hold on
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
h.Label.Interpreter='latex';h.Label.String='$\epsilon_{23} $';
caxis([-1 1].*max(abs(get(gca,'CLim'))))

subplot(3,1,3)
shz.plotShearZoneEdges(einv(2*shz.N+1:end))
hold on
rcv.plotPatch(nan(rcv.N,1)), shading flat
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
h.Label.Interpreter='latex';h.Label.String='$\epsilon_{33} $';
caxis([-1 1].*max(abs(get(gca,'CLim'))))

set(findall(gcf,'-property','FontSize'),'FontSize',15)
