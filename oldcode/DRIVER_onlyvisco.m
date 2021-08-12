% Testing a viscoelastic earthquake cycle model in planestrain (2d)
% Post-seismic relaxation with continuous loading superimposed
% Rishav Mallick, 2019, Earth Observatory of Singapore

clear
addpath ~/Dropbox/scripts/eqphysics/ODESolving/
addpath ~/Dropbox/scripts/topotoolbox/colormaps/
addpath ~/Dropbox/scripts/unicycle/matlab/
addpath ~/Dropbox/scripts/utils/
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

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%               S H E A R   Z O N E S                  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


% load shear zone
shz=unicycle.geometry.shearZoneReceiver('./faults/lowercrust',unicycle.greens.shearZone16(G,nu));
% shz.epsilonPlate = 1;

% uniform viscosity
shz.etaM= zeros(shz.N,1) + 1e18/1e6; % MPa s
% shz.etaM = 10.^(17 + (7)./(max(abs(shz.xc(:,3))) - min(abs(shz.xc(:,3)))).*(abs(shz.xc(:,3)) - min(abs(shz.xc(:,3)))))/1e6;
% shz.etaM = 10.^(19 - (2)./(max(abs(shz.xc(:,3))) - min(abs(shz.xc(:,3)))).*(abs(shz.xc(:,3)) - min(abs(shz.xc(:,3)))))/1e6;

% create slab as zone of high viscosity
% Plotshz = shz.xc(:,2)>1000e3;
Plotshz = (shz.xc(:,2)<-2*shz.xc(:,3) - 40e3) & (shz.xc(:,2)>-2*shz.xc(:,3) - 40e3);
% shz.etaM((shz.xc(:,2)<-1.9*shz.xc(:,3) - 40e3) & (shz.xc(:,2)>-1.9*shz.xc(:,3) - 120e3)) = (1e30/1e6);
shz.etaM(Plotshz) = (1e24/1e6);

% force one side to be low viscosity
% shz.etaM((shz.xc(:,2)<-1.9*shz.xc(:,3) - 40e3)) = (1e22/1e6);

shz.sigma0 = zeros(shz.N,1);
shz.sigma0(Plotshz) = 1;
shz.m = shz.m.*0 + 1;

figure(11),clf
subplot(211)
shz.plotShearZoneEdges((shz.etaM*1e6))
hold on
shz.plotShearZoneEdges
rcv.plotPatch
view(90,0)
axis tight equal, set(gca,'ColorScale','log')
colorbar
colormap jet
% caxis(minmax(log10(shz.etaM*1e6)))

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

%% total duration of simulation

tend = 10*s2y;

%% Initialize State Vector
tot_moment = rcv.N.*10;
coseismic_slip = zeros(rcv.N,1);
nslip = 2;
coseismic_slip(1:nslip) = tot_moment/(nslip);
coseismic_slip(1)
Y0=zeros(evl.shz.N*evl.shz.dgf,1);

% Initial Condition
Y0(4:evl.shz.dgf:evl.shz.dgf*evl.shz.N) = evl.KL{2,1}*coseismic_slip;
Y0(5:evl.shz.dgf:evl.shz.dgf*evl.shz.N) = evl.KL{2,3}*coseismic_slip;
Y0(6:evl.shz.dgf:evl.shz.dgf*evl.shz.N) = evl.KL{2,6}*coseismic_slip;

%% initialize the function handle with
% set constitutive parameters
yp=@(t,y) postseismicode_onlyvisco2d(t,y,evl);
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

e22 = ((Y(:,1:evl.shz.dgf:end)));
e23 = ((Y(:,2:evl.shz.dgf:end)));
e33 = ((Y(:,3:evl.shz.dgf:end)));
s22 = ((Y(:,4:evl.shz.dgf:end)));
s23 = ((Y(:,5:evl.shz.dgf:end)));
s33 = ((Y(:,6:evl.shz.dgf:end)));

emax = sqrt((0.5*(e22-e33)).^2 + (e23).^2);
smax = sqrt((0.5*(s22-s33)).^2 + (s23).^2);
% emax = abs(e22+e33);

ej2 = abs(sqrt(e22.*e33 - e23.^2));

Yp = zeros(size(Y));
for i = 1:length(t)
    Yp(i,:) = yp(t(i),Y(i,:)');
end
e22dot = Yp(:,1:evl.shz.dgf:end);
e23dot = Yp(:,2:evl.shz.dgf:end);
e33dot = Yp(:,3:evl.shz.dgf:end);
s22dot = Yp(:,4:evl.shz.dgf:end);
s23dot = Yp(:,5:evl.shz.dgf:end);
s33dot = Yp(:,6:evl.shz.dgf:end);
sdotmax = sqrt((0.5*(s22dot-s33dot)).^2 + (s23dot).^2);
edotmax = sqrt((0.5*(e22dot-e33dot)).^2 + (e23dot).^2);
ej2dot = sqrt(e22dot.*e33dot - e23dot.^2);
%% PLOT megathrust behaviour
down = 1; %temporal downslampling
% Plotshz = (true(1,shz.N));
% Plotshz = (shz.xc(:,2)<-1.9*shz.xc(:,3) - 30e3) & (shz.xc(:,2)>-1.9*shz.xc(:,3) - 40e3);
% Plotshz = abs(shz.xc(:,2))<=10e3;

% Plotshz(shz.xc(:,2)<-1.9*shz.xc(:,3) - 100e3 & shz.xc(:,2)>-1.9*shz.xc(:,3) - 120e3) = 1;
% shzindex = find(Plotshz==1);

x2pos = 60e3;
shzindex = find(abs(shz.xc(:,2)-x2pos) == min(abs(shz.xc(:,2)-x2pos)));
% shzindex = 1:shz.N;

% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
subplot(2,1,1)
pcolor(t(1:down:end)./3.15e7,shz.xc(shzindex,3)./1e3,log10(edotmax(1:down:end,shzindex)'))
shading flat
cb=colorbar;
cb.Location = 'northoutside';
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\log_{10}\dot{\epsilon}_{max}$';
cb.Label.FontSize = 15;
grid on
xlabel('Time (yrs)','Interpreter','latex')
ylabel('$\zeta_d$ (km)','Interpreter','latex','Fontsize',12)
set(gca,'FontSize',12,'XScale','linear')
caxis([-16 -12])

x2pos = -50e3;
shzindex = find(abs(shz.xc(:,2)-x2pos) == min(abs(shz.xc(:,2)-x2pos)));
% shzindex = [1:shz.N]';

subplot(2,1,2)
pcolor(t(1:down:end)./3.15e7,shz.xc(shzindex,3)./1e3,log10(edotmax(1:down:end,shzindex)'))
shading flat
cb=colorbar;
cb.Location = 'northoutside';
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\log_{10}\dot{\epsilon}_{max}$';
cb.Label.FontSize = 15;
grid on
xlabel('Time (yrs)','Interpreter','latex')
ylabel('$\zeta_d$ (km)','Interpreter','latex','Fontsize',12)
set(gca,'FontSize',12,'XScale','linear')
colormap(hot)
caxis([-16 -12])

% return
%% plot strain at snapshots
figure(4),clf

% tvec = [0.1,1,5,10,13];
tvec = [1,5,10,25,100];

for i = 1:length(tvec)
    tind = find(abs(t-tvec(i)*3.15e7) == min(abs(t-tvec(i)*3.15e7)));

    toplot = log10(emax);
    toplot(:,Plotshz) = nan;
    % strains
    subplot(length(tvec),2,2*(i-1) + 1)
    shz.plotShearZoneEdges((toplot(tind,:)))
    hold on
    %shz.plotShearZoneEdges
    rcv.plotPatch(nan(rcv.N,1))
    view(90,0)
    axis tight equal,grid off, box on
    xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
    h=colorbar;h.FontSize = 12;h.Location='eastoutside';
    set(gca,'Color','none');
    h.Label.Interpreter='latex';h.Label.String='$\log_{10} \epsilon_{max} $';
    title(['Time = ' num2str(t(tind)./s2y) ' yrs'])
    %caxis([1e-2 1].*max(abs(get(gca,'CLim'))))
    %caxis([1 1].*max((get(gca,'CLim'))) + [-2 0])
    caxis([-5.5 -3])
    %set(gca,'ColorScale','log')
    
    toplot = log10(edotmax);
    toplot(:,Plotshz) = nan;
    %strain rates
    subplot(length(tvec),2,2*i)
    shz.plotShearZoneEdges(((toplot(tind,:))))
    hold on
    rcv.plotPatch(nan(rcv.N,1))
    view(90,0)
    axis tight equal,grid off, box on
    xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
    h=colorbar;h.FontSize = 12;h.Location='eastoutside';%set(gca,'ColorScale','log');
    h.Label.Interpreter='latex';h.Label.String='$\log_{10} \dot{\epsilon}_{max} $';
    %caxis([1e-2 1].*max(abs(get(gca,'CLim'))))
    %caxis([1 1].*max((get(gca,'CLim'))) + [-2 0])
    caxis([-15 -12])
    title(['Time = ' num2str(t(tind)./s2y) ' yrs'])
    %set(gca,'ColorScale','log')
    set(gca,'Color','none');
end

% colormap(jet(20))
cspec = (magmacolor(9));
cspec(length(cspec(:,1))+1,:) = [1 1 1];
colormap(flipud(cspec))
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% plot of strain and strain rate
cspec = jet(length(shzindex));

figure(2),clf
for i = 1:length(shzindex)
    subplot(2,1,1)
    plot(t./s2y,emax(:,shzindex(i)),'-','LineWidth',1,'Color',cspec(i,:))
    axis tight, hold on
    set(gca,'YScale','log','XScale','log','FontSize',15)
    grid on
    xlim(10.^[-3 2])
    xlabel('Time (yrs)')
    ylabel('\epsilon_{max}')
    
    subplot(2,1,2)
    plot(t./s2y,edotmax(:,shzindex(i)),'-','LineWidth',1,'Color',cspec(i,:))
    axis tight, grid on, hold on
    set(gca,'YScale','log','XScale','log','FontSize',15)
    xlim(10.^[-3 2])
    xlabel('Time (yrs)')
    ylabel('$\dot{\epsilon}_{max}$','Interpreter','latex')
end

%% do simple matrix inversion
Aeq = [1.*eye(shz.N) eye(shz.N).*0 1.*eye(shz.N)];
beq = zeros(shz.N,1);
% Aeq = [];
% beq = [];
% if length(find(Plotshz==1))>0
%     index = find(Plotshz==1);
%     for i = 1:length(index)
%         dummy = zeros(1,shz.N);
%         dummy(index(i)) = 1;
%         Aeq = [Aeq;...
%             dummy zeros(1,2*shz.N);...
%             zeros(1,shz.N) dummy zeros(1,shz.N)];
%         beq = [beq;0;0];
%     end
% end
Gstress = [evl.LL{1,1} evl.LL{3,1} evl.LL{6,1};evl.LL{1,3} evl.LL{3,3} evl.LL{6,3};evl.LL{1,6} evl.LL{3,6} evl.LL{6,6}];
costress = [evl.KL{2,1}*coseismic_slip;evl.KL{2,3}*coseismic_slip;evl.KL{2,6}*coseismic_slip];
einv = lsqlin(Gstress,...
    -costress,[],[],...
    Aeq,beq);

% einv = lsqlin([evl.LL{4,4} evl.LL{5,4} evl.LL{6,4};evl.LL{4,5} evl.LL{5,5} evl.LL{6,5};evl.LL{4,6} evl.LL{5,6} evl.LL{6,6}],...
%     [evl.KL{2,4}*coseismic_slip;evl.KL{2,5}*coseismic_slip;evl.KL{2,6}*coseismic_slip],[],[],...
%     Aeq,beq);

figure(15),clf
subplot(2,1,1)
shz.plotShearZoneEdges(einv(1:shz.N))
% shz.plotShearZoneEdges(einv(2*shz.N+1:3*shz.N))
% shz.plotShearZoneEdges(costress(1:shz.N))
hold on
rcv.plotPatch(nan(rcv.N,1)), shading faceted
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
h.Label.Interpreter='latex';h.Label.String='$\epsilon_{22}$ or -$\epsilon_{33} $';
caxis([-1 1].*max(abs(get(gca,'CLim'))))
colormap bluewhitered(20)

subplot(2,1,2)
shz.plotShearZoneEdges(einv(shz.N+1:2*shz.N))
% shz.plotShearZoneEdges(costress(shz.N+1:2*shz.N))
hold on
rcv.plotPatch(nan(rcv.N,1)), shading faceted
view(90,0)
axis tight equal,grid on, box on
xlabel('X (m)'), ylabel('Y (m)'), zlabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
h.Label.Interpreter='latex';h.Label.String='$\epsilon_{23} $';
caxis([-1 1].*max(abs(get(gca,'CLim'))))

set(findall(gcf,'-property','FontSize'),'FontSize',15)
%% GPS deformation

gps = unicycle.manifold.gpsReceiver('./gps/gpsnetwork.dat',evl,3);

ushz = [gps.LO{1} gps.LO{3} gps.LO{6}]*[e22';e23';e33'];
vshz = [gps.LO{1} gps.LO{3} gps.LO{6}]*[e22dot';e23dot';e33dot'];

ufinal = [gps.LO{1} gps.LO{3} gps.LO{6}]*[einv(1:shz.N);einv(shz.N+1:2*shz.N);-einv(1:shz.N)];
tvec = [0:1:10];

cspec = jet(length(tvec));

figure(3),clf
for i = 1:length(tvec)
    tplot = find(abs(t-tvec(i)*3.15e7) == min(abs(t-tvec(i)*3.15e7)));    
    
    subplot(2,2,1)
    plot(gps.x(:,2)./1e3,1e3*ushz(2:3:end,tplot),'Color',cspec(i,:),'LineWidth',2), hold on
    plot(gps.x(:,2)./1e3,1e3*ufinal(2:3:end),'k-','LineWidth',2)
    plot(rcv.x(1,2)./1e3.*[1 1],get(gca,'YLim'),'k-')
    plot(rcv.x(end,2)./1e3.*[1 1],get(gca,'YLim'),'k-')
    axis tight
    ylabel('u (mm)')
    title('Horizontal Component')
    %ylim([-1 1].*max(abs(get(gca,'YLim'))))
    set(gca,'FontSize',15,'Color','none')
    xlim([-500 500])
    
    subplot(2,2,2)
    plot(gps.x(:,2)./1e3,1000*ushz(3:3:end,tplot),'Color',cspec(i,:),'LineWidth',2), hold on
    plot(gps.x(:,2)./1e3,1e3*ufinal(2:3:end),'k-','LineWidth',2)
    plot(rcv.x(1,2)./1e3.*[1 1],get(gca,'YLim'),'k-')
    plot(rcv.x(end,2)./1e3.*[1 1],get(gca,'YLim'),'k-')
    axis tight
    ylabel('u (mm)')
    title('Vertical Component')
    %ylim([-1 1].*max(abs(get(gca,'YLim'))))
    xlabel('Distance (km)')
    set(gca,'FontSize',15,'Color','none')
    xlim([-500 500])
    
    subplot(2,2,3)
    plot(gps.x(:,2)./1e3,s2y*1e3*vshz(2:3:end,tplot),'Color',cspec(i,:),'LineWidth',2), hold on
    plot(rcv.x(1,2)./1e3.*[1 1],get(gca,'YLim'),'k-')
    plot(rcv.x(end,2)./1e3.*[1 1],get(gca,'YLim'),'k-')
    axis tight
    ylabel('V (mm/yr)')
    title('Horizontal Component')
    %ylim([-1 1].*max(abs(get(gca,'YLim'))))
    set(gca,'FontSize',15,'Color','none')
    xlim([-500 500])
    
    subplot(2,2,4)
    plot(gps.x(:,2)./1e3,s2y*1e3*vshz(3:3:end,tplot),'Color',cspec(i,:),'LineWidth',2), hold on
    plot(rcv.x(1,2)./1e3.*[1 1],get(gca,'YLim'),'k-')
    plot(rcv.x(end,2)./1e3.*[1 1],get(gca,'YLim'),'k-')
    axis tight
    ylabel('V (mm/yr)')
    title('Vertical Component')
    %ylim([-1 1].*max(abs(get(gca,'YLim'))))
    xlabel('Distance (km)')
    set(gca,'FontSize',15,'Color','none')
    xlim([-500 500])
end
cb=colorbar; cb.Label.String = 'Time (yrs)'; cb.Location = 'east';
colormap(cspec)
caxis([min(tvec) max(tvec)])
