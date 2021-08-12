clear
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*
addpath ~/Dropbox/scripts/utils/
s2y=60*60*24*365;

x2extent = 100e3;
Nx2 = 20;

x3extent = 20e3;
Nx3 = 3;
x3shift = 50e3;
ngps = 200;

% ox = linspace(-x2extent,x2extent,ngps);
ox = linspace(-100,100,ngps).*1e3;

shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);

%% load kernels for displacement and stress
tic
evl = compute_shzkernels(shz,ox);
toc
%% plot geometry
dummy = zeros(shz.N,1);
dummy(26) = 1e-3;

toplot = evl.LL2323*dummy;

figure(1),clf
subplot(3,1,[2 3])
plotshz(shz,toplot)
axis tight equal
colormap bluewhitered
cb=colorbar;cb.Location = 'northoutside';

subplot(3,1,1)
plot(ox,evl.G23_2*dummy,'b-','LineWidth',1)
hold on
plot(ox,evl.G23_3*dummy,'r-','LineWidth',1)
axis tight
legend('horizontal','vertical')
% return
%% viscoelastic simulation

% uniform viscosity
shz.etaM = zeros(shz.N,1) + 1e17/1e6; % MPa s
shz.m = 3;
shz.coupling = 1;
evl.dgf = 6;
%% total duration of simulation
tend = 100*s2y;

% Initialize State Vector
% initial_edot = zeros(shz.N,1);
% index = 48;
% initial_edot(index) = 1e-3;
initial_strain = dummy;

Y0=zeros(shz.N*evl.dgf,1);

% Initial Condition
Y0(4:evl.dgf:evl.dgf*shz.N) = evl.LL3322*initial_strain;
Y0(5:evl.dgf:evl.dgf*shz.N) = evl.LL3323*initial_strain;
Y0(6:evl.dgf:evl.dgf*shz.N) = evl.LL3333*initial_strain;

%% initialize the function handle with
% set constitutive parameters
yp=@(t,y) ode_post_visco2d(t,y,shz,evl);
tic
% Solve the system
% options=odeset('Refine',1,'RelTol',1e-4,'AbsTol',1e-4,'InitialStep',1e-3,'MaxStep',1e7,'ODir','./odeout');
options=odeset('Refine',1,'RelTol',1e-10,'AbsTol',1e-7,'InitialStep',1e-6,'MaxStep',1e6);
tvec = [0,tend];
% [Y0,Yp0] = decic(yp,0,Y0,1,Y0,0);
[t,Y]=ode45(yp,tvec,Y0,options);
% t=t';
% Y=Y';
disp('Finished running ODE solver')
toc

%% Extract important quantities
e22 = ((Y(:,1:evl.dgf:end)));
e23 = ((Y(:,2:evl.dgf:end)));
e33 = ((Y(:,3:evl.dgf:end)));
s22 = ((Y(:,4:evl.dgf:end)));
s23 = ((Y(:,5:evl.dgf:end)));
s33 = ((Y(:,6:evl.dgf:end)));

emax = sqrt((0.5*(e22-e33)).^2 + (e23).^2);
smax = sqrt((0.5*(s22-s33)).^2 + (s23).^2);
% emax = abs(e22+e33);

ej2 = abs(sqrt(e22.*e33 - e23.^2));

Yp = zeros(size(Y));
for i = 1:length(t)
    Yp(i,:) = yp(t(i),Y(i,:)');
end
e22dot = Yp(:,1:evl.dgf:end);
e23dot = Yp(:,2:evl.dgf:end);
e33dot = Yp(:,3:evl.dgf:end);
s22dot = Yp(:,4:evl.dgf:end);
s23dot = Yp(:,5:evl.dgf:end);
s33dot = Yp(:,6:evl.dgf:end);
sdotmax = sqrt((0.5*(s22dot-s33dot)).^2 + (s23dot).^2);
edotmax = sqrt((0.5*(e22dot-e33dot)).^2 + (e23dot).^2);
ej2dot = sqrt(e22dot.*e33dot - e23dot.^2);



%% plot strain at snapshots
figure(4),clf

tvec = [0.1,1,5,10,100];
% tvec = [1,10,30,100];

for i = 1:length(tvec)
    tind = find(abs(t-tvec(i)*3.15e7) == min(abs(t-tvec(i)*3.15e7)));

    toplot = (emax);
    % strains
    subplot(length(tvec),2,2*(i-1) + 1)
    plotshz(shz,(toplot(tind,:)))
    hold on
    axis tight equal,grid off, box on
    xlabel('X (m)'), ylabel('Z (m)')
    h=colorbar;h.FontSize = 12;h.Location='eastoutside';
    set(gca,'Color','none');
    h.Label.Interpreter='latex';h.Label.String='$\epsilon_{max} $';
    title(['Time = ' num2str(t(tind)./s2y) ' yrs'])
    %caxis([1e-2 1].*max(abs(get(gca,'CLim'))))
    %caxis([1 1].*max((get(gca,'CLim'))) + [-2 0])
    
%     caxis([-5.5 -3])
    set(gca,'ColorScale','log')
    
    toplot = (edotmax);
    %strain rates
    subplot(length(tvec),2,2*i)
    plotshz(shz,((toplot(tind,:))))
    axis tight equal,grid off, box on
    xlabel('X (m)'), ylabel('Z (m)')
    h=colorbar;h.FontSize = 12;h.Location='eastoutside';%set(gca,'ColorScale','log');
    h.Label.Interpreter='latex';h.Label.String='$\dot{\epsilon}_{max} $';
    %caxis([1e-2 1].*max(abs(get(gca,'CLim'))))
    %caxis([1 1].*max((get(gca,'CLim'))) + [-2 0])
    
%     caxis([-15 -12])
    title(['Time = ' num2str(t(tind)./s2y) ' yrs'])
    set(gca,'ColorScale','log')
    set(gca,'Color','none');
end

colormap(jet(20))
% cspec = (magmacolor(9));

set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% plot of strain and strain rate
x2pos = 0;
shzindex = find(abs(shz.x2-x2pos) == min(abs(shz.x2-x2pos)));
shzindex = [1:shz.N];
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
%% matrix solution for strains
Aeq = [1.*eye(shz.N) eye(shz.N).*0 1.*eye(shz.N)];
beq = zeros(shz.N,1);
tauinit = -[evl.LL3322;evl.LL3323;evl.LL3333]*initial_strain;

einv = lsqlin(...
    [evl.LL2222 evl.LL2322 evl.LL3322;...
     evl.LL2223 evl.LL2323 evl.LL3323;...
     evl.LL2233 evl.LL2333 evl.LL3333],...
    tauinit...
    ,[],[],...
    Aeq,beq);

einvmax = sqrt((0.5*(einv(1:shz.N)-einv(2*shz.N+1:end))).^2 + (einv(shz.N+1:2*shz.N)).^2);
% einvmax = sqrt(einv(1:shz.N).^2 + einv(2*shz.N+1:end).^2 + 2.*(einv(shz.N+1:2*shz.N)).^2);

% figure(5),clf
% subplot(311)
% plotshz(shz,tauinit(1:shz.N))
% axis tight equal
% set(gca,'CLim',max(abs(get(gca,'CLim'))).*[-1 1])
% cb = colorbar;
% 
% subplot(312)
% plotshz(shz,tauinit(shz.N+1:2*shz.N))
% axis tight equal
% set(gca,'CLim',max(abs(get(gca,'CLim'))).*[-1 1])
% cb = colorbar;
% 
% subplot(313)
% plotshz(shz,tauinit(2*shz.N+1:end))
% axis tight equal
% set(gca,'CLim',max(abs(get(gca,'CLim'))).*[-1 1])
% cb = colorbar;
% colormap bluewhitered

figure(6),clf
plotshz(shz,einvmax), shading flat
hold on
axis tight equal,grid off, box on
xlabel('X (m)'), ylabel('Z (m)')
h=colorbar;h.FontSize = 12;h.Location='eastoutside';
set(gca,'Color','none');
h.Label.String='\epsilon_{max}';
% caxis([-5.5 -3])
set(gca,'ColorScale','log')

colormap(jet)
%% surface displacement
uz = ([evl.G22_3 evl.G23_3 evl.G33_3]*[e22';e23';e33'])';
ux = ([evl.G22_2 evl.G23_2 evl.G33_2]*[e22 e23 e33]')';

uzinv = [evl.G22_3 evl.G23_3 evl.G33_3]*einv;
uxinv = [evl.G22_2 evl.G23_2 evl.G33_2]*einv;
% uz = ([evl.G22_3 evl.G23_3 evl.G33_3]*[0*initial_strain;0*initial_strain;1.*initial_strain])';
% ux = ([evl.G22_2 evl.G23_2 evl.G33_2]*[0*initial_strain;0*initial_strain;1.*initial_strain])';

figure(3),clf
subplot(211)
plot(ox./1e3,uz(end,:),'k-','LineWidth',2), hold on
plot(ox./1e3,uzinv,'r-','LineWidth',2)
plot(ox./1e3,0.*ox,'k--')
axis tight
plot(shz.x2(find(dummy~=0))./1e3.*[1 1],get(gca,'YLim'),'k--')
ylabel('u_z')

subplot(212)
plot(ox./1e3,ux(end,:),'k-','LineWidth',2), hold on
plot(ox./1e3,uxinv,'r-','LineWidth',2)
plot(ox./1e3,0.*ox,'k--')
axis tight
plot(shz.x2(find(dummy~=0))./1e3.*[1 1],get(gca,'YLim'),'k--')
ylabel('u_x')















