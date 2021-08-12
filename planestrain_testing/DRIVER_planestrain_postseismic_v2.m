clear
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*
addpath ~/Dropbox/scripts/utils/
s2y=60*60*24*365;

x2extent = 20e3;
Nx2 = 20;

x3extent = 10e3;
Nx3 = 10;
x3shift = 10e3;
ngps = 200;

ox = linspace(-x2extent,x2extent,ngps);

shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift);

%% load kernels for displacement and stress
tic
evl = compute_shzkernels(shz,ox);
toc
%% plot geometry
dummy = zeros(shz.N,1);
dummy(32) = 1e-3;

toplot = evl.LL3322*dummy;

figure(1),clf
subplot(3,1,[2 3])
plotshz(shz,toplot)
axis tight equal
colormap bluewhitered
cb=colorbar;cb.Location = 'northoutside';

subplot(3,1,1)
plot(ox,evl.G33_2*dummy,'b-','LineWidth',1)
hold on
plot(ox,evl.G33_3*dummy,'r-','LineWidth',1)
axis tight

drawnow
%% viscoelastic simulation

% uniform viscosity
shz.etaM = zeros(shz.N,1) + 1e18/1e6; % MPa s
shz.eps0 = 1e-14.*ones(shz.N,1);
shz.coupling = 1;
evl.dgf = 3;
%% total duration of simulation
tend = 100*s2y;

% Initialize State Vector
initial_edot = dummy;

Y0=zeros(shz.N*evl.dgf,1);

% Initial Condition
Y0(1:evl.dgf:evl.dgf*shz.N) = ((evl.LL2222*initial_edot)./shz.etaM)./shz.eps0;
Y0(2:evl.dgf:evl.dgf*shz.N) = ((evl.LL2223*initial_edot)./shz.etaM)./shz.eps0;
Y0(3:evl.dgf:evl.dgf*shz.N) = ((evl.LL3322*initial_edot)./shz.etaM)./shz.eps0;

%% initialize the function handle with
% set constitutive parameters
yp=@(t,y) ode_post_visco2d_v2(t,y,shz,evl);
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

% extract important quantities
e22dot = ((Y(:,1:evl.dgf:end)));
e23dot = ((Y(:,2:evl.dgf:end)));
e33dot = ((Y(:,3:evl.dgf:end)));
edotmax = sqrt((0.5*(e22dot-e33dot)).^2 + (e23dot).^2);

%% plot strain rates
% based on constant x2
% x2pos = 10e3;
% shzindex = find(abs(shz.x2-x2pos) == min(abs(shz.x2-x2pos)));

%based on constant x3
x3pos = 19e3;
shzindex = find(abs(shz.x3-x3pos) == min(abs(shz.x3-x3pos)));
cspec = jet(length(shzindex));

toplot = 1.*e22dot + 0.*e23dot + 1.*e33dot;
% toplot = edotmax;

figure(2),clf
for i = 1:length(shzindex)
    plot(t./s2y,toplot(:,shzindex(i)),'-','LineWidth',2,'Color',cspec(i,:))
    axis tight, grid on, hold on
    set(gca,'YScale','linear','XScale','log','FontSize',15)
    xlim(10.^[-3 2])
    xlabel('Time (yrs)')
    ylabel('$\dot{\epsilon}$','Interpreter','latex','FontSize',20)
end
