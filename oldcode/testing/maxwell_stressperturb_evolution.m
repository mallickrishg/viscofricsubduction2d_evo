clear

properties.G = 0.3;
properties.eta = 4;
properties.vpl = 2;
properties.m = 1;

Y0 = [10,0]';% stress,strain
tend = 1e2;

yp=@(t,y) maxwell(t,y,properties);
tic

options=odeset('Refine',1,'RelTol',1e-7,'AbsTol',1e-7,'InitialStep',1e-6,'MaxStep',1e3);
tvec = [0,tend];
[t,Y]=ode45(yp,tvec,Y0,options);

figure(1),clf
subplot(211)
plot(t,Y(:,1),'LineWidth',2)
grid on, axis tight
ylabel('Stress')
set(gca,'YScale','linear','XScale','linear')
% ylim([0 10])

subplot(212)
plot(t,Y(:,1)/properties.eta,'k-','LineWidth',2)
grid on, axis tight
ylabel('Strain Rate')
set(gca,'YScale','log','XScale','linear')
% ylim([0 10])

yyaxis right
plot(t,Y(:,2),'-','LineWidth',2)
grid on, axis tight
ylabel('Strain')
xlabel('t')
% set(gca,'YScale','log')

% figure(2),clf
% plot(Y(:,1)/properties.eta,Y(:,1),'LineWidth',2)
% axis tight
% xlabel('Strain Rate')
% ylabel('Stress')