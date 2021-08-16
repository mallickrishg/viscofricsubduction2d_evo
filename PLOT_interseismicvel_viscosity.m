% compute interseismic slip rates and the surface velocities due to
% different viscosity structures in a simplified 2-d subduction zone 
% Rishav Mallick, EOS, 2021
clear

faultparams = [];
faultparams.dip = 10;
faultparams.Tplate = 30e3;
faultparams.Vpl = 1e-9;

shzparams = [];
shzparams.eta_arc = 1e15;
shzparams.eta_oc = 1e13;
etavec = [5e12,1e13,3e13,5e13,7e13,1e14,3e14,5e14,1e15,5e15];
% etavec = [1e14,1e15];

Trecur = 200;%in years
ncycles = 200;

Nruns = length(etavec);

%% run models and make plots
figure(1),clf
figure(2),clf
cspec = jet(Nruns);
p=[];
lgdvec={};
for i = 1:Nruns
    shzparamsval = shzparams;
    shzparamsval.eta_oc = etavec(i);
    [rcv,shz,src,results] = RUN_simpleSZ_imposedcycles(faultparams,shzparamsval,Trecur,ncycles);
    t = results.t;
    V = results.V;
    e12d = results.e12d;
    Vpl = mean(rcv.Vpl);
    %% displacement kernels
    ng = 400;
    ox = linspace(-500e3,700e3,ng)';
    
    GF_d = compute_displacementkernels([ox 0.*ox 0.*ox],rcv,shz,src);
    
    % surface velocity
    gps = [];
    gps.rcv.vh = (GF_d.rcv.Gh*V')';
    gps.rcv.vz = (GF_d.rcv.Gz*V')';
    gps.shz.vh = (GF_d.shz.Gh*e12d')';
    gps.shz.vz = (GF_d.shz.Gz*e12d')';
    gps.src.vh = Vpl.*(GF_d.src.Gh*src.Vpl)';
    gps.src.vz = Vpl.*(GF_d.src.Gz*src.Vpl)';
    
    %% plot surface velocities
    figure(1)
    subplot(211)
    plot(ox./1e3,(gps.src.vh + gps.shz.vh(end,:) + gps.rcv.vh(end,:))./Vpl,'-','LineWidth',2,'Color',cspec(i,:)); 
    hold on
    
    plot(rcv.xc(rcv.a<rcv.b,1)./1e3,0.*rcv.Vpl(rcv.a<rcv.b),'k-','LineWidth',3)
    axis tight, grid on
    ylim([-0.1 1.1])
    plot(max(rcv.xc(:,1)./1e3).*[1 1],[-1 2],'k--','Linewidth',2)
    plot(0.*[1 1],[-1 2],'k-','Linewidth',2)
    xlabel('x_2 (km)')
    ylabel('v_h/v_{pl}')
    set(gca,'Fontsize',20,'LineWidth',2)
    
    subplot(212)
    p(i) = plot(ox./1e3,(gps.src.vz + gps.shz.vz(end,:) + gps.rcv.vz(end,:))./Vpl,'-','LineWidth',2,'Color',cspec(i,:)); 
    dummy = log10(etavec(i));
    lgdvec{i} = ['\eta = ' num2str(round(10^(dummy-floor(dummy)))) 'x10^{' num2str(floor(dummy)) '} Pa-s'];
    hold on
    axis tight, grid on
    ylim([-1 1].*0.3)
    plot(max(rcv.xc(:,1)./1e3).*[1 1],[-1 1],'k--','Linewidth',2)
    plot(0.*[1 1],[-1 1],'k-','Linewidth',2)
    xlabel('x_2 (km)')
    ylabel('v_z/v_{pl}')
    set(gca,'Fontsize',20,'Linewidth',2)
    
    %% plot sliprate on fault and bottom shear zone
    figure(2)
    subplot(121)
    index = shz.Vpl>0;
    plot([rcv.xc(:,1);shz.xc(index,1)]./1e3,[V(end,:),e12d(end,index)]./Vpl,'r-','LineWidth',2,'Color',cspec(i,:)), hold on
    axis tight, grid on
    ylim([0 1])
    xlabel('x_2 (km)')
    ylabel('v/v_{pl}')
    set(gca,'Fontsize',20,'Linewidth',2)
    
    subplot(122)
    index = shz.Vpl<0;
    [xsort,isort] = sort(shz.xc(index,1)./1e3);
    toplot = e12d(end,index)./Vpl;
    plot(xsort,toplot(isort),'-','LineWidth',2,'Color',cspec(i,:)), hold on
    axis tight, grid on
    ylim([-1.1 0])
    xlabel('x_2 (km)')
    ylabel('v/v_{pl}')
    set(gca,'Fontsize',20,'Linewidth',2)
end
figure(1),subplot(212)
legend(p,lgdvec)
set(legend,'Box','off','Location','northwest','Fontsize',20)
print('Figures/interseismicvel_viscosity','-djpeg','-r300')

figure(2)
legend(lgdvec)
print('Figures/interseismicsliprates_viscosity','-djpeg','-r300')
