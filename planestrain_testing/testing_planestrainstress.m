clear
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*
addpath ~/Dropbox/scripts/utils/


nx2 = 70;
nx3 = 70;
x2 = linspace(-25e3,25e3,nx2);
x3 = linspace(-150e3,-100e3,nx3);

[X2,X3] = meshgrid(x2,x3);
% X2 = shz.xc(:,2);
% X3 = shz.xc(:,3);

% index = 209;

% q2 = shz.x(index,2);
% q3 = -shz.x(index,3);
% T = shz.L(index);
% W = shz.W(index);
q2 = 0;
q3 = 120e3;
T = 15e3;
W = 5e3;

phi = 90;
epsv22p = 1;
epsv23p = 0;
epsv33p = 1;
G = 30e3;
nu = 0.25;

[s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
    X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);

[u2,u3] = unicycle.greens.computeDisplacementPlaneStrainShearZone(...
    X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);

figure(13),clf

if epsv22p+epsv33p ~= 0
    toplot = (0.5.*s22 + 0.5.*s33)./G/(epsv22p+epsv23p+epsv33p);
    %toplot = (s23)./G/(epsv22p+epsv23p+epsv33p);
else
    toplot = (0.5.*s22 + 0.5.*s33)./G/abs(epsv22p+epsv23p);
end

% scatter(X2(:)./1e3,X3(:)./1e3,50,toplot,'filled')
imagesc(x2./1e3,x3./1e3,reshape(toplot,nx3,nx2)),hold on
% contour(x2./1e3,x3./1e3,reshape(toplot,nx3,nx2),50,'LineWidth',2), hold on
alpha 0.5
quiver(X2(:)./1e3,X3(:)./1e3,u2(:),-u3(:),1,'k-','LineWidth',1)
% quiver(X2(:),X3(:),u3(:),u2(:),'k-','LineWidth',1)

plot([q2-T/2,q2+T/2,q2+T/2,q2-T/2,q2-T/2]./1e3,-[q3,q3,q3+W,q3+W,q3]./1e3,'k-','LineWidth',2)

axis tight equal, shading flat
cb=colorbar; cb.Location = 'northoutside';
cb.Label.String = 'Normalized Pressure';
box on, grid on
% caxis([-1 1].*0.005)
colormap bluewhitered(200)
set(gca,'FontSize',15,'YDir','normal')
% print('Pressure_volumesource','-djpeg','-r100')
%%  figure(14),clf
% plot(x2,u2,'r-','LineWidth',2)
% hold on
% plot(x2,u3,'b-','LineWidth',2)
% plot(x2,x2.*0,'k--')
% axis tight
% ylim([-1 1].*max(abs(get(gca,'YLim'))))
% plot([0 0],get(gca,'YLim'),'k--')
% set(gca,'FontSize',15)
% legend('u2','u3')
figure(14),clf
quiver(X2(:),X3(:),u2(:),-u3(:),'k-')
axis tight equal
%% individual components
q2 = 0;
q3 = 120e3;
T = 3e3;
W = 3e3;

% for eps22 only
epsv22p = 1;
epsv23p = 0;
epsv33p = 0;
G = 30e3;
nu = 0.25;

[s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
    X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);

[u2,u3] = unicycle.greens.computeDisplacementPlaneStrainShearZone(...
    X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);

figure(10),clf
for i = 1:3
    subplot(1,3,i)
    
    if i == 1
        toplot = (0.5.*s22 + 0.5.*s33)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized Pressure';
    elseif i == 2
        toplot = (s22)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized \sigma_{22}';
    else
        toplot = (s33)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized \sigma_{33}';
    end
    imagesc(x2./1e3,x3./1e3,reshape(toplot,nx3,nx2)),hold on
    alpha 0.5
    quiver(X2(:)./1e3,X3(:)./1e3,u2(:),-u3(:),1,'k-','LineWidth',1)    
    plot([q2-T/2,q2+T/2,q2+T/2,q2-T/2,q2-T/2]./1e3,-[q3,q3,q3+W,q3+W,q3]./1e3,'k-','LineWidth',2)
    cb=colorbar; cb.Location = 'northoutside';
    cb.Label.String = cbstring;
    axis tight equal, shading flat
    box on, grid on
    caxis([-1 1].*0.5)
    xlabel('X_2 (km)')
    ylabel('X_3 (km)')
    colormap bluewhitered(40)
    set(gca,'FontSize',12,'YDir','normal')
end
% print('Stresscomponents_eps22','-djpeg','-r100')

% for eps33 only
epsv22p = 0;
epsv33p = 1;
G = 30e3;
nu = 0.25;

[s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
    X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);

[u2,u3] = unicycle.greens.computeDisplacementPlaneStrainShearZone(...
    X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);

figure(11),clf
for i = 1:3
    subplot(1,3,i)
    
    if i == 1
        toplot = (0.5.*s22 + 0.5.*s33)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized Pressure';
    elseif i == 2
        toplot = (s22)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized \sigma_{22}';
    else
        toplot = (s33)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized \sigma_{33}';
    end
    imagesc(x2./1e3,x3./1e3,reshape(toplot,nx3,nx2)),hold on
    alpha 0.5
    quiver(X2(:)./1e3,X3(:)./1e3,u2(:),-u3(:),1,'k-','LineWidth',1)    
    plot([q2-T/2,q2+T/2,q2+T/2,q2-T/2,q2-T/2]./1e3,-[q3,q3,q3+W,q3+W,q3]./1e3,'k-','LineWidth',2)
    cb=colorbar; cb.Location = 'northoutside';
    cb.Label.String = cbstring;
    axis tight equal, shading flat
    box on, grid on
    caxis([-1 1].*0.5)
    xlabel('X_2 (km)')
    ylabel('X_3 (km)')
    colormap bluewhitered(40)
    set(gca,'FontSize',12,'YDir','normal')
end
% print('Stresscomponents_eps33','-djpeg','-r100')

%% combined volume source
% for eps22 only
epsv22p = 1;
epsv23p = 0;
epsv33p = 1;
G = 30e3;
nu = 0.25;

[s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
    X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);

[u2,u3] = unicycle.greens.computeDisplacementPlaneStrainShearZone(...
    X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);

figure(12),clf
for i = 1:3
    subplot(1,3,i)
    
    if i == 1
        toplot = (0.5.*s22 + 0.5.*s33)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized Pressure';
    elseif i == 2
        toplot = (s22)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized \sigma_{22}';
    else
        toplot = (s33)./G/(epsv22p+epsv23p+epsv33p);
        cbstring = 'Normalized \sigma_{33}';
    end
    imagesc(x2./1e3,x3./1e3,reshape(toplot,nx3,nx2)),hold on
    alpha 0.5
    quiver(X2(:)./1e3,X3(:)./1e3,u2(:),-u3(:),1,'k-','LineWidth',1)    
    plot([q2-T/2,q2+T/2,q2+T/2,q2-T/2,q2-T/2]./1e3,-[q3,q3,q3+W,q3+W,q3]./1e3,'k-','LineWidth',2)
    cb=colorbar; cb.Location = 'northoutside';
    cb.Label.String = cbstring;
    axis tight equal, shading flat
    box on, grid on
    caxis([-1 1].*0.5)
    xlabel('X_2 (km)')
    ylabel('X_3 (km)')
    colormap bluewhitered(40)
    set(gca,'FontSize',12,'YDir','normal')
end
% print('Stresscomponents_epskk','-djpeg','-r100')


