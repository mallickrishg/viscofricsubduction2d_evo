clear
addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*
addpath ~/Dropbox/scripts/utils/

nx2 = 40;
nx3 = 40;
x2 = linspace(-5e3,5e3,nx2);
x3 = linspace(-10e3,-0e3,nx3);

[X2,X3] = meshgrid(x2,x3);

q2 = 0;
q3 = 3.5e3;
T = 1e3;
W = 1e3;

% shear zone deformation
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


% compute tractions on a fault
dip = 0;% i have checked this only for 0<dip<90

dhat = [-cosd(dip) sind(dip)];% dipvector
nhat = [sind(dip) cosd(dip)];% normalvector

Tvec = [s22.*nhat(1)+s23.*nhat(2) s23.*nhat(1)+s33.*nhat(2)]./G;
sn = Tvec(:,1).*nhat(1) + Tvec(:,2).*nhat(2);% fault normal traction
sd = Tvec(:,1).*dhat(1) + Tvec(:,2).*dhat(2);% fault parallel traction

toplot = sn;

figure(13),clf
imagesc(x2./1e3,x3./1e3,reshape(toplot,nx3,nx2)),hold on
% contour(x2./1e3,x3./1e3,reshape(toplot,nx3,nx2),50,'LineWidth',2), hold on
alpha 0.5
quiver(X2(:)./1e3,X3(:)./1e3,u2(:),-u3(:),1,'k-','LineWidth',1)

plot([q2-T/2,q2+T/2,q2+T/2,q2-T/2,q2-T/2]./1e3,-[q3,q3,q3+W,q3+W,q3]./1e3,'k-','LineWidth',2)
axis tight equal, shading flat
cb=colorbar; cb.Location = 'northoutside';
cb.Label.String = '\sigma_n (positive is unclamping)';
box on, grid on
caxis([-1 1].*.25)
colormap bluewhitered(40)
set(gca,'FontSize',15,'YDir','normal')
title(['\delta (fault dip) = ' num2str(dip) 'º'])

plot([-1 1].*2.5,[1 1].*-5.1,'k--','LineWidth',2)
print('Faulttractions_volumesource_0_deeper','-djpeg','-r100')
% print('Faulttractions_volumesource','-djpeg','-r100')
