%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
%                   Function of Time                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%% PLOT megathrust behaviour
down = 50; %temporal downslampling
Plotshz =(false(1,shz.N));
% shzindex = [1,2,3,41,42,43,82,83,123,124,163,164,165,204,205,206,245,246,257,287,288,289,328,329,330,...
%     371,372,373];
% Plotshz(shzindex) = 1;
Plotshz(shz.xc(:,2)<-1.9*shz.xc(:,3) - 100e3 & shz.xc(:,2)>-1.9*shz.xc(:,3) - 120e3) = 1;
shzindex = find(Plotshz==1);

% plot Time and Time-step evolution
figure(1);clf;set(gcf,'name','Time Evolution')
subplot(3,2,1)

pcolor(t(1:down:end)./3.15e7,dipdist/1e3,log10(V(1:down:end,:)')), shading flat, hold on

box on
set(gca,'YDir','reverse');h=colorbar('Location','NorthOutside');
title(h,'log_{10} V'),ylabel('\zeta_d (km)');
colormap(jet(40))
caxis([-13 0])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,[3 5])
pcolor(t(1:down:end)./3.15e7,shz.xc(shzindex,3)./1e3,log10(edotmax(1:down:end,Plotshz)')), shading flat, hold on
box on
set(gca,'YDir','normal');h=colorbar('Location','NorthOutside');
title(h,'log_{10} d\epsilon'),ylabel('Depth (km)');
colormap(jet(40))
caxis([-16 -13])
set(gca,'FontSize',15)
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                Function of Time Steps                %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
subplot(3,2,2);cla;
imagesc(1:down:length(t),dipdist/1e3,log10(V(1:down:end,:)')), shading flat, hold on

set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');h.Label.String='log_{10}V (m/s)';
ylabel('\zeta_d (km)');
caxis([-13 0])
set(gca,'XTickLabel','','FontSize',15)

subplot(3,2,[4 6])
pcolor((1:down:length(t))./3.15e7,shz.xc(shzindex,3)./1e3,log10(edotmax(1:down:end,Plotshz)')), shading flat, hold on
box on
set(gca,'YDir','normal');h=colorbar('Location','NorthOutside');
title(h,'log_{10} d\epsilon'),ylabel('Depth (km)');
caxis([-16 -13])
set(gca,'FontSize',15)

%% PLOT shear zone results

figure(2),clf
% plot movie
for i = 2.58e4
    subplot(2,1,1)
    plot(dipdist./1e3,V(i,:)./rcv.Vpl(1),'k','LineWidth',1)
    axis tight, grid on
    %set(gca,'YScale','log')
    ylim([1e-11 1])
    title(['Time = ' num2str(t(i)./3.15e7) ' yrs'])
    
    subplot(2,1,2)
    shz.plotShearZoneEdges(log10(edotmax(i,:))), hold on
    shz.plotShearZoneEdges
    view(89,0)
    axis tight equal, grid on
    xlabel('East (m)'),ylabel('North (m)'),zlabel('Z (m)')
    colormap jet
    colorbar
    %caxis([-15 -12])
    
end

%% PLOT shear zone results
if false
    down = 200;
    figure(2),clf
    % plot movie
    for i = 2:down:length(t)
        subplot(2,1,1)
        plot(dipdist./1e3,V(i,:),'k','LineWidth',1)
        axis tight, grid on
        set(gca,'YScale','log')
        ylim([1e-12 1])
        
        subplot(2,1,2)
        shz.plotShearZoneEdges(log10(edotmax(i,:))), hold on
        shz.plotShearZoneEdges
        view(89,0)
        axis tight equal, grid on
        xlabel('East (m)'),ylabel('North (m)'),zlabel('Z (m)')
        colormap jet
        colorbar
        caxis([-20 -10])
        
        pause(0.01)
    end
end
%%
Plotshz =(false(1,shz.N));
shzindex = [1,2,3,41,42,43,82,83,123,124,163,164,165,204,205,206,245,246,257,287,288,289,328,329,330,...
    371,372,373];
% Plotshz(shzindex) = 1;
Plotshz(shz.xc(:,2)<-1.9*shz.xc(:,3) - 100e3 & shz.xc(:,2)>-1.9*shz.xc(:,3) - 120e3) = 1;

figure(1),clf
subplot(2,1,1)
shz.plotShearZoneEdges(log10(shz.etaM)), hold on
shz.plotShearZoneEdges
view(90,0)
axis tight equal, grid on
xlabel('East (m)'),ylabel('North (m)'),zlabel('Z (m)')

subplot(2,1,2)
shz.plotShearZoneEdges(double(Plotshz)), hold on
shz.plotShearZoneEdges
view(90,0)
axis tight equal, grid on
xlabel('East (m)'),ylabel('North (m)'),zlabel('Z (m)')
