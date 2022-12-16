function plotpatch2d(rcv,toplot)

km2m = 1e3;

xtop = rcv.x(:,1);
xbot = rcv.x(:,1) + rcv.W.*cosd(rcv.dip);

ztop = -rcv.x(:,2);
zbot = -rcv.x(:,2) + rcv.W.*sind(rcv.dip);

if nargin==2
    for i = 1:rcv.N
        cline([xtop(i) xbot(i)]./km2m,-[ztop(i) zbot(i)]./km2m,toplot(i).*[1 1]), hold on
    end
else
    for i = 1:rcv.N
        plot([xtop(i) xbot(i)]./km2m,-[ztop(i) zbot(i)]./km2m,'k-','LineWidth',2), hold on
    end
end

end