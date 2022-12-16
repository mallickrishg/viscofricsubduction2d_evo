function [Ktau,Ksigma] = computetractionkernels(src,rcv)
% traction kernel computation
% Rishav Mallick, 2020

Ktau = zeros(rcv.N,src.N);
Ksigma = Ktau;
x = rcv.xc(:,1);
z = rcv.xc(:,2);

for i = 1:src.N
    m = [src.x(i,2) src.x(i,1) src.W(i) src.dip(i) 1];
    [Sxx,Sxz,Szz] = fgeom.EdgeStress(m,x,z,src.earthModel.nu,src.earthModel.G);
    
    t=[Sxx.*rcv.nv(:,1)+Sxz.*rcv.nv(:,2), ...
        Sxz.*rcv.nv(:,1)+Szz.*rcv.nv(:,2)];
    
    Ktau(:,i) = rcv.dv(:,1).*t(:,1) + rcv.dv(:,2).*t(:,2);
    Ksigma(:,i) = rcv.nv(:,1).*t(:,1) + rcv.nv(:,2).*t(:,2);
end



end