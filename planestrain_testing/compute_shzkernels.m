function evl = compute_shzkernels(shz,ox)
% compute stress interaction kernel for shear zones in plane strain
% compute displacement kernels at the surface

G = 30e3;
nu = 0.25;

% evluate stresses at these points (shear zone centers)
X2 = shz.x2;
X3 = shz.x3;

% always use vertical shear zones
phi = 90;

for i = 1:shz.N
    q2 = shz.x2(i);
    q3 = -shz.x3(i) - shz.W(i)/2;
    
    % for phi=90, T - x2; W - x3
    T = shz.T(i);
    W = shz.W(i);
    
    % eps22
    epsv22p = 1;
    epsv23p = 0;
    epsv33p = 0;
    [s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
        X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);
    evl.LL2222(:,i) = s22;
    evl.LL2223(:,i) = s23;
    evl.LL2233(:,i) = s33;

    [u2,u3] = unicycle.greens.computeDisplacementPlaneStrainShearZone(...
        ox(:),-0.*ox(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);
    evl.G22_2(:,i) = u2;
    evl.G22_3(:,i) = -u3;
    
    % eps23
    epsv22p = 0;
    epsv23p = 1;
    epsv33p = 0;
    [s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
        X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);
    evl.LL2322(:,i) = s22;
    evl.LL2323(:,i) = s23;
    evl.LL2333(:,i) = s33;
    
    [u2,u3] = unicycle.greens.computeDisplacementPlaneStrainShearZone(...
        ox(:),-0.*ox(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);
    evl.G23_2(:,i) = u2;
    evl.G23_3(:,i) = -u3;
    
    % eps33
    epsv22p = 0;
    epsv23p = 0;
    epsv33p = 1;
    [s22,s23,s33]=unicycle.greens.computeStressPlaneStrainShearZone( ...
        X2(:),-X3(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);
    evl.LL3322(:,i) = s22;
    evl.LL3323(:,i) = s23;
    evl.LL3333(:,i) = s33;
    
    [u2,u3] = unicycle.greens.computeDisplacementPlaneStrainShearZone(...
        ox(:),-0.*ox(:),q2,q3,T,W,phi,epsv22p,epsv23p,epsv33p,G,nu);
    evl.G33_2(:,i) = u2;
    evl.G33_3(:,i) = -u3;
end



end