function shz = create_shzmesh(x2extent,Nx2,x3extent,Nx3,x3shift)
% specify x2extent and x3extent (& x3shift) as positive values

N = Nx2*Nx3;
shz.N = N;

x2 = linspace(-x2extent,x2extent,Nx2);
T = abs(x2(1)-x2(2));

x3 = linspace(-(x3extent+x3shift),-x3shift,Nx3);
W = abs(x3(1)-x3(2));

[X2,X3] = meshgrid(x2,x3);
shz.x2 = X2(:);
shz.x3 = X3(:);
shz.T = T.*ones(N,1);
shz.W = W.*ones(N,1);

end