function Yp=maxwell(~,Y,properties)

G = properties.G;
eta = properties.eta;
m = properties.m;

edot = Y(1).^m/eta;

Yp = zeros(size(Y));
Yp(1) = -G*edot + G*properties.vpl;
Yp(2) = edot;

end