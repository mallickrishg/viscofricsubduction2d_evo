function Yp=ode_post_visco2d(~,Y,shz,evl)
% ODE function: call as follows
%   [t,y]=ode45(yp,[0 200],y0,options);
%
% where o is a class that contains all the relevent kernels and geometry
% information (evl class)
%        |      epsilon_ij   | (3 strains)
%        \      sigma_ij     / (3 stresses)
%
%
% for shear zones we use a maxwell rheology
%
%   d epsilon_ij          sigma_ij
%   ------------      =   --------
%      dt                  etaM
%
% etaM is the viscosity of the maxwell body
% 

% Initiate state derivative
Yp=zeros(size(Y));  

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%
%                 S H E A R   Z O N E S
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% isotropic stress
p=(Y(4:evl.dgf:end)+Y(6:evl.dgf:end))/2;

% deviatoric stress components
s22p=Y(4:evl.dgf:end)-p; 
s23p=Y(5:evl.dgf:end);   
s33p=Y(6:evl.dgf:end)-p; 

% Maxwell strain rates d epsilon / dt = sigma' / etaM
smax = sqrt((0.5*(s22p-s33p)).^2 + (s23p).^2);
e22d=(s22p).*(smax.^(shz.m-1))./shz.etaM; 
e23d=(s23p).*(smax.^(shz.m-1))./shz.etaM; 
e33d=(s33p).*(smax.^(shz.m-1))./shz.etaM; 

emaxd = sqrt((0.5*(e22d-e33d)).^2 + (e23d).^2);
% emaxd = sqrt(e22d.^2 + e33d.^2 + 2.*(e23d).^2);

if emaxd > 0
    g22d = (e22d);%./o.shz.epsilonPlate);
    g23d = (e23d);%./o.shz.epsilonPlate);
    g33d = (e33d);%./o.shz.epsilonPlate);
else
    g22d = 0;
    g23d = 0;
    g33d = 0;
end

% strainrate goes in the time-derivative
Yp(1:evl.dgf:end)=g22d;
Yp(2:evl.dgf:end)=g23d;
Yp(3:evl.dgf:end)=g33d;

coupling = shz.coupling;

taudot22 =     1.*evl.LL2222*g22d    + coupling.*evl.LL2322*g23d + coupling.*evl.LL3322*g33d;
taudot23 = coupling.*evl.LL2223*g22d + 1.*evl.LL2323*g23d        + coupling.*evl.LL3322*g33d;
taudot33 = coupling.*evl.LL2233*g22d + coupling.*evl.LL2333*g23d + 1.*evl.LL3333*g33d;

Yp(4:evl.dgf:end)=taudot22;
Yp(5:evl.dgf:end)=taudot23;
Yp(6:evl.dgf:end)=taudot33;


end 