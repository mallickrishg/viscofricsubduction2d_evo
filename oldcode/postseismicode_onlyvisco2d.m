function Yp=postseismicode_onlyvisco2d(~,Y,o)
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
p=1.*(Y(4:o.shz.dgf:end)+Y(6:o.shz.dgf:end))/2;

% deviatoric stress components
s22p=Y(4:o.shz.dgf:end)-p; 
s23p=Y(5:o.shz.dgf:end);   
s33p=Y(6:o.shz.dgf:end)-p; 

% Maxwell strain rates d epsilon / dt = sigma' / etaM
m = o.shz.m;
e22d=(s22p.*(s22p.^(m-1)))./o.shz.etaM; 
e23d=(s23p.*(s23p.^(m-1)))./o.shz.etaM; 
e33d=(s33p.*(s33p.^(m-1)))./o.shz.etaM; 

% zerostrain = 1e-18;
e22d(o.shz.sigma0==1) = 0;
e23d(o.shz.sigma0==1) = 0;
e33d(o.shz.sigma0==1) = 0;

g22d = (e22d);%./o.shz.epsilonPlate);
g23d = (e23d);%./o.shz.epsilonPlate);
g33d = (e33d);%./o.shz.epsilonPlate);

% strainrate goes in the time-derivative
Yp(1:o.shz.dgf:end)=g22d;
Yp(2:o.shz.dgf:end)=g23d;
Yp(3:o.shz.dgf:end)=g33d;

coupling = 1;

taudot22 =     1.*o.LL{1,1}*g22d    + coupling.*o.LL{3,1}*g23d + coupling.*o.LL{6,1}*g33d;
taudot23 = coupling.*o.LL{1,3}*g22d + 1.*o.LL{3,3}*g23d        + coupling.*o.LL{6,3}*g33d;
taudot33 = coupling.*o.LL{1,6}*g22d + coupling.*o.LL{3,6}*g23d + 1.*o.LL{6,6}*g33d;

Yp(4:o.shz.dgf:end)=taudot22;
Yp(5:o.shz.dgf:end)=taudot23;
Yp(6:o.shz.dgf:end)=taudot33;


end 