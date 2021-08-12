function Yp=ode_post_visco2d_v2(~,Y,shz,evl)
% ODE function: call as follows
%   [t,y]=ode45(yp,[0 200],y0,options);
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

% Maxwell strain rates d epsilon / dt 
ekkd = 0.*(Y(1:evl.dgf:end) + Y(3:evl.dgf:end))/2;

e22d=Y(1:evl.dgf:end)-ekkd; 
e23d=Y(2:evl.dgf:end); 
e33d=Y(3:evl.dgf:end)-ekkd; 

coupling = shz.coupling;

tau22d = (1.*evl.LL2222*e22d    + coupling.*evl.LL2322*e23d + coupling.*evl.LL3322*e33d);
tau23d = (coupling.*evl.LL2223*e22d + 1.*evl.LL2323*e23d        + coupling.*evl.LL3322*e33d);
tau33d = (coupling.*evl.LL2233*e22d + coupling.*evl.LL2333*e23d + 1.*evl.LL3333*e33d);
p = (tau22d + tau33d)/2;

e22dd = (tau22d-p)./shz.etaM;
e23dd = tau23d./shz.etaM;
e33dd = (tau33d-p)./shz.etaM;

% strainrate goes in the time-derivative
% Yp(1:evl.dgf:end)= sign(e22dd).*abs(e22dd);
% Yp(2:evl.dgf:end)= sign(e23dd).*abs(e23dd);
% Yp(3:evl.dgf:end)= sign(e33dd).*abs(e33dd);

Yp(1:evl.dgf:end)= (e22dd);
Yp(2:evl.dgf:end)= (e23dd);
Yp(3:evl.dgf:end)= (e33dd);
end 