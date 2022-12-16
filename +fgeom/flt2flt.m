function [flt]=flt2flt(xo,W,dip,wo,alphaw)
% [flt]=flt2flt(xo,L,W,strike,dip,lo,wo,alphal,alphaw)
%
% function flt2flt subsamples a fault patch into smaller segments
% of length and width starting from lo and wo, respectively
% increasing geometrically with down-dip distance with increment 
% alphal and alphaw (alphal>1 for increase).
%
% input:
%   xo     origin position vector [x1 (horizontal);x3 (down)]
%   W      total width
%   dip    dip angle in degrees
%   wo     approximative initial width of output segments
%   alphaw geometric factor for width increase
%
% output:
%   flt    list of output segments in the format
%          x1,x3,width,dip
%

% create wi
Wc=W;
k=0;
w=0;
while Wc>0
    Wt=wo*alphaw^k;
    if Wt > Wc/2
        Wt = Wc;
    end
    wn=min([Wt,Wc]);
    w=[w; wn];
    k=k+1;
    Wc=Wc-wn;
end
Nw=k;

% strike and dip direction normal vectors
Dv=[cosd(dip); sind(dip)];

flt=[];

% loop in dip direction
k=0;
for j=1:Nw
        
    x=xo + sum(w(1:j))*Dv;
    k=k+1;
    
    flt = [flt; [x',w(j+1),dip]];

end

end




        