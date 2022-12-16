function [rcv,shz,src] = create_flt_shz(earthModel,x0z0,dip,Fwidth,Mf,Vwidth,Mv,Vwidthbot,Mvbot,Tpl)

import fgeom.*

y2i = x0z0(1);
y3i = x0z0(2);

%% megathrust fault file
patchfname = 'megathrust2d.seg';

% width of patch segments
w = Fwidth/Mf;

fileID = fopen(patchfname,'w');
fprintf(fileID,'%s\n',...
    '# n  Vpl,  x   z    Width   Dip   W0   qW');
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %d\n',...
    1, 1, y2i, y3i, Fwidth, dip, w, 1.0);
fclose(fileID);

rcv = fgeom.receiver(patchfname,earthModel);


%% viscous shear zone (as a zero-width approximation)
y2hinge = -Tpl*tand(dip/2);

% width of patch segments
w = Vwidth/Mv;

patchvname = 'viscousfault2d.seg';

fileID = fopen(patchvname,'w');
fprintf(fileID,'%s\n',...
    '# n  Vpl x2   x3    Width   Dip W0   qW');

% Three segments

% segment 1 - below the fault
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %d\n',...
    1, 1, y2i+Fwidth*cosd(dip), y3i+Fwidth*sind(dip),Vwidth, dip, w, 1.0);

w = Vwidthbot/Mvbot;
% segment 2 - Tpl below the fault at dip
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %d\n',...
    1, -1, y2hinge, y3i+Tpl, Vwidthbot,dip, w, 1.0);

% segment 3 - Tpl below the fault (flat)
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %d\n',...
    1, -1,y2hinge-Vwidthbot, y3i+Tpl,Vwidthbot, 0, w, 1.0);

fclose(fileID);

shz = fgeom.receiver(patchvname,earthModel);

%% ESPM source
% width of patch segments
w = 100000e3;

patchsrcname = 'espmfault2d.seg';

fileID = fopen(patchsrcname,'w');
fprintf(fileID,'%s\n',...
    '# n  Vpl    x2   x3    Width   Dip  W0  qW');

% Three segments

% segment 1 - below the fault
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %d\n',...
    1, 1,y2i+(Vwidth+Fwidth)*cosd(dip), y3i+(Vwidth+Fwidth)*sind(dip),w, dip, w, 1.0);

% segment 2 - Tpl below the fault at dip
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %d\n',...
    1, -1,y2hinge+(Vwidthbot)*cosd(dip), y3i+Tpl+(Vwidthbot)*sind(dip),w+Fwidth, dip, w+Fwidth, 1.0);

% segment 3 - Tpl below the fault (flat)
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %d\n',...
    1, -1,y2hinge-Vwidthbot-w, y3i+Tpl, w, 0, w, 1.0);

fclose(fileID);

src = fgeom.receiver(patchsrcname,earthModel);


end