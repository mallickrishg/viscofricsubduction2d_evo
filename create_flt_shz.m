function [rcv,shz,src] = create_flt_shz(earthModel,y2i,y3i,dip,Fwidth,Mf,Vwidth,Mv,Vwidthbot,Mvbot,Tpl)

addpath ~/Dropbox/scripts/unicycle/matlab/
import unicycle.*

% dip = 10;
fault_length = 500e4;


%% megathrust fault file
patchfname = 'megathrust2d.seg';

% width of patch segments
w = Fwidth/Mf;

fileID = fopen(patchfname,'w');
fprintf(fileID,'%s\n',...
    '# n  Vpl    x1     x2   x3   Length   Width   Strike  Dip  Rake    L0    W0    qL  qW');
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %d %d\n',...
    1, 1,-fault_length/2, y2i, y3i,fault_length, Fwidth, 0, dip, 90, fault_length, w, 1.0, 1.0);
fclose(fileID);

rcv = unicycle.geometry.receiver(patchfname,earthModel);


%% viscous shear zone (as a zero-width approximation)

% width of patch segments
w = Vwidth/Mv;

patchvname = 'viscousfault2d.seg';

fileID = fopen(patchvname,'w');
fprintf(fileID,'%s\n',...
    '# n  Vpl    x1     x2   x3   Length   Width   Strike  Dip  Rake    L0    W0    qL  qW');

% Three segments

% segment 1 - below the fault
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %d %d\n',...
    1, 1,-fault_length/2, y2i+Fwidth*cosd(dip), y3i+Fwidth*sind(dip),fault_length, Vwidth, 0, dip, 90, fault_length, w, 1.0, 1.0);

w = Vwidthbot/Mvbot;
% segment 2 - Tpl below the fault at dip
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %d %d\n',...
    1, -1,-fault_length/2, y2i, y3i+Tpl,fault_length, Vwidthbot, 0, dip, 90, fault_length, w, 1.0, 1.0);

% segment 3 - Tpl below the fault (flat)
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %d %d\n',...
    1, -1,-fault_length/2, y2i-Vwidthbot, y3i+Tpl,fault_length, Vwidthbot, 0, 0*dip, 90, fault_length, w, 1.0, 1.0);

fclose(fileID);

shz = unicycle.geometry.receiver(patchvname,earthModel);

%% ESPM source
% width of patch segments
w = 8000e3;

patchsrcname = 'espmfault2d.seg';

fileID = fopen(patchsrcname,'w');
fprintf(fileID,'%s\n',...
    '# n  Vpl    x1     x2   x3   Length   Width   Strike  Dip  Rake    L0    W0    qL  qW');

% Three segments

% segment 1 - below the fault
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %d %d\n',...
    1, 1,-fault_length/2, y2i+(Vwidth+Fwidth)*cosd(dip), y3i+(Vwidth+Fwidth)*sind(dip),fault_length, w, 0, dip, 90, fault_length, w, 1.0, 1.0);

% segment 2 - Tpl below the fault at dip
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %d %d\n',...
    1, -1,-fault_length/2, y2i+(Vwidthbot)*cosd(dip), y3i+Tpl+(Vwidthbot)*sind(dip),fault_length, w+Fwidth, 0, dip, 90, fault_length, w+Fwidth, 1.0, 1.0);

% segment 3 - Tpl below the fault (flat)
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %.18f %d %d\n',...
    1, -1,-fault_length/2, y2i-Vwidthbot-w, y3i+Tpl,fault_length, w, 0, 0*dip, 90, fault_length, w, 1.0, 1.0);

fclose(fileID);

src = unicycle.geometry.receiver(patchsrcname,earthModel);


end