% script to make GPS transects
clear

Ng = 2000;



xg = 0;%linspace(-200e3,200e3,Ng);
yg = linspace(-500e3,1000e3,Ng);
[Xg,Yg] = meshgrid(xg,yg);
x1 = Yg(:);
x2 = Xg(:);

N = length(x1);

z = 0;
ntransect = length(x1);

T = [];

for i = 1:ntransect
%     T = [T; [((i-1)*N+1):(i*N)]' [((i-1)*N+1):(i*N)]' x1(i).*ones(size(x2)) x2 z.*ones(size(x2))];
    
    T = [T; i i x1(i) x2(i) z];
end

writetable(table(T),'./gps/gpsnetwork.dat','Delimiter','\t','WriteVariableNames',0)