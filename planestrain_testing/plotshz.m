function plotshz(shz,val)

dx = shz.T(:)/2;
dy = shz.W(:)/2;

x = [shz.x2(:)-dx,shz.x2(:)-dx,shz.x2(:)+dx,shz.x2(:)+dx,shz.x2(:)-dx]';
y = [shz.x3(:)-dy,shz.x3(:)+dy,shz.x3(:)+dy,shz.x3(:)-dy,shz.x3(:)-dy]';

patch(x,y,val)

end