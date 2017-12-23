load fort.199
x = fort(1,:);
y = fort(2,:);
xmin = -0.0001
xmax = 1.0001
ymin = -0.0001
ymax = 1.0001
figure(1)
axis([xmin xmax ymin ymax])
hold on
a1 = plotgraf(x',y',201,0);
a = plotgraf(x',y',200,1);
title(date)
