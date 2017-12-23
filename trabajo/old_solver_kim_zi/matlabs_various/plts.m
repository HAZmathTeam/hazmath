x=0:.1:1;
y=x;
[X,Y]=meshgrid(x,y);
%Z=4*(X.*X-Y) - Y.*Y;
%W=(2.*X-Y).*exp(8*(X-Y));
%W=10*(2.*X-Y).*sin(8*(X-Y));
Z=-exp(Y) + 2*X;
W=3*Y.*sin(pi*X);
mesh(X,Y,Z)
pause
surf(X,Y,Z)
title('Beta1 = 4(x^2 - y) - y^2')
pause
mesh(X,Y,W)
pause
surf(X,Y,W)
%title('Beta2 = (2x - y) e^{8(x-y)}')
title('Beta2 = 10(2x - y) sin 8(x-y)')
pause
mesh(X,Y,Z-W)
pause
surf(X,Y,Z-W)
pause
close all
