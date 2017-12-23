%n = 65
x = [0:1/(n-1):1];
y = x;
load fort.100
z = reshape(fort', [n n]);
figure(1)
mesh(x,y,z)
figure(2)
contour(x,y,z)

