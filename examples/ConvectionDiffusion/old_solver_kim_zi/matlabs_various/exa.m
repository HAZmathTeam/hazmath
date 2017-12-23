clear all
load fort.199
x = fort(1,:);
y = fort(2,:);
h = 1/(sqrt(length(x))-1);
exac = zeros(length(x),1);
for i = 1 : length(x)
exac(i) = (x(i)+x(i)*x(i));       
exac(i) = exac(i)*(y(i)+y(i)*y(i));
end
%exac = 100*exac;
load fort.100
sol = fort';
err = exac-sol;
[norm(err,2)*h,norm(err,inf)]
