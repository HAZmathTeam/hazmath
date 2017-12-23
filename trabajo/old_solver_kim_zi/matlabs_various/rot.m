%
% rot.m -- rotate the 3-D graph horizontally by 30 degrees or
%          steps of 30 degrees
%

function d=rot(n)

if nargin == 0
	angle = 20;
else
	angle = n * 20;
end

[az,el] = view;

if angle == 0  % want to go back to default
	az = -37.5;
else
	az = rem(az+angle, 360);
end

view(az, el)
