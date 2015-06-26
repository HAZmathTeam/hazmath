%%%  Plot solutions
clear all;
close all;

A = load('mat.dat');
b = load('rhs.dat');

p = load('coords.dat');

A = spconvert(A);

u = A\b;

%% read mesh from folder file


nnodes = size(p,1);
mydim = size(p,2);
myh = 1/(sqrt(nnodes)-1);

dim = [int2str(mydim),'D'];


    nov = 1;
    ndof = nnodes;

solblock = ndof;

start = 1;
vend = nnodes;
        
xmin = min(p(:,1));
ymin = min(p(:,2));
xmax = max(p(:,1));
ymax = max(p(:,2));

h = (xmax-xmin)/(sqrt(nnodes)-1);
spacing = round(1/(2*h));
spacing = 2;

    velocity = [];
    vtrue = [];
    for (i=1:nov)
        vstart = 1;
        vend = nnodes;
        velocity = [velocity u(vstart:vend)];
        
    end
    
    NumNod = nnodes; kdgof = mydim;
    
    p1nodes = ((sqrt(nnodes)-1)/2)^2;
    
    
    
    for(j=1:nnodes)
        
        n = round((p(j,1)-xmin)/h + 1);
        m = round((p(j,2)-ymin)/h + 1);
        x(j) = p(j,1);
        y(j) = p(j,2);
        X(m,n) = p(j,1);
        Y(m,n) = p(j,2);
        U(m,n) = velocity(j,1);
        Utrue(m,n) = sin(pi*y(j))*(sin(pi*x(j)));
  
    
    end
%         
        subplot(2,2,1)
        contour(X,Y,U)
        colorbar
        %axis([-1 1 -1 1 -1 1]);
        subplot(2,2,2)
        contour(X,Y,Utrue)
        colorbar
        % axis([-1 1 -1 1 -1 1]);
        title('u0')
