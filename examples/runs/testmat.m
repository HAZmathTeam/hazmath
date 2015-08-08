%%%  Plot solutions
clear all;
close all;

A = load('mat.dat');
b = load('rhs.dat');

p = load('coords.dat');

A = spconvert(A);

um = A\b;

u = load('sol.dat');

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
    vm = [];
    for (i=1:nov)
        vstart = 1;
        vend = nnodes;
        velocity = [velocity u(vstart:vend)];
        vm = [vm um(vstart:vend)];
        
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
        Um(m,n) = vm(j,1);
        Utrue(m,n) = sin(pi*y(j))*(sin(pi*x(j)));
  
    
    end
%         
        figure(1)
        subplot(2,2,1)
        contour(X,Y,U)
        colorbar
        title('HAZMAT')
        %axis([-1 1 -1 1 -1 1]);
        subplot(2,2,3)
        contour(X,Y,Utrue-U)
        colorbar
        % axis([-1 1 -1 1 -1 1]);
        title('Error with HAZMAT')
        subplot(2,2,2)
        contour(X,Y,Um)
        colorbar
        % axis([-1 1 -1 1 -1 1]);
        title('Matlab u')
        subplot(2,2,4)
        contour(X,Y,Utrue-Um)
        colorbar
        % axis([-1 1 -1 1 -1 1]);
        title('Error with Matlab')
        
        figure(2)
        subplot(2,2,1)
        mesh(X,Y,U)
        title('HAZMAT')
        %axis([-1 1 -1 1 -1 1]);
        subplot(2,2,3)
        mesh(X,Y,Utrue-U)
        % axis([-1 1 -1 1 -1 1]);
        title('Error with HAZMAT')
        subplot(2,2,2)
        mesh(X,Y,Um)
        % axis([-1 1 -1 1 -1 1]);
        title('Matlab u')
        subplot(2,2,4)
        mesh(X,Y,Utrue-Um)
        % axis([-1 1 -1 1 -1 1]);
        title('Error with Matlab')
