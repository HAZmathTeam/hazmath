xy=[...
    1, 1.5  ,3;
    2, 2.5  ,3;
    3, 2    ,1.5;
    4, 0.75 ,0.5;
    5, 2    ,3.75;
    6, 1    ,3.75;
    7, 2    ,5;
    8, 0.5, 2; 
     9, 1.5,2.1;
     10,0.5,3;
     11, 3.1,1.9...
    ];

seg=[...
        4,3;
	8,9;
	 10,1;
	 3,11;
	 9,2;
	 9,1;
	 1,5;
	 6,7;
	 6,1;
	 7,5;
	 10,8;
	 2,11;
	 9,3;
	 1,9;
	 8,4];
nv=size(xy,1);
ne = size(seg,1);
o=ones(ne,1);
A=sparse(seg(:,1),seg(:,2),o,nv,nv);
[i,j,o]=find(A+transpose(A));
ne1 = length(i);
o=ones(ne1,1);
A=sparse(i,j,o,nv,nv);
dim=2;
k=9;
I=speye(nv,nv);
ek=I(:,k);
%%a=(A+I)*ek; 
%%a=(A+I)*a-a;
a=A*(A+I)*ek;
jj=find(a==dim);
jj=jj(find(jj!=k))
hold on;
gplot(A,xy(:,2:3)); 
plot(xy(:,2),xy(:,3),'ro'); 
for i=1:nv
	set(text(xy(i,2)+0.01,xy(i,3)+0.01,int2str(i)),'FontSize',24);
end
axis off
hold off 
