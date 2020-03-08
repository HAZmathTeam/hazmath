%% R is the radius of the circle. the mesh size is about
%% h=R*2^(-level-1)
R=1; level=0;
[t,xy]=disk_mesh(R,level);
bcodes = bndrypts(t);
%% pick a filename
fname0=['circle_level',int2str(level),'.haz'];
stat0=m_hazw(fname0,t,xy,bcodes)
%%% just to check whether we really found the boundary points
bp=find(bcodes); hold on ; plot(xy(bp,1),xy(bp,2),'ro',xy(bp,1),xy(bp,2),'r*'); hold off;
