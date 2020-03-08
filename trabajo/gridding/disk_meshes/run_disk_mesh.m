%% R is the radius of the circle. the mesh size is about
%% h=R*2^(-level-1)
R=1; level=2;
[t,xy]=disk_mesh(R,level);
stat0=m_hazw(t,xy,'zzz.haz');

%%[VT,EV,ET,TT] = bndrypts(t);
