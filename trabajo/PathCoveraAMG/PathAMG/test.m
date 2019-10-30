%autorun

%generate graph on regular mesh
%[AG,xy,bpts,ipts]=get_mesh_graph(1);
%AG = A;
%translate graph laplacian AG to graph G

%%%%%%%%video output%%%%%%%%%%
%  writerObj = VideoWriter('cover.avi');
%  writerObj = VideoWriter('C:\Users\Junyuan Lin\Box Sync\Grad School Study\Research\candidacy exam\cover', 'MPEG-4');
%  witerObj = VideoWriter('cover.avi','Uncompressed AVI');
%  writerObj.FrameRate = 4;
%  writerObj.Quality= 100;
% open(writerObj)
% v = VideoWriter('cover.avi');
% v = VideoWriter('C:\Users\Junyuan Lin\Box Sync\Grad School Study\Research\candidacy exam\cover', 'MPEG-4');
% open(v);

G = graph(AG,'omitself');
G.Edges.Weight = -G.Edges.Weight;
n = numnodes(G);
%plot the graph
pl = plot(G,'EdgeLabel',G.Edges.Weight);
hold on

%place holder for the path cover
cover = {};

%set-up
set = [];

while ~isempty(find(AG))
%find a path of G
[ path,path_fwd,path_bkwd ] = pathfinder( AG );
cover{end+1} = path; 

%pause(2)
for i = 1:length(path_fwd)-1
    highlight(pl,[path_fwd(i) path_fwd(i+1)],'EdgeColor','r','LineWidth',5);
    %pause(1.5)
end
for i = 2:length(path_bkwd)-1
    highlight(pl,[path_bkwd(i) path_bkwd(i+1)],'EdgeColor','r','LineWidth',5);
    %pause(1.5)
end
set = union(set,path);

%get rid of the used nodes on the path
del_nodes = path;
for i = 1:length(del_nodes)
    AG(del_nodes(i),:) = zeros(1,n);
    AG(:,del_nodes(i)) = zeros(n,1);
end
% % regenerate the graph
% G = graph(AG,'omitself');
% G.Edges.Weight = -G.Edges.Weight;
%frame = getframe;
%writeVideo(writerObj,frame);
%close(v);
end


%close(writerObj);

isolated_pts = setdiff(1:n,set)

cover


