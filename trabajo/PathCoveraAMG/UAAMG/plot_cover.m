function plot_cover( G, xcoord, ycoord, cover )
% plot path covers of graph
%
% @ Xiaozhe Hu, Tufts University

pic = plot(G, 'XData', xcoord, 'YData', ycoord);
num_cover = length(cover);

for i=1:num_cover
    %highlight(pic,cover{i},'NodeColor','r','EdgeColor','r','LineWidth',3);
    highlight(pic,cover{i}, 'MarkerSize', 6, 'NodeColor','r');
    pause;
end


end

