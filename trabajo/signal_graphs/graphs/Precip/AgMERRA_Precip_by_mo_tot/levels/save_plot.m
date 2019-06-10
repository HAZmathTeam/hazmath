files = dir;
for i = 1:length(files)
  filename = files(i).name;
  [filepath, name, ext] = fileparts(filename);
  if strcmp(ext, '.data')
    data = dlmread(filename, '\t', 1, 0);
    log_th = log2(data(:,1));
    rel_err1 = data(:,3) ./ data(:,2);
    rel_err2 = data(:,4) ./ data(:,2);
    fig = figure;
    plot(log_th, rel_err1, '-o', log_th, rel_err2, '-o');
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    set(gcf,'color','w');
    frame = getframe(fig);
    imwrite(frame2im(frame), strcat('images/',name,'_levels.png'), 'png');
  end
end
