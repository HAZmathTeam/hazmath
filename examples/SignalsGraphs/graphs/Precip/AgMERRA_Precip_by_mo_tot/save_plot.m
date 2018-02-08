files = dir;
for i = 1:length(files)
  filename = files(i).name;
  if length(filename) >= 3 & filename(1:3) == 'agp'
    data = dlmread(filename);
    im = imagesc(data);
    saveas(im, strcat('images/',filename,'.png'), 'png');
  end
end
