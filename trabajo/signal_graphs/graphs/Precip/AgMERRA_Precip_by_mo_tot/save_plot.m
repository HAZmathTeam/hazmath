files = dir;
for i = 1:length(files)
  filename = files(i).name;
  if length(filename) >= 3 & filename(1:3) == 'agp'
    data = dlmread(filename);
    imwrite(ind2rgb(im2uint8(mat2gray(data)), parula(256)), strcat('images/',filename,'.png'), 'png');
  end
end
