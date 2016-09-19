%% Color histogram

% choose file
[FileName,PathName,FilterIndex] = uigetfile('*.tif')
%% change directory
cd(PathName);

imageNames = dir(fullfile(PathName,'*.TIF'));
imageNames = {imageNames.name}';

%%
if ~exist('ii')
    ii = 1; 
else ii = ii + 1; 
end

for ii = 2:length(imageNames)

flr = imread(strcat(PathName, imageNames{ii}));
ii
%
close all
%imshow(flr,'InitialMagnification', 50);
% Do this if you want to find the pollen using image segmentation
I = (flr);
imshow(I), title('original image');

% detect edges
% http://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
[~, threshold] = edge(I, 'sobel');
fudgeFactor = 0.5;
BWs = edge(I,'sobel', threshold * fudgeFactor);
imshow(BWs), title('binary gradient mask');

se90 = strel('disk', 10);

BWsdil = imdilate(BWs, [se90]);
imshow(BWsdil), title('dilated gradient mask'); 


BWdfill = imfill(BWsdil, 'holes');
imshow(BWdfill);
title('binary image with filled holes');

% remove border blobs
BWnobord = imclearborder(BWdfill, 4);
imshow(BWnobord), title('cleared border image');

% erode images
seD = strel('disk',10);
BWfinal = imerode(BWnobord,seD);
%BWfinal = imerode(BWfinal,seD);
%BWfinal = imerode(BWfinal,seD);
imshow(BWfinal), title('segmented image');


subplot(1,2,1); 
imshow(flr); 
subplot(1,2,2); 
imshow(BWfinal); 

%
figure;
rgb = flr;
mmme = imshow(rgb);
im = BWfinal;
hold on;
him=imshow(im);
set(him,'AlphaData',0.4);

%
bw = BWfinal >0;


% get region properties
s = regionprops(bw,'Centroid','MajorAxisLength','MinorAxisLength', 'Eccentricity', 'Area')
centroids = cat(1, s.Centroid);

imshow(bw)
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off

% manually set cutoff points, based on region props
% histogram([s.Eccentricity])
bigEnough = [s.MajorAxisLength] > 10 & [s.MajorAxisLength] < 50 ...
    & [s.MinorAxisLength] > 10 & [s.Eccentricity] < 0.8;
% show centroids
% 
% imshow(bw)
% hold on
% plot(centroids(bigEnough,1),centroids(bigEnough,2), 'b*')
% hold off

% keep specific regions (single pollen grains)
labeledImage = bwlabel(bw, 8);
imshow(labeledImage);


imshow(labeledImage, []);  % Show the gray scale image.
title('Labeled Image, from bwlabel()');

keeperIndexes = find(bigEnough);
keeperBlobsImage = ismember(labeledImage, keeperIndexes); 

newLabeledImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
% Now we're done.  We have a labeled image of blobs that meet our specified criteria.
% imshow(newLabeledImage , []);
imshow(newLabeledImage);
title('"Keeper" blobs');
hold on;
him=imshow(flr > 0);
set(him,'AlphaData',0.3);

%
foo = [s.Area]; 
foo(keeperIndexes);

plot(centroids(keeperIndexes,1),centroids(keeperIndexes,2), 'b*')

% mean area of individual pollen grains
median(foo(keeperIndexes))

% save binary image used for analysis
imwrite(newLabeledImage, strcat(PathName ,'BinaryThresh_', imageNames{ii}));


% save measurements
fid=fopen('PollenSizes.csv','at+');
fprintf(fid,'%s,',char(imageNames{ii})); % filename
fprintf(fid,'%s,',num2str(median(foo(keeperIndexes)))); % median pollen area (pixels)
fprintf(fid,'%s,',num2str(length(keeperIndexes))); % number of pollen grains used in calculation
fprintf(fid,'\n'); % newline
fclose(fid);

end

%%
imshow(newLabeledImage)
croppedImage = flr .* uint8(newLabeledImage > 0);
%%
imshow(croppedImage);





%%
%Split into RGB Channels
Red = croppedImage(:,:,1);
Green = croppedImage(:,:,2);
Blue = croppedImage(:,:,3);

%Get histValues for each channel
[yRed, x] = imhist(Red);
[yGreen, x] = imhist(Green);
[yBlue, x] = imhist(Blue);

%Plot them together in one plot
plot(x, yRed, 'Red', x, yGreen, 'Green', x, yBlue, 'Blue');
xlim([5,255])
%%
subplot(3,1,1)
plot(x, yRed, 'Red')
xlim([5,255])

subplot(3,1,2)
plot(x, yGreen, 'Green')
xlim([5,255])

subplot(3,1,3)
plot(x, yBlue, 'Blue')
xlim([5,255])