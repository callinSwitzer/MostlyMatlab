%% Color histogram

% choose file
[FileName,PathName,FilterIndex] = uigetfile('*.tif')
%% change directory
cd(PathName);

%%
flr1 = imread(strcat(PathName, FileName));
%%
close all
imshow(flr1,'InitialMagnification', 50);






%%
%%
%%
%%
%%
%%
%%
%% Do this if you want to select regions by hand
close all
% select region of interest
[BW, xi, yi] = roipoly(flr);


%%

croppedImage = flr .* repmat(uint8(BW),[1,1,3]);
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

% you can look at color schemes here: http://www.colorschemer.com/online.html


%%
%%
%%
%%
%%


%% Do this if you want to find the flowers using image segmentation
I = rgb2gray(flr);
imshow(I), title('original image');

%% detect edges
% http://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
[~, threshold] = edge(I, 'sobel', 0.016);
fudgeFactor = 0.5;
BWs = edge(I,'sobel', threshold * fudgeFactor);
imshow(BWs), title('binary gradient mask');

se90 = strel('line', 10, 90);
se0 = strel('line', 10, 0);

BWsdil = imdilate(BWs, [se90 se0]);
imshow(BWsdil), title('dilated gradient mask');


BWdfill = imfill(BWsdil, 'holes');
imshow(BWdfill);
title('binary image with filled holes');

% remove border blobs
BWnobord = imclearborder(BWdfill, 4);
imshow(BWnobord), title('cleared border image');

% erode images
seD = strel('diamond',8);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
BWfinal = imerode(BWfinal,seD);
imshow(BWfinal), title('segmented image');
%%
bw = BWfinal >0;


%% get region properties
s = regionprops(bw,'Centroid','MajorAxisLength','MinorAxisLength', 'Eccentricity')
centroids = cat(1, s.Centroid);

imshow(bw)
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off

%% manually set cutoff points, based on region props
histogram([s.Eccentricity])
bigEnough = [s.MajorAxisLength] > 200 & [s.MajorAxisLength] < 5000 ...
    & [s.MinorAxisLength] > 500 & [s.Eccentricity] < 0.5;
%% show centroids

imshow(bw)
hold on
plot(centroids(bigEnough,1),centroids(bigEnough,2), 'b*')
hold off

%% fill  small regions
labeledImage = bwlabel(bw, 8);
imshow(labeledImage);


imshow(labeledImage, []);  % Show the gray scale image.
title('Labeled Image, from bwlabel()');

keeperIndexes = find(bigEnough);
keeperBlobsImage = ismember(labeledImage, keeperIndexes); 

newLabeledImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
% Now we're done.  We have a labeled image of blobs that meet our specified criteria.
imshow(newLabeledImage , []);
title('"Keeper" blobs');

%%
imshow(newLabeledImage)
croppedImage = flr .* repmat(uint8(newLabeledImage > 0),[1,1,3]);
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