%% Color histogram

% choose file
[FileName,PathName,FilterIndex] = uigetfile({'*.tif'},'File Selector');

flr = imread(strcat(PathName, FileName));

close all
imshow(flr,'InitialMagnification', 50);

%%
cd(PathName);
flList = dir('*.tif');
fls = {flList.name};
length(fls)

%%
% find circles
for kk = 1:length(fls)

    %char(fls(kk))
    flr = imread(char(fls(kk)));
    imshow(flr); 
    [centers, radii, metric] = imfindcircles(flr,[10,20] ); 

    % get top 25% stongest circles
    ncirc  =  round(0.25 *length(metric), 0); 
    centersStrong5 = centers(1:ncirc,:);
    radiiStrong5 = radii(1:ncirc);
    metricStrong5 = metric(1:ncirc);
    viscircles(centersStrong5, radiiStrong5,'EdgeColor','b');
    
    
    % save circle images so I can see them later
    try
        B(:,:,1)=flr;
        B(:,:,2)=flr;
        B(:,:,3)=flr;
        for i=1:length(radiiStrong5)
            try
            theta=0:1:360;
            r=round(centersStrong5(i,1) + radiiStrong5(i)*sin(theta));
            c=round(centersStrong5(i,2) + radiiStrong5(i)*cos(theta));
              for j=1:length(r)
                  B(c(j),r(j),:)=[255 0 0];
              end
            catch
            disp('error drawing circles'); 
            end
        end
    catch
        disp('no image saved'); 
    end
    imwrite(B, strcat('circles_', char(fls(kk))));
    % save data to file
    fid=fopen('myFile.csv','at+');
    fprintf(fid,'%s,',char(fls(kk))); % filename
    fprintf(fid,'%s,',median(radiiStrong5)); % median of circle radii
    fprintf(fid,'%s,', num2str(length(radiiStrong5))); 
    fprintf(fid,'\n');
    fclose(fid);
end


%% Do this if you want to find the flowers using image segmentation

for kk = 1:length(fls)
    % read in image
    try
        flr = imread(char(fls(kk)));
    
    
        I = (flr);
        imshow(I), title('original image');

        % detect edges
        % http://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
        [~, threshold] = edge(I, 'sobel');
        fudgeFactor = 0.5;
        BWs = edge(I,'sobel', threshold * fudgeFactor);
        imshow(BWs,'InitialMagnification', 50), title('binary gradient mask');

        se90 = strel('line', 5, 90);
        se0 = strel('line', 5, 0);

        BWsdil = imdilate(BWs, [se90 se0]);
        imshow(BWsdil), title('dilated gradient mask');


        BWdfill = imfill(BWsdil, 'holes');
        imshow(BWdfill);
        title('binary image with filled holes');

        holes = BWdfill & ~BWsdil;
        bigholes = bwareaopen(holes, 1000000);
        %imshow(bigholes)
        smallholes = holes & ~bigholes;
        new = BWsdil | smallholes;
        imshow(new);
        %
        %new = BWdfill; 

        % remove border blobs
        BWnobord = imclearborder(new, 4);
        imshow(BWnobord), title('cleared border image');

        % erode images
        seD = strel('diamond',1);
        BWfinal = imerode(BWnobord,seD);
        BWfinal = imerode(BWfinal,seD);
        BWfinal = imerode(BWfinal,seD);
        imshow(BWfinal), title('segmented image');
        %
        bw = BWfinal >0;


        % get region properties
        s = regionprops(bw,'Centroid','MajorAxisLength','MinorAxisLength', 'Eccentricity', 'EquivDiameter','FilledArea');
        centroids = cat(1, s.Centroid);

        %imshow(bw)
        hold on
        plot(centroids(:,1),centroids(:,2), 'b*')
        a = [1:length([s.MajorAxisLength])]'; b = num2str(a); c = cellstr(b);
        text(centroids(:,1),centroids(:,2),  c);
        hold off

        % manually set cutoff points, based on region props
        %histogram([s.Eccentricity])
        bigEnough = [s.MajorAxisLength] > 200 & [s.MajorAxisLength] < 800 ...
            & [s.MinorAxisLength] > 100 & [s.Eccentricity] < 0.5 & [s.MajorAxisLength] ./ [s.MinorAxisLength] < 3;

        filledArea = bigEnough.*[s.FilledArea];

        [sortedX,sortingIndices] = sort(filledArea,'descend');

        %imshow(bw)
%         hold on 
%         plot(centroids(sortingIndices(1:3),1),centroids(sortingIndices(1:3),2), 'b*')
%         a = [1:length([s.MajorAxisLength])]'; b = num2str(a); c = cellstr(b);
%         text(centroids(sortingIndices(1:3),1),centroids(sortingIndices(1:3),2),  c(sortingIndices(1:3)));
%         hold off

        % remove other blobs
        labeledImage = bwlabel(bw, 8);
        %imshow(labeledImage);


        %imshow(labeledImage, []);  % Show the gray scale image.
        title('Labeled Image, from bwlabel()');

        keeperIndexes = sortingIndices(1:3);
        keeperBlobsImage = ismember(labeledImage, keeperIndexes); 

        newLabeledImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
        % Now we're done.  We have a labeled image of blobs that meet our specified criteria.
        %imshow(newLabeledImage , []);
        title('"Keeper" blobs');

        %
        %%imshow(newLabeledImage)
        %croppedImage = flr .* repmat(uint8(newLabeledImage > 0),[1,1,3]);

        % draw mask over center of flower
        foobar = s(sortingIndices);
        dd = [foobar.MajorAxisLength];
        equivDiams =  cat(1, s.EquivDiameter);
        foo = centroids(sortingIndices, :);
        bar = equivDiams(sortingIndices, :); 

        % Next create the ellipse in the image.
        centerX = foo(1,1);
        centerY = foo(1,2);
        radiusX = .25 * dd(1);

        centerX2 = foo(2,1);
        centerY2 = foo(2,2);
        radiusX2 = .25 * dd(2);

        I = newLabeledImage;
        pos1 = [centerX centerY radiusX];
        pos2 = [centerX2 centerY2 radiusX2];
        pos3 = [foo(3,1) foo(3,2) 0.25*dd(3)]; 

        rrr = insertShape(I, 'FilledCircle', pos1,  'Color' ,'Black','Opacity', 1 );
        rrr = insertShape(rrr, 'FilledCircle', pos2,  'Color' ,'Black','Opacity', 1 );
        rrr = insertShape(rrr, 'FilledCircle', pos3,  'Color' ,'Black','Opacity', 1 );
        %imshow(rrr, 'InitialMagnification', 8); 


        croppedImage = flr .* repmat(uint8(rrr(:, :, 1) > 0),[1,1,3]);
        %imshow(croppedImage); 


        %
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

        %
        vec = vertcat(yRed, yBlue, yGreen); 

        fid=fopen('myFile.csv','at+');

        x = vec;
        [rows,cols]=size(x); % row and col are backwards
        fprintf(fid,'%s,',char(fls(kk))); % filename
        fprintf(fid,'%s,',num2str(dd(1))); % major axis length for the three flowers
        fprintf(fid,'%s,',num2str(dd(2)));
        fprintf(fid,'%s,',num2str(dd(3)));
        for i=1:rows
              fprintf(fid,'%s,',num2str(x(i)));  % values from histogram
        end
        fprintf(fid,'\n');
        fclose(fid);

        % save cropped image
        imwrite(croppedImage,strcat('cropped_', char(fls(kk))));

    catch 
         warning('Problem using function.');
    end
    strcat(num2str(kk), ' of ', num2str(length(fls)))
end


%% Debugging zone

I = rgb2gray(flr);
imshow(I), title('original image');

% detect edges
% http://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
[~, threshold] = edge(I, 'sobel');
fudgeFactor = 0.5;
BWs = edge(I,'sobel', threshold * fudgeFactor);
%imshow(BWs,'InitialMagnification', 8), title('binary gradient mask');

se90 = strel('line', 10, 90);
se0 = strel('line', 10, 0);

BWsdil = imdilate(BWs, [se90 se0]);
%imshow(BWsdil), title('dilated gradient mask');


BWdfill = imfill(BWsdil, 'holes');
%imshow(BWdfill);
title('binary image with filled holes');

holes = BWdfill & ~BWsdil;
bigholes = bwareaopen(holes, 1000000);
%imshow(bigholes)
smallholes = holes & ~bigholes;
new = BWsdil | smallholes;
%imshow(new);
%
%new = BWdfill; 

% remove border blobs
BWnobord = imclearborder(new, 4);
%imshow(BWnobord), title('cleared border image');

% erode images
seD = strel('diamond',8);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
BWfinal = imerode(BWfinal,seD);
imshow(BWfinal), title('segmented image');
%
bw = BWfinal >0;


% get region properties
s = regionprops(bw,'Centroid','MajorAxisLength','MinorAxisLength', 'Eccentricity', 'EquivDiameter','FilledArea');
centroids = cat(1, s.Centroid);

imshow(bw)
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
a = [1:length([s.MajorAxisLength])]'; b = num2str(a); c = cellstr(b);
text(centroids(:,1),centroids(:,2),  c);
hold off

% manually set cutoff points, based on region props
%histogram([s.Eccentricity])
bigEnough = [s.MajorAxisLength] > 200 & [s.MajorAxisLength] < 800 ...
    & [s.MinorAxisLength] > 100 & [s.Eccentricity] < 0.5 & [s.MajorAxisLength] ./ [s.MinorAxisLength] < 3;

filledArea = bigEnough.*[s.FilledArea];

[sortedX,sortingIndices] = sort(filledArea,'descend');

%imshow(bw)
%         hold on 
%         plot(centroids(sortingIndices(1:3),1),centroids(sortingIndices(1:3),2), 'b*')
%         a = [1:length([s.MajorAxisLength])]'; b = num2str(a); c = cellstr(b);
%         text(centroids(sortingIndices(1:3),1),centroids(sortingIndices(1:3),2),  c(sortingIndices(1:3)));
%         hold off

% remove other blobs
labeledImage = bwlabel(bw, 8);
imshow(labeledImage);


%imshow(labeledImage, []);  % Show the gray scale image.
title('Labeled Image, from bwlabel()');

keeperIndexes = sortingIndices(1:3);
keeperBlobsImage = ismember(labeledImage, keeperIndexes); 

newLabeledImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
% Now we're done.  We have a labeled image of blobs that meet our specified criteria.
imshow(newLabeledImage , []);
title('"Keeper" blobs');

%
%%imshow(newLabeledImage)
%croppedImage = flr .* repmat(uint8(newLabeledImage > 0),[1,1,3]);

% draw mask over center of flower
foobar = s(sortingIndices);
dd = [foobar.MajorAxisLength];
equivDiams =  cat(1, s.EquivDiameter);
foo = centroids(sortingIndices, :);
bar = equivDiams(sortingIndices, :); 

% Next create the ellipse in the image.
centerX = foo(1,1);
centerY = foo(1,2);
radiusX = .25 * dd(1);

centerX2 = foo(2,1);
centerY2 = foo(2,2);
radiusX2 = .25 * dd(2);

I = newLabeledImage;
pos1 = [centerX centerY radiusX];
pos2 = [centerX2 centerY2 radiusX2];
pos3 = [foo(3,1) foo(3,2) 0.25*dd(3)]; 

rrr = insertShape(I, 'FilledCircle', pos1,  'Color' ,'Black','Opacity', 1 );
rrr = insertShape(rrr, 'FilledCircle', pos2,  'Color' ,'Black','Opacity', 1 );
rrr = insertShape(rrr, 'FilledCircle', pos3,  'Color' ,'Black','Opacity', 1 );
%imshow(rrr, 'InitialMagnification', 8); 


croppedImage = flr .* repmat(uint8(rrr(:, :, 1) > 0),[1,1,3]);
imshow(croppedImage); 


%
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

%
vec = vertcat(yRed, yBlue, yGreen); 

fid=fopen('myFile.csv','at+');

x = vec;
[rows,cols]=size(x); % row and col are backwards
fprintf(fid,'%s,',char(fls(kk))); % filename
fprintf(fid,'%s,',num2str(dd(1))); % major axis length for the three flowers
fprintf(fid,'%s,',num2str(dd(2)));
fprintf(fid,'%s,',num2str(dd(3)));
for i=1:rows
      fprintf(fid,'%s,',num2str(x(i)));  % values from histogram
end
fprintf(fid,'\n');
fclose(fid);

% save cropped image
imwrite(croppedImage,strcat('cropped_', char(fls(kk))));