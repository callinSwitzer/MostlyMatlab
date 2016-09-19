[filename pathname]  = uigetfile('*.avi', 'Select video');
mov = VideoReader(strcat(pathname, '/',filename));
cd(pathname);
%%
[tagfile tagpath] = uigetfile('*.csv', 'select taglist');
taglist = csvread(strcat(tagpath, '/', tagfile), 1,0);
taglist = taglist(:,1);
%%
[trackedFile trackedPath] = uigetfile('*.mat', 'select matlab file with tracked data');
load(strcat(trackedPath, '/', trackedFile));
trackingData(trackingData == 0) = NaN;
trackingData = trackingData./2;
for i = 1:size(trackingData,2)
    trackingData(:,i,1) = fixShortNanGaps(trackingData(:,i,1),10);
end

for i = 1:size(trackingData,2)
    trackingData(:,i,2) = fixShortNanGaps(trackingData(:,i,2),10);
end
frameInds = round(linspace(1,7200,61));

%% Load file with manually digitized behavior, loads as "behaviorData" object
uiopen('load');

%% Load list of behaviors
[trackedFile trackedPath] = uigetfile('*.csv', 'select behavior list');

behaviorList = readtable([trackedPath trackedFile]);

%%
t = tabulate(reshape(behaviorData, numel(behaviorData), 1));
t = array2table(t);
t.behChar = behaviorList.Behavior;
t = t(:,[4 2 3]);
t.Properties.VariableNames = {'Behavior', 'Count', 'Percent'};
%%

videoIndex = 60;
beeID = 8;

currentFrames = frameInds(videoIndex):frameInds(videoIndex+1);
currentFrames = currentFrames(55:64);
%for j = 1:
for i = 1:numel(currentFrames)
    imshow(imadjust(rgb2gray(read(mov, currentFrames(i)))));
    text(500,500, behaviorList.Behavior(behaviorData(videoIndex,beeID)), 'Color', 'G');
    hold on;
    plot(trackingData(currentFrames(i),beeID,1), trackingData(currentFrames(i),beeID,2), 'r.', 'MarkerSize', 30);
    pause(0.1);
    hold off
end

%end
%%
%Write app data to global
setappdata(0,'mov', mov);
setappdata(0,'trackingData', trackingData);
setappdata(0, 'taglist', taglist);
setappdata(0, 'behaviorData', behaviorData);
setappdata(0, 'frameInds', frameInds);
videoIndex = 1;
setappdata(0, 'videoIndex', videoIndex);

setappdata(0,'currentBee',1);
currentFrames = frameInds(videoIndex):frameInds(videoIndex+1);
setappdata(0,'currentFrames', currentFrames);
currentBeesIndex  = find(sum(~isnan(trackingData(currentFrames(55:64),:,1))) > 0);
setappdata(0, 'currentBeesIndex',currentBeesIndex);
setappdata(0,'numberOfBeesInVideo', numel(currentBeesIndex));

%Pre-load video frames
currentVideoFrames = nan(mov.Height/2, mov.Width/2, 60);
h = waitbar(0, 'Loading new video...');
for i = 1:60
    im = read(mov, (currentFrames(i+29)));
    im = rgb2gray(im);
    imr = imresize(im,0.5);
    currentVideoFrames(:,:,i) = imr;
    waitbar(i/60);
end
close(h);
setappdata(0,'currentVideoFrames', currentVideoFrames);


%Cumbersome, horrible way to extract and print individual data
behaviorData = getappdata(0,'behaviorData');
videoIndex = getappdata(0,'videoIndex');
currentBee = getappdata(0,'currentBee');
currentBeesIndex = getappdata(0,'currentBeesIndex');
numBees = getappdata(0,'numberOfBeesInVideo');
taglist = getappdata(0,'taglist');
set(handles.text1, 'String', strcat('Bee_',num2str(taglist(currentBeesIndex(currentBee))), ':_', num2str(currentBee), '_of_', num2str(numBees)))
set(handles.text2, 'String', strcat('Video_', num2str(videoIndex), '_of_', num2str(60)));
set(handles.text3, 'String', strcat('recorded behavior:_', num2str(behaviorData(videoIndex, currentBeesIndex(currentBee)))));