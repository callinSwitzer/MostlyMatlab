%%Load in compiled data (i.e. 'allData.mat' file)
%cd(uigetdir());
fps = 2;
%% Remove any extra frames
nestPos = preNest(1:7200,:,:);
postNest = postNest(1:7200,:,:);

%%
preDiffVel = abs(diff(preNest(:,:,1:2)));
preVels = sqrt(preDiffVel(:,:,1).^2 + preDiffVel(:,:,2).^2);
preVels = preVels*fps; %Correct for frame rate

postDiffVel = abs(diff(postNest(:,:,1:2)));
postVels = sqrt(postDiffVel(:,:,1).^2 + postDiffVel(:,:,2).^2);
postVels = postVels*fps; %Correct for frame rate

%%
activeThresh = 10^-3.9;

%% Work in progress - classify behaviors
onNestBin = nan(size(preNest,1), size(preNest,2));


for i = 1:numel(tags)
    %%
    curPos = nan(size(preNest,1), 3);
    curPos(:,1) = preNest(:,i,1);
    curPos(:,2) = preNest(:,i,2);
    
    for j = 1:size(preNest,1)
        %%
        if ~isnan(curPos(j,1))
            dists = sqrt((broodPre(:,1) - curPos(j,1)).^2+(broodPre(:,2)-curPos(j,2)).^2);
            ind = dists < 0.015;
            if sum(ind) > 0 %If the bee is close to some aspect of the hive
                curPos(j,3) = 1;
                
            else
                curPos(j,3) = 0;
            end
        end
        
    end
    onNestBin(:,i) = curPos(:,3);
end

%% Visualize nest location over time

binPic = onNestBin;

binPic(isnan(onNestBin)) =-1;
binPic(onNestBin == 0) = 1;
binPic(onNestBin == 1) = 0;

colormap summer;
imagesc(binPic);

%%
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
onNestSubset = nan(size(behaviorData,1), size(behaviorData,2));



for i = 1:size(onNestSubset,1)
    
    %%
    currentFrames = frameInds(i):frameInds(i+1);
    currentFrames = currentFrames(55:64);
    
    %% loop across bees
    for j = 1:size(onNestSubset,2)
        %%
        onNestSubset(i,j) = nanmedian(onNestBin(currentFrames,j));
        
    end
    
end


%% tabulate behavior by nest binary
t = tabulate(behaviorData(onNestSubset == 1));
if size(t,1) < 20
    t(20,1) = 0
end
t = array2table(t);
t.behChar = behaviorList.Behavior;
t = t(:,[4 2 3]);
t.Properties.VariableNames = {'Behavior', 'Count', 'Percent'};

onNestPort = t;

t = tabulate(behaviorData(onNestSubset == 0));
if size(t,1) < 20
    t(20,1) = 0
end
t = array2table(t);
t.behChar = behaviorList.Behavior;
t = t(:,[4 2 3]);
t.Properties.VariableNames = {'Behavior', 'Count', 'Percent'};
offNestPort = t;


%% Make single table

totPort = onNestPort;
totPort.onNestPerc = onNestPort.Percent;
totPort.offNestPerc = offNestPort.Percent;
totPort = totPort(:,[1 4 5]);

%Remove behaviors that don't represent at least 5% in either off or on nest
totPort = totPort((totPort.onNestPerc > 3 | totPort.offNestPerc > 3),:);

diffs = totPort.onNestPerc./totPort.offNestPerc;
[x i] = sort(diffs);
totPort = totPort(i,:);
% for i = 1:size(totPort,1)
%     s = totPort.onNestPerc(i) + totPort.offNestPerc(i);
%     totPort.onNestPerc(i) = totPort.onNestPerc(i)/s;
%     totPort.offNestPerc(i) = totPort.offNestPerc(i)/s;
% end
figure(1);
bar([totPort.offNestPerc totPort.onNestPerc]);

% Create total sums for portions in "on-nest behaviors" and "off-nest"
% behaviors
bp = nan(2,2);
bp(1,1) = sum(totPort.offNestPerc(1:5));
bp(1,2) = sum(totPort.offNestPerc(7:10));
bp(2,1) = sum(totPort.onNestPerc(1:5));
bp(2,2) = sum(totPort.onNestPerc(7:10));

%Normalize rows
figure(2);
bp(1,:) = bp(1,:)./sum(bp(1,:));
bp(2,:) = bp(2,:)./sum(bp(2,:));
bar(bp,'stacked');
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


