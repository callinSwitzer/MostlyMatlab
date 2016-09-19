%%Load in compiled data (i.e. 'allData.mat' file)
%cd(uigetdir());
fps = 2;
%% Remove any extra frames
preNest = preNest(1:7200,:,:);
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
framesActivePre = sum(preVels > activeThresh);
velocityFramesPre = sum(~isnan(preVels));
framesDetectedPre = sum(~isnan(preNest(:,:,1)));
portionActivePre = framesActivePre./framesDetectedPre;

framesActivePost = sum(postVels > activeThresh);
velocityFramesPost = sum(~isnan(postVels));
framesDetectedPost = sum(~isnan(postNest(:,:,1)));
portionActivePost = framesActivePost./framesDetectedPost;


activeVelsPre = preVels;
activeVelsPre(activeVelsPre < activeThresh) = NaN;

activeVelsPost = postVels;
activeVelsPost(activeVelsPost < activeThresh) = NaN;


activeVelPre = nanmean(activeVelsPre); %Take mean velocity when bee is active
activeVelPost = nanmean(activeVelsPost); %Take mean velocity when active

meanVelPre = nanmean(preVels);
meanVelPost = nanmean(postVels);

medianVelPre = nanmedian(preVels);
medianVelPost = nanmedian(postVels);
%% Calculate average hive position
preMean = nanmean(nanmean(preNest));
preMean = [preMean(:,:,1) preMean(:,:,2)];
distx = preNest(:,:,1) - preMean(1);
disty = preNest(:,:,2) - preMean(2);
totDist = sqrt(distx.^2 + disty.^2);
socDistPre = nanmean(totDist);

postMean = nanmean(nanmean(postNest));
postMean = [postMean(:,:,1) postMean(:,:,2)];
distx = postNest(:,:,1) - postMean(1);
disty = postNest(:,:,2) - postMean(2);
totDist = sqrt(distx.^2 + disty.^2);
socDistPost = nanmean(totDist);
%% Alternative Social Distance - mean istantaneous distance to all other bees
% socDistPre = nanmean(nanmean(distMatPre,3));
% socDistPost = nanmean(nanmean(distMatPost,3));
% clear preMean;
% clear postMean;
%% Calculate distance from the queen

distx = bsxfun(@minus, preNest(:,:,1), preNest(:,1,1));
disty = bsxfun(@minus, preNest(:,:,2), preNest(:,1,2));
totDist = sqrt(distx.^2 + disty.^2);
queenDistPre = nanmean(totDist);
queenMeanPre(1) = nanmean(preNest(:,1,1));
queenMeanPre(2) = nanmean(preNest(:,1,2));

distx = bsxfun(@minus, postNest(:,:,1), postNest(:,1,1));
disty = bsxfun(@minus, postNest(:,:,2), postNest(:,1,2));
totDist = sqrt(distx.^2 + disty.^2);
queenDistPost = nanmean(totDist);
queenMeanPost(1) = nanmean(postNest(:,1,1));
queenMeanPost(2) = nanmean(postNest(:,1,2));

%% Look at social interactions
distThresh = 0.01; %Distance threshold from another bee to be considered touching, in m
intMatPre = distMatPre < distThresh & ~isnan(distMatPre);
intMatPost = distMatPost < distThresh & ~isnan(distMatPost);

numInteractionsPre = sum(sum(intMatPre,3),1)/2; %How many total interactions for each bee pre treatment??
numInteractionsPost = sum(sum(intMatPost,3),1)/2; %How many total interactions post treatment?

%This metric seems problematic, some unreasonably high outliers (one bee
%connected to 32 others per frame)
interactionRatePre = numInteractionsPre./framesDetectedPre*fps;
interactionRatePost = numInteractionsPost./framesDetectedPost*fps;


qualThresh = 100; %Need to be found in at least 200 frames to have your interaction rate tested

interactionRatePre(framesDetectedPre < qualThresh) = NaN;
interactionRatePost(framesDetectedPost < qualThresh) = NaN;
%% Look at foraging data
%Work in progress
timeInt = 120/(3600*24); %time interval to count unique foraging bouts in hours

nectarPre = preFeederData(:,:,1) < 600;
nectarBoutsPre = nan(numel(tags),1);
for i = 1:numel(tags)
    dat = nectarPre(:,i);
    times = preFeederTimes(dat);
    if numel(times) >0
        nectarBoutsPre(i) = sum(diff(times) > timeInt) +1;
    end
end

pollenPre = preFeederData(:,:,1) > 600;
pollenBoutsPre = nan(numel(tags),1);
for i = 1:numel(tags)
    dat = pollenPre(:,i);
    times = preFeederTimes(dat);
    if numel(times) >0
        pollenBoutsPre(i) = sum(diff(times) > timeInt) +1;
    end
end

totalPre = ~isnan(preFeederData(:,:,1));
totalBoutsPre = nan(numel(tags),1);
for i = 1:numel(tags)
    dat = totalPre(:,i);
    times = preFeederTimes(dat);
    if numel(times) >0
        totalBoutsPre(i) = sum(diff(times) > timeInt) +1;
    end
end

nectarPost = postFeederData(:,:,1) < 600;
nectarBoutsPost = nan(numel(tags),1);
for i = 1:numel(tags)
    dat = nectarPost(:,i);
    times = postFeederTimes(dat);
    if numel(times) >0
        nectarBoutsPost(i) = sum(diff(times) > timeInt) +1;
    end
end

pollenPost = postFeederData(:,:,1) > 600;
pollenBoutsPost = nan(numel(tags),1);
for i = 1:numel(tags)
    dat = pollenPost(:,i);
    times = postFeederTimes(dat);
    if numel(times) >0
        pollenBoutsPost(i) = sum(diff(times) > timeInt) +1;
    end
end

totalPost = ~isnan(postFeederData(:,:,1));
totalBoutsPost = nan(numel(tags),1);
for i = 1:numel(tags)
    dat = totalPost(:,i);
    times = preFeederTimes(dat);
    if numel(times) >0
        totalBoutsPost(i) = sum(diff(times) > timeInt) +1;
    end
end

nectarTimePre = sum(nectarPre);
nectarTimePost = sum(nectarPost);
pollenTimePre = sum(pollenPre);
pollenTimePost = sum(pollenPost);

pollenBoutsPre(isnan(pollenBoutsPre)) = 0;
pollenBoutsPost(isnan(pollenBoutsPost)) = 0;
nectarBoutsPre(isnan(nectarBoutsPre)) = 0;
nectarBoutsPost(isnan(nectarBoutsPost)) = 0;
totalBoutsPre(isnan(totalBoutsPre)) = 0;
totalBoutsPost(isnan(totalBoutsPost)) = 0;

%% Look at correlation of space use between and after treatment
% Calculate individual spatial pdfs
xrange = -0.02:.025:.27; %x scale: -2:27 cm in 1 cm intervals
yrange = -.02:.025:.21; %yscale: -2:22 in 1 cm intervals
spatialPrefsPre = nan(numel(yrange)-1, numel(xrange) - 1,numel(tags));

for j = 1:numel(tags)
    %%
    
    if framesDetectedPre(j) > 200 %If there are more than 100 counts for today, include, otherwise leave blank
        out = count2D(preNest(:,j,1), preNest(:,j,2),xrange,yrange);
        outNorm = out./sum(sum(out));
        
        spatialPrefsPre(:,:,j) = outNorm;
        
    end
end

spatialPrefsPost = nan(numel(yrange)-1, numel(xrange) - 1,numel(tags));
for j = 1:numel(tags)
    %%
    
    if framesDetectedPost(j) > 200 %If there are more than 100 counts for today, include, otherwise leave blank
        out = count2D(postNest(:,j,1), postNest(:,j,2),xrange,yrange);
        outNorm = out./sum(sum(out));
        spatialPrefsPost(:,:,j) = outNorm;
        
    end
end

% %visual check
% for i = 1:10
%    subplot(2,1,1)
%    imagesc(spatialPrefsPre(:,:,i));
%    subplot(2,1,2);
%    plot(preNest(:,i,1), preNest(:,i,2));
%    xlim([min(xrange) max(xrange)]);
%    ylim([min(yrange) max(yrange)]);
%    pause(0.5);
%
% end

%Get correlation between pre and post spatial preferences
prePostCorr = nan(numel(tags),1);
for i = 1:numel(tags)
    
    prePostCorr(i) = corr2(spatialPrefsPre(:,:,i), spatialPrefsPost(:,:,i));
end

%% Work in progress - classify behaviors
framesOnNestPre = nan(numel(tags),1);
framesOnNestPost = nan(numel(tags),1);
nestDistThresh = 0.015;
vis = 0;
for i = 1:numel(tags)
    %%
    curPos = nan(size(preNest,1), 3);
    curPos(:,1) = preNest(:,i,1);
    curPos(:,2) = preNest(:,i,2);
    
    for j = 1:size(preNest,1)
        %%
        dists = sqrt((broodPre(:,1) - curPos(j,1)).^2+(broodPre(:,2)-curPos(j,2)).^2);
        ind = dists < nestDistThresh; %distance threshold for detecting "on nest"
      if sum(ind) > 0 %If the bee is close to some aspect of the hive
          curPos(j,3) = 1;
      end
          
        if vis == 1
        plot(curPos(j,1), curPos(j,2), '.', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
        hold on;
        plotbroodTrans(broodPre, 200, 0.3);
        plot(broodPre(ind,1), broodPre(ind,2), 'ro','MarkerSize', 20);
        xlim([0 0.25]);
        ylim([0 0.2]);
        title(num2str(curPos(j,3)));
        drawnow;
        hold off;
        end
        
    end
    framesOnNestPre(i) = nansum(curPos(:,3));
end


%Repeat for post
for i = 1:numel(tags)
    %%
    curPos = nan(size(postNest,1), 3);
    curPos(:,1) = postNest(:,i,1);
    curPos(:,2) = postNest(:,i,2);
    
    for j = 1:size(postNest,1)
        %%
        dists = sqrt((broodPost(:,1) - curPos(j,1)).^2+(broodPost(:,2)-curPos(j,2)).^2);
        ind = dists < nestDistThresh;
      if sum(ind) > 0 %If the bee is close to some aspect of the hive
          curPos(j,3) = 1;
      end
          
        if vis == 1
        plot(curPos(j,1), curPos(j,2), '.', 'MarkerSize', 20, 'MarkerFaceColor', 'b');
        hold on;
        plotbroodTrans(broodPost, 200, 0.3);
        plot(broodPost(ind,1), broodPost(ind,2), 'ro','MarkerSize', 20);
        xlim([0 0.25]);
        ylim([0 0.2]);
        title(num2str(curPos(j,3)));
        drawnow;
        hold off;
        end
        
    end
    framesOnNestPost(i) = nansum(curPos(:,3));
end

%% Generate social network
weightedNetworkPost = sum(distMatPost < 0.01,3);

% I think this should say "distMatPre"
weightedNetworkPre = sum(distMatPre < 0.01, 3);
%%
% writetable(array2table(weightedNetworkPre), 'weightedInteractionNetworkPre.csv')
% writetable(array2table(weightedNetworkPost), 'weightedInteractionNetworkPost.csv')

%% Output simplified csv
masterData = array2table([tags orTagTreat framesActivePre' framesActivePost' framesDetectedPre' framesDetectedPost' activeVelPre' activeVelPost' meanVelPre' meanVelPost' socDistPre' socDistPost' queenDistPre' queenDistPost' interactionRatePre' interactionRatePost' nectarBoutsPost nectarBoutsPre pollenBoutsPost pollenBoutsPre framesOnNestPre framesOnNestPost velocityFramesPre' velocityFramesPost' medianVelPre' medianVelPost' totalBoutsPre totalBoutsPost]);
masterData.Properties.VariableNames = {'tags', 'treatment', 'framesActivePre', 'framesActivePost', 'framesDetectedPre', 'framesDetectedPost','activeVelPre', 'activeVelPost', 'meanVelPre', 'meanVelPost', 'socDistPre', 'socDistPost', 'queenDistPre', 'queenDistPost', 'interactionRatePre', 'interactionRatePost', 'nectarBoutsPost','nectarBoutsPre', 'pollenBoutsPost', 'pollenBoutsPre', 'framesOnNestPre', 'framesOnNestPost', 'numVelocityFramesPre', 'numVelocityFramesPost', 'medianVelPre', 'medianVelPost', 'totalForageBoutsPre', 'totalForageBoutsPost'};
writetable(masterData,'masterDataForR.csv');