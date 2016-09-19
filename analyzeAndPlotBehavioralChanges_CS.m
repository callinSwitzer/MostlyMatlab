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
distThresh = 0.015; %Distance threshold from another bee to be considered touching, in m
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
vis = 0;
for i = 1:numel(tags)
    %%
    curPos = nan(size(preNest,1), 3);
    curPos(:,1) = preNest(:,i,1);
    curPos(:,2) = preNest(:,i,2);
    
    for j = 1:size(preNest,1)
        %%
        dists = sqrt((broodPre(:,1) - curPos(j,1)).^2+(broodPre(:,2)-curPos(j,2)).^2);
        ind = dists < 0.0075;
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
        ind = dists < 0.0075;
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
writetable(array2table(weightedNetworkPre), 'weightedInteractionNetworkPre.csv')
writetable(array2table(weightedNetworkPost), 'weightedInteractionNetworkPost.csv')
%% Make plots!

% Plot all bee positions
figure(1);
plot(preNest(:,:,1), preNest(:,:,2));
hold on;
plot(preNest(:,1,1), preNest(:,1,2), 'b', 'LineWidth', 1.5);
plot(queenMeanPre(1), queenMeanPre(2), 'b.', 'MarkerSize', 50);
plot(preMean(1), preMean(2), 'g.', 'MarkerSize',50);

%change in velocity
figure(2);
currVar = medianVelPost - medianVelPre;
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
text([1 2 3 4]-.1, repmat(.004,4,1), {'untreated', 'control', 'low', 'high'})
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Change in velocity (m/s)');
hold off;

%Change in portion of time active
figure(3);
currVar = portionActivePost - portionActivePre;
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Change in portion of time active (m/s)');

%Change in distance from hive center
figure(4);
currVar = socDistPost - socDistPre;
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Mean distance from other bees (m)');

%Change in distance from queen
figure(5);
currVar = queenDistPost - queenDistPre;
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Change in distance from queen (m)');

%Change in interaction rate
figure(6);
currVar = interactionRatePost - interactionRatePre;
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Change in interaction rate (m)');

figure(7);
plot(socDistPre,interactionRatePre,'bo');
hold on;
plot(socDistPost,interactionRatePost,'r.', 'MarkerSize', 20);
text(.07,2.4, 'Pre-Treatment', 'Color', 'b');
text(.03,.5, 'Post-Treatment', 'Color', 'r');
xlabel('Mean distance to nest center (m)');
ylabel('Interaction frequency (interactions per second)');
hold off;


figure(8);
currVar = (nectarBoutsPost + pollenBoutsPost) - (nectarBoutsPre + pollenBoutsPre);
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Times found on feeder (nectar or pollen)');

% figure(8);
% currVar = nectarTimePost - nectarTimePre;
% boxplot(currVar, orTagTreat);
% hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
% hline(0, 'k', 'No change');
% xlabel('treatment');
% ylabel('Times found on pollen');
figure(9);
currVar = prePostCorr;
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Correlation in spatial position pre- and post- treatment');
%%
figure(10);
currVar = framesOnNestPost./framesDetectedPost'- framesOnNestPre./framesDetectedPre';
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('portion of time on nest');
%% Plot bee positions coded by treatment
figure(11);
xl = [0 0.25];
yl = [0 0.2];
lw = 0.5;
bs = 100;
alph = 0.2;
contInd = find(orTagTreat == 1);
lowInd = find(orTagTreat == 2);
highInd = find(orTagTreat == 3);

subplot(3,2,1);
plot(preNest(:,contInd,1), preNest(:,contInd,2), 'LineWidth', lw, 'Color', 'b');
hold on;
xlim(xl);
ylim(yl);
axis equal
title('control');
plotbroodGreyTrans(broodPre,bs, alph);

subplot(3,2,3);
plot(preNest(:,lowInd,1), preNest(:,lowInd,2), 'LineWidth', lw, 'Color','g');
hold on
xlim(xl);
ylim(yl);
axis equal
title('low');
plotbroodGreyTrans(broodPre, bs, alph);

subplot(3,2,5);
plot(preNest(:,highInd,1), preNest(:,highInd,2), 'LineWidth', lw, 'Color', 'r');
hold on
xlim(xl);
ylim(yl);
axis equal
title('high');
plotbroodGreyTrans(broodPre, bs, alph);

subplot(3,2,2);
plot(postNest(:,contInd,1), postNest(:,contInd,2), 'LineWidth', lw, 'Color', 'b');
hold on
xlim(xl);
ylim(yl);
axis equal
title('control');
plotbroodGreyTrans(broodPost, bs, alph);

subplot(3,2,4);
plot(postNest(:,lowInd,1), postNest(:,lowInd,2), 'LineWidth', lw, 'Color', 'g');
hold on
xlim(xl);
ylim(yl);
axis equal
title('low');
plotbroodGreyTrans(broodPost, bs, alph);

subplot(3,2,6);
plot(postNest(:,highInd,1), postNest(:,highInd,2), 'LineWidth', lw, 'Color', 'r');
hold on
xlim(xl);
ylim(yl);
axis equal
title('high');
plotbroodGreyTrans(broodPost, bs, alph);

%%
figure(12);
hist(log10(reshape(preVels, numel(preVels), 1)), 500);
vline(-3.9);
xlabel('log10(nest velocity)');
ylabel('count');


%% Plot foraging data pre and post
figure(13);
xl = [0 1000];
yl = [0 1000];
ms = 20;
contInd = find(orTagTreat == 1);
lowInd = find(orTagTreat == 2);
highInd = find(orTagTreat == 3);

subplot(3,2,1);
plot(preFeederData(:,contInd,1), preFeederData(:,contInd,2), '.','MarkerSize', ms);
xlim(xl);
ylim(yl);
axis equal
title('control');
vline(600);

subplot(3,2,3);
plot(preFeederData(:,lowInd,1), preFeederData(:,lowInd,2), '.','MarkerSize', ms);
xlim(xl);
ylim(yl);
axis equal
title('low');
vline(600);

subplot(3,2,5);
plot(preFeederData(:,highInd,1), preFeederData(:,highInd,2), '.','MarkerSize', ms);
xlim(xl);
ylim(yl);
axis equal
title('high');
vline(600);

subplot(3,2,2);
plot(postFeederData(:,contInd,1), postFeederData(:,contInd,2), '.','MarkerSize', ms);
xlim(xl);
ylim(yl);
axis equal
title('control');
vline(600);

subplot(3,2,4);
plot(postFeederData(:,lowInd,1), postFeederData(:,lowInd,2), '.','MarkerSize', ms);
xlim(xl);
ylim(yl);
axis equal
title('low');
vline(600);

subplot(3,2,6);
plot(postFeederData(:,highInd,1), postFeederData(:,highInd,2), '.','MarkerSize', ms);
xlim(xl);
ylim(yl);
axis equal
title('high');
vline(600);


%% Separate nectar and pollen bouts
figure(14);
currVar = nectarBoutsPost - nectarBoutsPre;
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Times found on feeder (nectar or pollen)');

figure(15);
currVar = pollenBoutsPost - pollenBoutsPre;
boxplot(currVar, orTagTreat);
hold on; scatter(1+orTagTreat, currVar, 70,'r', 'filled', 'MarkerFaceAlpha', 0.3);
hline(0, 'k', 'No change');
xlabel('treatment');
ylabel('Times found on feeder (nectar or pollen)');

%% Plot indiviudal patterns, by hive
figure(16);
subplot(2,1,1);
plot(preNest(:,orTagTreat == 1,1), preNest(:,orTagTreat == 1,2), 'b');
hold on;
plot(preNest(:,orTagTreat == 2,1), preNest(:,orTagTreat == 2,2), 'g');
plot(preNest(:,orTagTreat == 3,1), preNest(:,orTagTreat == 3,2), 'r');
plot(preNest(:,orTagTreat == 1,1), preNest(:,orTagTreat == 1,2), 'b.');
plot(preNest(:,orTagTreat == 2,1), preNest(:,orTagTreat == 2,2), 'g.');
plot(preNest(:,orTagTreat == 3,1), preNest(:,orTagTreat == 3,2), 'r.');
plotbroodGreyTrans(broodPre,200, 0.8);
axis equal

subplot(2,1,2);
plot(postNest(:,orTagTreat == 1,1), postNest(:,orTagTreat == 1,2), 'b');
hold on;
plot(postNest(:,orTagTreat == 2,1), postNest(:,orTagTreat == 2,2), 'g');
plot(postNest(:,orTagTreat == 3,1), postNest(:,orTagTreat == 3,2), 'r');
plot(postNest(:,orTagTreat == 1,1), postNest(:,orTagTreat == 1,2), 'b.');
plot(postNest(:,orTagTreat == 2,1), postNest(:,orTagTreat == 2,2), 'g.');
plot(postNest(:,orTagTreat == 3,1), postNest(:,orTagTreat == 3,2), 'r.');
plotbroodGreyTrans(broodPre,200, 0.8);
axis equal

%% similar to above, but subplot instead of overlay each treatment
figure(17);
xl = [0 0.25];
yl = [0 0.2];
sz = 8;
alph = 0.01;
broodAlph = 0.1
subplot(2,3,1);
X = preNest(:,orTagTreat == 1,1);
X = reshape(X, numel(X), 1);
Y = preNest(:,orTagTreat == 1,2);
Y = reshape(Y,numel(Y),1);

%scatter(reshape(X,numel(X), 1), reshape(Y, numel(Y),1), sz,'b', 'filled', 'MarkerFaceAlpha', alph);
%[bandwidth density X Y] = kde2d(data, 2^5);

hold on;
%plot(preNest(:,orTagTreat == 1,1), preNest(:,orTagTreat == 1,2), 'b');
scatter(reshape(X,numel(X), 1), reshape(Y, numel(Y),1), sz,'g', 'filled', 'MarkerFaceAlpha', alph);
plotbroodGreyTrans(broodPre,60, broodAlph);
xlim(xl);
ylim(yl);
axis equal
hold off

subplot(2,3,2);
X = preNest(:,orTagTreat == 2,1);
Y = preNest(:,orTagTreat == 2,2);
scatter(reshape(X,numel(X), 1), reshape(Y, numel(Y),1), sz,'g', 'filled', 'MarkerFaceAlpha', alph);
hold on;
%plot(preNest(:,orTagTreat == 1,1), preNest(:,orTagTreat == 1,2), 'b');
plotbroodGreyTrans(broodPre,60, broodAlph);
xlim(xl);
ylim(yl);
axis equal
hold off

subplot(2,3,3);
X = preNest(:,orTagTreat == 3,1);
Y = preNest(:,orTagTreat == 3,2);
scatter(reshape(X,numel(X), 1), reshape(Y, numel(Y),1), sz,'r', 'filled', 'MarkerFaceAlpha', alph);
hold on;
%plot(preNest(:,orTagTreat == 1,1), preNest(:,orTagTreat == 1,2), 'b');
plotbroodGreyTrans(broodPre,60, broodAlph);
xlim(xl);
ylim(yl);
axis equal
hold off

subplot(2,3,4);
X = postNest(:,orTagTreat == 1,1);
Y = postNest(:,orTagTreat == 1,2);
scatter(reshape(X,numel(X), 1), reshape(Y, numel(Y),1), sz,'b', 'filled', 'MarkerFaceAlpha', alph);
hold on;
%plot(postNest(:,orTagTreat == 1,1), postNest(:,orTagTreat == 1,2), 'b');
plotbroodGreyTrans(broodPost,60, broodAlph);
xlim(xl);
ylim(yl);
axis equal
hold off

subplot(2,3,5);
X = postNest(:,orTagTreat == 2,1);
Y = postNest(:,orTagTreat == 2,2);
scatter(reshape(X,numel(X), 1), reshape(Y, numel(Y),1), sz,'g', 'filled', 'MarkerFaceAlpha', alph);
hold on;
%plot(postNest(:,orTagTreat == 1,1), postNest(:,orTagTreat == 1,2), 'b');
plotbroodGreyTrans(broodPost,60, broodAlph);
xlim(xl);
ylim(yl);
axis equal
hold off

subplot(2,3,6);
X = postNest(:,orTagTreat == 3,1);
Y = postNest(:,orTagTreat == 3,2);
scatter(reshape(X,numel(X), 1), reshape(Y, numel(Y),1), sz,'r', 'filled', 'MarkerFaceAlpha', alph);
hold on;
%plot(postNest(:,orTagTreat == 1,1), postNest(:,orTagTreat == 1,2), 'b');
plotbroodGreyTrans(broodPost,60, broodAlph);
xlim(xl);
ylim(yl);
axis equal
axis equal
hold off
%% Output simplified csv
masterData = array2table([tags orTagTreat framesActivePre' framesActivePost' framesDetectedPre' framesDetectedPost' activeVelPre' activeVelPost' meanVelPre' meanVelPost' socDistPre' socDistPost' queenDistPre' queenDistPost' interactionRatePre' interactionRatePost' nectarBoutsPost nectarBoutsPre pollenBoutsPost pollenBoutsPre framesOnNestPre framesOnNestPost velocityFramesPre' velocityFramesPost' medianVelPre' medianVelPost' totalBoutsPre totalBoutsPost]);
masterData.Properties.VariableNames = {'tags', 'treatment', 'framesActivePre', 'framesActivePost', 'framesDetectedPre', 'framesDetectedPost','activeVelPre', 'activeVelPost', 'meanVelPre', 'meanVelPost', 'socDistPre', 'socDistPost', 'queenDistPre', 'queenDistPost', 'interactionRatePre', 'interactionRatePost', 'nectarBoutsPost','nectarBoutsPre', 'pollenBoutsPost', 'pollenBoutsPre', 'framesOnNestPre', 'framesOnNestPost', 'numVelocityFramesPre', 'numVelocityFramesPost', 'medianVelPre', 'medianVelPost', 'totalForageBoutsPre', 'totalForageBoutsPost'};
writetable(masterData,'masterDataForRCS.csv');