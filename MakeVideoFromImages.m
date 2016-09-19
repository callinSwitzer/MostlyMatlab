% Make video from images
% Callin Switzer
% 4/26/2016

[FileName,PathName,FilterIndex] = uigetfile
%% find image file names

imageNames = dir(fullfile(PathName,'*.tif'));
imageNames = {imageNames.name}';
length(imageNames)


%% Write video
outputVideo = VideoWriter(fullfile(PathName,'Rhodo_cadis_video.mp4'),'MPEG-4');
outputVideo.FrameRate = 30; % original vid at 0.5 fps
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(fullfile(PathName,imageNames{ii}));
   writeVideo(outputVideo,img)
   ii
end

close(outputVideo)
