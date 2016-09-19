%% Make movie into mp4 instead of avi to put into keynote
%% choose a file
[FileName,PathName,FilterIndex] = uigetfile
%%
%PathName = '/Users/callinswitzer/Desktop/ImidaclopridVids/'; 

%FileName = 'btag_2015-10-16-121953-0001.avi';

% where you want to save it
%file = '/Users/callinswitzer/Desktop/btag0.mp4' 
foo = strsplit(FileName, '.')

file = char(strcat(PathName, foo(1), '.mp4'));

vidObj = VideoReader(strcat(PathName, FileName));

vidObj.currentTime = 0; % start 381 seconds in

% get info
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

writer = VideoWriter(file, 'MPEG-4');

writer.FrameRate = vidObj.FrameRate;

open(writer);
% Read and write each frame.


% read in frames -- do 1000 frames
k = 1;
while hasFrame(vidObj)
    s(k).cdata = readFrame(vidObj);
    %
%     positions = [450 550; 450 565];
%     test_str = cell({'Concord Field Station' ,'Callin Switzer'});
%     
%     RGB = insertText(s(k).cdata, positions,...
%         test_str,'BoxOpacity',0, 'TextColor', 'white',...
%         'FontSize', 15);
%     
%     imshow(imadjust(RGB,[;],[], .9)); % adjust gamma
   
    writeVideo(writer,s(k).cdata);
    k = k+1
    
end
close(writer);

