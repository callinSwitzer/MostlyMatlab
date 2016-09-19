%% Make movie into mp4 instead of avi to put into keynote
%% choose a file
[FileName,PathName,FilterIndex] = uigetfile
%%
%PathName = '/Users/callinswitzer/Desktop/ImidaclopridVids/'; 

%FileName = 'btag_2015-10-16-121953-0000.avi';

% where you want to save it
%file = '/Users/callinswitzer/Desktop/btag0.mp4' 
foo = strsplit(FileName, '.')

file = char(strcat(PathName, foo(1), '.mp4'));

mm = VideoReader(strcat(PathName, FileName));

% Get Some info on the video 
NumberOfFrames  = mm.NumberOfFrames
Width           = mm.Width
Height          = mm.Height

% For speed, preallocate the array 
frames = uint8(zeros(Height, Width, NumberOfFrames));

% Load in the frames into the preallocated array
for kk=1:mm.NumberOfFrames
%for kk=2000:2600
    tmp = read(mm,kk);              
    frames(:,:,kk) = tmp(1:Height, 1:Width);  
    kk
end

% Save the video as .mp4
video = VideoWriter(file,'MPEG-4');
video.FrameRate = 30;  % playback speed
video.Quality = 100;     % quality


open(video);

cmap = gray(2);
display('saving video')
for kk=1:mm.NumberOfFrames
    tmp = frames(:,:,kk); 
    % adjust gamma
    writeVideo(video, imresize(tmp,   0.9)  ); % had to resize, b/c the video was too big
    kk
end

close(video)


