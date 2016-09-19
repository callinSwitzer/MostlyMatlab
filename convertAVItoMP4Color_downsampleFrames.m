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

% get info
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

writer = VideoWriter(file, 'MPEG-4');

writer.FrameRate = vidObj.FrameRate;


open(writer);
% Read and write each frame.


% read in frames
k = 1;
while hasFrame(vidObj)
    s(k).cdata = readFrame(vidObj);
    if mod(k, 10) == 0 % one frame in 10 is written (making new 100 fps)
        writeVideo(writer,imadjust(s(k).cdata,[;],[], .6));
    end
    k = k+1
end
close(writer);

%% file list for tagging
file_list = {
    'btag_2015-10-16-121953-0000.avi', 
    'btag_2015-10-16-121953-0001.avi', 
    'btag_2015-10-16-121953-0002.avi', 
    'btag_2015-10-16-121953-0003.avi',
    'btag_2015-10-16-121953-0004.avi',
    'btag_2015-10-16-121953-0005.avi',
    'btag_2015-10-16-121953-0006.avi',
    'btag_2015-10-16-121953-0007.avi',
    'btag_2015-10-16-121953-0008.avi',
    };

file = '/Users/callinswitzer/Desktop/ImidaclopridVids/tagging.mp4'

%%
file_list = {
'fc2_save_2015-10-13-162506-0000.avi',
'fc2_save_2015-10-13-162506-0001.avi',
'fc2_save_2015-10-13-162506-0002.avi',
'fc2_save_2015-10-13-162506-0003.avi',
'fc2_save_2015-10-13-162506-0004.avi',
'fc2_save_2015-10-13-162506-0005.avi',
'fc2_save_2015-10-13-162506-0006.avi',
'fc2_save_2015-10-13-162506-0007.avi',
'fc2_save_2015-10-13-162506-0008.avi', 
'fc2_save_2015-10-13-162506-0009.avi',
};

file = '/Users/callinswitzer/Desktop/ImidaclopridVids/starving.mp4'

%%

file_list = {
'fc2_save_2015-10-13-160716-0000.avi',
'fc2_save_2015-10-13-160716-0001.avi',
'fc2_save_2015-10-13-160716-0002.avi',
'fc2_save_2015-10-13-160716-0003.avi',
};

file = '/Users/callinswitzer/Desktop/ImidaclopridVids/feeding1.mp4'
%%

% loop to combine a bunch of files and convert to .mp4
PathName = '/Users/callinswitzer/Desktop/ImidaclopridVids/';

writer = VideoWriter(file, 'MPEG-4');

writer.FrameRate = 30;

open(writer);

for ii=1:length(file_list)
    
    FileName = char(file_list(ii));
    % where you want to save it
    %file = '/Users/callinswitzer/Desktop/btag0.mp4' 
    foo = strsplit(FileName, '.')

    %file = char(strcat(PathName, foo(1), '.mp4'));

    vidObj = VideoReader(strcat(PathName, FileName));

    % get info
    vidHeight = vidObj.Height;
    vidWidth = vidObj.Width;

    s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
        'colormap',[]);

    
    % Read and write each frame.


    % read in frames
    k = 1;
    while hasFrame(vidObj)
        s(k).cdata = readFrame(vidObj);
        writeVideo(writer,s(k).cdata);
        k = k+1
    end
    
end
close(writer);
