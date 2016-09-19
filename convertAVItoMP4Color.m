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
    writeVideo(writer,imadjust(s(k).cdata,[;],[], .9));
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
'26July_2016-07-26-083042-0000.avi',
'26July_2016-07-26-083042-0001.avi',
'26July_2016-07-26-083042-0002.avi',
'26July_2016-07-26-083042-0003.avi',
'26July_2016-07-26-083042-0004.avi',
'26July_2016-07-26-083042-0005.avi',
'26July_2016-07-26-083042-0006.avi',
'26July_2016-07-26-083042-0007.avi',
'26July_2016-07-26-083042-0008.avi',
'26July_2016-07-26-083042-0009.avi',
'26July_2016-07-26-083042-0010.avi',
'26July_2016-07-26-083042-0011.avi',
'26July_2016-07-26-083042-0012.avi',
'26July_2016-07-26-083042-0013.avi',
'26July_2016-07-26-083042-0014.avi',
'26July_2016-07-26-083042-0015.avi',
'26July_2016-07-26-083042-0016.avi',


};

file = '\Users\Combes4\Desktop\Yellow.mp4'
%%
%%

file_list = {
'29July2016_2016-07-29-073219-0000.avi',
'29July2016_2016-07-29-073219-0001.avi',
'29July2016_2016-07-29-073219-0002.avi',
'29July2016_2016-07-29-073219-0003.avi',
'29July2016_2016-07-29-073219-0004.avi',
'29July2016_2016-07-29-073219-0005.avi',
'29July2016_2016-07-29-073219-0006.avi',
'29July2016_2016-07-29-073219-0007.avi',
'29July2016_2016-07-29-073219-0008.avi',
'29July2016_2016-07-29-073219-0009.avi',
'29July2016_2016-07-29-073219-0010.avi',
'29July2016_2016-07-29-073219-0011.avi',
'29July2016_2016-07-29-073219-0012.avi',
'29July2016_2016-07-29-073219-0013.avi',
'29July2016_2016-07-29-073219-0014.avi',
'29July2016_2016-07-29-073219-0015.avi',
'29July2016_2016-07-29-073219-0016.avi',
'29July2016_2016-07-29-073219-0017.avi'

};

file = '\Users\Combes4\Desktop\Turquoise_29July.mp4'
%%



% loop to combine a bunch of files and convert to .mp4
PathName = 'C:\Users\Combes4\Desktop\Turquoise_29July_TestTrain\Turquoise_29July2016_train_thresh0.5\Turquoise_train_UncompressedVideos_29July2016\';
cd(PathName);

writer = VideoWriter(file, 'MPEG-4');

writer.FrameRate = 20;

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
